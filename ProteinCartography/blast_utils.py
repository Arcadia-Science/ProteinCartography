import os
import pathlib
import re
import shutil
import subprocess
import time
from typing import Optional

import requests

# If this env variable is set to the path of a canned blastp results file, `run_blast` copies that
# file to its output path instead of querying the remote NCBI server.
# This exists because the real remote blastp search queues on NCBI's public servers and takes
# 10-40 minutes, which makes end-to-end verification of the pipeline impractical.
# It is intended for testing only (see the `smoke-test` make target); never set it in production.
#
# Note: an env variable is used (rather than a CLI flag on `run_blast.py`) because the snakemake
# rule environments inherit their env variables from the process that calls snakemake,
# so no changes to the Snakefile or the pipeline config are needed to activate the stub.
STUB_RESULTS_FILEPATH_ENV_VAR = "PROTEINCARTOGRAPHY_BLAST_STUB_RESULTS_FILEPATH"

# The env variable used to override the contact email address sent to NCBI (see `resolve_email`).
EMAIL_ENV_VAR = "PROTEINCARTOGRAPHY_BLAST_EMAIL"

# The URL of NCBI's BLAST URL API.
# See https://ncbi.github.io/blast-cloud/dev/api.html for the documentation of this API.
NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

# NCBI asks that automated clients identify themselves with a tool name and a contact email address
# so that they can get in touch before blocking a misbehaving client.
NCBI_TOOL_NAME = "proteincartography"

# The default contact email address sent to NCBI. This is the ProteinCartography repository's
# public issue tracker rather than any individual's address, so that NCBI has somewhere to
# report a problem without a personal address being baked into the source.
# Override it with the `blast_email` config parameter or the env variable named above.
DEFAULT_EMAIL = "proteincartography@arcadiascience.com"

# NCBI asks that the status of a search is polled no more than once every 60 seconds.
MIN_POLLING_INTERVAL_SECONDS = 60

# The timeout applied to each individual HTTP request. This bounds how long a single request can
# hang; the overall search is bounded separately by the caller-supplied `timeout_seconds`.
HTTP_TIMEOUT_SECONDS = 60

# The fields returned by NCBI's `Tabular` format, in the order in which they are returned.
# These are the same twelve fields that the `blastp` CLI emits for a bare `-outfmt 6`.
NCBI_TABULAR_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

# The value written for the requested output fields that NCBI's URL API does not return
# (see `tabular_row_to_output_fields`). `blastp` itself writes 'N/A' for fields it cannot populate.
UNAVAILABLE_FIELD_VALUE = "N/A"


class BlastApiError(Exception):
    """
    Raised when a search submitted to NCBI's BLAST URL API cannot be completed.
    """


def blast_call_failed(result: subprocess.CompletedProcess) -> bool:
    """
    Determine whether a call to blastp failed.

    The exit code alone is not sufficient: when the remote server refuses to queue a request,
    `blastp -remote` writes an error to stderr, writes an empty results file, and nevertheless
    exits with a status of zero. For example:

        Error: [blastp] bad_request: Could not queue request: DB operation failed.

    Checking only the exit code therefore treats these failures as successes, which prevents
    the word-size backoff in `run_blast.py` from ever being attempted and defers the failure
    to `extract_blast_hits.py`, where it appears as a misleading "no hits were returned" error.

    Note: `run_blast` no longer calls the `blastp` CLI, but it still reports its result as a
    `CompletedProcess` whose stderr is prefixed with "Error:" when the search failed, so that
    this check (and the retry loop in `run_blast.py` that depends on it) is unchanged.

    Note: an empty results file is deliberately *not* treated as a failure here, because a query
    that legitimately has no hits also produces one.

    Args:
        result (subprocess.CompletedProcess): the result of the attempted blastp search.

    Returns:
        True if the call failed.
    """
    if result.returncode != 0:
        return True

    stderr = result.stderr
    if isinstance(stderr, bytes):
        stderr = stderr.decode(errors="replace")

    return "Error:" in (stderr or "")


def run_blast_stub(stub_results_filepath: str, out: str):
    """
    Stand in for a blastp search by copying a canned blastp results file to `out`.

    Note: the canned file is copied verbatim, so it is not filtered by, or consistent with,
    the query sequence or any of the blastp parameters that `run_blast` is called with.

    Args:
        stub_results_filepath (str): path of the canned blastresults.tsv file to copy.
        out (str): path of the destination blastresults.tsv file.
    """
    stub_results_filepath = pathlib.Path(stub_results_filepath)
    if not stub_results_filepath.exists():
        raise FileNotFoundError(
            f"The blast stub results file '{stub_results_filepath}' does not exist "
            f"(it was set by the '{STUB_RESULTS_FILEPATH_ENV_VAR}' env variable)."
        )

    print(f"Using the stubbed blast results in '{stub_results_filepath}' instead of calling blastp")
    shutil.copy(stub_results_filepath, out)
    return subprocess.CompletedProcess(args=[], returncode=0, stdout=b"", stderr=b"")


def resolve_email(email: Optional[str] = None) -> str:
    """
    Determine the contact email address to send to NCBI.

    The env variable takes precedence over the value from the pipeline config, so that the address
    can be set without editing the config file (which matters when the config file is shared).

    Args:
        email (str): the address from the pipeline config, if any.

    Returns:
        The address to send to NCBI.
    """
    return os.environ.get(EMAIL_ENV_VAR) or email or DEFAULT_EMAIL


def output_field_names(outfmt: str) -> list[str]:
    """
    Parse the names of the output fields from a blastp `-outfmt` format string.

    The leading '6' that selects the tabular format is not a field name, so it is dropped.

    Args:
        outfmt (str): a blastp `-outfmt` format string (usually `constants.BLAST_OUTFMT`).

    Returns:
        The names of the requested output fields, in order.
    """
    return [name for name in outfmt.split(" ") if name and name != "6"]


def parse_accession(sseqid: str) -> str:
    """
    Parse the subject accession (with its version) out of a subject sequence id.

    NCBI sequence ids are pipe-delimited and put the accession in the second field, whether they
    have two fields ('ref|XP_052610122.1|') or three ('sp|P60709|ACTB_HUMAN'). Ids without any
    pipes are already bare accessions.

    Args:
        sseqid (str): a subject sequence id, as returned in NCBI's tabular output.

    Returns:
        The accession, including its version suffix if the id had one.
    """
    if "|" not in sseqid:
        return sseqid

    fields = sseqid.split("|")
    return fields[1] if len(fields) > 1 and fields[1] else sseqid


def strip_accession_version(saccver: str) -> str:
    """
    Strip the trailing version suffix from an accession ('XP_052610122.1' -> 'XP_052610122').
    """
    return re.sub(r"\.\d+$", "", saccver)


def tabular_row_to_output_fields(row: list[str], field_names: list[str]) -> list[str]:
    """
    Convert one row of NCBI's tabular output into the fields requested by the pipeline.

    NCBI's URL API does not accept a custom `-outfmt` field list; it only returns the twelve
    standard tabular fields. Of the five additional fields the pipeline requests, 'sacc' and
    'saccver' are derived from the subject sequence id, and the remaining three ('sgi', 'staxids'
    and 'scomnames') are not available from this API and are written as 'N/A'.

    Note: only 'sacc' is read downstream (by `extract_blast_hits.py`), so the unavailable fields
    do not affect the pipeline's results; they are emitted only to keep the column layout
    identical to the one produced by `blastp -outfmt 6`.

    Args:
        row (list): the tab-separated values of one row of NCBI's tabular output.
        field_names (list): the names of the fields to return, in order.

    Returns:
        The values of the requested fields, in order.
    """
    if len(row) < len(NCBI_TABULAR_FIELDS):
        raise BlastApiError(
            f"NCBI returned a tabular row with {len(row)} fields, "
            f"but at least {len(NCBI_TABULAR_FIELDS)} were expected: {row!r}"
        )

    values = dict(zip(NCBI_TABULAR_FIELDS, row))
    values["saccver"] = parse_accession(values["sseqid"])
    values["sacc"] = strip_accession_version(values["saccver"])

    return [values.get(field_name, UNAVAILABLE_FIELD_VALUE) for field_name in field_names]


def format_results(tabular_results: str, outfmt: str) -> str:
    """
    Convert NCBI's tabular output into the exact column layout that the pipeline expects.

    NCBI's tabular output interleaves comment lines (which start with '#') with the rows of hits;
    the comment lines are dropped, because `blastp -outfmt 6` does not emit them.

    Args:
        tabular_results (str): the body of NCBI's `FORMAT_TYPE=Tabular` response.
        outfmt (str): the blastp `-outfmt` format string naming the fields to write.

    Returns:
        The results as a tab-separated string with one hit per line.
    """
    field_names = output_field_names(outfmt)

    lines = []
    for line in tabular_results.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        row = line.rstrip("\n").split("\t")
        lines.append("\t".join(tabular_row_to_output_fields(row, field_names)))

    return "".join(f"{line}\n" for line in lines)


def parse_qblast_info(response_text: str, key: str) -> Optional[str]:
    """
    Parse a value out of the 'QBlastInfo' block that NCBI embeds in its responses.

    The block is formatted as, for example:

        <!--QBlastInfoBegin
            RID = M8K4ARMH013
            RTOE = 27
        QBlastInfoEnd-->

    Args:
        response_text (str): the body of a response from NCBI's BLAST URL API.
        key (str): the name of the value to parse (for example 'RID' or 'Status').

    Returns:
        The value, or None if it is not present.
    """
    match = re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)\s*$", response_text, re.MULTILINE)
    return match.group(1) if match else None


def summarize_response(response_text: str, max_length: int = 500) -> str:
    """
    Condense a response body into a single line, so that it can be included in an error message.

    NCBI reports some failures (notably "Temporarily unavailable") as a full HTML page, so the
    body is stripped of markup and truncated rather than being quoted verbatim.
    """
    text = re.sub(r"<[^>]+>", " ", response_text)
    text = " ".join(text.split())
    if len(text) > max_length:
        text = f"{text[:max_length]}..."
    return text


def request(params: dict) -> str:
    """
    Make a request to NCBI's BLAST URL API and return the body of the response.

    POST is used rather than GET because the query sequence can easily exceed the length of a URL.

    Args:
        params (dict): the query parameters to send.

    Returns:
        The body of the response.

    Raises:
        BlastApiError: if the request failed or NCBI returned an error status.
    """
    params = dict(params, tool=NCBI_TOOL_NAME)
    try:
        response = requests.post(NCBI_BLAST_URL, data=params, timeout=HTTP_TIMEOUT_SECONDS)
    except requests.RequestException as exception:
        raise BlastApiError(f"The request to NCBI failed: {exception}") from exception

    if response.status_code != 200:
        raise BlastApiError(
            f"NCBI returned a status of {response.status_code} "
            f"for a '{params.get('CMD')}' request: {summarize_response(response.text)}"
        )

    return response.text


def submit_search(
    query_sequence: str,
    max_target_seqs: int,
    word_size: int,
    evalue: float,
    email: str,
    database: str = "nr",
) -> tuple[str, int]:
    """
    Submit a blastp search to NCBI and return the id of the queued search.

    Args:
        query_sequence (str): the contents of the query FASTA file.
        max_target_seqs (int): the maximum number of hits to return.
        word_size (int): the length of the exact matches used to seed alignments.
        evalue (float): the maximum e-value of the hits to return.
        email (str): the contact email address to send to NCBI.
        database (str): the name of the NCBI database to search.

    Returns:
        A tuple of the request id (RID) and NCBI's estimate, in seconds, of how long the search
        will take (RTOE).

    Raises:
        BlastApiError: if NCBI did not queue the search.
    """
    response_text = request(
        {
            "CMD": "Put",
            "PROGRAM": "blastp",
            "DATABASE": database,
            "QUERY": query_sequence,
            "HITLIST_SIZE": max_target_seqs,
            "WORD_SIZE": word_size,
            "EXPECT": evalue,
            "email": email,
        }
    )

    request_id = parse_qblast_info(response_text, "RID")
    if request_id is None:
        raise BlastApiError(
            "NCBI did not return a request id (RID), so the search was not queued. "
            f"NCBI responded with: {summarize_response(response_text)}"
        )

    # NCBI's time estimate is advisory and is occasionally missing or unparseable,
    # in which case the first poll happens after the minimum polling interval instead.
    raw_time_estimate = parse_qblast_info(response_text, "RTOE")
    try:
        time_estimate = int(raw_time_estimate)
    except (TypeError, ValueError):
        time_estimate = MIN_POLLING_INTERVAL_SECONDS

    return request_id, time_estimate


def poll_search_status(request_id: str, email: str) -> str:
    """
    Ask NCBI for the status of a queued search.

    Args:
        request_id (str): the id (RID) of the search.
        email (str): the contact email address to send to NCBI.

    Returns:
        The status, which is one of 'WAITING', 'READY', 'FAILED' or 'UNKNOWN'.

    Raises:
        BlastApiError: if NCBI did not return a status at all.
    """
    response_text = request(
        {
            "CMD": "Get",
            "RID": request_id,
            "FORMAT_OBJECT": "SearchInfo",
            "email": email,
        }
    )

    status = parse_qblast_info(response_text, "Status")
    if status is None:
        raise BlastApiError(
            f"NCBI did not report a status for the search '{request_id}'. "
            f"NCBI responded with: {summarize_response(response_text)}"
        )

    return status


def wait_for_search(request_id: str, time_estimate: int, timeout_seconds: float, email: str):
    """
    Poll NCBI until a queued search is ready, or until the timeout is exceeded.

    The first poll happens after NCBI's own time estimate (RTOE) has elapsed, and subsequent polls
    happen no more often than once every 60 seconds, as NCBI asks. Both waits are shortened when
    less than that much of the timeout remains, so that the timeout is honored exactly.

    Args:
        request_id (str): the id (RID) of the search.
        time_estimate (int): NCBI's estimate (RTOE), in seconds, of how long the search will take.
        timeout_seconds (float): how long to wait, in total, before giving up.
        email (str): the contact email address to send to NCBI.

    Raises:
        BlastApiError: if the search failed, expired, or did not finish before the timeout.
    """
    started_at = time.monotonic()
    deadline = started_at + timeout_seconds
    delay = max(time_estimate, 0)

    while True:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            break

        time.sleep(min(delay, remaining))
        if time.monotonic() >= deadline:
            break

        status = poll_search_status(request_id, email)
        elapsed = round(time.monotonic() - started_at)
        print(f"The blast search '{request_id}' is {status} after {elapsed}s")

        if status == "READY":
            return
        if status == "FAILED":
            raise BlastApiError(
                f"NCBI reported that the blast search '{request_id}' failed. "
                "This usually indicates a problem on NCBI's side; if it persists, "
                "report the request id to blast-help@ncbi.nlm.nih.gov."
            )
        if status == "UNKNOWN":
            raise BlastApiError(
                f"NCBI no longer recognizes the blast search '{request_id}'. "
                "Request ids expire after 24 hours, and are also reported as unknown "
                "when the search was never queued."
            )
        if status != "WAITING":
            raise BlastApiError(
                f"NCBI reported the unexpected status '{status}' "
                f"for the blast search '{request_id}'."
            )

        delay = MIN_POLLING_INTERVAL_SECONDS

    raise BlastApiError(
        f"The blast search '{request_id}' did not finish within the timeout of "
        f"{round(timeout_seconds)}s. NCBI estimated that it would take {time_estimate}s "
        f"({round(time_estimate / 60)} minutes) when the search was queued. "
        "NCBI's public queue is shared, so this usually means the queue is busy rather than "
        "that anything is wrong with the query; either retry later or raise the "
        "'blast_timeout_seconds' config parameter."
    )


def fetch_results(request_id: str, max_target_seqs: int, email: str) -> str:
    """
    Retrieve the results of a completed search in NCBI's tabular format.

    Args:
        request_id (str): the id (RID) of the search.
        max_target_seqs (int): the maximum number of hits to return.
        email (str): the contact email address to send to NCBI.

    Returns:
        The body of NCBI's tabular response.

    Raises:
        BlastApiError: if NCBI returned something other than tabular results.
    """
    response_text = request(
        {
            "CMD": "Get",
            "RID": request_id,
            "FORMAT_TYPE": "Tabular",
            "ALIGNMENTS": max_target_seqs,
            "DESCRIPTIONS": max_target_seqs,
            "HITLIST_SIZE": max_target_seqs,
            "email": email,
        }
    )

    # A tabular response always carries a comment block; anything else is an error page.
    # (A search with no hits still returns the comment block, and must not be treated as an error.)
    if not any(line.startswith("#") for line in response_text.splitlines()):
        raise BlastApiError(
            f"NCBI did not return tabular results for the blast search '{request_id}'. "
            f"NCBI responded with: {summarize_response(response_text)}"
        )

    return response_text


def run_blast(
    query: str,
    out: str,
    max_target_seqs: int,
    outfmt: str,
    word_size: int,
    evalue: float,
    timeout_seconds: float,
    email: Optional[str] = None,
):
    """
    Run a blastp search against NCBI's `nr` database using NCBI's BLAST URL API.

    This uses the URL API directly rather than `blastp -remote`, because `blastp -remote` polls
    NCBI's queue indefinitely with no way for the caller to bound how long it waits: when NCBI's
    public queue is busy, it blocks for hours while producing no output and consuming no CPU.
    Driving the same API directly lets the pipeline give up after `timeout_seconds` and report
    why it gave up.

    Note: the argument names used here correspond exactly to (a subset of the) blastp CLI arguments

    Note: if the env variable named by `STUB_RESULTS_FILEPATH_ENV_VAR` is set,
    this function copies a canned results file instead of searching (see `run_blast_stub`).

    Note: the result is reported as a `CompletedProcess` for the benefit of `blast_call_failed`
    and the retry loop in `run_blast.py`, both of which predate this function calling an API
    rather than a subprocess.

    Args:
        query (str): path of input peptide FASTA file.
        out (str): path of destination blastresults.tsv file.
        max_target_seqs (int): maximum number of hits to return
        outfmt (str): the field layout to write, in the notation of blastp '-outfmt'
        word_size (int): passed to NCBI as 'WORD_SIZE'
        evalue (float): passed to NCBI as 'EXPECT'
        timeout_seconds (float): how long to wait for NCBI to complete the search before failing.
        email (str): the contact email address to send to NCBI (see `resolve_email`).

    Returns:
        A `CompletedProcess` whose return code is zero if the search succeeded.
    """
    stub_results_filepath = os.environ.get(STUB_RESULTS_FILEPATH_ENV_VAR)
    if stub_results_filepath:
        return run_blast_stub(stub_results_filepath, out)

    email = resolve_email(email)

    try:
        with open(query) as file:
            query_sequence = file.read()

        request_id, time_estimate = submit_search(
            query_sequence=query_sequence,
            max_target_seqs=max_target_seqs,
            word_size=word_size,
            evalue=evalue,
            email=email,
        )
        print(
            f"NCBI queued the blast search as '{request_id}' "
            f"and estimated that it will take {time_estimate}s"
        )

        wait_for_search(
            request_id=request_id,
            time_estimate=time_estimate,
            timeout_seconds=timeout_seconds,
            email=email,
        )
        tabular_results = fetch_results(
            request_id=request_id, max_target_seqs=max_target_seqs, email=email
        )
    except (BlastApiError, OSError) as exception:
        return subprocess.CompletedProcess(
            args=[], returncode=1, stdout=b"", stderr=f"Error: [blastp] {exception}"
        )

    with open(out, "w") as file:
        file.write(format_results(tabular_results, outfmt))

    return subprocess.CompletedProcess(args=[], returncode=0, stdout=b"", stderr="")
