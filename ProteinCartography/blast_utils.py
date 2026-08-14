import json
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

# The format in which the results are retrieved from NCBI.
#
# `Tabular` would need no conversion, but NCBI does not actually deliver it: fetching a completed
# search with `FORMAT_TYPE=Tabular` returns a 62-byte body containing the status block and no hits,
# with or without `FORMAT_OBJECT=Alignment`. Probing one completed search (RID 7TJ2P2DR014) gave:
#
#     Tabular  : 62 bytes (nothing but the status block)
#     Text     : 157,935 bytes
#     XML2_S   : 511,934 bytes
#     JSON2_S  : 422,458 bytes
#
# JSON2_S is used because it is the machine-readable format that carries the most of what the
# pipeline's output fields need, including the subject accession and taxonomy id.
NCBI_RESULTS_FORMAT = "JSON2_S"

# The value written for the output fields that cannot be populated from NCBI's response
# (see `hsp_to_output_fields`). `blastp` itself writes 'N/A' for fields it cannot populate.
UNAVAILABLE_FIELD_VALUE = "N/A"

# The value written for the 'sgi' field. NCBI has retired GI numbers, and `blastp -outfmt 6`
# writes a literal 0 for them rather than 'N/A' (as the committed blastresults.tsv artifact shows).
RETIRED_GI_VALUE = "0"


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


def format_evalue(evalue: float) -> str:
    """
    Format an e-value the way `blastp` formats it in its tabular output.

    The thresholds and precisions reproduce NCBI's own `ScoreAndEvalueToBuffers`, so that the
    results file looks like one written by `blastp -outfmt 6` rather than like Python's repr
    of a float (which would render, for example, 0.0 as '0' and 1e-05 as '1e-05').
    """
    evalue = float(evalue)
    if evalue < 1.0e-180:
        return "0.0"
    if evalue < 1.0e-99:
        return f"{evalue:2.0e}".strip()
    if evalue < 0.0009:
        return f"{evalue:3.0e}".strip()
    if evalue < 0.1:
        return f"{evalue:4.3f}".strip()
    if evalue < 1.0:
        return f"{evalue:3.2f}".strip()
    if evalue < 10.0:
        return f"{evalue:2.1f}".strip()
    return f"{evalue:5.0f}".strip()


def format_bitscore(bitscore: float) -> str:
    """
    Format a bit score the way `blastp` formats it in its tabular output.

    As with `format_evalue`, the thresholds reproduce NCBI's own formatting: bit scores above
    99.9 are written as integers, which is why the committed blastresults.tsv artifact contains
    '787' rather than '787.334'.
    """
    bitscore = float(bitscore)
    if bitscore > 9999:
        return f"{bitscore:4.3e}".strip()
    if bitscore > 99.9:
        return f"{bitscore:.0f}"
    return f"{bitscore:.1f}"


def count_gap_openings(hsp: dict) -> int:
    """
    Count the gap openings in an alignment, which is what the 'gapopen' output field reports.

    NCBI's JSON reports `gaps`, the total number of gapped *positions*, which is not the same
    number: a single five-residue gap is five gapped positions but one gap opening. The openings
    are therefore counted directly, as the number of runs of '-' in the two aligned sequences.

    Args:
        hsp (dict): one HSP from NCBI's JSON response.

    Returns:
        The number of gap openings.
    """
    return sum(len(re.findall(r"-+", hsp.get(key) or "")) for key in ("qseq", "hseq"))


def query_seqid(search: dict) -> str:
    """
    Determine the query sequence id that `blastp` would report in its 'qseqid' output field.

    NCBI assigns its own internal id (for example 'Query_5706473') to a sequence submitted
    through the URL API and keeps the original FASTA defline in `query_title`, so the id the
    pipeline expects is the first whitespace-delimited token of the title rather than `query_id`.

    Args:
        search (dict): the `search` object from NCBI's JSON response.

    Returns:
        The query sequence id.
    """
    query_title = search.get("query_title")
    if query_title and query_title.split():
        return query_title.split()[0]
    return search.get("query_id") or UNAVAILABLE_FIELD_VALUE


def hsp_to_output_fields(qseqid: str, description: dict, hsp: dict) -> dict:
    """
    Convert one HSP of one hit into the output fields that `blastp -outfmt 6` would write.

    Most fields are read straight out of NCBI's JSON. The ones that are not:

      - 'pident' is computed as the fraction of identities over the alignment length, and
        'mismatch' as the alignment length less the identities and the gapped positions,
        which are the definitions `blastp` uses.
      - 'gapopen' is counted from the aligned sequences (see `count_gap_openings`), because
        NCBI's JSON reports total gapped positions rather than gap openings.
      - 'sgi' is written as 0, matching `blastp`, because NCBI has retired GI numbers.
      - 'scomnames' is written as 'N/A'. NCBI's JSON carries `sciname`, but that is a scientific
        name and this field is defined as the subject's *common* name, so filling it with the
        scientific name would mislabel the data. (The committed artifact from a real
        `blastp -remote` run also has 'N/A' throughout this column.)

    Note: only 'sacc' is read downstream, by `extract_blast_hits.py`. Every other field is
    written solely to keep the column layout identical to `blastp -outfmt 6`, so the
    approximations above cannot affect the pipeline's results.

    Args:
        qseqid (str): the query sequence id (see `query_seqid`).
        description (dict): the first description of the hit (see `format_results`).
        hsp (dict): the HSP to convert.

    Returns:
        The output fields, keyed by the names used in `constants.BLAST_OUTPUT_FIELDS`.
    """
    align_len = hsp.get("align_len") or 0
    identity = hsp.get("identity") or 0
    gaps = hsp.get("gaps") or 0

    sseqid = description.get("id") or UNAVAILABLE_FIELD_VALUE
    saccver = parse_accession(sseqid)
    sacc = description.get("accession") or strip_accession_version(saccver)
    taxid = description.get("taxid")

    return {
        "qseqid": qseqid,
        "sseqid": sseqid,
        "pident": f"{(100.0 * identity / align_len) if align_len else 0.0:.3f}",
        "length": str(align_len),
        "mismatch": str(align_len - identity - gaps),
        "gapopen": str(count_gap_openings(hsp)),
        "qstart": str(hsp.get("query_from", "")),
        "qend": str(hsp.get("query_to", "")),
        "sstart": str(hsp.get("hit_from", "")),
        "send": str(hsp.get("hit_to", "")),
        "evalue": format_evalue(hsp.get("evalue") or 0),
        "bitscore": format_bitscore(hsp.get("bit_score") or 0),
        "sacc": sacc,
        "saccver": saccver,
        "sgi": RETIRED_GI_VALUE,
        "staxids": str(taxid) if taxid is not None else UNAVAILABLE_FIELD_VALUE,
        "scomnames": UNAVAILABLE_FIELD_VALUE,
    }


def iter_searches(results: dict):
    """
    Iterate over the per-query search results in NCBI's JSON response.

    The response holds one report per query, so a multi-sequence FASTA file produces several.

    Args:
        results (dict): the parsed body of NCBI's `FORMAT_TYPE=JSON2_S` response.

    Yields:
        The `search` object of each report.

    Raises:
        BlastApiError: if the response is not a BLAST JSON report at all.
    """
    reports = results.get("BlastOutput2")
    if reports is None:
        raise BlastApiError(
            "NCBI returned JSON that is not a blast report (it has no 'BlastOutput2' key)."
        )

    if isinstance(reports, dict):
        reports = [reports]

    for report in reports:
        search = report.get("report", {}).get("results", {}).get("search")
        if search is not None:
            yield search


def format_results(results: dict, outfmt: str) -> str:
    """
    Convert NCBI's JSON response into the exact column layout that the pipeline expects.

    One row is written per HSP, and only the *first* description of each hit is used. This
    matches `blastp -outfmt 6`: when `nr` merges identical sequences, the hit carries one
    description per redundant accession, and the tabular 'sacc' field reports only the
    representative one (reporting all of them is what the separate 'sallacc' field is for).
    This was checked against the committed results file from a real `blastp -remote` run: each
    of its rows corresponds to one hit, including a hit whose sequence is shared by 1053
    accessions, and in every case the accession it reports is the first description's.

    Args:
        results (dict): the parsed body of NCBI's `FORMAT_TYPE=JSON2_S` response.
        outfmt (str): the blastp `-outfmt` format string naming the fields to write.

    Returns:
        The results as a tab-separated string with one HSP per line.
    """
    field_names = output_field_names(outfmt)

    lines = []
    for search in iter_searches(results):
        qseqid = query_seqid(search)
        for hit in search.get("hits") or []:
            descriptions = hit.get("description") or [{}]
            for hsp in hit.get("hsps") or []:
                values = hsp_to_output_fields(qseqid, descriptions[0], hsp)
                lines.append(
                    "\t".join(
                        values.get(field_name, UNAVAILABLE_FIELD_VALUE)
                        for field_name in field_names
                    )
                )

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

        # The status is always checked after waiting, even when the wait used up the rest of the
        # budget. NCBI's estimate is routinely larger than the whole timeout, and it is only an
        # estimate: a search whose estimate was hours may well be ready. Returning without ever
        # asking would report a timeout for a search that had already finished.
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
                "Request ids expire after 36 hours, and are also reported as unknown "
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


def fetch_results(request_id: str, max_target_seqs: int, email: str) -> dict:
    """
    Retrieve the results of a completed search (see `NCBI_RESULTS_FORMAT` for the format used).

    Args:
        request_id (str): the id (RID) of the search.
        max_target_seqs (int): the maximum number of hits to return.
        email (str): the contact email address to send to NCBI.

    Returns:
        The parsed body of NCBI's response.

    Raises:
        BlastApiError: if NCBI returned something other than the expected results.
    """
    response_text = request(
        {
            "CMD": "Get",
            "RID": request_id,
            "FORMAT_TYPE": NCBI_RESULTS_FORMAT,
            "ALIGNMENTS": max_target_seqs,
            "DESCRIPTIONS": max_target_seqs,
            "HITLIST_SIZE": max_target_seqs,
            "email": email,
        }
    )

    # An error page, or the status-block-only body that NCBI returns for some format types,
    # is not JSON and must not be mistaken for a search that simply had no hits.
    try:
        return json.loads(response_text)
    except ValueError as exception:
        raise BlastApiError(
            f"NCBI did not return {NCBI_RESULTS_FORMAT} results for the blast search "
            f"'{request_id}' ({exception}). "
            f"NCBI responded with: {summarize_response(response_text)}"
        ) from exception


def run_blast(
    query: str,
    out: str,
    max_target_seqs: int,
    outfmt: str,
    word_size: int,
    evalue: float,
    timeout_seconds: float,
    email: Optional[str] = None,
    database: str = "nr",
):
    """
    Run a blastp search against an NCBI database using NCBI's BLAST URL API.

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
        database (str): the name of the NCBI database to search. Changing this changes which
            homologs are found, so it is a change to the results and not only to the runtime.

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
            database=database,
        )
        # The request id is printed as soon as it exists, because NCBI keeps the results of a
        # search for 36 hours: if the pipeline gives up, or the process dies, the search can
        # still be retrieved from this id rather than being resubmitted to the queue.
        print(
            f"NCBI queued the blast search of '{database}' as '{request_id}' "
            f"and estimated that it will take {time_estimate}s"
        )

        wait_for_search(
            request_id=request_id,
            time_estimate=time_estimate,
            timeout_seconds=timeout_seconds,
            email=email,
        )
        results = fetch_results(request_id=request_id, max_target_seqs=max_target_seqs, email=email)
        formatted_results = format_results(results, outfmt)
    except (BlastApiError, OSError) as exception:
        return subprocess.CompletedProcess(
            args=[], returncode=1, stdout=b"", stderr=f"Error: [blastp] {exception}"
        )

    with open(out, "w") as file:
        file.write(formatted_results)

    return subprocess.CompletedProcess(args=[], returncode=0, stdout=b"", stderr="")
