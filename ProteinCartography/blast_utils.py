import os
import pathlib
import shutil
import subprocess

# If this env variable is set to the path of a canned blastp results file, `run_blast` copies that
# file to its output path instead of calling `blastp` against the remote NCBI server.
# This exists because the real remote blastp call queues on NCBI's public servers and takes
# 10-40 minutes, which makes end-to-end verification of the pipeline impractical.
# It is intended for testing only (see the `smoke-test` make target); never set it in production.
#
# Note: an env variable is used (rather than a CLI flag on `run_blast.py`) because the snakemake
# rule environments inherit their env variables from the process that calls snakemake,
# so no changes to the Snakefile or the pipeline config are needed to activate the stub.
STUB_RESULTS_FILEPATH_ENV_VAR = "PROTEINCARTOGRAPHY_BLAST_STUB_RESULTS_FILEPATH"


def blast_call_failed(result: subprocess.CompletedProcess) -> bool:
    """
    Determine whether a call to `blastp` failed.

    The exit code alone is not sufficient: when the remote server refuses to queue a request,
    `blastp -remote` writes an error to stderr, writes an empty results file, and nevertheless
    exits with a status of zero. For example:

        Error: [blastp] bad_request: Could not queue request: DB operation failed.

    Checking only the exit code therefore treats these failures as successes, which prevents
    the word-size backoff in `run_blast.py` from ever being attempted and defers the failure
    to `extract_blast_hits.py`, where it appears as a misleading "no hits were returned" error.

    Note: an empty results file is deliberately *not* treated as a failure here, because a query
    that legitimately has no hits also produces one.

    Args:
        result (subprocess.CompletedProcess): the result of the call to `blastp`.

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
    Stand in for a call to `blastp` by copying a canned blastp results file to `out`.

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


def run_blast(
    query: str,
    out: str,
    max_target_seqs: int,
    outfmt: str,
    word_size: int,
    evalue: float,
):
    """
    Call `blastp` against the remote server

    Note: the argument names used here correspond exactly to (a subset of the) blastp CLI arguments

    Note: if the env variable named by `STUB_RESULTS_FILEPATH_ENV_VAR` is set,
    this function copies a canned results file instead of calling `blastp` (see `run_blast_stub`).

    Args:
        query (str): path of input peptide FASTA file.
        out (str): path of destination blastresults.tsv file.
        max_target_seqs (int): maximum number of hits to return
        outfmt (str): passed to blastp '-outfmt'
        word_size (str): passed to blastp '-word_size'
        evalue (str): passed to blastp '-evalue'
    """
    stub_results_filepath = os.environ.get(STUB_RESULTS_FILEPATH_ENV_VAR)
    if stub_results_filepath:
        return run_blast_stub(stub_results_filepath, out)

    database = "nr"
    result = subprocess.run(
        " ".join(
            [
                "blastp",
                "-remote",
                "-db",
                database,
                "-query",
                query,
                "-out",
                out,
                "-max_target_seqs",
                str(max_target_seqs),
                "-outfmt",
                f"'{outfmt}'",
                "-word_size",
                str(word_size),
                "-evalue",
                str(evalue),
            ]
        ),
        capture_output=True,
        shell=True,
    )
    return result
