import subprocess


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

    Note: this function is in its own module, rather than in `run_blast.py`,
    so that it can be mocked by calling `tests.mocks.mock_run_blast` at the top of `run_blast.py`

    Args:
        query (str): path of input peptide FASTA file.
        out (str): path of destination blastresults.tsv file.
        max_target_seqs (int): maximum number of hits to return
        outfmt (str): passed to blastp '-outfmt'
        word_size (str): passed to blastp '-word_size'
        evalue (str): passed to blastp '-evalue'
    """
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
    )
    return result
