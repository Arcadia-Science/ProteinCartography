# default column names for a Foldseek run in this pipeline
import enum
import os

FOLDSEEK_COLUMN_NAMES = [
    "query",
    "target",
    "fident",
    "alnlen",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "prob",
    "evalue",
    "bits",
    "qcov",
    "tcov",
    "qlan",
    "taln",
    "coord",
    "tseq",
    "taxid",
    "taxname",
]

FOLDSEEK_OUT_COLUMN_NAMES = ["protid", "fident", "prob", "evalue"]


class UniProtService(enum.Enum):
    """
    The services that can be used to query UniProt
    """

    # the REST APIs provided by UniProt itself
    REST = "rest"

    # the APIs provided by the `bioservices` package
    BIOSERVICES = "bioservices"


BLAST_OUTPUT_FIELDS = [
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
    "sacc",
    "saccver",
    "sgi",
    "staxids",
    "scomnames",
]

# the prepended '6' is how we specify that the output format is tabular
BLAST_OUTFMT = " ".join(["6"] + BLAST_OUTPUT_FIELDS)


def blast_soft_fail_enabled() -> bool:
    """True only for an explicit PC_BLAST_SOFT_FAIL opt-in. Default is hard-fail.

    Shared by ``run_blast`` and ``extract_blast_hits`` so an empty BLAST TSV
    cannot succeed one step and abort the next.
    """
    return os.environ.get("PC_BLAST_SOFT_FAIL", "0").strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }
