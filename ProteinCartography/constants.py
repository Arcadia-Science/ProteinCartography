# default column names for a Foldseek run in this pipeline
import enum

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


class UniProtServices(enum.Enum):
    """
    The services that can be used to query UniProt
    """

    # the REST APIs provided by UniProt itself
    REST = "rest"

    # the APIs provided by the `bioservices` package
    BIOSERVICES = "bioservices"
