#!/usr/bin/env python
import argparse
import os
import sys
import warnings

# depends on api_utils.py
from api_utils import session_with_retry
from Bio import SeqIO
from requests.packages.urllib3.exceptions import InsecureRequestWarning

### NOTES
# ESMFold API example from website:
#
# ```
# data=KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL # noqa: E501
# curl -X POST --data $data https://api.esmatlas.com/foldSequence/v1/pdb/
# ```

# only import these functions when using import *
__all__ = ["post_esmfold_apiquery", "esmfold_apiquery"]

# set acceptable fasta format suffixes
FASTA_FORMATS = ["fa", "fna", "fasta", "faa", "ffa"]


# parse command line arguments
def parse_args():
    # Set command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=(
            "Name of input file. Must be a single-entry peptide FASTA file "
            "(ends with .fa, .fna, .fasta)."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        help='Name of output file. If not provided, replaces ".fasta" from input with ".pdb".',
    )
    args = parser.parse_args()
    return args


def post_esmfold_apiquery(fasta: str):
    """
    Posts a query to the ESMfold API.

    Args:
        fasta (str): string of valid amino acids for query.
    """
    # Here we're making the request with verify=False to disable SSL
    # This may be a security concern, but is (hopefully)
    # temporary until the ESM Atlas SSL certificates are fixed
    result = session_with_retry().post(
        "https://api.esmatlas.com/foldSequence/v1/pdb/", data=fasta, verify=False
    )
    if result.status_code == 200:
        return result.text
    else:
        print(f"Error: {result.status_code}")
        return None


def esmfold_apiquery(input_file: str, output_file=None):
    """
    Takes an input peptide FASTA file with one entry
    and submits it for folding using the ESMfold API.
    Creates a PDB file with the resulting output.

    Args:
        input_file (str): path to input FASTA file
            (must end with .fa, .fna, or .fasta; not case-sensitive).
        output_file (str): path to the destination PDB file.
    """
    # Check to make sure input file has a correct FASTA suffix
    if not any([fmt in input_file.rsplit(".", 1)[1].lower() for fmt in FASTA_FORMATS]):
        sys.exit(f"Input expects a FASTA file ({FASTA_FORMATS})")

    # Makes sure that the input file exists
    if not os.path.exists(input_file):
        sys.exit(f"File {input_file} not found.")

    # if output file is not provided, generate a name for output file
    if output_file is None:
        output_filepath = input_file.replace(input_file.rsplit(".", 1)[1], "pdb")
    else:
        output_filepath = output_file

    # parse input_file and assign it to the variable records
    records = list(SeqIO.parse(input_file, "fasta"))

    if len(records) > 1:
        raise Exception(
            "This script expects a single FASTA entry in the input file. Please try again."
        )

    record = records[0]
    fasta = str(record.seq)

    prot_len = len(fasta)
    if prot_len > 400:
        print(
            f"The input protein is {prot_len} AA long.\n"
            f"ESMFold API query only allows proteins up to 400 AA long.\n"
            f"Try using ColabFold instead.\n"
            f"Skipping..."
        )
        return

    # submit a new job via the API
    result = post_esmfold_apiquery(fasta)

    if result is not None:
        with open(output_filepath, "w+") as file:
            file.write(result)


# run this if called from the interpreter
def main():
    # parse args
    args = parse_args()

    input_file = args.input
    output_file = args.output

    # Ignore warnings when we make the ESMFold API request
    # We're disabling SSL which `requests` will otherwise warn us about
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=InsecureRequestWarning)
        esmfold_apiquery(input_file, output_file)


# check if called from interpreter
if __name__ == "__main__":
    main()
