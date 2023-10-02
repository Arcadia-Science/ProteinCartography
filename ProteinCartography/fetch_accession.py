#!/usr/bin/env python
import argparse
import os
from pathlib import Path
import subprocess

from api_utils import UniProtWithExpBackoff, USER_AGENT_HEADER


# only import these functions when using import *
__all__ = ["fetch_fasta", "fetch_pdb"]


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--accession",
        required=True,
        nargs="+",
        help="UniprotKB accession of target.",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output directory for resulting files."
    )
    parser.add_argument(
        "-f",
        "--format",
        nargs="+",
        default=["fasta", "pdb"],
        help='Formats to acquire.\nIf "fasta", requests the FASTA sequence using bioservices UniProt.\nIf "pdb", downloads the pdb from AlphaFold.\nCan accept multiple arguments.',
    )
    args = parser.parse_args()
    return args


def fetch_fasta(accession: str, output_dir: str):
    """
    Fetches a FASTA file from Uniprot, given an accession. Places the file in the output_dir.

    Args:
        accession (str): a valid UniprotKB accession.
        output_dir (str): path to the output directory. File will be saved as "{output_dir}/{accession}.fasta".
    """
    u = UniProtWithExpBackoff()
    output_path = Path(output_dir) / (accession + ".fasta")

    if not os.path.exists(output_path):
        res = u.retrieve(accession, frmt="fasta")
        with open(output_path, "w+") as f:
            f.write(res)


def fetch_pdb(accession: str, output_dir: str):
    """
    Fetches a PDB file from AlphaFold, given an accession. Places the file in the output_dir.

    Args:
        accession (str): a valid UniprotKB accession.
        output_dir (str): path to the output directory. File will be saved as "{output_dir}/{accession}.pdb".
    """
    output_path = Path(output_dir) / (accession + ".pdb")
    source = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb".format(accession)

    if not os.path.exists(output_path):
        user_agent = USER_AGENT_HEADER["User-Agent"]
        subprocess.run(
            ["curl", "-JLo", output_path, source, "--user-agent", f"'{user_agent}'"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


# run this if called from the interpreter
def main():
    args = parse_args()
    accessions = args.accession
    output_dir = args.output

    formats = args.format

    for accession in accessions:
        if "fasta" in formats:
            fetch_fasta(accession, output_dir)
        if "pdb" in formats:
            fetch_pdb(accession, output_dir)


# check if called from interpreter
if __name__ == "__main__":
    main()
