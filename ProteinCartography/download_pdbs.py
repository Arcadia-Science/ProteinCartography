#!/usr/bin/env python
import argparse
import concurrent.futures
from pathlib import Path

import api_utils
import fetch_accession
import tqdm
from ratelimiter import RateLimiter


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input file path of a .txt file with one accession per line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output directory in which to save the PDB files.",
    )
    parser.add_argument(
        "-M",
        "--max-structures",
        type=int,
        required=False,
        help="Maximum number of PDB files to download.",
    )
    args = parser.parse_args()
    return args


def download_pdbs(input_file: str, output_dir: str, maximum=None):
    """
    Download PDBs for the accessions listed in `input_file` from AlphaFold.

    Args:
        input_file (str): path to an text file containing one accession per line.
        output_dir (str): path to output directory in which to save the PDB files.
        maximum (int): maximum number of accessions to download. If None, downloads all.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    with open(input_file) as file:
        accessions = file.read().splitlines()

    if maximum is not None:
        accessions = accessions[:maximum]

    session = api_utils.session_with_retry()
    rate_limiter = RateLimiter(max_calls=100, period=1)
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
        futures_to_accessions = {}
        for accession in accessions:
            future = executor.submit(
                rate_limiter(fetch_accession.fetch_pdb),
                accession=accession,
                output_dir=output_dir,
                session=session,
            )
            futures_to_accessions[future] = accession

        for future in tqdm.tqdm(
            concurrent.futures.as_completed(futures_to_accessions),
            total=len(futures_to_accessions),
            desc="Downloading PDBs from AlphaFold",
        ):
            try:
                future.result()
            except Exception as exception:
                print(f"Error fetching PDB '{futures_to_accessions[future]}': {exception}")


def main():
    args = parse_args()
    download_pdbs(input_file=args.input, output_dir=args.output, maximum=args.max_structures)


if __name__ == "__main__":
    main()
