#!/usr/bin/env python
import argparse
import os
import re

import numpy as np
import pandas as pd
from api_utils import UniProtWithExpBackoff, session_with_retry
from constants import UniProtService
from tests import mocks

# if necessary, mock the `uniprot.search` method (used by `query_uniprot`)
# see comments in `tests.mocks` for more details
if os.environ.get("PROTEINCARTOGRAPHY_SHOULD_USE_MOCKS") == "true":
    mocks.mock_bioservices_uniprot_search()

# only import these functions when using import *
__all__ = ["query_uniprot"]

# global fields for querying using REST API
REQUIRED_FIELDS_DICT = {
    "Entry": "accession",
    "Entry Name": "id",
    "Protein names": "protein_name",
    "Gene Names (primary)": "gene_primary",
    "Annotation": "annotation_score",
    "Organism": "organism_name",
    "Taxonomic lineage": "lineage",
    "Length": "length",
    "Fragment": "fragment",
    "Sequence": "sequence",
}
OTHER_FIELDS_DICT = {
    "Reviewed": "reviewed",
    "Gene Names": "gene_names",
    "Protein existence": "protein_existence",
    "Sequence version": "sequence_version",
    "RefSeq": "xref_refseq",
    "GeneID": "xref_geneid",
    "EMBL": "xref_embl",
    "AlphaFoldDB": "xref_alphafolddb",
    "PDB": "xref_pdb",
    "Pfam": "xref_pfam",
    "InterPro": "xref_interpro",
}
DEFAULT_FIELDS_DICT = REQUIRED_FIELDS_DICT | OTHER_FIELDS_DICT
REQUIRED_FIELDS = list(REQUIRED_FIELDS_DICT.values())
DEFAULT_FIELDS = list(DEFAULT_FIELDS_DICT.values())


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="path of input merged hits file")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="path of destination uniprot_features.tsv file",
    )
    parser.add_argument(
        "-s", "--service", default=UniProtService.BIOSERVICES.value, help="how to fetch metadata"
    )
    parser.add_argument(
        "-a",
        "--additional-fields",
        nargs="*",
        help="additional non-default fields to fetch from uniprot if using REST",
    )
    args = parser.parse_args()

    return args


def load_existing_data(temp_filepath: str) -> tuple[set, bool]:
    """
    Load any existing results from previous incomplete runs, to remove them from the query list.

    Args:
        temp_filepath (str): Path to the temporary file containing previously downloaded entries.

    Returns:
        tuple:
        - set: Set of existing entries (accessions) found in the file.
        - bool: True if header was written (file exists and has data), False otherwise.
    """

    existing_data = set()
    header_written = False

    if not os.path.exists(temp_filepath):
        print(f"{temp_filepath} does not exist. Skipping existing data loading.")
        return existing_data, header_written

    try:
        existing_df = pd.read_csv(temp_filepath, sep="\t", usecols=["Entry"])
        if "Entry" not in existing_df.columns:
            print(f"No 'Entry' column found in {temp_filepath}. Skipping existing data loading.")
            return existing_data, header_written

        existing_data = set(existing_df["Entry"].values)
        print(f"Loaded {len(existing_data)} entries from {temp_filepath}")
        header_written = True

    except pd.errors.EmptyDataError:
        print(f"{temp_filepath} is empty. Skipping existing data loading.")

    return existing_data, header_written


def query_uniprot(
    input_file: str,
    output_file: str,
    batch_size=100,
    sub_batch_size=100,
    fmt="tsv",
    fields=DEFAULT_FIELDS,
    service=UniProtService.BIOSERVICES,
):
    """
    Takes an input list of accessions and gets the full information set from Uniprot
    for those proteins.

    Args:
        input_file (str): path to a text file of accessions (with one accession per line)
        output_file (str): path of destination tsv file with all uniprot features
        batch_size (int): number of entries to query per batch.
        sub_batch_size (int): number of entries to pull per page of batch.
        fmt (str): output suffix format (default 'tsv')
        fields (list): list of UniProt fields to retrieve.
        service (str): which API to use: 'rest' or 'bioservices'.

    Returns:
        a pandas.DataFrame of the resulting features
    """

    temp_filepath = output_file + ".temp"

    if os.path.exists(output_file):
        print("Output file already exists at this location. Aborting.")
        return

    existing_data, header_written = load_existing_data(temp_filepath)

    # Define regular expression pattern to extract the next URL from the response headers
    re_next_link = re.compile(r'<(.+)>; rel="next"')

    session = session_with_retry()

    def get_next_link(headers):
        # Extract the next URL from the "Link" header if present
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    fields_string = ",".join(fields)

    def get_batch(query_string: str):
        # Generator function to fetch data in batches
        if service == UniProtService.REST:
            batch_url = f"https://rest.uniprot.org/uniprotkb/search?query={query_string}&format={fmt}&fields={fields_string}&size={sub_batch_size}"
            while batch_url:
                response = session.get(batch_url)
                response.raise_for_status()
                yield response.text
                batch_url = get_next_link(response.headers)
        elif service == UniProtService.BIOSERVICES:
            uniprot = UniProtWithExpBackoff()
            results = uniprot.search(
                query_string, columns=fields_string, size=sub_batch_size, progress=False
            )
            # bioservices doesn't directly return the number of results.
            yield results
        else:
            raise ValueError(f"Unknown service {service}")

    with open(input_file) as q:
        query_accessions = [line.rstrip("\n") for line in q.readlines()]
        # Remove accessions that have already been downloaded, keeping the order.
        query_accessions = [e for e in query_accessions if e not in existing_data]

    # Split the query accessions into batches
    accession_batches = [
        query_accessions[i : i + batch_size] for i in range(0, len(query_accessions), batch_size)
    ]

    # Process each batch separately
    for i, accession_batch in enumerate(accession_batches):
        print(f">> Starting batch {i + 1} of {len(accession_batches)}")

        # Construct the query string for the batch
        query_string = " OR ".join(f"accession:{accession}" for accession in accession_batch)
        query_string = f"({query_string})"

        # Construct the URL with the constructed query string
        progress = 0
        accessions_in_batch = len(accession_batch)

        with open(temp_filepath, "a") as temp_file:
            for batch in get_batch(query_string):
                lines = batch.splitlines()

                if not header_written:
                    print_lines = lines
                    header_written = True
                else:
                    print_lines = lines[1:]

                for line in print_lines:
                    print(line, file=temp_file)

                progress += len(lines) - 1
                print(f"downloaded {progress} / {accessions_in_batch} hits for batch {i + 1}")

    df = pd.read_csv(temp_filepath, sep="\t")
    df.insert(0, "protid", df["Entry"].values)

    def lineage_string_splitter(lineage_string):
        if lineage_string is np.nan:
            return np.nan
        else:
            return [rank.split(" (")[0] for rank in lineage_string.split(", ")]

    try:
        df["Lineage"] = df["Taxonomic lineage"].apply(lineage_string_splitter)
    except Exception:
        pass

    df.to_csv(output_file, sep="\t", index=None)
    os.remove(temp_filepath)

    return df


# run this if called from the interpreter
def main():
    args = parse_args()

    input_file = args.input
    output_file = args.output
    service = UniProtService(args.service)
    additional_fields = args.additional_fields

    fields = DEFAULT_FIELDS
    if additional_fields is not None:
        fields += additional_fields

    query_uniprot(input_file, output_file, fields=fields, service=service)


# check if called from interpreter
if __name__ == "__main__":
    main()
