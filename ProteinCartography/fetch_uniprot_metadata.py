#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import os
import re

from api_utils import session_with_retry, UniProtWithExpBackoff


# only import these functions when using import *
__all__ = ["query_uniprot_bioservices", "query_uniprot_rest"]

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
    parser.add_argument(
        "-i", "--input", required=True, help="path of input merged hits file"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="path of destination uniprot_features.tsv file",
    )
    parser.add_argument("-s", "--service", default="rest", help="how to fetch mapping")
    parser.add_argument(
        "-a",
        "--additional-fields",
        nargs="?",
        default=[],
        help="additional non-default fields to fetch from uniprot if using REST",
    )
    args = parser.parse_args()

    return args


#############################################
## Query using bioservices (not preferred) ##
#############################################


def query_uniprot_bioservices(
    input_file: str, output_file: str, save=True
) -> pd.DataFrame:
    """
    Takes an input list of accessions and gets the full information set from Uniprot for those proteins.

    Args:
        input_file (str): path of input list text file where each accession is on a new line
        output_file (str): path of destination tsv file with all uniprot features

    Returns:
        a pandas.DataFrame of the resulting features
    """

    # open and read the file
    with open(input_file, "r") as f:
        id_list = [i.rstrip("\n") for i in f.readlines()]

    # perform ID mapping using bioservices UniProt
    # should probably do this differently in the future because it often results in weird memory leaks
    u = UniProtWithExpBackoff()
    results = u.mapping("UniProtKB_AC-ID", "UniProtKB", query=id_list, progress=False)

    # read the results as a normalized json
    results_df = pd.json_normalize(results["results"])
    # remove the "to" prefix for tidier columns
    results_df.columns = results_df.columns.str.removeprefix("to.")
    # add a protid column for later merging
    results_df.insert(0, "protid", results_df["primaryAccession"])

    # save to file if needed
    if save:
        results_df.to_csv(output_file, index=None, sep="\t")

    return results_df


######################################
## Query using REST API (preferred) ##
######################################


def query_uniprot_rest(
    query_list: str,
    output_file: str,
    batch_size=300,
    sub_batch_size=300,
    fmt="tsv",
    fields=DEFAULT_FIELDS,
):
    """
    Takes an input list of accessions and gets the full information set from Uniprot for those proteins.

    Args:
        input_file (str): path of input list text file where each accession is on a new line
        output_file (str): path of destination tsv file with all uniprot features
        batch_size (int): number of entries to query per batch.
        sub_batch_size (int): number of entries to pull per page of batch.
        fmt (str): output suffix format (default 'tsv')
        fields (list): list of UniProt fields to retrieve.

    Returns:
        a pandas.DataFrame of the resulting features
    """

    temp_file = output_file + ".temp"

    if os.path.exists(output_file):
        print("Output file already exists at this location. Aborting.")
        return

    header_written = False  # Flag to check if header has been written

    # Load any existing results from previous incomplete runs, to remove them from the query list.
    existing_data = set()
    if os.path.exists(temp_file):
        existing_df = pd.read_csv(temp_file, sep="\t", usecols=["Entry"])
        existing_data = set(existing_df["Entry"].values)
        print(f"Loaded {len(existing_data)} entries from {temp_file}")
        header_written = True

    # Define regular expression pattern to extract the next URL from the response headers
    re_next_link = re.compile(r'<(.+)>; rel="next"')

    session = session_with_retry()

    def get_next_link(headers):
        # Extract the next URL from the "Link" header if present
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(batch_url):
        # Generator function to fetch data in batches
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)

    fields_string = ",".join(fields)

    with open(query_list, "r") as q:
        query_accessions = [line.rstrip("\n") for line in q.readlines()]
        # Remove accessions that have already been downloaded, keeping the order.
        query_accessions = [e for e in query_accessions if e not in existing_data]

    # Split the query accessions into batches
    accession_batches = [
        query_accessions[i : i + batch_size]
        for i in range(0, len(query_accessions), batch_size)
    ]

    # Process each batch separately
    for i, accession_batch in enumerate(accession_batches):
        print(f">> Starting batch {i + 1} of {len(accession_batches)}")

        # Construct the query string for the batch
        query_string = " OR ".join(
            f"accession:{accession}" for accession in accession_batch
        )
        query_string = f"({query_string})"

        # Construct the URL with the constructed query string
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query_string}&format={fmt}&fields={fields_string}&size={sub_batch_size}"
        progress = 0

        with open(temp_file, "a") as f:
            for batch, total in get_batch(url):
                lines = batch.text.splitlines()

                if not header_written:
                    print_lines = lines
                    header_written = True
                else:
                    print_lines = lines[1:]

                for line in print_lines:
                    print(line, file=f)

                progress += len(lines) - 1
                print(f"downloaded {progress} / {total} hits for batch {i + 1}")

    df = pd.read_csv(temp_file, sep="\t")

    lineage_string_splitter = (
        lambda lineage_string: [
            rank.split(" (")[0] for rank in lineage_string.split(", ")
        ]
        if lineage_string is not np.nan
        else np.nan
    )

    df.insert(0, "protid", df["Entry"].values)

    try:
        df["Lineage"] = df["Taxonomic lineage"].apply(lineage_string_splitter)
    except:
        pass

    df.to_csv(output_file, sep="\t", index=None)
    os.remove(temp_file)

    return df


# run this if called from the interpreter
def main():
    args = parse_args()

    input_file = args.input
    output_file = args.output
    service = args.service
    additional_fields = args.additional_fields

    if additional_fields is not None:
        fields = DEFAULT_FIELDS + additional_fields
    else:
        fields = DEFAULT_FIELDS

    if service == "bioservices":
        query_uniprot_bioservices(input_file, output_file)
    else:
        query_uniprot_rest(input_file, output_file, fields=fields)


# check if called from interpreter
if __name__ == "__main__":
    main()
