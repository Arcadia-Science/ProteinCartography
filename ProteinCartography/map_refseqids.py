#!/usr/bin/env python
import argparse
import os
import sys
from time import sleep

import pandas as pd
from api_utils import (
    UniProtWithExpBackoff,
    session_with_retry,
)
from tests import api_mocks

# if necessary, mock the `uniprot.mapping` method (used by `map_refseqids_bioservices`)
# see comments in `tests.api_mocks` for more details
if os.environ.get("PROTEINCARTOGRAPHY_WAS_CALLED_BY_PYTEST") == "true":
    api_mocks.mock_bioservices_uniprot_mapping()

# only import these functions when using import *
__all__ = ["map_refseqids_bioservices", "map_refseqids_rest"]

# check through these default databases
DEFAULT_DBS = ["EMBL-GenBank-DDBJ_CDS", "RefSeq_Protein"]

# id mapping link
UNIPROT_IDMAPPING_API = "https://rest.uniprot.org/idmapping"

# requests constants
REQUESTS_LIMIT = 10
REQUESTS_SLEEP_TIME = 30


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="path to input .txt file containing one accession per line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="path to .txt file where uniquely-mapped Uniprot accessions will be printed.",
    )
    parser.add_argument(
        "-d",
        "--databases",
        nargs="+",
        default=DEFAULT_DBS,
        help=f"which databases to use for mapping. defaults to {DEFAULT_DBS}",
    )
    parser.add_argument("-s", "--service", default="rest", help="how to fetch mapping")
    args = parser.parse_args()
    return args


# takes a list of IDs and maps them to Uniprot using bioservices
# might make a more generalizable version of this and put it somewhere else
def map_refseqids_bioservices(
    input_file: str, output_file: str, query_dbs: list, return_full=False
):
    """
    Takes an input .txt file of accessions and maps to UniProt accessions.

    Args:
        input_file (str): path to input .txt file containing one accession per line.
        output_file (str): path to destination .txt file.
        query_dbs (list): list of valid databases to query using the Uniprot ID mapping API.
            Each database will be queried individually.
            The results are compiled and unique results are printed to output_file.
    """

    # make object that references UniProt database
    uniprot = UniProtWithExpBackoff()

    # open the input file to extract ids
    with open(input_file) as f:
        ids = f.read().splitlines()

    # limit the number of ids to prevent the uniprot mapping API from timing out
    max_num_ids = 100000
    ids = ids[:max_num_ids]

    # make an empty collector dataframe for mapping
    dummy_df = pd.DataFrame()

    # for each query database, map
    for i, db in enumerate(query_dbs):
        # uniprot.mapping returns a gross json file
        results = uniprot.mapping(db, "UniProtKB", query=",".join(ids))

        # pandas can normalize the json and make it more tractable
        results_df = pd.json_normalize(results["results"])

        # if there are no results, move on
        if len(results_df) == 0:
            continue

        # if it's the first database, replace it with the dummy dataframe
        if i == 0:
            dummy_df = results_df
        # otherwise append to the dataframe
        else:
            dummy_df = pd.concat([dummy_df, results_df], axis=0)

    # extract just the unique Uniprot accessions
    hits = dummy_df["to.primaryAccession"].unique()

    # save those accessions to a .txt file
    with open(output_file, "w+") as f:
        f.writelines(hit + "\n" for hit in hits)

    if return_full:
        return dummy_df


# Example curl POST request
# ```
# % curl --request POST 'https://rest.uniprot.org/idmapping/run' \
#   --form 'ids="P21802,P12345"' \
#   --form 'from="UniProtKB_AC-ID"' \
#   --form 'to="UniRef90"'
# ```


def map_refseqids_rest(input_file: str, output_file: str, query_dbs: list, return_full=False):
    """
    Takes an input .txt file of accessions and maps to UniProt accessions.

    Args:
        input_file (str): path to input .txt file containing one accession per line.
        output_file (str): path to destination .txt file.
        query_dbs (list): list of valid databases to query using the Uniprot ID mapping API.
            Each database will be queried individually.
            The results are compiled and unique results are printed to output_file.
        return_full (bool): whether to return all of the results as a dataframe
    """
    # open the input file to extract ids
    with open(input_file) as f:
        input_lines = f.read().splitlines()
        input_ids = list(set(input_lines))
        input_string = ",".join(input_ids)

    dummy_df = pd.DataFrame()

    for i, db in enumerate(query_dbs):
        ticket = (
            session_with_retry()
            .post(
                f"{UNIPROT_IDMAPPING_API}/run",
                {"ids": input_string, "from": db, "to": "UniProtKB"},
            )
            .json()
        )

        # poll until the job was successful or failed
        repeat = True
        tries = 0
        while repeat and tries < REQUESTS_LIMIT:
            status = (
                session_with_retry()
                .get(
                    f'{UNIPROT_IDMAPPING_API}/status/{ticket["jobId"]}',
                )
                .json()
            )

            # wait a short time between poll requests
            sleep(REQUESTS_SLEEP_TIME)
            tries += 1
            repeat = "results" not in status

        if tries == 10:
            sys.exit(f"The ticket failed to complete after {tries * REQUESTS_SLEEP_TIME} seconds.")

        results = (
            session_with_retry().get(f'{UNIPROT_IDMAPPING_API}/stream/{ticket["jobId"]}').json()
        )
        results_df = pd.DataFrame(results["results"])

        # if there are no results, move on
        if len(results_df) == 0:
            continue

        # if it's the first database, replace it with the dummy dataframe
        if i == 0:
            dummy_df = results_df
        # otherwise append to the dataframe
        else:
            dummy_df = pd.concat([dummy_df, results_df], axis=0)

    # extract just the unique Uniprot accessions
    hits = dummy_df["to"].unique()

    # save those accessions to a .txt file
    with open(output_file, "w+") as f:
        f.writelines(hit + "\n" for hit in hits)

    if return_full:
        return dummy_df


# run this if called from the interpreter
def main():
    # parse arguments
    args = parse_args()

    # collect arguments individually
    input_file = args.input
    output_file = args.output
    query_dbs = args.databases
    service = args.service

    if service == "bioservices":
        map_refseqids_bioservices(input_file, output_file, query_dbs)
    else:
        map_refseqids_rest(input_file, output_file, query_dbs)


# check if called from interpreter
if __name__ == "__main__":
    main()
