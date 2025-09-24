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
from constants import UniProtService
from tests import mocks

# If necessary, mock the `uniprot.mapping` method (used by `map_refseqids_bioservices`).
# See comments in `tests.mocks` for more details.
if os.environ.get("PROTEINCARTOGRAPHY_WAS_CALLED_BY_PYTEST") == "true":
    mocks.mock_bioservices_uniprot_mapping()

__all__ = ["map_refseqids_bioservices", "map_refseqids_rest"]

# The default databases to query.
DEFAULT_DBS = ["EMBL-GenBank-DDBJ_CDS", "RefSeq_Protein"]

UNIPROT_IDMAPPING_API_URL = "https://rest.uniprot.org/idmapping"

# Time in seconds to wait between job status checks.
SLEEP_TIME_BETWEEN_JOB_STATUS_CHECKS = 30

# Allow a generous number of status checks because sometimes the UniProt ID Mapping API is slow.
MAX_NUM_JOB_STATUS_CHECKS = 100


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="File path to the input text file containing one accession per line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="File path to the output text file of uniquely-mapped Uniprot accessions.",
    )
    parser.add_argument(
        "-d",
        "--databases",
        nargs="+",
        default=DEFAULT_DBS,
        help=f"The databases to use for mapping. Defaults to {DEFAULT_DBS}",
    )
    parser.add_argument(
        "-s", "--service", default=UniProtService.REST.value, help="How to fetch the mapping."
    )
    args = parser.parse_args()
    return args


def map_refseqids_bioservices(
    input_file: str, output_file: str, query_dbs: list, return_full=False
):
    """
    Takes an input .txt file of accessions and maps to UniProt accessions using `bioservices`.

    Args:
        input_file (str): path to input .txt file containing one accession per line.
        output_file (str): path to destination .txt file.
        query_dbs (list): list of valid databases to query using the Uniprot ID mapping API.
            Each database will be queried individually.
            The results are compiled and unique results are printed to output_file.
    """

    uniprot = UniProtWithExpBackoff()

    with open(input_file) as f:
        ids = f.read().splitlines()

    # Limit the number of ids to prevent the uniprot mapping API from timing out.
    max_num_ids = 100000
    ids = ids[:max_num_ids]

    all_results_df = pd.DataFrame()

    for query_db in query_dbs:
        results = uniprot.mapping(query_db, "UniProtKB", query=",".join(ids))
        results_df = pd.json_normalize(results["results"])

        # If there are no results, skip to the next database.
        if len(results_df) == 0:
            continue

        all_results_df = pd.concat([all_results_df, results_df], axis=0)

    # Extract just the unique Uniprot accessions.
    hits = all_results_df["to.primaryAccession"].unique()

    # Save those accessions to a .txt file.
    with open(output_file, "w+") as f:
        f.writelines(hit + "\n" for hit in hits)

    if return_full:
        return all_results_df


def map_refseqids_rest(input_file: str, output_file: str, query_dbs: list, return_full=False):
    """
    Takes an input .txt file of accessions and maps to UniProt accessions
    using the Uniprot ID mapping REST API.

    Args:
        input_file (str): path to input .txt file containing one accession per line.
        output_file (str): path to destination .txt file.
        query_dbs (list): list of valid databases to query using the Uniprot ID mapping API.
            Each database will be queried individually.
            The results are compiled and unique results are printed to output_file.
        return_full (bool): whether to return all of the results as a dataframe

    For reference, here is an example curl POST request using the REST API:
    ```
    curl --request POST 'https://rest.uniprot.org/idmapping/run' \
      --form 'ids="P21802,P12345"' \
      --form 'from="UniProtKB_AC-ID"' \
      --form 'to="UniRef90"'
    ```
    """
    with open(input_file) as f:
        input_lines = f.read().splitlines()
        input_ids = list(set(input_lines))
        input_string = ",".join(input_ids)

    all_results_df = pd.DataFrame()

    for query_db in query_dbs:
        print(f"Mapping database '{query_db}' to UniProtKB")

        submission_response = (
            session_with_retry()
            .post(
                f"{UNIPROT_IDMAPPING_API_URL}/run",
                {"ids": input_string, "from": query_db, "to": "UniProtKB"},
            )
            .json()
        )

        print(f"Submission response for mapping request from {query_db} to UniProtKB:")
        print(submission_response)

        job_id = submission_response["jobId"]

        # Poll until the job is successful or we time out.
        num_tries = 0
        while num_tries < MAX_NUM_JOB_STATUS_CHECKS:
            status_response = (
                session_with_retry().get(f"{UNIPROT_IDMAPPING_API_URL}/status/{job_id}").json()
            )
            num_tries += 1

            # The status response is supposed to include a "jobStatus" key, but anecdotally,
            # this is only true while the job is running. Once the job is complete,
            # there is instead a "results" key.
            job_is_complete = (
                status_response.get("results") is not None
                or status_response.get("jobStatus").lower() == "finished"
            )

            if job_is_complete:
                print(f"Job is complete after {num_tries * SLEEP_TIME_BETWEEN_JOB_STATUS_CHECKS}s")
                break

            sleep(SLEEP_TIME_BETWEEN_JOB_STATUS_CHECKS)

        if num_tries == MAX_NUM_JOB_STATUS_CHECKS:
            sys.exit(
                f"The mapping request failed to complete after "
                f"{num_tries * SLEEP_TIME_BETWEEN_JOB_STATUS_CHECKS}s."
            )

        results_response = (
            session_with_retry().get(f"{UNIPROT_IDMAPPING_API_URL}/stream/{job_id}").json()
        )

        results_df = pd.DataFrame(results_response["results"])
        all_results_df = pd.concat([all_results_df, results_df], axis=0)

    # Extract just the unique Uniprot accessions and save them to a .txt file.
    hits = all_results_df["to"].unique()
    with open(output_file, "w+") as f:
        f.writelines(hit + "\n" for hit in hits)

    if return_full:
        return all_results_df


def main():
    args = parse_args()

    service = UniProtService(args.service)
    if service == UniProtService.BIOSERVICES:
        map_refseqids_bioservices(
            input_file=args.input, output_file=args.output, query_dbs=args.databases
        )
    elif service == UniProtService.REST:
        map_refseqids_rest(input_file=args.input, output_file=args.output, query_dbs=args.databases)


if __name__ == "__main__":
    main()
