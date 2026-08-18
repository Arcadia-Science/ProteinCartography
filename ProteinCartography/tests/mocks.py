import re
import shutil
import subprocess
from unittest import mock

import file_utils
import requests

# TODO (KC): eliminate the hard-coded dataset name ("actin"),
# which needs to match the dataset name used in the `stage_inputs` pytest fixture.
ARTIFACTS_DIRPATH = (
    file_utils.find_repo_dirpath()
    / "ProteinCartography"
    / "tests"
    / "integration-test-artifacts"
    / "search-mode"
    / "actin"
)

API_RESPONSE_ARTIFACTS_DIRPATH = ARTIFACTS_DIRPATH / "api_response_content"


def mock_run_blast():
    """
    Mock the `run_blast.run_blast` method to prevent it from calling `blastp`
    by copying a mock blastresults.tsv file to the output directory.
    """

    result = mock.Mock(spec=subprocess.CompletedProcess)
    result.returncode = 0

    def side_effect(out=None, **_):
        shutil.copy(ARTIFACTS_DIRPATH / "output" / "P60709.blastresults.tsv", out)
        return result

    patch = mock.patch("blast_utils.run_blast", side_effect=side_effect)
    patch.start()


# Distinct domain-path neighborhoods used by the two-domain Snakemake test.
# These accessions have AlphaFold artifacts. Protein-path mocks still return the
# full actin hit list, so a parent-PDB search that only relabels outputs cannot
# produce this partitioned union.
DOMAIN_SEARCH_HITS = {
    "d01": ["A0A286Q506"],
    "d02": ["Q6QAQ1"],
}


def maybe_write_per_domain_hits(output_file: str) -> bool:
    """If ``output_file`` is a domain-path hit list, write that domain's neighborhood.

    Returns True when the caller should skip normal extraction.
    """
    name = output_file.replace("\\", "/")
    for suffix, hits in DOMAIN_SEARCH_HITS.items():
        token = f"__{suffix}."
        if token in name or name.endswith(f"__{suffix}"):
            from pathlib import Path

            path = Path(output_file)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("".join(hit + "\n" for hit in hits))
            return True
    return False


def mock_requests_session_request():
    """
    Mock the `request` method of `requests.Session` to return mock responses
    (constructed in `mock_response` below) instead of making real API calls.
    """
    patch = mock.patch("requests.Session.request", side_effect=mock_response)
    # we start the patch but never stop it, because we want to mock all requests
    # for the duration of the current python process
    patch.start()


def mock_bioservices_uniprot_search():
    """
    Mock the `search` method of `bioservices.UniProt` to prevent it from making real API calls.
    Note: because `UniProt.search` queries the UniProtKB REST API, this patch can reuse the response
    constructed in `mock_uniprotkb_rest_api_responses`.
    """
    patch = mock.patch(
        "api_utils.UniProtWithExpBackoff.search",
        return_value=mock_uniprotkb_rest_api_responses().text,
    )
    patch.start()


def mock_bioservices_uniprot_mapping():
    """
    Mock the `mapping` method of `bioservices.UniProt` to prevent it from making real API calls.
    Note: because `UniProt.mapping` queries the UniProt ID-Mapping API,
    this patch can reuse the response constructed in `mock_uniprot_id_mapping_api_responses`.
    """
    response = mock_uniprot_id_mapping_api_responses("GET", "stream/0")

    # restructure the response to match the format returned by uniprot.mapping
    response_json = response.json()
    for row in response_json["results"]:
        row["to"] = {"primaryAccession": row["to"]}

    patch = mock.patch("api_utils.UniProtWithExpBackoff.mapping", return_value=response_json)
    patch.start()


def mock_response(method, url, **_):
    """
    Return a mock response for a given method and url.
    """

    print(f"Mocking the response to: {method} {url}")

    # Requests to the UniProt ID Mapping REST API.
    if url.startswith("https://rest.uniprot.org/idmapping"):
        return mock_uniprot_id_mapping_api_responses(method, url)

    # Requests to the Foldseek API.
    elif url.startswith("https://search.foldseek.com/api"):
        return mock_foldseek_api_responses(method, url)

    # Requests to the UniProtKB REST API (search *or* /uniprotkb/accessions).
    elif url.startswith("https://rest.uniprot.org/uniprotkb/search") or url.startswith(
        "https://rest.uniprot.org/uniprotkb/accessions"
    ):
        return mock_uniprotkb_rest_api_responses()

    # AlphaFold prediction API (isoform-correct pdbUrl, then files download).
    elif url.startswith("https://alphafold.ebi.ac.uk/api/prediction"):
        return mock_alphafold_prediction_api_responses(url)

    # Requests to the alphafold files API.
    elif url.startswith("https://alphafold.ebi.ac.uk/files"):
        return mock_alphafold_files_api_responses(url)

    # TED domain summaries (query gate and domain-path hit assignment).
    elif "ted.cathdb.info" in url:
        return mock_ted_api_responses(url)

    else:
        raise ValueError(f"Unexpected url: {url}")


def mock_uniprot_id_mapping_api_responses(method, url):
    """
    Mock the responses to calls made to the UniProt ID Mapping REST API
    (these are made by `map_refseqids.map_refseqids_rest`).
    """
    mock_response = mock.Mock(spec=requests.Response)
    mock_response.status_code = 200

    job_id = "0"
    payload = {}

    # The initial POST request.
    if method == "POST" and url.endswith("run"):
        payload = {"jobId": job_id}

    # The subsequent polling requests.
    # ('success' is defined by the presence of a 'results' key in the response.)
    elif url.endswith(f"status/{job_id}"):
        payload = {"results": "some-results"}

    # The request to get the results.
    # Note: this payload is manually aggregated from the results of real API calls
    # to both of the default databases ("EMBL-GenBank-DDBJ_CDS" and "RefSeq_Protein").
    elif url.endswith(f"stream/{job_id}"):
        payload = {
            "results": [
                {"from": "XP_007129366", "to": "A0A2Y9FRR4"},
                {"from": "NP_001009784", "to": "P60713"},
                {"from": "NP_001009784", "to": "D7RIF5"},
                {"from": "XP_026547382", "to": "A0A6J1VWC1"},
                {"from": "KAF0882893", "to": "A0A6G1B5T4"},
                {"from": "NWI03924", "to": "A0A850ZFV5"},
                {"from": "TEA41296", "to": "A0A484H1H1"},
                {"from": "AAS55927", "to": "Q6QAQ1"},
                {"from": "KAF6447643", "to": "A0A7J8FIQ0"},
                {"from": "KAF6447654", "to": "A0A7J8FIS9"},
                {"from": "RLW01512", "to": "A0A3L8SFX2"},
                {"from": "KAF6480625", "to": "A0A7J8I8Z0"},
                {"from": "BAD96645", "to": "Q53GK6"},
                {"from": "KAF6081813", "to": "A0A833YM98"},
            ]
        }

    else:
        raise ValueError(f"Unexpected url: {url}")

    mock_response.json.return_value = payload
    return mock_response


def mock_foldseek_api_responses(method, url):
    """
    Mock the responses to calls made to the Foldseek API (by `foldseek_apiquery`).

    Note: this makes no attempt to parse the POST request body;
    it just returns the response from a real API call.
    """

    mock_response = mock.Mock(spec=requests.Response)
    mock_response.status_code = 200

    # The initial POST request.
    job_id = "0"
    if method == "POST" and url.endswith("api/ticket"):
        mock_response.json.return_value = {"id": job_id, "status": "COMPLETE"}

    # The polling request.
    # (This request is made even if the response to the initial POST request is "COMPLETE".)
    elif url.endswith(f"api/ticket/{job_id}"):
        mock_response.json.return_value = {"id": job_id, "status": "COMPLETE"}

    # The request to get the results.
    elif url.endswith(f"api/result/download/{job_id}"):
        with open(
            API_RESPONSE_ARTIFACTS_DIRPATH / "search.foldseek.com_api_result_download", "rb"
        ) as file:
            content = file.read()

        # This allows the response content to be streamed.
        mock_response.iter_content.return_value = [content]

    else:
        raise ValueError(f"Unexpected url: {url}")

    return mock_response


def mock_ted_api_responses(url):
    """TED summary API. Default: one domain (gate off). P99999: two domains (gate on)."""
    mock_response = mock.Mock(spec=requests.Response)
    accession = url.rstrip("/").split("/")[-1].split("?")[0]
    if accession == "P99999":
        mock_response.status_code = 200
        mock_response.ok = True
        mock_response.json.return_value = {
            "data": [
                {
                    "ted_id": "AF-P99999-F1-model_v4_TED01",
                    "uniprot_acc": "P99999",
                    "chopping": "1-80",
                    "nres_domain": 80,
                    "cath_label": "3.40.50.300",
                },
                {
                    "ted_id": "AF-P99999-F1-model_v4_TED02",
                    "uniprot_acc": "P99999",
                    "chopping": "81-160",
                    "nres_domain": 80,
                    "cath_label": "1.10.10.10",
                },
            ],
            "count": 2,
        }
        return mock_response

    mock_response.status_code = 200
    mock_response.ok = True
    mock_response.json.return_value = {
        "data": [
            {
                "ted_id": f"AF-{accession}-F1-model_v4_TED01",
                "uniprot_acc": accession,
                "chopping": "1-100",
                "nres_domain": 100,
                "cath_label": "1.10.10.10",
            }
        ],
        "count": 1,
    }
    return mock_response


def mock_uniprotkb_rest_api_responses():
    """
    Mock the response to calls made to the UniProtKB REST API
    (these are made by `fetch_uniprot_metadata.query_uniprot`).

    Note: this makes no attempt to parse the query string; it just returns a manually curated
    response from a real API call.
    """

    mock_response = mock.Mock(spec=requests.Response)
    mock_response.status_code = 200
    mock_response.ok = True

    # Define an empty header, specifically one without a 'Link' key
    # to prevent `fetch_uniprot_metadata.query_uniprot` from requesting a second batch of results.
    mock_response.headers = {}

    with open(
        API_RESPONSE_ARTIFACTS_DIRPATH / "rest.uniprot.org_uniprotkb_search", encoding="utf-8"
    ) as file:
        mock_response.text = file.read()

    return mock_response


def mock_alphafold_prediction_api_responses(url):
    """Return a pdbUrl that the files mock can serve (canonical F1 v6 path)."""
    accession = url.rstrip("/").split("/")[-1].split("?")[0]
    mock_response = mock.Mock(spec=requests.Response)
    mock_response.status_code = 200
    mock_response.ok = True
    mock_response.json.return_value = [
        {"pdbUrl": f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.pdb"}
    ]
    return mock_response


def mock_alphafold_files_api_responses(url):
    """
    Mock the response to calls made to the AlphaFold API to download PDB files.
    """

    # Parse the accession from the url. Isoform ids look like AF-Q9Y6V0-3-F1-...
    result = re.findall(r"AF-([A-Z0-9]+)(?:-\d+)?-F1-model_v[46]\.pdb", url)
    if not result:
        raise ValueError(f"Unexpected url: {url}")
    accession = result[0]

    mock_response = mock.Mock(spec=requests.Response)
    mock_response.status_code = 200
    mock_response.ok = True

    # The AlphaFold URL requests v6 files, but the test artifacts contain v4 files.
    # Because we can safely use the "old" v4 files for testing purposes,
    # we point to the v4 files here.
    artifact_filepath = (
        API_RESPONSE_ARTIFACTS_DIRPATH / f"alphafold.ebi.ac.uk_files_AF-{accession}-F1-model_v4.pdb"
    )
    if not artifact_filepath.exists():
        raise ValueError(
            f"No artifact found for the AlphaFold PDB file for {accession} at {artifact_filepath})"
        )

    with open(artifact_filepath, encoding="utf-8") as file:
        mock_response.text = file.read()

    return mock_response
