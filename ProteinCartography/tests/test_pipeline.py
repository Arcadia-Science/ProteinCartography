import os
import pathlib
import shutil

import pandas as pd
import pytest
import snakemake
import yaml


def _load_config(filepath):
    """
    Convenience function to load a yaml config file
    """
    with open(filepath) as file:
        config = yaml.safe_load(file)
    return config


@pytest.fixture
def config_filepath(tmp_path):
    """
    Generate a config file for testing the pipeline
    """
    config = {
        "input_dir": str(tmp_path / "input"),
        "output_dir": str(tmp_path / "output"),
        "analysis_name": "test",
        "foldseek_databases": ["afdb50", "afdb-swissprot", "afdb-proteome"],
        "plotting_modes": ["pca_umap"],
        "max_blasthits": 10,
        "max_foldseekhits": 10,
        "max_structures": 10,
        "taxon_focus": "euk",
        "uniprot_additional_fields": [],
        "min_length": 0,
        "max_length": 0,
        "resources": {"mem_mb": 16 * 1000},
        "cores": 8,
        "max-jobs-per-second": 1,
        "max-status-checks-per-second": 10,
        "local-cores": 10,
    }

    filepath = tmp_path / "config.yaml"
    with open(filepath, "w") as file:
        yaml.dump(config, file)

    return str(filepath)


@pytest.fixture
def stage_inputs(integration_test_artifacts_dirpath, config_filepath):
    """
    Create the input directory and input files for the pipeline
    """

    config = _load_config(config_filepath)

    # for now, hard-code the dataset name
    dataset_name = "actin"
    shutil.copytree(
        integration_test_artifacts_dirpath / dataset_name / "input", config["input_dir"]
    )

    # copy the blastresults.csv file to the input directory
    # so it can be used by the `mock_run_blast` rule to mock the output of the `run_blast` rule
    shutil.copy(
        integration_test_artifacts_dirpath / dataset_name / "output" / "P60709.blastresults.tsv",
        config["input_dir"],
    )


@pytest.fixture
def snakefile_filepath(repo_dirpath):
    """
    Generate a snakefile with the `run_blast` rule overridden by a `mock_run_blast` rule
    to eliminate the call to `blastp`

    TODO (KC): to avoid this super ugly hack, wrap the call to blastp in a python method
    that can be patched
    """

    # copy the snakefile to a temporary file
    # note: snakemake requires that the `envs/` dir be relative to the dir containing the snakefile,
    # so we cannot write the modified snakefile to the tmp_path directory but instead must write it
    # to the repo directory and manually remove it after the test has run
    filepath = repo_dirpath / "Snakefile_tmp"
    shutil.copy(repo_dirpath / "Snakefile", filepath)

    # define a rule to copy the blastresults.tsv file from the input to the output directory
    # and prioritize it over the `run_blast` rule to prevent the latter from being called
    mock_run_blast_rule = """
rule mock_run_blast:
    input:
        cached_blast_results=input_dir / "{protid}.blastresults.tsv"
    output:
        blastresults=output_dir / blastresults_dir / "{protid}.blastresults.tsv"
    shell:
        "cp {input.cached_blast_results} {output.blastresults}"
ruleorder: mock_run_blast > run_blast
    """

    with open(repo_dirpath / "Snakefile_tmp", "a") as file:
        file.write(mock_run_blast_rule)

    yield filepath

    os.remove(filepath)


@pytest.fixture
def set_env_variables():
    """
    Set the env variable used to mock the API responses
    made by the python scripts that are called by the snakemake rules

    Note: this works because the rule environments inherit their env variables from the environment
    in which `snakemake` was called (which is this pytest python process)
    """

    # the names of the proteincartography-specific env variables
    was_called_by_pytest = "PROTEINCARTOGRAPHY_WAS_CALLED_BY_PYTEST"
    should_log_api_requests = "PROTEINCARTOGRAPHY_SHOULD_LOG_API_REQUESTS"

    os.environ[was_called_by_pytest] = "true"

    # don't log API requests during the tests
    should_log_api_requests_value = os.environ.pop(should_log_api_requests, None)

    yield

    os.environ.pop(was_called_by_pytest)

    # as a convenience, restore the logging env variable to its original value
    if should_log_api_requests_value is not None:
        os.environ[should_log_api_requests] = should_log_api_requests_value


@pytest.mark.usefixtures("stage_inputs")
@pytest.mark.usefixtures("set_env_variables")
def test_pipeline_with_mocked_api_calls(snakefile_filepath, config_filepath):
    """
    Run the pipeline with the test config file, the temporary snakefile,
    and mocked API calls
    """

    snakemake.snakemake(
        snakefile=snakefile_filepath,
        configfiles=[config_filepath],
        use_conda=True,
        cores=8,
        verbose=True,
    )

    config = _load_config(config_filepath)
    output_dirpath = pathlib.Path(config["output_dir"])

    # check (some of) the expected output files
    expected_output_filepaths = [
        output_dirpath / "clusteringresults" / f"{config['analysis_name']}{appendix}"
        for appendix in [
            "_leiden_similarity.html",
            "_strucluster_similarity.html",
            "_semantic_analysis.pdf",
            "_semantic_analysis.html",
        ]
    ]
    for filepath in expected_output_filepaths:
        # TODO (KC): check that the content of the files looks correct
        # (not sure we can do a literal comparison because of timestamps, umap stochasticity, etc.)
        assert filepath.exists()

    # check that the shape of the all-by-all similarity matrix is correct:
    # there should be 11 structures clustered by foldseek
    # (the 10 determined by the `max_structures` config param, plus the input structure),
    # so the dataframe should have 11 rows and 12 columns (since the first column is the index)
    similarity_matrix_filepath = (
        output_dirpath / "clusteringresults" / "all_by_all_tmscore_pivoted.tsv"
    )
    similarity_matrix = pd.read_csv(similarity_matrix_filepath, sep="\t")
    assert similarity_matrix.shape == (11, 12)
