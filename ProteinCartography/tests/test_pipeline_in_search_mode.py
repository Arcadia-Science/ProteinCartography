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
        "mode": "search",
        "analysis_name": "test",
        "input_dir": str(tmp_path / "input"),
        "output_dir": str(tmp_path / "output"),
        "plotting_modes": ["pca_umap"],
        "max_blast_hits": 10,
        "max_foldseek_hits": 10,
        "max_structures": 10,
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
        integration_test_artifacts_dirpath / "search-mode" / dataset_name / "input",
        config["input_dir"],
    )


@pytest.fixture
def set_env_variables(pytestconfig):
    """
    Set the env variable used to mock the API responses
    made by the python scripts that are called by the snakemake rules

    Note: this works because the rule environments inherit their env variables from the environment
    in which `snakemake` was called (which is this pytest python process)
    """

    # the names of the proteincartography-specific env variables
    should_use_mocks = "PROTEINCARTOGRAPHY_SHOULD_USE_MOCKS"
    should_log_api_requests = "PROTEINCARTOGRAPHY_SHOULD_LOG_API_REQUESTS"

    if not pytestconfig.getoption("no_mocks"):
        os.environ[should_use_mocks] = "true"

    # don't log API requests during the tests
    should_log_api_requests_value = os.environ.pop(should_log_api_requests, None)

    yield

    os.environ.pop(should_use_mocks, None)

    # as a convenience, restore the logging env variable to its original value
    if should_log_api_requests_value is not None:
        os.environ[should_log_api_requests] = should_log_api_requests_value


@pytest.mark.usefixtures("stage_inputs")
@pytest.mark.usefixtures("set_env_variables")
def test_pipeline_in_search_mode_with_mocked_api_calls(repo_dirpath, config_filepath):
    """
    Run the pipeline in "search" mode with the test config file, the temporary snakefile,
    and mocked API calls
    """

    snakemake.snakemake(
        snakefile=(repo_dirpath / "Snakefile"),
        configfiles=[config_filepath],
        use_conda=True,
        cores=8,
        verbose=True,
    )

    config = _load_config(config_filepath)
    output_dirpath = pathlib.Path(config["output_dir"])

    # check (some of) the expected output files
    expected_output_filepaths = [
        output_dirpath / "final_results" / f"{config['analysis_name']}_{appendix}"
        for appendix in [
            "leiden_similarity.html",
            "strucluster_similarity.html",
            "semantic_analysis.pdf",
            "semantic_analysis.html",
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
        output_dirpath / "foldseek_clustering_results" / "all_by_all_tmscore_pivoted.tsv"
    )
    similarity_matrix = pd.read_csv(similarity_matrix_filepath, sep="\t")
    assert similarity_matrix.shape == (11, 12)
