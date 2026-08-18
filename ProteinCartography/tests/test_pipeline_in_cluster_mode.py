import os
import pathlib
import shutil

import pandas as pd
import pytest
import snakemake
import yaml


@pytest.fixture
def config_filepath(tmp_path):
    """
    Generate a config file for testing the pipeline in "cluster" mode.
    """
    config = {
        "mode": "cluster",
        "analysis_name": "test",
        "input_dir": str(tmp_path / "input"),
        "output_dir": str(tmp_path / "output"),
        "plotting_modes": ["pca_umap"],
        "features_file": "uniprot_features.tsv",
        "key_protids": ["P60709"],
    }

    filepath = tmp_path / "config.yaml"
    with open(filepath, "w") as file:
        yaml.dump(config, file)

    return str(filepath)


@pytest.fixture
def stage_inputs(integration_test_artifacts_dirpath, config_filepath):
    """
    Create the input directory and the PDB files for the pipeline in "cluster" mode.
    """
    with open(config_filepath) as file:
        config = yaml.safe_load(file)

    # For now, hard-code the dataset name.
    dataset_name = "actin"
    shutil.copytree(
        integration_test_artifacts_dirpath / "cluster-mode" / dataset_name / "input",
        config["input_dir"],
    )


@pytest.fixture
def set_env_variables(pytestconfig):
    """Mock TED (and other HTTP) so domain_map auto does not call live APIs."""
    should_use_mocks = "PROTEINCARTOGRAPHY_SHOULD_USE_MOCKS"
    if not pytestconfig.getoption("no_mocks"):
        os.environ[should_use_mocks] = "true"
    yield
    os.environ.pop(should_use_mocks, None)


@pytest.mark.usefixtures("stage_inputs")
@pytest.mark.usefixtures("set_env_variables")
def test_pipeline_in_cluster_mode(repo_dirpath, config_filepath):
    """
    Run the pipeline in "cluster" mode with the test config file.
    """

    snakemake.snakemake(
        snakefile=(repo_dirpath / "Snakefile"),
        configfiles=[config_filepath],
        use_conda=True,
        cores=8,
        verbose=True,
    )

    with open(config_filepath) as file:
        config = yaml.safe_load(file)

    input_dirpath = pathlib.Path(config["input_dir"])
    output_dirpath = pathlib.Path(config["output_dir"])

    # Check (some of) the expected output files.
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
        assert filepath.exists()

    # Check that the shape of the all-by-all similarity matrix is correct.
    num_structures = len(list(input_dirpath.glob("*.pdb")))
    similarity_matrix_filepath = (
        output_dirpath / "foldseek_clustering_results" / "all_by_all_tmscore_pivoted.tsv"
    )
    similarity_matrix = pd.read_csv(similarity_matrix_filepath, sep="\t")

    # The matrix should have one row and one column per structure (plus one column for the index).
    assert similarity_matrix.shape == (num_structures, num_structures + 1)

    # Actin cluster inputs are single-domain under the TED mock; the domain DAG is a no-op.
    domain_html_name = f"{config['analysis_name']}_leiden_similarity_domain.html"
    domain_html = output_dirpath / "final_results" / domain_html_name
    assert not domain_html.exists()
