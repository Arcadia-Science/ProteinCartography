import pathlib
import shutil

import pandas as pd
import pytest
import snakemake
import yaml


@pytest.fixture
def config_filepath(tmp_path):
    """
    Generate a config file for testing the pipeline in "cluster" mode
    """
    config = {
        "input_dir": str(tmp_path / "input"),
        "output_dir": str(tmp_path / "output"),
        "analysis_name": "test",
        "plotting_modes": ["pca_umap"],
        "taxon_focus": "euk",
    }

    filepath = tmp_path / "config.yaml"
    with open(filepath, "w") as file:
        yaml.dump(config, file)

    return str(filepath)


@pytest.fixture
def stage_inputs(integration_test_artifacts_dirpath, config_filepath):
    """
    Create the input directory and the PDB files for the pipeline in "cluster" mode
    """
    with open(config_filepath) as file:
        config = yaml.safe_load(file)

    # for now, hard-code the dataset name
    dataset_name = "actin"
    shutil.copytree(
        integration_test_artifacts_dirpath / "cluster-mode" / dataset_name / "input",
        config["input_dir"],
    )


@pytest.fixture
def snakefile_filepath(repo_dirpath):
    """
    The path to the cluster-mode Snakefile
    """
    return repo_dirpath / "Snakefile_ff"


@pytest.mark.usefixtures("stage_inputs")
def test_pipeline_in_cluster_mode(snakefile_filepath, config_filepath):
    """
    Run the pipeline in "cluster" mode with the test config file
    """

    snakemake.snakemake(
        snakefile=snakefile_filepath,
        configfiles=[config_filepath],
        use_conda=True,
        cores=8,
        verbose=True,
    )

    with open(config_filepath) as file:
        config = yaml.safe_load(file)

    input_dirpath = pathlib.Path(config["input_dir"])
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
        assert filepath.exists()

    # check that the shape of the all-by-all similarity matrix is correct
    num_structures = len(list(input_dirpath.glob("*.pdb")))
    similarity_matrix_filepath = (
        output_dirpath / "foldseek_clustering_results" / "all_by_all_tmscore_pivoted.tsv"
    )
    similarity_matrix = pd.read_csv(similarity_matrix_filepath, sep="\t")

    # matrix should have one row and one column per structure (plus one column for the index)
    assert similarity_matrix.shape == (num_structures, num_structures + 1)
