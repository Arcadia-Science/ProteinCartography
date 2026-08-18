"""Mocked Snakemake run of the parallel domain path for a two-domain query."""

from __future__ import annotations
import os
import pathlib

import pandas as pd
import pytest
import snakemake
import yaml


def _ca_pdb(nres: int) -> str:
    lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
        f"{float(i):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 90.00           C  "
        for i in range(1, nres + 1)
    ]
    return "\n".join(lines) + "\nEND\n"


def _load_config(filepath):
    with open(filepath) as file:
        return yaml.safe_load(file)


@pytest.fixture
def config_filepath(tmp_path):
    config = {
        "mode": "search",
        "analysis_name": "testdomain",
        "input_dir": str(tmp_path / "input"),
        "output_dir": str(tmp_path / "output"),
        "plotting_modes": ["pca_umap"],
        "max_blast_hits": 10,
        "max_foldseek_hits": 10,
        "max_structures": 10,
        "domain_map": "auto",
    }
    filepath = tmp_path / "config.yaml"
    with open(filepath, "w") as file:
        yaml.dump(config, file)
    return str(filepath)


@pytest.fixture
def stage_inputs(config_filepath):
    config = _load_config(config_filepath)
    input_dir = pathlib.Path(config["input_dir"])
    input_dir.mkdir(parents=True)
    (input_dir / "P99999.fasta").write_text(">P99999\n" + "A" * 160 + "\n")
    (input_dir / "P99999.pdb").write_text(_ca_pdb(160))


@pytest.fixture
def set_env_variables(pytestconfig):
    should_use_mocks = "PROTEINCARTOGRAPHY_SHOULD_USE_MOCKS"
    should_log_api_requests = "PROTEINCARTOGRAPHY_SHOULD_LOG_API_REQUESTS"
    if not pytestconfig.getoption("no_mocks"):
        os.environ[should_use_mocks] = "true"
    should_log_api_requests_value = os.environ.pop(should_log_api_requests, None)
    yield
    os.environ.pop(should_use_mocks, None)
    if should_log_api_requests_value is not None:
        os.environ[should_log_api_requests] = should_log_api_requests_value


@pytest.mark.usefixtures("stage_inputs")
@pytest.mark.usefixtures("set_env_variables")
def test_pipeline_in_search_mode_domain_map(repo_dirpath, config_filepath):
    snakemake.snakemake(
        snakefile=(repo_dirpath / "Snakefile"),
        configfiles=[config_filepath],
        use_conda=True,
        cores=8,
        verbose=True,
    )

    config = _load_config(config_filepath)
    output_dirpath = pathlib.Path(config["output_dir"])
    name = config["analysis_name"]

    protein_html = output_dirpath / "final_results" / f"{name}_leiden_similarity.html"
    domain_html = output_dirpath / "final_results" / f"{name}_leiden_similarity_domain.html"
    assert protein_html.exists()
    assert domain_html.exists()

    domain_matrix = (
        output_dirpath / "foldseek_clustering_results_domain" / "all_by_all_tmscore_pivoted.tsv"
    )
    assert domain_matrix.exists()
    df = pd.read_csv(domain_matrix, sep="\t")
    assert "protid" in df.columns
    assert df["protid"].astype(str).str.contains("__d").all()
    assert "P99999__d01" in set(df["protid"].astype(str))
    assert "P99999__d02" in set(df["protid"].astype(str))

    features = pd.read_csv(
        output_dirpath / "final_results" / f"{name}_aggregated_features_domain.tsv", sep="\t"
    )
    tm_cols = [c for c in features.columns if c.startswith("TMscore_v_")]
    assert tm_cols
    assert all(c.replace("TMscore_v_", "").startswith("P99999__d") for c in tm_cols)
    assert "TMscore_v_P99999" not in tm_cols
    assert not any(c.startswith("fident_v_") for c in features.columns)

    query_fastas = list(
        (output_dirpath / "domain_path" / "query_structures").glob("P99999__d*.fasta")
    )
    assert {p.name for p in query_fastas} == {"P99999__d01.fasta", "P99999__d02.fasta"}

    domain_hits = {
        line
        for line in (output_dirpath / "domain_path" / "features" / "aggregated_hits.txt")
        .read_text()
        .splitlines()
        if line
    }
    assert domain_hits == {"A0A286Q506", "Q6QAQ1"}
    found = (output_dirpath / "domain_path" / "features" / "found_by.tsv").read_text()
    assert "A0A286Q506\tP99999__d01" in found
    assert "Q6QAQ1\tP99999__d02" in found
    protein_hits_file = output_dirpath / "protein_features" / "aggregated_hits.txt"
    if protein_hits_file.is_file():
        protein_hits = {line for line in protein_hits_file.read_text().splitlines() if line}
        assert domain_hits != protein_hits
