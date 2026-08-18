"""Search-correctness: domain-path hit union is per-query-domain, not full-length search."""

from __future__ import annotations
from pathlib import Path

import aggregate_domain_hits
import domain_utils as du


def test_union_of_per_domain_hits_includes_both_neighborhoods(tmp_path: Path):
    d01 = tmp_path / "P99999__d01.foldseek_hits.txt"
    d02 = tmp_path / "P99999__d02.foldseek_hits.txt"
    d01.write_text("A1\n")
    d02.write_text("B1\n")
    out = tmp_path / "aggregated_hits.txt"
    found = tmp_path / "found_by.tsv"
    accessions = aggregate_domain_hits.aggregate_domain_hits(
        [str(d01), str(d02)], str(out), str(found)
    )
    assert accessions == ["A1", "B1"]
    text = found.read_text()
    assert "A1\tP99999__d01" in text
    assert "B1\tP99999__d02" in text


def test_protein_path_hits_are_a_separate_list(tmp_path: Path):
    """Full-length search of Q returning only A1 must not be used as the domain union."""
    protein_hits = {"A1"}
    d01 = tmp_path / "P99999__d01.blast_hits.uniprot.txt"
    d02 = tmp_path / "P99999__d02.blast_hits.uniprot.txt"
    d01.write_text("A1\n")
    d02.write_text("B1\n")
    domain_hits = set(
        aggregate_domain_hits.aggregate_domain_hits(
            [str(d01), str(d02)],
            str(tmp_path / "agg.txt"),
            str(tmp_path / "found.tsv"),
        )
    )
    assert protein_hits == {"A1"}
    assert domain_hits == {"A1", "B1"}
    assert "B1" not in protein_hits


def test_domain_search_inputs_are_cropped_domain_files_not_parent(tmp_path: Path):
    """Per-domain search must use Q__dNN files; searching Q.pdb then splitting would fail this."""
    parent_fasta = tmp_path / "Q.fasta"
    parent_fasta.write_text(">Q\n" + "A" * 80 + "C" * 80 + "\n")
    d01 = tmp_path / "Q__d01.fasta"
    d02 = tmp_path / "Q__d02.fasta"
    du.crop_fasta_file(parent_fasta, "1-80", d01, "Q__d01")
    du.crop_fasta_file(parent_fasta, "81-160", d02, "Q__d02")
    _, seq01 = du.parse_fasta(d01.read_text())
    _, seq02 = du.parse_fasta(d02.read_text())
    _, parent_seq = du.parse_fasta(parent_fasta.read_text())
    assert seq01 != parent_seq
    assert seq02 != parent_seq
    assert seq01 != seq02
    assert d01.name == "Q__d01.fasta"
    assert d02.name == "Q__d02.fasta"


def test_domain_id_from_hit_filename():
    assert (
        aggregate_domain_hits._domain_id_from_hit_file("out/P99999__d01.blast_hits.uniprot.txt")
        == "P99999__d01"
    )
    parent, index = du.parse_domain_id("P99999__d02")
    assert parent == "P99999" and index == 2


def test_maybe_write_per_domain_hits_partitions_neighborhoods(tmp_path: Path):
    from tests import mocks

    d01 = tmp_path / "P99999__d01.foldseek_hits.txt"
    d02 = tmp_path / "P99999__d02.blast_hits.uniprot.txt"
    parent = tmp_path / "P99999.foldseek_hits.txt"
    assert mocks.maybe_write_per_domain_hits(str(d01))
    assert d01.read_text() == "A0A286Q506\n"
    assert mocks.maybe_write_per_domain_hits(str(d02))
    assert d02.read_text() == "Q6QAQ1\n"
    assert not mocks.maybe_write_per_domain_hits(str(parent))
    assert not parent.exists()
