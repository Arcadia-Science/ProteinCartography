"""Unit tests for domain identity, chopping, FASTA/PDB crop, and the multi-domain gate."""

from __future__ import annotations
import pathlib

import config_utils
import domain_utils as du
import pytest

P00698_CHOPPING = "22-141"


def _ca_atom(serial: int, resseq: int, x: float = 0.0) -> str:
    return (
        f"ATOM  {serial:5d}  CA  ALA A{resseq:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 90.00           C  "
    )


def _pdb_for_residues(resseqs: list[int]) -> str:
    lines = [_ca_atom(i + 1, resseq, float(i)) for i, resseq in enumerate(resseqs)]
    lines.append("END")
    return "\n".join(lines) + "\n"


def test_parse_chopping_continuous():
    assert du.parse_chopping("22-141") == [(22, 141)]


def test_parse_chopping_discontinuous():
    assert du.parse_chopping("22-50_80-120") == [(22, 50), (80, 120)]


def test_parse_chopping_rejects_garbage():
    with pytest.raises(du.DomainChoppingError):
        du.parse_chopping("")
    with pytest.raises(du.DomainChoppingError):
        du.parse_chopping("abc")
    with pytest.raises(du.DomainChoppingError):
        du.parse_chopping("10-5")


def test_domain_id_round_trip():
    assert du.make_domain_id("P12345", 1) == "P12345__d01"
    assert du.parse_domain_id("P12345__d01") == ("P12345", 1)
    assert du.make_domain_id("P12345", 10) == "P12345__d10"
    assert du.parse_domain_id("P12345__d10") == ("P12345", 10)


def test_slice_fasta_one_based_inclusive():
    sequence = "A" * 141
    sliced = du.slice_fasta_sequence(sequence, P00698_CHOPPING)
    assert sliced == sequence[21:141]
    assert len(sliced) == 120


def test_write_fasta_uses_domain_id_header(tmp_path: pathlib.Path):
    out = tmp_path / "P00698__d01.fasta"
    du.write_fasta(out, "P00698__d01", "ACDE")
    header, seq = du.parse_fasta(out.read_text())
    assert header == "P00698__d01"
    assert seq == "ACDE"


def test_crop_pdb_by_resseq_not_atom_serial():
    pdb = _pdb_for_residues(list(range(1, 201)))
    cropped, n_ca = du.crop_pdb_text(pdb, P00698_CHOPPING)
    resseqs = []
    for line in cropped.splitlines():
        if line.startswith("ATOM"):
            resseqs.append(int(line[22:26]))
    assert min(resseqs) == 22
    assert max(resseqs) == 141
    assert n_ca == 120


def test_crop_pdb_discontinuous_omits_gap():
    pdb = _pdb_for_residues(list(range(1, 121)))
    cropped, n_ca = du.crop_pdb_text(pdb, "22-50_80-120")
    resseqs = {int(line[22:26]) for line in cropped.splitlines() if line.startswith("ATOM")}
    assert 51 not in resseqs
    assert 79 not in resseqs
    assert 22 in resseqs and 50 in resseqs
    assert 80 in resseqs and 120 in resseqs
    assert n_ca == (50 - 22 + 1) + (120 - 80 + 1)


def test_crop_pdb_uses_resseq_when_file_starts_at_20():
    pdb = _pdb_for_residues(list(range(20, 50)))
    cropped, n_ca = du.crop_pdb_text(pdb, "22-30")
    resseqs = [int(line[22:26]) for line in cropped.splitlines() if line.startswith("ATOM")]
    assert resseqs == list(range(22, 31))
    assert n_ca == 9


def test_min_domain_length_drops_short_span():
    row = du.domain_row("P1", 1, "1-10", "ted", nres_domain=10)
    kept = du.filter_domain_rows([row], min_domain_length=30)
    assert kept == []


def test_is_multidomain_gate():
    two = [
        du.domain_row("P1", 1, "1-80", "ted"),
        du.domain_row("P1", 2, "81-160", "ted"),
    ]
    assert du.is_multidomain(len(du.filter_domain_rows(two, 30)))
    one = [du.domain_row("P1", 1, "22-141", "ted")]
    assert not du.is_multidomain(len(du.filter_domain_rows(one, 30)))
    mixed = [
        du.domain_row("P1", 1, "1-80", "ted"),
        du.domain_row("P1", 2, "81-90", "ted", nres_domain=10),
    ]
    kept = du.filter_domain_rows(mixed, 30)
    assert not du.is_multidomain(len(kept))
    assert kept[0]["protid"] == "P1__d01"


def test_looks_like_uniprot():
    assert du.looks_like_uniprot_accession("P60709")
    assert du.looks_like_uniprot_accession("P99999")
    assert du.looks_like_uniprot_accession("A0A2Y9FRR4")
    assert not du.looks_like_uniprot_accession("QTEST")


def test_discontinuous_chopping_is_one_domain_instance():
    row = du.domain_row("P12345", 1, "22-50_80-120", "ted")
    assert row["protid"] == "P12345__d01"
    assert row["chopping"] == "22-50_80-120"
    assert row["nres_domain"] == (50 - 22 + 1) + (120 - 80 + 1)


def test_domain_map_config_values():
    assert config_utils._get_domain_map({}) == "auto"
    assert config_utils._get_domain_map({"domain_map": "OFF"}) == "off"
    with pytest.raises(config_utils.ProteinCartographyInputError):
        config_utils._get_domain_map({"domain_map": "domain-only"})


def test_missing_user_domains_file_raises(tmp_path: pathlib.Path):
    with pytest.raises(config_utils.ProteinCartographyInputError, match="user_domains_file"):
        config_utils._get_user_domains_file(
            {"input_dir": str(tmp_path), "user_domains_file": "missing.tsv"}
        )
