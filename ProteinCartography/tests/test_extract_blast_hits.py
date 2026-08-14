"""Unit tests for BLAST hit extraction."""

from __future__ import annotations
import pathlib

import constants
import extract_blast_hits
import pytest


def test_extract_blast_hits_empty_hard_fails_by_default(tmp_path: pathlib.Path):
    infile = tmp_path / "blast.tsv"
    outfile = tmp_path / "hits.txt"
    infile.write_text("")
    with pytest.raises(Exception, match="No hits were returned"):
        extract_blast_hits.extract_blast_hits(
            str(infile),
            str(outfile),
            column_names=["qseqid", "sacc"],
        )
    assert outfile.read_text() == ""


def test_extract_blast_hits_empty_soft_fail(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setenv("PC_BLAST_SOFT_FAIL", "1")
    infile = tmp_path / "blast.tsv"
    outfile = tmp_path / "hits.txt"
    infile.write_text("")
    extract_blast_hits.extract_blast_hits(
        str(infile),
        str(outfile),
        column_names=["qseqid", "sacc"],
    )
    assert outfile.read_text() == ""


def test_extract_blast_hits_reads_sacc(tmp_path: pathlib.Path):
    infile = tmp_path / "blast.tsv"
    outfile = tmp_path / "hits.txt"
    infile.write_text("q1\tNP_1\nq1\tNP_2\nq1\tNP_1\n")
    extract_blast_hits.extract_blast_hits(
        str(infile),
        str(outfile),
        column_names=["qseqid", "sacc"],
    )
    hits = outfile.read_text().splitlines()
    assert hits == ["NP_1", "NP_2"]


def test_soft_fail_flag_is_opt_in_only(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("PC_BLAST_SOFT_FAIL", raising=False)
    assert constants.blast_soft_fail_enabled() is False
    monkeypatch.setenv("PC_BLAST_SOFT_FAIL", "1")
    assert constants.blast_soft_fail_enabled() is True
    monkeypatch.setenv("PC_BLAST_SOFT_FAIL", "maybe")
    assert constants.blast_soft_fail_enabled() is False
