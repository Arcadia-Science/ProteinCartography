"""Unit tests for empty BLAST hit extraction."""

from __future__ import annotations
import pathlib

import extract_blast_hits


def test_extract_blast_hits_empty_input_writes_empty_file(tmp_path: pathlib.Path):
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
