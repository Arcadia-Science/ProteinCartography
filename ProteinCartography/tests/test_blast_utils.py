"""Unit tests for blast_utils WWW XML → TSV conversion."""

from __future__ import annotations
import pathlib

import blast_utils


def test_xml_to_blast_tsv_writes_sacc(tmp_path: pathlib.Path):
    xml = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_query-def>P07338</BlastOutput_query-def>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_id>ref|NP_001234.1|</Hit_id>
          <Hit_accession>NP_001234</Hit_accession>
          <Hit_def>example protein</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>100</Hsp_bit-score>
              <Hsp_evalue>1e-20</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>50</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>50</Hsp_hit-to>
              <Hsp_identity>45</Hsp_identity>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>50</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""
    out = tmp_path / "blast.tsv"
    n = blast_utils._xml_to_blast_tsv(xml, str(out))
    assert n == 1
    line = out.read_text().strip().split("\t")
    assert line[0] == "P07338"
    assert "NP_001234" in line  # sacc column


def test_xml_to_blast_tsv_empty_without_blast_output(tmp_path: pathlib.Path):
    out = tmp_path / "empty.tsv"
    n = blast_utils._xml_to_blast_tsv("<html>waiting</html>", str(out))
    assert n == 0
    assert out.read_text() == ""
