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


# ------------------------------------------------------------------------------------------------
# Detecting a `blastp -remote` call that failed but exited zero.
# ------------------------------------------------------------------------------------------------


class _Completed:
    def __init__(self, returncode=0, stderr=b""):
        self.returncode = returncode
        self.stderr = stderr
        self.stdout = b""


def test_a_successful_call_did_not_fail():
    assert not blast_utils.blast_call_failed(_Completed())


def test_a_nonzero_exit_code_is_a_failure():
    assert blast_utils.blast_call_failed(_Completed(returncode=1))


def test_an_error_on_stderr_is_a_failure_despite_a_zero_exit_code():
    """
    The case this exists for: when the remote server refuses to queue the request, `blastp
    -remote` reports the error, writes an empty results file, and still exits zero. Branching on
    the exit code alone treats that as a success, so the word-size backoff never runs and the
    failure surfaces later as a misleading "no hits were returned".
    """
    stderr = (
        b"Error: [blastp] bad_request: Could not queue request: DB operation failed."
        b"(ERR:-99  )((Severe Error) DB Put Request error: sp_NewRequestEx failed)\n"
    )
    assert blast_utils.blast_call_failed(_Completed(returncode=0, stderr=stderr))


def test_a_warning_on_stderr_is_not_a_failure():
    stderr = b"Warning: [blastp] Examining 5 or more matches is recommended\n"
    assert not blast_utils.blast_call_failed(_Completed(returncode=0, stderr=stderr))


def test_stderr_may_be_bytes_a_string_or_none():
    assert blast_utils.blast_call_failed(_Completed(stderr="Error: something went wrong"))
    assert not blast_utils.blast_call_failed(_Completed(stderr="all good"))
    assert not blast_utils.blast_call_failed(_Completed(stderr=None))
