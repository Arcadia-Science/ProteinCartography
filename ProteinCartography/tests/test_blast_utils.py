"""Unit tests for blast_utils WWW XML → TSV conversion."""

from __future__ import annotations
import pathlib

import blast_utils
import pytest


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


# ------------------------------------------------------------------------------------------------
# Waiting on a queued search.
#
# The clock and NCBI are both faked, so these exercise multi-hour timeouts instantly and without
# network access. `_run_blast_www` measures elapsed time with `time.time`, so that is what the
# fake clock replaces.
# ------------------------------------------------------------------------------------------------


PUT_RESPONSE = """<html><body>
<!--QBlastInfoBegin
    RID = M8K4ARMH013
    RTOE = {rtoe}
QBlastInfoEnd
-->
</body></html>
"""

READY_XML = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_query-def>P60709</BlastOutput_query-def>
  <BlastOutput_iterations><Iteration><Iteration_hits></Iteration_hits></Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


def waiting_response():
    return "<html><body><pre>QBlastInfoBegin\n\tStatus=WAITING\nQBlastInfoEnd</pre></body></html>"


class FakeClock:
    """
    A clock that only advances when something sleeps, so a 2-hour wait costs no real time.
    """

    def __init__(self):
        self.now = 1000.0
        self.sleeps = []

    def time(self):
        return self.now

    def monotonic(self):
        return self.now

    def sleep(self, seconds):
        self.sleeps.append(seconds)
        self.now += seconds


class FakeNcbi:
    """
    Stand in for NCBI's Blast.cgi, dispatching on the CMD in the request like NCBI does.

    `statuses` are returned to successive polls; the last one repeats forever so that a test can
    describe a search that never finishes without listing a status per poll. An `Exception` in the
    list is raised instead, to describe a poll that does not reach NCBI at all.
    """

    def __init__(self, clock, rtoe=27, statuses=("READY",)):
        self.clock = clock
        self.rtoe = rtoe
        self.statuses = list(statuses)
        self.poll_times = []
        self.submitted_at = None

    def urlopen(self, req, timeout=None):
        url = req.full_url if hasattr(req, "full_url") else str(req)
        if "CMD=Get" in url:
            self.poll_times.append(self.clock.time())
            status = self.statuses.pop(0) if len(self.statuses) > 1 else self.statuses[0]
            if isinstance(status, Exception):
                raise status
            body = READY_XML if status == "READY" else waiting_response()
            if status not in {"READY", "WAITING"}:
                body = f"<html><pre>Status={status}</pre></html>"
            return _FakeResponse(body)
        self.submitted_at = self.clock.time()
        return _FakeResponse(PUT_RESPONSE.format(rtoe=self.rtoe))


class _FakeResponse:
    def __init__(self, text):
        self._text = text

    def read(self):
        return self._text.encode()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


@pytest.fixture
def clock(monkeypatch):
    fake = FakeClock()
    monkeypatch.setattr(blast_utils.time, "time", fake.time)
    monkeypatch.setattr(blast_utils.time, "monotonic", fake.monotonic)
    monkeypatch.setattr(blast_utils.time, "sleep", fake.sleep)
    return fake


@pytest.fixture
def ncbi(monkeypatch, clock):
    def install(**kwargs):
        fake = FakeNcbi(clock, **kwargs)
        monkeypatch.setattr(blast_utils.urllib.request, "urlopen", fake.urlopen)
        return fake

    return install


@pytest.fixture
def query_filepath(tmp_path):
    path = tmp_path / "P60709.fasta"
    path.write_text(">sp|P60709|ACTB_HUMAN\nMDDDIAALVVDNGSGMCKAGF\n")
    return path


def run_www(query_filepath, out, **kwargs):
    kwargs.setdefault("timeout_seconds", 600)
    return blast_utils._run_blast_www(
        str(query_filepath), str(out), max_target_seqs=3000, word_size=5, evalue=1.0, **kwargs
    )


def test_an_estimate_longer_than_the_timeout_does_not_extend_it(
    clock, ncbi, query_filepath, tmp_path
):
    """
    The configured timeout is a ceiling, not a floor.

    NCBI's estimate is routinely hours and tracks how busy its shared queue is rather than this
    search. Treating it as a reason to wait longer than configured -- which
    `max(timeout, min(rtoe + 300, 7200))` did -- means the setting cannot bound anything, which
    is the failure mode `blastp -remote` had.
    """
    ncbi(rtoe=9168, statuses=["WAITING"])

    result = run_www(query_filepath, tmp_path / "out.tsv", timeout_seconds=600)

    assert result.returncode == 124
    assert sum(clock.sleeps) <= 600
    assert "TIMEOUT after 600s" in result.stderr.decode()


def test_a_search_that_is_ready_despite_a_long_estimate_succeeds(
    clock, ncbi, query_filepath, tmp_path
):
    """
    A search is polled at least once even when the estimate exceeds the whole budget, so a search
    that turned out to be ready is not reported as a timeout.
    """
    ncbi(rtoe=9168, statuses=["READY"])

    result = run_www(query_filepath, tmp_path / "out.tsv", timeout_seconds=600)

    assert result.returncode == 0


def test_polls_that_never_reach_ncbi_end_the_search(clock, ncbi, query_filepath, tmp_path):
    """
    A poll that fails is retried, but not forever: an unreachable NCBI is bounded by
    MAX_CONSECUTIVE_POLL_FAILURES rather than only by the timeout.
    """
    fake = ncbi(rtoe=1, statuses=[OSError("connection refused")])

    result = run_www(query_filepath, tmp_path / "out.tsv", timeout_seconds=100_000)

    assert result.returncode == 1
    assert len(fake.poll_times) == blast_utils.MAX_CONSECUTIVE_POLL_FAILURES
    assert "consecutive polls" in result.stderr.decode()


def test_a_transient_poll_failure_is_retried(clock, ncbi, query_filepath, tmp_path):
    """
    A single refused connection does not discard the wait already spent.
    """
    ncbi(rtoe=1, statuses=[OSError("connection refused"), "READY"])

    result = run_www(query_filepath, tmp_path / "out.tsv", timeout_seconds=100_000)

    assert result.returncode == 0


def test_an_unrecognized_response_is_not_treated_as_ready(clock, ncbi, query_filepath, tmp_path):
    """
    An HTML error page from NCBI used to fall through to READY, which wrote an empty results file
    and reported success. It must not be mistaken for a finished search.
    """
    ncbi(rtoe=1, statuses=["<html><h1>Temporarily unavailable</h1></html>"])

    out = tmp_path / "out.tsv"
    result = run_www(query_filepath, out, timeout_seconds=100_000)

    assert result.returncode != 0
    assert not (out.exists() and out.read_text() == "" and result.returncode == 0)


def test_a_long_wait_reports_progress(capsys, clock, ncbi, query_filepath, tmp_path):
    """
    A long wait is not silent, so a queued search is distinguishable from a hung process on an
    unattended worker.
    """
    ncbi(rtoe=9168, statuses=["WAITING"])

    run_www(query_filepath, tmp_path / "out.tsv", timeout_seconds=600)

    reports = [
        line for line in capsys.readouterr().out.splitlines() if "still waiting on NCBI" in line
    ]
    assert len(reports) > 1
