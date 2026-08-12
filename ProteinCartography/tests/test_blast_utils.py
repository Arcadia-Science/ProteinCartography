"""
Tests for the blast search that `run_blast.py` runs against NCBI, and for detecting when it fails.

These are unit tests: every HTTP request to NCBI is mocked and the clock is faked, so unlike the
other tests in this directory they do not need network access, conda environments, or real
API responses, and they run without ever sleeping.
"""

import subprocess

import blast_utils
import constants
import pytest


def make_result(returncode=0, stderr=b""):
    """
    Construct a `CompletedProcess` that stands in for the result of a call to `blastp`.
    """
    return subprocess.CompletedProcess(args=[], returncode=returncode, stdout=b"", stderr=stderr)


def test_successful_call_did_not_fail():
    assert not blast_utils.blast_call_failed(make_result())


def test_nonzero_exit_code_is_a_failure():
    assert blast_utils.blast_call_failed(make_result(returncode=1))


def test_error_on_stderr_is_a_failure_despite_a_zero_exit_code():
    """
    Tests the case that motivates checking stderr at all: when the remote server refuses to
    queue the request, `blastp -remote` reports the error but still exits with a status of zero.
    """
    stderr = (
        b"Error: [blastp] bad_request: Could not queue request: DB operation failed."
        b"(ERR:-99  )((Severe Error) DB Put Request error: sp_NewRequestEx failed)\n"
    )
    assert blast_utils.blast_call_failed(make_result(returncode=0, stderr=stderr))


def test_warning_on_stderr_is_not_a_failure():
    """
    Tests that stderr output which is not an error does not cause a successful call
    to be treated as a failure.
    """
    stderr = b"Warning: [blastp] Examining 5 or more matches is recommended\n"
    assert not blast_utils.blast_call_failed(make_result(returncode=0, stderr=stderr))


def test_stderr_may_be_a_string_or_none():
    """
    Tests that stderr is handled whether it is bytes, str, or absent, since the stub
    and the real call do not construct it the same way.
    """
    assert blast_utils.blast_call_failed(make_result(stderr="Error: something went wrong"))
    assert not blast_utils.blast_call_failed(make_result(stderr="all good"))
    assert not blast_utils.blast_call_failed(make_result(stderr=None))


def test_the_stub_does_not_look_like_a_failure(tmp_path):
    """
    Tests that the result returned by the blast stub is not mistaken for a failed call.
    """
    stub_results_filepath = tmp_path / "stub.tsv"
    stub_results_filepath.write_text("qseqid\tsseqid\n")
    out_filepath = tmp_path / "out.tsv"

    result = blast_utils.run_blast_stub(str(stub_results_filepath), str(out_filepath))

    assert not blast_utils.blast_call_failed(result)
    assert out_filepath.read_text() == "qseqid\tsseqid\n"


# ------------------------------------------------------------------------------------------------
# Canned responses from NCBI's BLAST URL API.
# These reproduce the shape of the real responses, including the HTML that NCBI wraps them in.
# ------------------------------------------------------------------------------------------------

REQUEST_ID = "M8K4ARMH013"

# The response to a `CMD=Put` request that was successfully queued.
PUT_RESPONSE = f"""<html><head><title>NCBI Blast</title></head><body>
<!--QBlastInfoBegin
    RID = {REQUEST_ID}
    RTOE = 27
QBlastInfoEnd
-->
<p>Your request has been submitted.</p>
</body></html>
"""

# The response NCBI returns when the BLAST CGI itself is not accepting requests.
# Note that this is returned in place of the QBlastInfo block, with a status of 200.
TEMPORARILY_UNAVAILABLE_RESPONSE = """<html><body>
<h1>Temporarily unavailable</h1>
<p>The BLAST service is temporarily unavailable. Please try again later.</p>
</body></html>
"""

# The tabular results used by the tests below. These are copied from a real `blastp -remote`
# results file, truncated to the twelve fields that NCBI's URL API returns.
TABULAR_RESPONSE = (
    "# blastp\n"
    "# Iteration: 0\n"
    "# Query: sp|P60709|ACTB_HUMAN Actin, cytoplasmic 1\n"
    f"# RID: {REQUEST_ID}\n"
    "# Database: nr\n"
    "# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, "
    "q. start, q. end, s. start, s. end, evalue, bit score\n"
    "# 2 hits found\n"
    "sp|P60709|ACTB_HUMAN\tref|XP_052610122.1|\t100.000\t375\t0\t0\t1\t375\t46\t420\t0.0\t787\n"
    "sp|P60709|ACTB_HUMAN\tgb|TEA41296.1|\t99.733\t375\t1\t0\t1\t375\t36\t410\t0.0\t785\n"
)

# A tabular response for a query that legitimately has no hits.
TABULAR_RESPONSE_WITH_NO_HITS = (
    "# blastp\n"
    "# Iteration: 0\n"
    "# Query: sp|P60709|ACTB_HUMAN Actin, cytoplasmic 1\n"
    f"# RID: {REQUEST_ID}\n"
    "# Database: nr\n"
    "# 0 hits found\n"
)


def search_info_response(status):
    """
    Build the response to a `CMD=Get&FORMAT_OBJECT=SearchInfo` request reporting `status`.
    """
    return (
        "<html><body><pre>\n"
        "QBlastInfoBegin\n"
        f"\tStatus={status}\n"
        "QBlastInfoEnd\n"
        "</pre></body></html>\n"
    )


class FakeResponse:
    """
    Stand in for a `requests.Response`.
    """

    def __init__(self, text, status_code=200):
        self.text = text
        self.status_code = status_code


class FakeNcbi:
    """
    Stand in for NCBI's BLAST URL API.

    Requests are dispatched on the `CMD` parameter in the same way NCBI dispatches them, and every
    request is recorded so that the tests can assert on what was sent.

    Args:
        put_response: the response to the `CMD=Put` request that submits the search.
        statuses: the statuses to report, in order, to successive polling requests.
        results_response: the response to the `CMD=Get` request that retrieves the results.
    """

    def __init__(self, put_response=PUT_RESPONSE, statuses=("READY",), results_response=None):
        self.put_response = self._as_response(put_response)
        self.statuses = list(statuses)
        self.results_response = self._as_response(
            TABULAR_RESPONSE if results_response is None else results_response
        )
        self.requests = []

    @staticmethod
    def _as_response(response):
        return response if isinstance(response, FakeResponse) else FakeResponse(response)

    def post(self, url, data=None, timeout=None):
        self.requests.append(data)

        if data["CMD"] == "Put":
            return self.put_response

        if data.get("FORMAT_OBJECT") == "SearchInfo":
            # Repeat the last status once the scripted ones are exhausted, so that a test can
            # describe a search that waits forever without listing a status for every poll.
            status = self.statuses.pop(0) if len(self.statuses) > 1 else self.statuses[0]
            return self._as_response(search_info_response(status))

        return self.results_response

    @property
    def polling_requests(self):
        return [
            request for request in self.requests if request.get("FORMAT_OBJECT") == "SearchInfo"
        ]


class FakeClock:
    """
    A monotonic clock that only advances when something sleeps.

    This lets the tests assert on how long `wait_for_search` waited, and lets them exercise
    multi-hour timeouts instantly.
    """

    def __init__(self):
        self.now = 1000.0
        self.sleeps = []

    def monotonic(self):
        return self.now

    def sleep(self, seconds):
        self.sleeps.append(seconds)
        self.now += seconds


@pytest.fixture
def clock(monkeypatch):
    fake_clock = FakeClock()
    monkeypatch.setattr(blast_utils.time, "monotonic", fake_clock.monotonic)
    monkeypatch.setattr(blast_utils.time, "sleep", fake_clock.sleep)
    return fake_clock


@pytest.fixture
def ncbi(monkeypatch):
    """
    Replace the HTTP requests made to NCBI with a `FakeNcbi` that the test can configure.
    """

    def install(**kwargs):
        fake_ncbi = FakeNcbi(**kwargs)
        monkeypatch.setattr(blast_utils.requests, "post", fake_ncbi.post)
        return fake_ncbi

    return install


@pytest.fixture
def query_filepath(tmp_path):
    filepath = tmp_path / "P60709.fasta"
    filepath.write_text(">sp|P60709|ACTB_HUMAN Actin, cytoplasmic 1\nMDDDIAALVVDNGSGMCKAGF\n")
    return filepath


def run_blast(query_filepath, out_filepath, timeout_seconds=1800, **kwargs):
    """
    Call `blast_utils.run_blast` with the parameters the pipeline uses by default.
    """
    return blast_utils.run_blast(
        query=str(query_filepath),
        out=str(out_filepath),
        max_target_seqs=3000,
        outfmt=constants.BLAST_OUTFMT,
        word_size=5,
        evalue=1.0,
        timeout_seconds=timeout_seconds,
        **kwargs,
    )


# ------------------------------------------------------------------------------------------------
# The happy path: submit, poll, fetch.
# ------------------------------------------------------------------------------------------------


def test_a_successful_search_writes_the_expected_columns(clock, ncbi, query_filepath, tmp_path):
    """
    Tests the whole submit-poll-fetch sequence, and checks that the results are written in the
    exact column layout that `constants.BLAST_OUTFMT` asks for, since `extract_blast_hits.py`
    reads the results back using that same layout.
    """
    ncbi()
    out_filepath = tmp_path / "results.tsv"

    result = run_blast(query_filepath, out_filepath)

    assert not blast_utils.blast_call_failed(result)

    rows = [line.split("\t") for line in out_filepath.read_text().splitlines()]
    assert len(rows) == 2
    for row in rows:
        assert len(row) == len(constants.BLAST_OUTPUT_FIELDS)

    assert rows[0] == [
        "sp|P60709|ACTB_HUMAN",
        "ref|XP_052610122.1|",
        "100.000",
        "375",
        "0",
        "0",
        "1",
        "375",
        "46",
        "420",
        "0.0",
        "787",
        # 'sacc' and 'saccver' are derived from the subject sequence id ...
        "XP_052610122",
        "XP_052610122.1",
        # ... and the remaining fields are not available from NCBI's URL API.
        "N/A",
        "N/A",
        "N/A",
    ]

    # `extract_blast_hits.py` reads only this column, so it is the one that has to be right.
    assert [row[constants.BLAST_OUTPUT_FIELDS.index("sacc")] for row in rows] == [
        "XP_052610122",
        "TEA41296",
    ]


def test_the_search_is_submitted_with_the_pipelines_parameters(
    clock, ncbi, query_filepath, tmp_path
):
    """
    Tests that the blastp arguments the pipeline is configured with are translated into the
    parameters of NCBI's URL API, and that the search identifies itself to NCBI.
    """
    fake_ncbi = ncbi()

    run_blast(query_filepath, tmp_path / "results.tsv", email="someone@example.com")

    submission = fake_ncbi.requests[0]
    assert submission["CMD"] == "Put"
    assert submission["PROGRAM"] == "blastp"
    assert submission["DATABASE"] == "nr"
    assert submission["HITLIST_SIZE"] == 3000
    assert submission["WORD_SIZE"] == 5
    assert submission["EXPECT"] == 1.0
    assert "ACTB_HUMAN" in submission["QUERY"]
    assert submission["tool"] == "proteincartography"
    assert submission["email"] == "someone@example.com"

    # Every subsequent request has to carry the request id NCBI handed out.
    for request in fake_ncbi.requests[1:]:
        assert request["RID"] == REQUEST_ID
        assert request["tool"] == "proteincartography"


def test_a_search_with_no_hits_is_not_a_failure(clock, ncbi, query_filepath, tmp_path):
    """
    Tests that a query which legitimately has no hits produces an empty results file rather than
    an error, which is the behavior `blast_call_failed` documents.
    """
    ncbi(results_response=TABULAR_RESPONSE_WITH_NO_HITS)
    out_filepath = tmp_path / "results.tsv"

    result = run_blast(query_filepath, out_filepath)

    assert not blast_utils.blast_call_failed(result)
    assert out_filepath.read_text() == ""


# ------------------------------------------------------------------------------------------------
# Polling.
# ------------------------------------------------------------------------------------------------


def test_a_waiting_search_is_polled_until_it_is_ready(clock, ncbi, query_filepath, tmp_path):
    """
    Tests that a search which is not ready immediately is polled until it is.
    """
    fake_ncbi = ncbi(statuses=["WAITING", "WAITING", "READY"])
    out_filepath = tmp_path / "results.tsv"

    result = run_blast(query_filepath, out_filepath)

    assert not blast_utils.blast_call_failed(result)
    assert len(fake_ncbi.polling_requests) == 3
    assert out_filepath.read_text() != ""


def test_polling_waits_for_the_estimate_and_then_polls_once_a_minute(
    clock, ncbi, query_filepath, tmp_path
):
    """
    Tests NCBI's polling etiquette: the first poll happens only after the estimated time to
    completion (RTOE) that NCBI returned, and subsequent polls are at least 60 seconds apart.
    """
    ncbi(statuses=["WAITING", "WAITING", "READY"])

    run_blast(query_filepath, tmp_path / "results.tsv")

    # The RTOE in the canned `CMD=Put` response is 27 seconds.
    assert clock.sleeps[0] == 27
    assert all(delay >= blast_utils.MIN_POLLING_INTERVAL_SECONDS for delay in clock.sleeps[1:])


# ------------------------------------------------------------------------------------------------
# Failures.
# ------------------------------------------------------------------------------------------------


def test_a_failed_search_is_reported_as_a_failure(clock, ncbi, query_filepath, tmp_path):
    """
    Tests the case where NCBI runs the search and reports that it failed.
    """
    ncbi(statuses=["FAILED"])

    result = run_blast(query_filepath, tmp_path / "results.tsv")

    assert blast_utils.blast_call_failed(result)
    assert "failed" in result.stderr
    assert REQUEST_ID in result.stderr


def test_an_unknown_request_id_is_reported_as_a_failure(clock, ncbi, query_filepath, tmp_path):
    """
    Tests the case where NCBI no longer recognizes the request id, which happens when the id has
    expired and also when the search was silently never queued.
    """
    ncbi(statuses=["UNKNOWN"])

    result = run_blast(query_filepath, tmp_path / "results.tsv")

    assert blast_utils.blast_call_failed(result)
    assert "no longer recognizes" in result.stderr


def test_an_unexpected_status_is_reported_as_a_failure(clock, ncbi, query_filepath, tmp_path):
    """
    Tests that a status NCBI has never been observed to return is treated as a failure rather
    than being mistaken for either success or a search still in progress.
    """
    ncbi(statuses=["SOMETHING_NEW"])

    result = run_blast(query_filepath, tmp_path / "results.tsv")

    assert blast_utils.blast_call_failed(result)
    assert "SOMETHING_NEW" in result.stderr


def test_a_response_without_a_request_id_is_a_failure(clock, ncbi, query_filepath, tmp_path):
    """
    Tests the case where NCBI accepts the submission but returns a body with no QBlastInfo block,
    which is the kind of failure that `blastp -remote` reports while still exiting with zero.
    """
    ncbi(put_response="<html><body><p>Something unexpected</p></body></html>")

    result = run_blast(query_filepath, tmp_path / "results.tsv")

    assert blast_utils.blast_call_failed(result)
    assert "did not return a request id" in result.stderr


def test_temporarily_unavailable_is_reported_as_a_failure(clock, ncbi, query_filepath, tmp_path):
    """
    Tests the response NCBI returns when its BLAST service is refusing requests. NCBI returns
    this with a status of 200, so it has to be detected from the body rather than the status code.
    """
    ncbi(put_response=TEMPORARILY_UNAVAILABLE_RESPONSE)

    result = run_blast(query_filepath, tmp_path / "results.tsv")

    assert blast_utils.blast_call_failed(result)
    # The message has to name the actual problem, not just say that there was no request id.
    assert "Temporarily unavailable" in result.stderr
    # The HTML tags are stripped so that the message stays readable on one line.
    assert "<h1>" not in result.stderr


def test_an_http_error_status_is_reported_as_a_failure(clock, ncbi, query_filepath, tmp_path):
    """
    Tests that an error status from NCBI is reported rather than being parsed as a response.
    """
    ncbi(put_response=FakeResponse("<html><body>Service Unavailable</body></html>", 503))

    result = run_blast(query_filepath, tmp_path / "results.tsv")

    assert blast_utils.blast_call_failed(result)
    assert "503" in result.stderr


def test_a_non_tabular_results_response_is_a_failure(clock, ncbi, query_filepath, tmp_path):
    """
    Tests the case where the search becomes ready but NCBI returns an error page instead of the
    results. Without this check the pipeline would write the error page to the results file.
    """
    ncbi(results_response="<html><body><p>cannot retrieve results</p></body></html>")
    out_filepath = tmp_path / "results.tsv"

    result = run_blast(query_filepath, out_filepath)

    assert blast_utils.blast_call_failed(result)
    assert "did not return tabular results" in result.stderr
    # The results file must not be written at all when the search failed.
    assert not out_filepath.exists()


def test_a_network_error_is_reported_as_a_failure(monkeypatch, clock, query_filepath, tmp_path):
    """
    Tests that a connection error or a timeout on an individual HTTP request is reported as a
    failed search rather than raising out of `run_blast`.
    """

    def raise_connection_error(*args, **kwargs):
        raise blast_utils.requests.ConnectionError("connection refused")

    monkeypatch.setattr(blast_utils.requests, "post", raise_connection_error)

    result = run_blast(query_filepath, tmp_path / "results.tsv")

    assert blast_utils.blast_call_failed(result)
    assert "The request to NCBI failed" in result.stderr


def test_a_missing_query_file_is_reported_as_a_failure(clock, ncbi, tmp_path):
    """
    Tests that a missing input file is reported the same way as any other failure,
    rather than raising a bare `FileNotFoundError` out of `run_blast`.
    """
    ncbi()

    result = run_blast(tmp_path / "does-not-exist.fasta", tmp_path / "results.tsv")

    assert blast_utils.blast_call_failed(result)


# ------------------------------------------------------------------------------------------------
# The timeout, which is the reason this code calls NCBI's URL API instead of `blastp -remote`.
# ------------------------------------------------------------------------------------------------


def test_a_search_that_never_finishes_times_out(clock, ncbi, query_filepath, tmp_path):
    """
    Tests that a search which stays in NCBI's queue is abandoned once the timeout is exceeded,
    rather than being polled forever as `blastp -remote` does.
    """
    fake_ncbi = ncbi(statuses=["WAITING"])
    out_filepath = tmp_path / "results.tsv"

    result = run_blast(query_filepath, out_filepath, timeout_seconds=300)

    assert blast_utils.blast_call_failed(result)
    assert not out_filepath.exists()

    # The timeout is honored exactly: no more time is spent waiting than the caller allowed.
    assert sum(clock.sleeps) == 300
    assert len(fake_ncbi.polling_requests) > 0


def test_the_timeout_message_reports_the_wait_and_the_estimate(
    clock, ncbi, query_filepath, tmp_path
):
    """
    Tests that the timeout says how long it waited and what NCBI's own estimate was, so that a
    user can tell the difference between a busy queue and a broken pipeline. The estimate used
    here is the ~2.5 hour one that NCBI returned while this failure was being investigated.
    """
    ncbi(put_response=PUT_RESPONSE.replace("RTOE = 27", "RTOE = 9168"), statuses=["WAITING"])

    result = run_blast(query_filepath, tmp_path / "results.tsv", timeout_seconds=1800)

    assert blast_utils.blast_call_failed(result)
    assert "1800s" in result.stderr
    assert "9168s" in result.stderr
    assert "blast_timeout_seconds" in result.stderr


def test_an_estimate_longer_than_the_timeout_does_not_overshoot_it(
    clock, ncbi, query_filepath, tmp_path
):
    """
    Tests that an estimate longer than the whole timeout does not cause the first wait to run
    past the deadline. NCBI's estimate is routinely far larger than the timeout the caller set,
    so this is the common case rather than an edge case.
    """
    ncbi(put_response=PUT_RESPONSE.replace("RTOE = 27", "RTOE = 9168"), statuses=["WAITING"])

    result = run_blast(query_filepath, tmp_path / "results.tsv", timeout_seconds=600)

    assert blast_utils.blast_call_failed(result)
    assert clock.sleeps == [600]


def test_a_missing_estimate_falls_back_to_the_polling_interval(
    clock, ncbi, query_filepath, tmp_path
):
    """
    Tests that a response without a parseable RTOE still polls, rather than crashing or polling
    immediately. NCBI's estimate is advisory and is occasionally absent.
    """
    ncbi(put_response=PUT_RESPONSE.replace("RTOE = 27", "RTOE = unknown"), statuses=["READY"])

    result = run_blast(query_filepath, tmp_path / "results.tsv")

    assert not blast_utils.blast_call_failed(result)
    assert clock.sleeps[0] == blast_utils.MIN_POLLING_INTERVAL_SECONDS


# ------------------------------------------------------------------------------------------------
# Parsing and formatting helpers.
# ------------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "sseqid, expected_saccver, expected_sacc",
    [
        # The two-field form used by RefSeq and GenBank subjects in `nr`.
        ("ref|XP_052610122.1|", "XP_052610122.1", "XP_052610122"),
        ("gb|TEA41296.1|", "TEA41296.1", "TEA41296"),
        ("emb|CAI5791831.1|", "CAI5791831.1", "CAI5791831"),
        # The three-field form used by SwissProt subjects, where the accession is still second.
        ("sp|P60709|ACTB_HUMAN", "P60709", "P60709"),
        # The bare form, which a blast database built with a recent makeblastdb can return.
        ("XP_052610122.1", "XP_052610122.1", "XP_052610122"),
    ],
)
def test_the_accession_is_parsed_from_the_subject_sequence_id(
    sseqid, expected_saccver, expected_sacc
):
    saccver = blast_utils.parse_accession(sseqid)
    assert saccver == expected_saccver
    assert blast_utils.strip_accession_version(saccver) == expected_sacc


def test_comment_lines_are_dropped_from_the_results():
    """
    Tests that NCBI's comment lines are stripped, since `blastp -outfmt 6` does not emit them
    and `extract_blast_hits.py` reads the file with no comment handling.
    """
    results = blast_utils.format_results(TABULAR_RESPONSE, constants.BLAST_OUTFMT)
    assert "#" not in results
    assert len(results.splitlines()) == 2


def test_a_truncated_tabular_row_is_an_error():
    """
    Tests that a row with fewer fields than expected is rejected rather than silently producing
    misaligned columns in the results file.
    """
    with pytest.raises(blast_utils.BlastApiError, match="tabular row"):
        blast_utils.format_results("# blastp\nqueryid\tsubjectid\t100.0\n", constants.BLAST_OUTFMT)


def test_a_custom_output_format_is_respected():
    """
    Tests that the fields written are the ones named by `outfmt`, so that the results stay
    consistent with whatever `extract_blast_hits.py` is told to read them back with.
    """
    results = blast_utils.format_results(TABULAR_RESPONSE, "6 sacc evalue")
    assert results.splitlines()[0] == "XP_052610122\t0.0"


def test_the_status_is_parsed_from_the_qblastinfo_block():
    assert blast_utils.parse_qblast_info(search_info_response("WAITING"), "Status") == "WAITING"
    assert blast_utils.parse_qblast_info(PUT_RESPONSE, "RID") == REQUEST_ID
    assert blast_utils.parse_qblast_info(PUT_RESPONSE, "RTOE") == "27"
    assert blast_utils.parse_qblast_info(PUT_RESPONSE, "Status") is None


# ------------------------------------------------------------------------------------------------
# The contact email address that NCBI asks automated clients to send.
# ------------------------------------------------------------------------------------------------


def test_the_email_falls_back_to_the_default(monkeypatch):
    monkeypatch.delenv(blast_utils.EMAIL_ENV_VAR, raising=False)
    assert blast_utils.resolve_email(None) == blast_utils.DEFAULT_EMAIL


def test_the_email_can_be_set_from_the_config(monkeypatch):
    monkeypatch.delenv(blast_utils.EMAIL_ENV_VAR, raising=False)
    assert blast_utils.resolve_email("config@example.com") == "config@example.com"


def test_the_env_variable_overrides_the_email_from_the_config(monkeypatch):
    monkeypatch.setenv(blast_utils.EMAIL_ENV_VAR, "env@example.com")
    assert blast_utils.resolve_email("config@example.com") == "env@example.com"


def test_no_personal_email_address_is_hardcoded():
    """
    Tests that the default address is a project address rather than an individual's, since it is
    sent to NCBI with every search made by every user of the pipeline.
    """
    assert blast_utils.DEFAULT_EMAIL.startswith("proteincartography@")


# ------------------------------------------------------------------------------------------------
# The stub used by the smoke test.
# ------------------------------------------------------------------------------------------------


def test_the_stub_short_circuits_the_call_to_ncbi(monkeypatch, query_filepath, tmp_path):
    """
    Tests that setting the stub env variable prevents any request from being made to NCBI,
    which is what lets `make smoke-test` run the whole pipeline without network access.
    """

    def fail_if_called(*args, **kwargs):
        raise AssertionError("run_blast made a request to NCBI even though the stub was set")

    monkeypatch.setattr(blast_utils.requests, "post", fail_if_called)

    stub_results_filepath = tmp_path / "stub.tsv"
    stub_results_filepath.write_text("qseqid\tsseqid\n")
    monkeypatch.setenv(blast_utils.STUB_RESULTS_FILEPATH_ENV_VAR, str(stub_results_filepath))

    out_filepath = tmp_path / "results.tsv"
    result = run_blast(query_filepath, out_filepath)

    assert not blast_utils.blast_call_failed(result)
    assert out_filepath.read_text() == "qseqid\tsseqid\n"
