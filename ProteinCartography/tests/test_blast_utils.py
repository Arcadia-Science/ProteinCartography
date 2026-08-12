"""
Tests for detecting failed calls to `blastp`.

These are unit tests and do not call `blastp`, so unlike the other tests in this
directory they do not need network access, conda environments, or mocked API responses.
"""

import subprocess

import blast_utils


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
