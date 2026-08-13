"""
Tests for the retry loop in `run_blast.py`, and in particular for its timeout budget.

The blast search itself is replaced with a stand-in, so these tests make no network calls and
never sleep; what they exercise is how the attempts share the one timeout the user configured.
"""

import subprocess
import sys

import pytest
import run_blast

# The arguments the `run_blast` snakemake rule passes, minus the ones each test varies.
BASE_ARGV = [
    "run_blast.py",
    "--query",
    "input/P60709.fasta",
    "--out",
    "output/P60709.blast_results.tsv",
    "--max_target_seqs",
    "3000",
    "--word_size",
    "5",
    "--word_size_backoff",
    "6",
    "--evalue",
    "1.0",
]


class FakeSearch:
    """
    Stand in for `blast_utils.run_blast`.

    Each call consumes the time the scripted outcome says it takes, so that the tests can check
    how the timeout budget is spent across attempts.

    Args:
        outcomes: one (duration, succeeded) pair per attempt. The last pair repeats.
    """

    def __init__(self, clock, outcomes):
        self.clock = clock
        self.outcomes = list(outcomes)
        self.calls = []

    def __call__(self, **kwargs):
        self.calls.append(kwargs)

        index = min(len(self.calls) - 1, len(self.outcomes) - 1)
        duration, succeeded = self.outcomes[index]

        # A search never runs for longer than the budget it was given.
        duration = min(duration, kwargs["timeout_seconds"])
        self.clock.advance(duration)

        if succeeded:
            return subprocess.CompletedProcess(args=[], returncode=0, stdout=b"", stderr="")
        return subprocess.CompletedProcess(
            args=[], returncode=1, stdout=b"", stderr="Error: [blastp] the search timed out"
        )

    @property
    def timeouts(self):
        """
        The budget handed to each attempt, rounded for readability.
        """
        return [round(call["timeout_seconds"]) for call in self.calls]


class FakeClock:
    def __init__(self):
        self.now = 1000.0

    def monotonic(self):
        return self.now

    def advance(self, seconds):
        self.now += seconds


@pytest.fixture
def clock(monkeypatch):
    fake_clock = FakeClock()
    monkeypatch.setattr(run_blast.time, "monotonic", fake_clock.monotonic)
    return fake_clock


@pytest.fixture
def run_main(monkeypatch, clock):
    """
    Run `run_blast.main` with a stand-in search, and return it along with the exit that resulted.
    """

    def run(outcomes, num_attempts=3, timeout_seconds=1800):
        search = FakeSearch(clock, outcomes)
        monkeypatch.setattr(run_blast.blast_utils, "run_blast", search)
        monkeypatch.setattr(
            sys,
            "argv",
            BASE_ARGV
            + [
                "--num_attempts",
                str(num_attempts),
                "--timeout_seconds",
                str(timeout_seconds),
            ],
        )

        with pytest.raises(SystemExit) as exit_info:
            run_blast.main()

        return search, exit_info.value

    return run


def test_a_search_that_succeeds_first_time_exits_zero(run_main):
    search, exit_info = run_main([(120, True)])

    assert exit_info.code == 0
    assert len(search.calls) == 1


def test_a_failed_search_is_retried_with_the_backoff_word_size(run_main):
    search, exit_info = run_main([(60, False), (60, True)])

    assert exit_info.code == 0
    assert [call["word_size"] for call in search.calls] == [5, 6]


def test_the_timeout_is_a_total_budget_rather_than_a_per_attempt_one(run_main):
    """
    Tests the property that makes the configured timeout mean what it says: each attempt is
    given only what is left of the budget, so retrying cannot extend the total wait.
    """
    search, exit_info = run_main([(600, False)], num_attempts=3, timeout_seconds=1800)

    # Three attempts of 600s each exactly consume the 1800s budget.
    assert search.timeouts == [1800, 1200, 600]
    assert exit_info.code != 0


def test_the_budget_stops_further_attempts_once_it_is_exhausted(run_main):
    """
    Tests that attempts stop when the budget runs out, even though attempts remain. Without this
    the worst case would be `num_attempts` times the configured timeout.
    """
    search, exit_info = run_main([(1800, False)], num_attempts=3, timeout_seconds=1800)

    # The first attempt used the whole budget, so the other two are never made.
    assert len(search.calls) == 1
    assert exit_info.code != 0


def test_a_retry_does_not_reset_the_budget(run_main):
    """
    Tests that the second attempt is given the remaining budget and not a fresh one.
    """
    search, _ = run_main([(1000, False)], num_attempts=3, timeout_seconds=1800)

    assert search.timeouts == [1800, 800]


def test_the_failure_message_reports_the_total_time_and_the_budget(run_main):
    _, exit_info = run_main([(600, False)], num_attempts=3, timeout_seconds=1800)

    message = str(exit_info.code)
    assert "3 attempt(s)" in message
    assert "1800s total budget" in message
    assert "blast_timeout_seconds" in message
    # The time actually spent is reported, so a user can see where the budget went.
    assert "1800s of the" in message


def test_an_exhausted_budget_is_reported_even_if_no_attempt_was_made(run_main):
    """
    Tests the edge case of a budget so small that it is gone before the first attempt.
    """
    search, exit_info = run_main([(10, False)], num_attempts=3, timeout_seconds=0)

    assert len(search.calls) == 0
    assert "0 attempt(s)" in str(exit_info.code)
