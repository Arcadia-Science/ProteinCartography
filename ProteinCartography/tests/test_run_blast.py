"""Unit tests for run_blast empty-result handling."""

from __future__ import annotations
import pathlib
import sys

import blast_utils
import pytest
import run_blast


def _argv(tmp_path: pathlib.Path, out: pathlib.Path) -> list[str]:
    query = tmp_path / "q.fa"
    query.write_text(">q\nM\n")
    return [
        "run_blast",
        "--query",
        str(query),
        "--out",
        str(out),
        "--max_target_seqs",
        "10",
        "--word_size",
        "3",
        "--word_size_backoff",
        "2",
        "--num_attempts",
        "1",
        "--evalue",
        "0.05",
    ]


def test_run_blast_empty_tsv_hard_fails_by_default(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
):
    out = tmp_path / "blast.tsv"

    def fake_run_blast(**kwargs):
        pathlib.Path(kwargs["out"]).write_text("")
        return blast_utils.BlastResult(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(run_blast.blast_utils, "run_blast", fake_run_blast)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, out))
    with pytest.raises(SystemExit) as exc:
        run_blast.main()
    assert exc.value.code != 0


def test_run_blast_empty_tsv_soft_fail_exits_zero(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
):
    out = tmp_path / "blast.tsv"
    monkeypatch.setenv("PC_BLAST_SOFT_FAIL", "1")

    def fake_run_blast(**kwargs):
        pathlib.Path(kwargs["out"]).write_text("")
        return blast_utils.BlastResult(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(run_blast.blast_utils, "run_blast", fake_run_blast)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, out))
    with pytest.raises(SystemExit) as exc:
        run_blast.main()
    assert exc.value.code == 0
