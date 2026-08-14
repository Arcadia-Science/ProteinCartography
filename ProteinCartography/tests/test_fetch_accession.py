"""Unit tests for AlphaFold PDB fetch (isoform URLs + no error-body writes)."""

from __future__ import annotations
import pathlib
import sys
from unittest.mock import MagicMock, Mock

sys.modules.setdefault("bioservices", MagicMock())

import fetch_accession  # noqa: E402


class _FakeResponse:
    def __init__(self, *, ok=True, status_code=200, text="", json_data=None):
        self.ok = ok
        self.status_code = status_code
        self.text = text
        self._json_data = json_data

    def json(self):
        return self._json_data


def test_fetch_pdb_uses_isoform_pdburl(tmp_path: pathlib.Path):
    pdb_body = "ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00 80.00\n"
    isoform_url = "https://alphafold.ebi.ac.uk/files/AF-Q9Y6V0-3-F1-model_v6.pdb"
    guessed_f1 = "https://alphafold.ebi.ac.uk/files/AF-Q9Y6V0-F1-model_v6.pdb"

    def fake_get(url, timeout=None):
        if "api/prediction/Q9Y6V0" in url:
            return _FakeResponse(json_data=[{"pdbUrl": isoform_url}])
        if url == isoform_url:
            return _FakeResponse(text=pdb_body)
        if url == guessed_f1:
            return _FakeResponse(ok=False, status_code=404, text="<Error>NoSuchKey</Error>")
        raise AssertionError(f"unexpected url {url}")

    session = Mock()
    session.get.side_effect = fake_get

    fetch_accession.fetch_pdb("Q9Y6V0", str(tmp_path), session=session)
    out = tmp_path / "Q9Y6V0.pdb"
    assert out.read_text() == pdb_body
    assert "<Error>" not in out.read_text()


def test_fetch_pdb_does_not_write_error_xml(tmp_path: pathlib.Path):
    error = "<Error>NoSuchKey</Error>"

    def fake_get(url, timeout=None):
        if "api/prediction" in url:
            return _FakeResponse(ok=False, status_code=404, text=error)
        return _FakeResponse(ok=False, status_code=404, text=error)

    session = Mock()
    session.get.side_effect = fake_get

    fetch_accession.fetch_pdb("Q9Y6V0", str(tmp_path), session=session)
    assert not (tmp_path / "Q9Y6V0.pdb").exists()


def test_looks_like_pdb_rejects_error_body():
    assert fetch_accession._looks_like_pdb("<Error>NoSuchKey</Error>") is False
    assert fetch_accession._looks_like_pdb("ATOM      1  N") is True
