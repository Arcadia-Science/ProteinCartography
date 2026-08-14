#!/usr/bin/env python
"""Fetch FASTA / AlphaFold PDB for a UniProt accession.

Prefer AlphaFold ``model_v6``, fall back to ``v4``, then the AFDB prediction
API ``pdbUrl``. Never persist HTML/XML error bodies as PDB files (retired AFDB
URLs previously wrote S3 ``NoSuchKey`` XML that poisoned downstream Foldseek).
"""

from __future__ import annotations
import argparse
import os
from pathlib import Path

from api_utils import UniProtWithExpBackoff, session_with_retry

__all__ = ["fetch_fasta", "fetch_pdb"]

AFDB_FILES = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_{ver}.pdb"
AFDB_PREDICTION = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
# Newest first; v4 kept as a short-lived fallback for mirrors that lag.
AFDB_VERSIONS = ("v6", "v4")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--accession",
        required=True,
        nargs="+",
        help="UniprotKB accession of target.",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output directory for resulting files."
    )
    parser.add_argument(
        "-f",
        "--format",
        nargs="+",
        default=["fasta", "pdb"],
        help=(
            "Formats to acquire.\n"
            'If "fasta", requests the FASTA sequence using bioservices UniProt.\n'
            'If "pdb", downloads the pdb from AlphaFold.\nCan accept multiple arguments.'
        ),
    )
    args = parser.parse_args()
    return args


def fetch_fasta(accession: str, output_dir: str):
    """
    Fetches a FASTA file from Uniprot, given an accession. Places the file in the output_dir.

    Args:
        accession (str): a valid UniprotKB accession.
        output_dir (str): path to the output directory.
        File will be saved as "{output_dir}/{accession}.fasta".
    """
    uniprot = UniProtWithExpBackoff()
    output_path = Path(output_dir) / (accession + ".fasta")

    if not os.path.exists(output_path):
        res = uniprot.retrieve(accession, frmt="fasta")
        with open(output_path, "w+") as file:
            file.write(res)


def _looks_like_pdb(text: str) -> bool:
    if not text or not text.strip():
        return False
    # AF / ESMFold PDBs always have ATOM/HETATM; error pages have <Error> or HTML.
    if "<Error>" in text or "<html" in text.lower():
        return False
    return "ATOM" in text or "HETATM" in text


def _pdb_from_response(result) -> str | None:
    """Return PDB text only for a successful HTTP body that looks like a PDB.

    AlphaFold file URLs 404 with a ~127-byte ``<Error>`` body for accessions
    whose real entry is an isoform (e.g. Q9Y6V0 → AF-Q9Y6V0-3-F1). Never treat
    that as a structure.
    """
    if result is None:
        return None
    text = getattr(result, "text", "") or ""
    if getattr(result, "ok", False) and _looks_like_pdb(text):
        return text
    return None


def _download_pdb_text(accession: str, session) -> str | None:
    # Prefer the prediction API: guessed AF-{acc}-F1 URLs 404 when the AFDB
    # entry uses an isoform id (Q9Y6V0 → AF-Q9Y6V0-3-F1-model_v6.pdb).
    try:
        meta = session.get(AFDB_PREDICTION.format(acc=accession), timeout=60)
        if getattr(meta, "ok", False):
            entries = meta.json()
            if isinstance(entries, list):
                for entry in entries:
                    pdb_url = (entry or {}).get("pdbUrl") or (entry or {}).get("pdbFile")
                    if not pdb_url:
                        continue
                    result = session.get(pdb_url, timeout=120)
                    text = _pdb_from_response(result)
                    if text:
                        return text
                    print(
                        f"AlphaFold pdbUrl miss for {accession}: "
                        f"status={getattr(result, 'status_code', '?')} "
                        f"bytes={len(getattr(result, 'text', '') or '')}"
                    )
    except Exception as exc:
        print(f"AlphaFold API fetch failed for {accession}: {exc}")

    for ver in AFDB_VERSIONS:
        url = AFDB_FILES.format(acc=accession, ver=ver)
        try:
            result = session.get(url, timeout=120)
        except Exception as exc:
            print(f"AlphaFold {ver} fetch failed for {accession}: {exc}")
            continue
        text = _pdb_from_response(result)
        if text:
            return text
        print(
            f"AlphaFold {ver} miss for {accession}: "
            f"status={getattr(result, 'status_code', '?')} "
            f"bytes={len(getattr(result, 'text', '') or '')}"
        )
    return None


def fetch_pdb(accession: str, output_dir: str, session=None):
    """
    Fetches a PDB file from AlphaFold, given an accession. Places the file in the output_dir.

    Args:
        accession (str): a valid UniprotKB accession.
        output_dir (str): path to the output directory.
            File will be saved as "{output_dir}/{accession}.pdb".
        session (requests.Session, optional): the requests session to use for the request.
    """
    output_path = Path(output_dir) / f"{accession}.pdb"

    if session is None:
        session = session_with_retry()

    if os.path.exists(output_path):
        # Re-fetch if a previous run left an AF error stub.
        try:
            existing = output_path.read_text(errors="replace")
        except OSError:
            existing = ""
        if _looks_like_pdb(existing):
            return
        try:
            output_path.unlink()
        except OSError:
            pass

    text = _download_pdb_text(accession, session)
    if text is None:
        print(f"No AlphaFold PDB for {accession}; skipping file write")
        return
    with open(output_path, "w") as file:
        file.write(text)


def main():
    args = parse_args()
    accessions = args.accession
    output_dir = args.output

    formats = args.format

    for accession in accessions:
        if "fasta" in formats:
            fetch_fasta(accession, output_dir)
        if "pdb" in formats:
            fetch_pdb(accession, output_dir)


if __name__ == "__main__":
    main()
