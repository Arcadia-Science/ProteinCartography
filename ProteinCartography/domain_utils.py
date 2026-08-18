"""Domain identity, TED chopping, FASTA/PDB cropping, and the multi-domain gate.

The protein pipeline does not import this module. Domain-map code uses it to
turn a parent protein into ``{parent}__d01``-style domain instances.
"""

from __future__ import annotations
import json
import re
from collections.abc import Iterable, Sequence
from pathlib import Path

DOMAIN_ID_SEP = "__d"
DOMAIN_ID_PATTERN = re.compile(r"^(?P<parent>.+)" + DOMAIN_ID_SEP + r"(?P<index>\d+)$")
CHOPPING_SEGMENT = re.compile(r"^(\d+)-(\d+)$")
# Canonical UniProtKB accessions (6 or 10 characters), optional isoform suffix.
UNIPROT_ACCESSION = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]"
    r"|A0A[A-Z0-9]{7})(?:-\d+)?$"
)

DOMAIN_FEATURES_COLUMNS = [
    "protid",
    "parent_protid",
    "domain_index",
    "chopping",
    "start",
    "end",
    "nres_domain",
    "assignment_source",
    "ted_id",
    "cath_label",
]

DEFAULT_MIN_DOMAIN_LENGTH = 30
TED_SUMMARY_URL = "https://ted.cathdb.info/api/v1/uniprot/summary/{acc}"
TED_PAGE_SIZE = 50
TED_MAX_429_RETRIES = 5


class DomainChoppingError(ValueError):
    """Raised when a TED/user chopping string cannot be parsed."""


def parse_chopping(chopping: str) -> list[tuple[int, int]]:
    """Parse a TED chopping into 1-based inclusive ``(start, end)`` segments.

    Continuous: ``\"22-141\"``. Discontinuous: ``\"22-50_80-120\"``.
    """
    if not chopping or not str(chopping).strip():
        raise DomainChoppingError("empty chopping")
    segments: list[tuple[int, int]] = []
    for part in str(chopping).strip().split("_"):
        match = CHOPPING_SEGMENT.match(part)
        if not match:
            raise DomainChoppingError(f"invalid chopping segment: {part!r} in {chopping!r}")
        start, end = int(match.group(1)), int(match.group(2))
        if start < 1 or end < start:
            raise DomainChoppingError(f"invalid residue range {start}-{end} in {chopping!r}")
        segments.append((start, end))
    return segments


def chopping_nres(segments: Sequence[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in segments)


def residue_in_segments(resseq: int, segments: Sequence[tuple[int, int]]) -> bool:
    return any(start <= resseq <= end for start, end in segments)


def make_domain_id(parent_protid: str, domain_index: int) -> str:
    if domain_index < 1:
        raise ValueError("domain_index is 1-based")
    return f"{parent_protid}{DOMAIN_ID_SEP}{domain_index:02d}"


def parse_domain_id(domain_id: str) -> tuple[str, int]:
    match = DOMAIN_ID_PATTERN.match(domain_id)
    if not match:
        raise ValueError(f"not a domain id: {domain_id!r}")
    return match.group("parent"), int(match.group("index"))


def looks_like_uniprot_accession(protid: str) -> bool:
    return bool(UNIPROT_ACCESSION.match(protid))


def canonical_uniprot_accession(protid: str) -> str:
    """Strip an isoform suffix so TED can be queried with the canonical acc."""
    if "-" in protid and looks_like_uniprot_accession(protid):
        return protid.rsplit("-", 1)[0]
    return protid


def is_multidomain(n_kept_domains: int) -> bool:
    return n_kept_domains >= 2


def keep_domain(nres: int, min_domain_length: int = DEFAULT_MIN_DOMAIN_LENGTH) -> bool:
    return nres >= min_domain_length


def domain_row(
    parent_protid: str,
    domain_index: int,
    chopping: str,
    assignment_source: str,
    ted_id: str = "",
    cath_label: str = "",
    nres_domain: int | None = None,
) -> dict:
    segments = parse_chopping(chopping)
    nres = chopping_nres(segments) if nres_domain is None else int(nres_domain)
    return {
        "protid": make_domain_id(parent_protid, domain_index),
        "parent_protid": parent_protid,
        "domain_index": domain_index,
        "chopping": chopping,
        "start": segments[0][0],
        "end": segments[-1][1],
        "nres_domain": nres,
        "assignment_source": assignment_source,
        "ted_id": ted_id or "",
        "cath_label": cath_label or "",
    }


def filter_domain_rows(
    rows: Iterable[dict], min_domain_length: int = DEFAULT_MIN_DOMAIN_LENGTH
) -> list[dict]:
    kept = [row for row in rows if keep_domain(int(row["nres_domain"]), min_domain_length)]
    # Re-number domain_index / protid after filtering so ids stay contiguous.
    renumbered = []
    by_parent: dict[str, list[dict]] = {}
    for row in kept:
        by_parent.setdefault(row["parent_protid"], []).append(row)
    for parent, parent_rows in by_parent.items():
        for i, row in enumerate(parent_rows, start=1):
            updated = dict(row)
            updated["domain_index"] = i
            updated["protid"] = make_domain_id(parent, i)
            renumbered.append(updated)
    return renumbered


def parse_fasta(text: str) -> tuple[str, str]:
    lines = [line.rstrip("\n") for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError("FASTA must start with a '>' header")
    header = lines[0][1:].split()[0]
    sequence = "".join(lines[1:]).replace(" ", "").replace("\r", "")
    return header, sequence


def slice_fasta_sequence(sequence: str, chopping: str) -> str:
    """Slice a protein sequence with 1-based inclusive UniProt numbering."""
    parts = []
    for start, end in parse_chopping(chopping):
        if end > len(sequence):
            raise ValueError(f"chopping {start}-{end} exceeds sequence length {len(sequence)}")
        parts.append(sequence[start - 1 : end])
    return "".join(parts)


def write_fasta(path: Path, seq_id: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    wrapped = "\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60))
    path.write_text(f">{seq_id}\n{wrapped}\n")


def _pdb_resseq(line: str) -> int | None:
    if len(line) < 26:
        return None
    if not line.startswith("ATOM") and not line.startswith("HETATM"):
        return None
    raw = line[22:26]
    try:
        return int(raw)
    except ValueError:
        return None


def _is_c_alpha(line: str) -> bool:
    if len(line) < 16:
        return False
    return line.startswith("ATOM") and line[12:16].strip() == "CA"


def crop_pdb_text(pdb_text: str, chopping: str) -> tuple[str, int]:
    """Keep ATOM/HETATM records whose residue sequence number is in ``chopping``.

    Returns ``(cropped_pdb, n_c_alpha)``. Residue numbers are PDB ``resSeq``
    fields (author/UniProt numbering), not atom serials or list indices.
    """
    segments = parse_chopping(chopping)
    kept_lines = []
    ca_resseqs = set()
    for line in pdb_text.splitlines():
        resseq = _pdb_resseq(line)
        if resseq is None:
            continue
        if residue_in_segments(resseq, segments):
            kept_lines.append(line.rstrip("\n"))
            if _is_c_alpha(line):
                ca_resseqs.add(resseq)
    kept_lines.append("END")
    return "\n".join(kept_lines) + "\n", len(ca_resseqs)


def crop_pdb_file(pdb_path: Path, chopping: str, output_path: Path) -> int:
    text = pdb_path.read_text(errors="replace")
    cropped, n_ca = crop_pdb_text(text, chopping)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(cropped)
    return n_ca


def crop_fasta_file(fasta_path: Path, chopping: str, output_path: Path, seq_id: str) -> str:
    _, sequence = parse_fasta(fasta_path.read_text())
    sliced = slice_fasta_sequence(sequence, chopping)
    write_fasta(output_path, seq_id, sliced)
    return sliced


def ted_cache_path(cache_dir: Path, accession: str) -> Path:
    return Path(cache_dir) / f"{accession}.ted.json"


def load_ted_cache(cache_dir: Path, accession: str) -> dict | None:
    path = ted_cache_path(cache_dir, accession)
    if not path.is_file():
        return None
    return json.loads(path.read_text())


def save_ted_cache(cache_dir: Path, accession: str, payload: dict) -> None:
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    ted_cache_path(cache_dir, accession).write_text(json.dumps(payload))


def rows_from_ted_payload(parent_protid: str, payload: dict) -> list[dict]:
    data = payload.get("data") or []
    rows = []
    for i, item in enumerate(data, start=1):
        chopping = item.get("chopping") or ""
        if not chopping:
            continue
        rows.append(
            domain_row(
                parent_protid=parent_protid,
                domain_index=i,
                chopping=chopping,
                assignment_source="ted",
                ted_id=item.get("ted_id") or "",
                cath_label=item.get("cath_label") or "",
                nres_domain=item.get("nres_domain"),
            )
        )
    return rows


def rows_from_user_tsv_records(records: Iterable[dict]) -> list[dict]:
    """Build domain rows from user TSV dicts with ``parent_protid`` and chopping/start-end."""
    grouped: dict[str, list[dict]] = {}
    for rec in records:
        parent = str(rec["parent_protid"]).strip()
        if rec.get("chopping"):
            chopping = str(rec["chopping"]).strip()
        else:
            chopping = f"{int(rec['start'])}-{int(rec['end'])}"
        grouped.setdefault(parent, []).append({"chopping": chopping, **rec})

    rows = []
    for parent, items in grouped.items():
        for i, item in enumerate(items, start=1):
            index = int(item["domain_index"]) if item.get("domain_index") else i
            rows.append(
                domain_row(
                    parent_protid=parent,
                    domain_index=index,
                    chopping=item["chopping"],
                    assignment_source="user",
                    ted_id=str(item.get("ted_id") or ""),
                    cath_label=str(item.get("cath_label") or ""),
                )
            )
    return rows
