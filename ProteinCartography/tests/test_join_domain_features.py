"""Join UniProt columns onto domain rows without replacing domain ids."""

from __future__ import annotations
from pathlib import Path

import join_domain_features


def test_join_keeps_domain_id_and_copies_parent_metadata(tmp_path: Path):
    domain = tmp_path / "domains.tsv"
    domain.write_text(
        "protid\tparent_protid\tchopping\tcath_label\tnres_domain\n"
        "P12345__d01\tP12345\t1-80\t3.40.50.300\t80\n"
    )
    uniprot = tmp_path / "uniprot.tsv"
    uniprot.write_text(
        "protid\tProtein names\tOrganism\tLength\nP12345\tExample protein\tHomo sapiens\t375\n"
    )
    found = tmp_path / "found.tsv"
    found.write_text("accession\tfound_by\nP12345\tP99999__d01\n")
    out = tmp_path / "joined.tsv"
    df = join_domain_features.join_domain_uniprot(
        str(domain), str(uniprot), str(out), found_by_file=str(found)
    )
    assert list(df["protid"]) == ["P12345__d01"]
    assert df.loc[0, "Protein names"] == "Example protein"
    assert df.loc[0, "found_by"] == "P99999__d01"
    assert int(df.loc[0, "Length"]) == 80
