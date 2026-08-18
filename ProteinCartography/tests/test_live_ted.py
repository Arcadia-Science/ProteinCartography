"""Live TED schema check. Not run in default CI: pytest -m live_ted."""

from __future__ import annotations

import pytest
import requests

TED_URL = "https://ted.cathdb.info/api/v1/uniprot/summary/P00698"


@pytest.mark.live_ted
def test_p00698_chopping_is_still_22_141():
    response = requests.get(TED_URL, timeout=60)
    response.raise_for_status()
    payload = response.json()
    choppings = [item["chopping"] for item in payload.get("data") or []]
    assert "22-141" in choppings
