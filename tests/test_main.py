"""
Smoke tests for primary modules and user-facing flows (mocked HTTP where needed).
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import process_reactivity_oxygen as pr
import visualize_pathways


class _FakeResp:
    def __init__(self, products: list):
        self.status_code = 200
        self.text = json.dumps({"result": products})

    def json(self):
        return json.loads(self.text)


def _mock_post_cco_chain(*args, **kwargs):
    payload = kwargs.get("json") or {}
    smis = (payload.get("smiles") or [""])[0]
    n = pr.normalize_smiles(smis)
    if n is None:
        return _FakeResp([])
    if n == pr.normalize_smiles("CCO"):
        return _FakeResp([{"smiles": "CC=O", "probability": 0.9, "score": 0.5}])
    if n == "CC=O":
        return _FakeResp([{"smiles": "O=C=O", "probability": 0.95, "score": 0.5}])
    return _FakeResp([])


def test_core_modules_import():
    """Ensure main library modules import (CLI entrypoint lives in process_reactivity_oxygen)."""
    assert pr.ASKCOS_ENDPOINT == "/api/forward/call-sync"
    assert callable(visualize_pathways.parse_pathway_segments)


def test_single_molecule_combustion_delta_g_numeric():
    """Methane (SMILES C): full oxidation balance yields a finite ΔG (kJ/mol) via thermo."""
    r = pr.normalize_smiles("C")
    assert r is not None
    prods = ["O=C=O", "[OH2]", "[OH2]"]
    bal = pr.balance_multi_product(r, prods)
    assert bal is not None
    a, b = bal
    norms: list[str] = []
    for p in prods:
        n = pr.normalize_smiles(p)
        assert n is not None
        norms.append(n)
    dg = pr.delta_g_multi_product(r, norms, a, b)
    assert isinstance(dg, float)
    assert dg == dg  # not NaN


def test_invalid_smiles_graceful_askcos_wrapper():
    """Invalid SMILES must not crash the ASKCOS request wrapper; no HTTP on validation failure."""
    out = pr.send_askcs_forward_synthesis(
        "%%%not_a_smiles%%%", base_url="http://127.0.0.1"
    )
    assert out["success"] is False
    assert out["products"] == []
    assert "Invalid" in (out.get("error") or "") or out.get("normalized_smiles") is None


def test_iterative_oxidation_mocked_http(monkeypatch, tmp_path: Path):
    """Iterative pathway expansion runs end-to-end with requests.post patched."""
    monkeypatch.setattr(pr.requests, "post", _mock_post_cco_chain)

    inp = tmp_path / "in.csv"
    out = tmp_path / "pathways"
    inp.write_text("smiles\nCCO\n", encoding="utf-8")

    pr.process_smiles_iterative(
        str(inp),
        str(out),
        save_interval=1,
        max_depth=8,
        prob_threshold=0.01,
        base_url="http://test",
    )

    outs = list(tmp_path.glob("pathways_*.csv"))
    assert len(outs) == 1
    df = pd.read_csv(outs[0])
    assert "cumulative_delta_g_kj_per_mol" in df.columns
    assert len(df) >= 1
