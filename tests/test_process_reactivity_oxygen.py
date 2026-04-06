"""
Tests for original (single-call) and iterative ASKCOS modes with mocked HTTP (no live server).
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

import process_reactivity_oxygen as pr


class FakeResponse:
    def __init__(self, products: list):
        self.status_code = 200
        body = {"result": products}
        self.text = json.dumps(body)

    def json(self):
        return json.loads(self.text)


def _mock_askcos_chain_cco_to_co2(*args, **kwargs):
    """CCO -> CC=O -> O=C=O (terminal); payload uses normalized SMILES in json['smiles'][0]."""
    payload = kwargs.get("json") or {}
    smis = (payload.get("smiles") or [""])[0]
    n = pr.normalize_smiles(smis)
    if n is None:
        return FakeResponse([])
    if n == pr.normalize_smiles("CCO"):
        return FakeResponse([{"smiles": "CC=O", "probability": 0.9, "score": 0.5}])
    if n == "CC=O":
        return FakeResponse([{"smiles": "O=C=O", "probability": 0.95, "score": 0.5}])
    return FakeResponse([])


def test_is_fully_oxidized_co2_water():
    assert pr.is_fully_oxidized("O=C=O") is True
    assert pr.is_fully_oxidized("[OH2]") is True
    assert pr.is_fully_oxidized("CCO") is False


def test_balance_o2_h2o_ethanol_to_acetaldehyde():
    bal = pr.balance_o2_h2o("CCO", "CC=O")
    assert bal is not None
    a, b = bal
    assert a >= 0 and b >= 0


def test_original_mode_writes_products(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(pr.requests, "post", _mock_askcos_chain_cco_to_co2)

    inp = tmp_path / "in.csv"
    out = tmp_path / "out"
    inp.write_text("smiles\nCCO\n")

    pr.process_smiles_for_reactivity(str(inp), str(out), save_interval=1, base_url="http://test")

    outs = list(tmp_path.glob("out_*.csv"))
    assert len(outs) == 1
    df = pd.read_csv(outs[0])
    assert len(df) >= 1
    assert "product_smiles" in df.columns
    row = df.iloc[0]
    assert row["normalized_smiles"] == pr.normalize_smiles("CCO")
    assert float(row["probability"]) == pytest.approx(0.9)


def test_iterative_mode_pathways_to_co2(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(pr.requests, "post", _mock_askcos_chain_cco_to_co2)

    inp = tmp_path / "in.csv"
    out = tmp_path / "pathways"
    inp.write_text("smiles\nCCO\n")

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
    assert list(df.columns) == [
        "original_smiles",
        "pathway",
        "probability",
        "cumulative_delta_g_kj_per_mol",
    ]
    assert len(df) >= 1
    last = df.iloc[0]
    assert last["original_smiles"] == "CCO"
    assert "O=C=O" in last["pathway"]
    assert float(last["probability"]) == pytest.approx(0.9 * 0.95)


def test_recursive_combustion_returns_empty_for_invalid(monkeypatch):
    monkeypatch.setattr(pr.requests, "post", _mock_askcos_chain_cco_to_co2)
    paths = pr.recursive_combustion(
        "not_a_smiles",
        "O=O",
        0,
        1.0,
        0.0,
        set(),
        max_depth=5,
        prob_threshold=0.01,
        max_products=12,
        base_url="http://test",
        combine_products_in_step=False,
    )
    assert paths == []


def test_balance_multi_dehydration_ethanol():
    """C2H5OH -> C2H4 + H2O; no O2/H2O on reactant side (a=b=0)."""
    bal = pr.balance_multi_product("CCO", ["C=C", "[OH2]"])
    assert bal is not None
    a, b = bal
    assert a == pytest.approx(0.0)
    assert b == pytest.approx(0.0)


def test_balance_multi_benzene_to_co2_water():
    """C6H6 + 7.5 O2 -> 6 CO2 + 3 H2O."""
    prods = ["O=C=O"] * 6 + ["[OH2]"] * 3
    bal = pr.balance_multi_product("c1ccccc1", prods)
    assert bal is not None
    a, b = bal
    assert a == pytest.approx(7.5)
    assert b == pytest.approx(0.0)
