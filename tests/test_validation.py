"""Tests for validation/benchmark.py (mocked ASKCOS / fast paths)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from validation import benchmark as bm


def _fake_process_smiles_iterative(inp: str, out_prefix: str, **kwargs) -> str:
    """Write a minimal iterative CSV so benchmarking logic runs without ASKCOS."""
    out_path = Path(out_prefix).parent / "mock_iterative_out.csv"
    pd.DataFrame(
        {
            "original_smiles": ["C"],
            "pathway": ["C -> O=C=O -> [OH2]"],
            "probability": [1.0],
            "cumulative_delta_g_kj_per_mol": [-123.45],
        }
    ).to_csv(out_path, index=False)
    return str(out_path)


def test_run_benchmark_smoke(monkeypatch):
    monkeypatch.setattr(
        bm.pr, "process_smiles_iterative", _fake_process_smiles_iterative
    )

    df = bm.run_benchmark(
        [
            {
                "smiles": "C",
                "name": "methane",
                "delta_g_comb_exp_kj_per_mol": -818.0,
            }
        ],
        askcos_url="http://unused.test",
        max_depth=3,
        prob_threshold=0.01,
        combine_products=False,
    )
    assert len(df) == 1
    assert df.iloc[0]["smiles"] == "C"
    pred = df.iloc[0]["delta_g_comb_pred_kj_per_mol"]
    assert pred is not None and float(pred) == pytest.approx(-123.45)


def test_plot_parity_writes_png(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "delta_g_comb_exp_kj_per_mol": [-100.0, -200.0],
            "delta_g_comb_pred_kj_per_mol": [-105.0, -195.0],
        }
    )
    png = tmp_path / "parity.png"
    bm.plot_parity(df, png)
    assert png.is_file()
    assert png.stat().st_size > 100
