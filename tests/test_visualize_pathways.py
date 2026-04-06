"""Pathway visualization (requires matplotlib, networkx)."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("matplotlib")
pytest.importorskip("networkx")

from visualize_pathways import render_pathway_png, visualize_from_csv


def test_render_pathway_png_writes_file(tmp_path: Path):
    out = tmp_path / "t.png"
    render_pathway_png("CCO -> CC=O -> O=C=O", 0.5, -100.0, out)
    assert out.is_file()
    assert out.stat().st_size > 100


def test_visualize_from_csv_mock(tmp_path: Path):
    csv = tmp_path / "p.csv"
    csv.write_text(
        "original_smiles,pathway,probability,cumulative_delta_g_kj_per_mol\n"
        'CCO,"CCO -> CC=O",0.9,-10.0\n'
    )
    outs = visualize_from_csv(csv, tmp_path)
    assert len(outs) == 1
    assert outs[0].suffix == ".png"
    assert outs[0].is_file()
