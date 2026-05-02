#!/usr/bin/env python3
"""
Compare iterative-pathway cumulative ΔG (ASKCOS + thermo steps) to experimental
standard Gibbs energies of combustion ΔG°comb (kJ/mol).

Requires a running ASKCOS instance (see repository README).
Experimental values below are **illustrative placeholders** — replace with
literature data for publication-quality benchmarking.
"""

from __future__ import annotations

import argparse
import math
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import process_reactivity_oxygen as pr  # noqa: E402


# ΔG°comb (kJ/mol), standard-state conventions as in typical chemistry tables;
# negative values = spontaneous combustion. Replace with peer-reviewed sources.
BENCHMARK_MOLECULES: List[Dict[str, Any]] = [
    {"smiles": "C", "name": "methane", "delta_g_comb_exp_kj_per_mol": -818.0},
    {"smiles": "CC", "name": "ethane", "delta_g_comb_exp_kj_per_mol": -1468.0},
    {"smiles": "CCC", "name": "propane", "delta_g_comb_exp_kj_per_mol": -2044.0},
    {"smiles": "CCCC", "name": "butane", "delta_g_comb_exp_kj_per_mol": -2658.0},
    {"smiles": "CO", "name": "methanol", "delta_g_comb_exp_kj_per_mol": -638.0},
    {"smiles": "CCO", "name": "ethanol", "delta_g_comb_exp_kj_per_mol": -1317.0},
    {"smiles": "CCCO", "name": "1-propanol", "delta_g_comb_exp_kj_per_mol": -1818.0},
    {"smiles": "C=C", "name": "ethylene", "delta_g_comb_exp_kj_per_mol": -1332.0},
    {"smiles": "CC(=O)C", "name": "acetone", "delta_g_comb_exp_kj_per_mol": -1690.0},
    {"smiles": "CC(=O)O", "name": "acetic acid", "delta_g_comb_exp_kj_per_mol": -874.0},
    {"smiles": "c1ccccc1", "name": "benzene", "delta_g_comb_exp_kj_per_mol": -3202.0},
    {"smiles": "Cc1ccccc1", "name": "toluene", "delta_g_comb_exp_kj_per_mol": -3770.0},
    {
        "smiles": "CCOC(C)=O",
        "name": "ethyl acetate",
        "delta_g_comb_exp_kj_per_mol": -2236.0,
    },
    {
        "smiles": "O=CC=O",
        "name": "glyoxal (demo)",
        "delta_g_comb_exp_kj_per_mol": -900.0,
    },
]


def _pathway_looks_like_full_combustion(pathway: str) -> bool:
    """Heuristic: pathway lists CO₂ and water ([OH2]) as products."""
    if not pathway or not isinstance(pathway, str):
        return False
    return "O=C=O" in pathway and "[OH2]" in pathway


def predicted_delta_g_from_iterative_csv(
    df: pd.DataFrame,
) -> Tuple[Optional[float], str]:
    """
    Pick cumulative ΔG (kJ/mol) from iterative output: prefer highest-probability
    pathway that appears to reach CO₂ + H₂O; otherwise best-effort highest probability.
    """
    if df.empty:
        return None, "empty_output"
    need = {"pathway", "probability", "cumulative_delta_g_kj_per_mol"}
    if not need.issubset(df.columns):
        return None, "missing_columns"

    complete = df[df["pathway"].map(_pathway_looks_like_full_combustion)]
    pick = complete if not complete.empty else df
    idx = pick["probability"].idxmax()
    row = pick.loc[idx]
    val = float(row["cumulative_delta_g_kj_per_mol"])
    note = "complete_pathway" if not complete.empty else "fallback_max_probability"
    return val, note


def run_benchmark(
    molecules: List[Dict[str, Any]],
    *,
    askcos_url: str,
    max_depth: int,
    prob_threshold: float,
    combine_products: bool,
) -> pd.DataFrame:
    rows_out: List[Dict[str, Any]] = []
    with tempfile.TemporaryDirectory(dir=REPO_ROOT) as tmp:
        tmp_path = Path(tmp)
        for entry in molecules:
            smiles = entry["smiles"]
            name = entry["name"]
            exp = float(entry["delta_g_comb_exp_kj_per_mol"])

            inp = tmp_path / f"in_{abs(hash(smiles))}.csv"
            out_prefix = tmp_path / f"out_{abs(hash(smiles))}"
            inp.write_text("smiles\n" + smiles + "\n", encoding="utf-8")

            csv_written = pr.process_smiles_iterative(
                str(inp),
                str(out_prefix),
                save_interval=1,
                max_depth=max_depth,
                prob_threshold=prob_threshold,
                base_url=askcos_url.rstrip("/"),
                combine_products_in_step=combine_products,
            )
            df = pd.read_csv(csv_written)
            pred, sel_note = predicted_delta_g_from_iterative_csv(df)
            err = (
                (pred - exp)
                if pred is not None and not math.isnan(pred)
                else float("nan")
            )

            rows_out.append(
                {
                    "smiles": smiles,
                    "name": name,
                    "delta_g_comb_exp_kj_per_mol": exp,
                    "delta_g_comb_pred_kj_per_mol": pred,
                    "error_kj_per_mol": err,
                    "selection": sel_note,
                }
            )
    return pd.DataFrame(rows_out)


def _mae_rmse(errors: List[float]) -> Tuple[float, float]:
    clean = [e for e in errors if not math.isnan(e)]
    if not clean:
        return float("nan"), float("nan")
    mae = sum(abs(x) for x in clean) / len(clean)
    rmse = math.sqrt(sum(x * x for x in clean) / len(clean))
    return mae, rmse


def plot_parity(df: pd.DataFrame, out_png: Path) -> None:
    exp = df["delta_g_comb_exp_kj_per_mol"].astype(float)
    pred = df["delta_g_comb_pred_kj_per_mol"].astype(float)
    mask = pred.notna()
    exp_p = exp[mask].values
    pred_p = pred[mask].values

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(exp_p, pred_p, alpha=0.75, edgecolors="k", linewidths=0.5)
    lims = [min(exp_p.min(), pred_p.min()), max(exp_p.max(), pred_p.max())]
    pad = 0.05 * (lims[1] - lims[0] + 1e-9)
    lo, hi = lims[0] - pad, lims[1] + pad
    ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="y = x")
    ax.set_xlabel(r"Experimental $\Delta G^\circ_{\mathrm{comb}}$ (kJ mol$^{-1}$)")
    ax.set_ylabel(r"Predicted cumulative $\Delta G$ (kJ mol$^{-1}$)")
    ax.set_title("Parity: iterative pathway model vs. experimental placeholders")
    ax.legend(loc="upper left")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark combustion ΔG vs. experimental placeholders."
    )
    parser.add_argument(
        "--askcos-url",
        default="http://127.0.0.1",
        help="ASKCOS base URL (default http://127.0.0.1).",
    )
    parser.add_argument("--max-depth", type=int, default=pr.MAX_DEPTH)
    parser.add_argument("--prob-threshold", type=float, default=pr.PROB_THRESHOLD)
    parser.add_argument(
        "--combine-products",
        action="store_true",
        help="Pass combine_products_in_step=True to iterative expansion.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only first N molecules (for quick smoke tests).",
    )
    parser.add_argument(
        "--output-table",
        type=Path,
        default=REPO_ROOT / "validation" / "benchmark_results.csv",
        help="Where to write the numeric results table.",
    )
    parser.add_argument(
        "--parity-plot",
        type=Path,
        default=REPO_ROOT / "validation" / "parity_plot.png",
        help="Where to save the parity plot PNG.",
    )
    args = parser.parse_args()

    mols = BENCHMARK_MOLECULES[: args.limit] if args.limit else BENCHMARK_MOLECULES

    print(
        f"Running benchmark on {len(mols)} molecules (ASKCOS at {args.askcos_url})...\n"
    )
    df = run_benchmark(
        mols,
        askcos_url=args.askcos_url,
        max_depth=args.max_depth,
        prob_threshold=args.prob_threshold,
        combine_products=args.combine_products,
    )

    pd.set_option("display.max_rows", None)
    pd.set_option("display.width", 120)
    print(df.to_string(index=False))

    errs = df["error_kj_per_mol"].tolist()
    mae, rmse = _mae_rmse([float(e) for e in errs])
    print(f"\nMAE (kJ/mol):  {mae:.2f}")
    print(f"RMSE (kJ/mol): {rmse:.2f}")

    args.output_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_table, index=False)
    print(f"\nWrote table: {args.output_table}")

    plot_parity(df, args.parity_plot)
    print(f"Wrote parity plot: {args.parity_plot}")


if __name__ == "__main__":
    main()
