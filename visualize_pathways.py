"""
Draw pathway diagrams from iterative CSV output (matplotlib + networkx).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


def parse_pathway_segments(pathway_str: str) -> List[str]:
    """Split pathway on arrows (whitespace-tolerant)."""
    if not pathway_str or not str(pathway_str).strip():
        return []
    return [p.strip() for p in re.split(r"\s*->\s*", str(pathway_str)) if p.strip()]


def pathway_to_graph(pathway_str: str) -> Tuple[nx.DiGraph, dict]:
    """Build a directed path graph; node keys n0, n1, ..."""
    segments = parse_pathway_segments(pathway_str)
    G = nx.DiGraph()
    labels: dict[str, str] = {}
    for i, seg in enumerate(segments):
        nid = f"n{i}"
        labels[nid] = seg if len(seg) <= 48 else seg[:45] + "…"
        G.add_node(nid)
        if i > 0:
            G.add_edge(f"n{i-1}", nid)
    return G, labels


def render_pathway_png(
    pathway_str: str,
    probability: float,
    cumulative_delta_g_kj_per_mol: float,
    out_path: Path,
    *,
    title: Optional[str] = None,
) -> None:
    """Render a single pathway as a horizontal directed graph; save PNG."""
    G, labels = pathway_to_graph(pathway_str)
    if G.number_of_nodes() == 0:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "Empty pathway", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    pos = nx.spring_layout(G, seed=42, k=2.0)
    # Prefer left-to-right for a path
    seg_list = parse_pathway_segments(pathway_str)
    if len(seg_list) > 1:
        pos = {f"n{i}": (i, 0) for i in range(len(seg_list))}

    fig, ax = plt.subplots(figsize=(max(8, len(seg_list) * 2.2), 3))
    nx.draw_networkx_nodes(G, pos, node_color="#e8f4fc", node_size=2800, ax=ax)
    nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)
    nx.draw_networkx_edges(G, pos, arrows=True, arrowsize=20, ax=ax)

    edge_labels = {}
    edges = list(G.edges())
    for i, (u, v) in enumerate(edges):
        edge_labels[(u, v)] = f"step {i+1}"

    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=7, ax=ax)

    ttl = title or "Pathway"
    fig.suptitle(
        f"{ttl}\nP={probability:.4g}  ΔG_cum={cumulative_delta_g_kj_per_mol:.2f} kJ/mol",
        fontsize=10,
    )
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def sanitize_filename(s: str) -> str:
    s = re.sub(r"[^\w\-.+]+", "_", s)[:80]
    return s or "pathway"


def visualize_from_csv(
    csv_path: Path,
    out_dir: Path,
    filter_smiles: Optional[str] = None,
    first_only: bool = False,
) -> List[Path]:
    """
    Read iterative CSV; write one PNG per row (or first row per original_smiles if first_only).
    Returns list of written paths.
    """
    df = pd.read_csv(csv_path)
    required = {"original_smiles", "pathway", "probability", "cumulative_delta_g_kj_per_mol"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required}")

    out_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    seen_orig: set[str] = set()
    row_idx = 0

    for _, row in df.iterrows():
        orig = str(row["original_smiles"])
        if filter_smiles and orig != filter_smiles:
            continue
        if first_only and orig in seen_orig:
            continue
        seen_orig.add(orig)

        pathway = str(row["pathway"])
        prob = float(row["probability"])
        dg = float(row["cumulative_delta_g_kj_per_mol"])
        stem = sanitize_filename(orig)
        out_path = out_dir / f"{stem}_tree_{row_idx}.png"
        row_idx += 1
        render_pathway_png(pathway, prob, dg, out_path, title=f"{orig}")
        written.append(out_path)

    return written


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Visualize iterative pathway CSV as network graphs.")
    p.add_argument("csv", type=Path, help="Pathways CSV from process_reactivity_oxygen.py --iterative")
    p.add_argument("-o", "--out-dir", type=Path, default=Path("."), help="Output directory for PNGs")
    p.add_argument("--smiles", default=None, help="Only plot rows with this original_smiles")
    p.add_argument("--first-only", action="store_true", help="Only first pathway per original_smiles")
    args = p.parse_args(argv)

    paths = visualize_from_csv(args.csv, args.out_dir, filter_smiles=args.smiles, first_only=args.first_only)
    for w in paths:
        print(f"Wrote {w}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
