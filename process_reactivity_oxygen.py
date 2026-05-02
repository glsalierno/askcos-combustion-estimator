"""
ASKCOS forward oxidation with optional iterative pathway expansion.

Original mode: one ASKCOS call per input SMILES, up to three products per row.
Iterative mode: recursive ASKCOS until fully oxidized, max depth, or probability floor;
                 optional stoichiometric O2/H2O balance per step and cumulative ΔG (thermo).
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import pandas as pd
import requests
from rdkit import Chem
from thermo import Chemical

# --- Iterative-mode configuration (override via CLI where available) ---
MAX_DEPTH = 5
PROB_THRESHOLD = 0.01
FULLY_OXIDIZED_SMILES = [
    "O=C=O",
    "[OH2]",
    "O=S(=O)=O",
]
MAX_PRODUCTS_PER_STEP = 12

# ASKCOS server (same default as original)
ASKCOS_BASE_URL = "http://0.0.0.0"
ASKCOS_ENDPOINT = "/api/forward/call-sync"


def validate_smiles(smiles: str) -> bool:
    """Validate if a SMILES string can be parsed by RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except Exception as e:
        print(f"Error validating SMILES {smiles}: {str(e)}")
        return False


def normalize_smiles(smiles: str) -> Optional[str]:
    """Normalize a SMILES string using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception as e:
        print(f"Error normalizing SMILES {smiles}: {str(e)}")
        return None


def _flatten_askcos_product_list(products: Any) -> List[Any]:
    """
    ASKCOS forward responses often wrap outcomes as [[{...}, {...}, ...]] (one inner list
    per query). Flatten one level so callers iterate over outcome dicts.
    """
    if not isinstance(products, list) or not products:
        return []
    flat: List[Any] = []
    for item in products:
        if isinstance(item, list):
            flat.extend(item)
        else:
            flat.append(item)
    return flat


def extract_product_info(product: Union[Dict[str, Any], str]) -> Dict[str, Any]:
    """Extract product information from ASKCOS outcome dicts or legacy shapes."""
    if isinstance(product, dict):
        smiles = product.get("smiles") or product.get("outcome") or ""
        prob = product.get("probability")
        if prob is None:
            prob = product.get("prob")
        mw = product.get("molecular_weight")
        if mw is None:
            mw = product.get("mol_wt")
        return {
            "product_smiles": (smiles or "").strip(),
            "probability": float(prob or 0.0),
            "score": float(product.get("score", 0.0) or 0.0),
            "molecular_weight": float(mw or 0.0),
        }
    return {
        "product_smiles": str(product),
        "probability": 0.0,
        "score": 0.0,
        "molecular_weight": 0.0,
    }


def _count_CHO(mol: Chem.Mol) -> Tuple[int, int, int]:
    c = h = o = 0
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z == 6:
            c += 1
        elif z == 1:
            h += 1
        elif z == 8:
            o += 1
    return c, h, o


def balance_o2_h2o(
    reactant_smiles: str, product_smiles: str
) -> Optional[Tuple[float, float]]:
    """
    Stoichiometry R + a O2 + b H2O -> P with atom balance (C/H/O only).
    Returns (a, b) or None if carbon count differs or balance is infeasible.
    """
    mr = Chem.MolFromSmiles(reactant_smiles)
    mp = Chem.MolFromSmiles(product_smiles)
    if mr is None or mp is None:
        return None
    cr, hr, or_ = _count_CHO(mr)
    cp, hp, op = _count_CHO(mp)
    if cr != cp:
        return None
    b = (hp - hr) / 2.0
    a = (op - or_ - b) / 2.0
    if a < -1e-6 or b < -1e-6:
        return None
    return (max(a, 0.0), max(b, 0.0))


def balance_multi_product(
    reactant_smiles: str,
    products: List[str],
) -> Optional[Tuple[float, float]]:
    """
    Stoichiometry R + a O2 + b H2O -> P1 + P2 + ... (C/H/O only).
    Carbon must be conserved (sum of C in products = reactant C).
    Returns (a, b) or None if infeasible.
    """
    if not products:
        return None
    if len(products) == 1:
        return balance_o2_h2o(reactant_smiles, products[0])
    mr = Chem.MolFromSmiles(reactant_smiles)
    if mr is None:
        return None
    cr, hr, or_ = _count_CHO(mr)
    cp = hp = op = 0
    for ps in products:
        m = Chem.MolFromSmiles(ps)
        if m is None:
            return None
        c, h, o = _count_CHO(m)
        cp += c
        hp += h
        op += o
    if cr != cp:
        return None
    b = (hp - hr) / 2.0
    a = (op - or_ - b) / 2.0
    if a < -1e-6 or b < -1e-6:
        return None
    return (max(a, 0.0), max(b, 0.0))


def delta_g_multi_product(
    reactant_norm: str,
    product_norms: List[str],
    o2_mol: float,
    h2o_mol: float,
) -> float:
    """ΔG (kJ/mol) for R + a O2 + b H2O -> sum(products)."""
    gr = _gf_j_per_mol(reactant_norm)
    if gr is None:
        return 0.0
    gp_sum = 0.0
    for p in product_norms:
        g = _gf_j_per_mol(p)
        if g is None:
            return 0.0
        gp_sum += g
    gh2o = _gf_water_j_per_mol()
    dg_j = gp_sum - gr - o2_mol * 0.0 - h2o_mol * gh2o
    return dg_j / 1000.0


def _gf_j_per_mol(smiles: str) -> Optional[float]:
    """Standard Gibbs formation energy (J/mol) from thermo; None if unavailable."""
    try:
        c = Chemical(smiles)
        g = getattr(c, "Gf", None)
        if g is not None:
            return float(g)
    except Exception:
        pass
    return None


def _gf_water_j_per_mol() -> float:
    g = _gf_j_per_mol("H2O")
    if g is not None:
        return g
    c = Chemical("water")
    return float(c.Gf)


def _canonical_set_from_patterns(patterns: List[str]) -> Set[str]:
    out: Set[str] = set()
    for p in patterns:
        n = normalize_smiles(p)
        if n:
            out.add(n)
    return out


_FULLY_OXIDIZED_CANONICAL = _canonical_set_from_patterns(FULLY_OXIDIZED_SMILES)


def is_fully_oxidized(smiles: str) -> bool:
    """True if normalized SMILES matches a terminal combustion product pattern."""
    n = normalize_smiles(smiles)
    if n is None:
        return False
    if n in _FULLY_OXIDIZED_CANONICAL:
        return True
    return False


def send_askcs_forward_synthesis(
    smiles: str,
    reagent: str = "O=O",
    max_retries: int = 3,
    retry_delay: int = 5,
    base_url: Optional[str] = None,
) -> Dict[str, Any]:
    """Send request to ASKCOS API with retry logic and log full response."""
    bu = (base_url or ASKCOS_BASE_URL).rstrip("/")
    endpoint = ASKCOS_ENDPOINT
    headers = {
        "accept": "application/json",
        "Content-Type": "application/json",
    }

    if not validate_smiles(smiles):
        return {
            "success": False,
            "error": "Invalid SMILES string",
            "original_smiles": smiles,
            "normalized_smiles": None,
            "products": [],
            "raw_response": None,
        }

    normalized_smiles = normalize_smiles(smiles)
    if normalized_smiles is None:
        return {
            "success": False,
            "error": "Failed to normalize SMILES",
            "original_smiles": smiles,
            "normalized_smiles": None,
            "products": [],
            "raw_response": None,
        }

    payload = {
        "smiles": [normalized_smiles],
        "reagents": reagent,
    }

    for attempt in range(max_retries):
        try:
            response = requests.post(
                f"{bu}{endpoint}",
                headers=headers,
                json=payload,
                timeout=30,
            )

            raw_response_text = response.text
            if response.status_code == 200:
                try:
                    result = response.json()
                    products = []
                    if isinstance(result, dict):
                        if "result" in result:
                            products = result["result"]
                        elif "products" in result:
                            products = result["products"]
                        else:
                            print(
                                f"[WARNING] Unexpected API response structure: {result}"
                            )
                    else:
                        print(f"[WARNING] API response is not a dict: {result}")
                    products = _flatten_askcos_product_list(products)
                    return {
                        "success": True,
                        "error": None,
                        "original_smiles": smiles,
                        "normalized_smiles": normalized_smiles,
                        "products": products,
                        "raw_response": raw_response_text,
                    }
                except json.JSONDecodeError as e:
                    if attempt < max_retries - 1:
                        time.sleep(retry_delay)
                        continue
                    return {
                        "success": False,
                        "error": f"Invalid JSON response: {str(e)}",
                        "original_smiles": smiles,
                        "normalized_smiles": normalized_smiles,
                        "products": [],
                        "raw_response": raw_response_text,
                    }
            else:
                if attempt < max_retries - 1:
                    time.sleep(retry_delay)
                    continue
                return {
                    "success": False,
                    "error": f"API error: {response.status_code}",
                    "original_smiles": smiles,
                    "normalized_smiles": normalized_smiles,
                    "products": [],
                    "raw_response": raw_response_text,
                }

        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
                continue
            return {
                "success": False,
                "error": f"Request failed: {str(e)}",
                "original_smiles": smiles,
                "normalized_smiles": normalized_smiles,
                "products": [],
                "raw_response": None,
            }

    return {
        "success": False,
        "error": "Max retries exceeded",
        "original_smiles": smiles,
        "normalized_smiles": normalized_smiles,
        "products": [],
        "raw_response": None,
    }


def recursive_combustion(
    smiles: str,
    reagent: str,
    depth: int,
    cum_prob: float,
    cum_delta_g_kj: float,
    visited: Set[str],
    *,
    max_depth: int,
    prob_threshold: float,
    max_products: int,
    base_url: str,
    combine_products_in_step: bool = False,
) -> List[Tuple[List[str], float, float]]:
    """
    Returns list of (pathway_segment_list, cumulative_probability, cumulative_delta_g_kj_per_mol).
    Each segment is a SMILES string or \"P1 + P2 + ...\" for a multi-product step.
    """
    n = normalize_smiles(smiles)
    if n is None:
        return []

    if depth >= max_depth or is_fully_oxidized(n):
        return [([n], cum_prob, cum_delta_g_kj)]

    if n in visited:
        return []

    result = send_askcs_forward_synthesis(smiles, reagent, base_url=base_url)
    if not result["success"] or not result["products"]:
        return [([n], cum_prob, cum_delta_g_kj)]

    products_sorted: List[Tuple[float, str]] = []
    for product in result["products"][:max_products]:
        inf = extract_product_info(product)
        ps = (inf.get("product_smiles") or "").strip()
        if not ps or not validate_smiles(ps):
            continue
        pr = float(inf.get("probability", 0.0) or 0.0)
        products_sorted.append((pr, ps))
    products_sorted.sort(key=lambda x: -x[0])

    valid: List[Tuple[float, str, str]] = []
    for pr, ps in products_sorted:
        if pr < prob_threshold:
            continue
        pn = normalize_smiles(ps)
        if pn is None:
            continue
        valid.append((pr, ps, pn))

    if not valid:
        return [([n], cum_prob, cum_delta_g_kj)]

    out: List[Tuple[List[str], float, float]] = []

    if combine_products_in_step:
        norms = [t[2] for t in valid]
        if any(x in visited for x in norms):
            return []
        joint_prob = 1.0
        for pr, _, _ in valid:
            joint_prob *= pr
        label = " + ".join(norms)
        bal = balance_multi_product(n, norms)
        if bal is not None:
            o2_mol, h2o_mol = bal
            step_dg = delta_g_multi_product(n, norms, o2_mol, h2o_mol)
        else:
            step_dg = 0.0
        new_prob = cum_prob * joint_prob
        new_cum = cum_delta_g_kj + step_dg
        organic = [x for x in norms if not is_fully_oxidized(x)]
        vis_next = visited | {n}
        if not organic:
            return [([n, label], new_prob, new_cum)]
        any_child = False
        for pi in organic:
            if pi in vis_next:
                continue
            ps_for_pi = ""
            for _, ps, pn in valid:
                if pn == pi:
                    ps_for_pi = ps
                    break
            if not ps_for_pi:
                continue
            sub = recursive_combustion(
                ps_for_pi,
                reagent,
                depth + 1,
                new_prob,
                new_cum,
                vis_next,
                max_depth=max_depth,
                prob_threshold=prob_threshold,
                max_products=max_products,
                base_url=base_url,
                combine_products_in_step=combine_products_in_step,
            )
            for path_rest, pprob, pdg in sub:
                out.append(([n, label] + path_rest[1:], pprob, pdg))
                any_child = True
        if not any_child:
            return [([n, label], new_prob, new_cum)]
        return out

    any_child = False
    for pr, ps, pn in valid:
        if pn in visited:
            continue
        bal = balance_multi_product(n, [pn])
        if bal is not None:
            o2_mol, h2o_mol = bal
            step_dg = delta_g_multi_product(n, [pn], o2_mol, h2o_mol)
        else:
            step_dg = 0.0
        new_prob = cum_prob * pr
        new_cum = cum_delta_g_kj + step_dg
        sub = recursive_combustion(
            ps,
            reagent,
            depth + 1,
            new_prob,
            new_cum,
            visited | {n},
            max_depth=max_depth,
            prob_threshold=prob_threshold,
            max_products=max_products,
            base_url=base_url,
            combine_products_in_step=combine_products_in_step,
        )
        for path_rest, pprob, pdg in sub:
            out.append(([n] + path_rest, pprob, pdg))
            any_child = True

    if not any_child:
        return [([n], cum_prob, cum_delta_g_kj)]
    return out


def process_smiles_for_reactivity(
    input_csv: str,
    output_csv: str,
    reagent: str = "O=O",
    save_interval: int = 10,
    base_url: str = ASKCOS_BASE_URL,
) -> None:
    """Process SMILES strings and save results to CSV (original single-call behavior)."""
    df = pd.read_csv(input_csv)

    if "smiles" not in df.columns:
        raise ValueError("CSV file must contain a 'smiles' column")

    results = []
    total = len(df)
    processed = 0
    successful = 0
    failed = 0

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    print(f"\nStarting processing at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total SMILES to process: {total}\n")

    for i, smiles in enumerate(df["smiles"], 1):
        print(f"Processing {i}/{total}: {smiles}")

        result = send_askcs_forward_synthesis(smiles, reagent, base_url=base_url)
        processed += 1

        if result["success"]:
            successful += 1
            products = result["products"]
            if products:
                for j, product in enumerate(products[:3], 1):
                    product_info = extract_product_info(product)
                    product_data = {
                        "original_smiles": result["original_smiles"],
                        "normalized_smiles": result["normalized_smiles"],
                        "product_rank": j,
                        **product_info,
                        "error": None,
                        "raw_response": result["raw_response"],
                    }
                    results.append(product_data)
            else:
                print(
                    f"[WARNING] No products found for SMILES: {smiles}. Raw response: {result['raw_response']}"
                )
                results.append(
                    {
                        "original_smiles": result["original_smiles"],
                        "normalized_smiles": result["normalized_smiles"],
                        "product_rank": 1,
                        "product_smiles": "",
                        "probability": 0.0,
                        "score": 0.0,
                        "molecular_weight": 0.0,
                        "error": "No products found",
                        "raw_response": result["raw_response"],
                    }
                )
        else:
            failed += 1
            print(
                f"[ERROR] Failed to process SMILES: {smiles}. Error: {result['error']}. Raw response: {result['raw_response']}"
            )
            results.append(
                {
                    "original_smiles": result["original_smiles"],
                    "normalized_smiles": result["normalized_smiles"],
                    "product_rank": 1,
                    "product_smiles": "",
                    "probability": 0.0,
                    "score": 0.0,
                    "molecular_weight": 0.0,
                    "error": result["error"],
                    "raw_response": result["raw_response"],
                }
            )

        if i % save_interval == 0 or i == total:
            df_results = pd.DataFrame(results)
            df_results.to_csv(f"{output_csv}_{timestamp}.csv", index=False)
            print("\nProgress update:")
            print(f"Processed: {processed}/{total} ({processed / total * 100:.1f}%)")
            print(f"Successful: {successful}")
            print(f"Failed: {failed}")
            print(f"Results saved to: {output_csv}_{timestamp}.csv\n")

    df_results = pd.DataFrame(results)
    df_results.to_csv(f"{output_csv}_{timestamp}.csv", index=False)

    print(f"\nProcessing complete at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("Final statistics:")
    print(f"Total processed: {processed}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Success rate: {successful / total * 100:.1f}%")
    print(f"Results saved to: {output_csv}_{timestamp}.csv")


def process_smiles_iterative(
    input_csv: str,
    output_csv: str,
    reagent: str = "O=O",
    save_interval: int = 1,
    *,
    max_depth: int = MAX_DEPTH,
    prob_threshold: float = PROB_THRESHOLD,
    max_products: int = MAX_PRODUCTS_PER_STEP,
    base_url: str = ASKCOS_BASE_URL,
    combine_products_in_step: bool = False,
) -> str:
    """
    For each input SMILES, expand pathways recursively; write one row per pathway.
    Columns: original_smiles, pathway, probability, cumulative_delta_g_kj_per_mol
    Returns path to the written CSV (timestamped filename).
    """
    df = pd.read_csv(input_csv)
    if "smiles" not in df.columns:
        raise ValueError("CSV file must contain a 'smiles' column")

    rows: List[Dict[str, Any]] = []
    total = len(df)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    print(
        f"\n[Iterative mode] max_depth={max_depth}, prob_threshold={prob_threshold}, "
        f"combine_products_in_step={combine_products_in_step}"
    )
    print(f"Starting at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}, N={total}\n")

    for i, smiles in enumerate(df["smiles"], 1):
        print(f"{i}/{total}: {smiles}")
        orig = str(smiles)
        paths = recursive_combustion(
            orig,
            reagent,
            0,
            1.0,
            0.0,
            set(),
            max_depth=max_depth,
            prob_threshold=prob_threshold,
            max_products=max_products,
            base_url=base_url,
            combine_products_in_step=combine_products_in_step,
        )
        if not paths:
            rows.append(
                {
                    "original_smiles": orig,
                    "pathway": orig,
                    "probability": 1.0,
                    "cumulative_delta_g_kj_per_mol": 0.0,
                }
            )
        else:
            for path, pprob, pdg in paths:
                rows.append(
                    {
                        "original_smiles": orig,
                        "pathway": " -> ".join(path),
                        "probability": pprob,
                        "cumulative_delta_g_kj_per_mol": round(pdg, 4),
                    }
                )

        if i % save_interval == 0 or i == total:
            pd.DataFrame(rows).to_csv(f"{output_csv}_{timestamp}.csv", index=False)
            print(f"  checkpoint -> {output_csv}_{timestamp}.csv ({len(rows)} rows)")

    out_path = f"{output_csv}_{timestamp}.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"\nDone. Wrote {len(rows)} rows to {out_path}")
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ASKCOS combustion estimator (original or iterative)."
    )
    parser.add_argument(
        "--iterative",
        action="store_true",
        help="Recursive ASKCOS expansion with pathway CSV (see README).",
    )
    parser.add_argument(
        "--input",
        "-i",
        default="DFP_smiles.csv",
        help="Input CSV with column 'smiles'.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="output_reactivity_oxygen",
        help="Output prefix (timestamp appended).",
    )
    parser.add_argument(
        "--reagent", default="O=O", help="ASKCOS reagent SMILES (default O=O)."
    )
    parser.add_argument(
        "--save-interval",
        type=int,
        default=10,
        help="Checkpoint every N rows (original mode).",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=MAX_DEPTH,
        help="Max oxidation steps (iterative).",
    )
    parser.add_argument(
        "--prob-threshold",
        type=float,
        default=PROB_THRESHOLD,
        help="Min branch probability (iterative).",
    )
    parser.add_argument(
        "--max-products",
        type=int,
        default=MAX_PRODUCTS_PER_STEP,
        help="Max ASKCOS products per step (iterative).",
    )
    parser.add_argument(
        "--askcos-url", default=ASKCOS_BASE_URL, help="ASKCOS server base URL."
    )
    parser.add_argument(
        "--combine-products",
        action="store_true",
        help="Iterative: treat all ASKCOS products in one response as co-products of a single step (balancing + joint probability).",
    )
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="After iterative run, render pathway PNGs via visualize_pathways.py (requires matplotlib, networkx).",
    )
    parser.add_argument(
        "--visualize-smiles",
        default=None,
        help="Only plot pathways for this original_smiles (optional; default: all rows).",
    )
    parser.add_argument(
        "--visualize-dir",
        type=str,
        default=".",
        help="Directory for pathway PNGs when using --visualize.",
    )
    args = parser.parse_args()

    askcos_base = args.askcos_url.rstrip("/")

    if args.iterative:
        csv_path = process_smiles_iterative(
            args.input,
            args.output,
            reagent=args.reagent,
            save_interval=max(1, args.save_interval),
            max_depth=args.max_depth,
            prob_threshold=args.prob_threshold,
            max_products=args.max_products,
            base_url=askcos_base,
            combine_products_in_step=args.combine_products,
        )
        if args.visualize:
            from visualize_pathways import visualize_from_csv

            visualize_from_csv(
                Path(csv_path),
                Path(args.visualize_dir),
                filter_smiles=args.visualize_smiles,
                first_only=False,
            )
    else:
        process_smiles_for_reactivity(
            args.input,
            args.output,
            reagent=args.reagent,
            save_interval=args.save_interval,
            base_url=askcos_base,
        )


if __name__ == "__main__":
    main()
