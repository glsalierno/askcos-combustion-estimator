#!/usr/bin/env python3
"""Emit a CSV with `smiles` column for batch testing (default 100 diverse inputs)."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from rdkit import Chem

# Diverse organic SMILES (validated below); expand pool then take first N unique canonical.
_RAW = [
    "CC",
    "CCC",
    "CCCC",
    "CCCCC",
    "CCCCCC",
    "CCCCCCC",
    "CCCCCCCC",
    "CC(C)C",
    "CC(C)(C)C",
    "C=C",
    "CC=CC",
    "C#CC",
    "C=CC=C",
    "c1ccccc1",
    "Cc1ccccc1",
    "CCc1ccccc1",
    "CO",
    "CCO",
    "CCCO",
    "CC(C)O",
    "CC(O)C",
    "OCc1ccccc1",
    "CC(=O)O",
    "CC(=O)OC",
    "CC(=O)C",
    "CCN(CC)CC",
    "CCOC",
    "COC",
    "CCCl",
    "CCBr",
    "CCI",
    "FC(F)(F)C",
    "CC(F)(F)F",
    "O=C(O)c1ccccc1",
    "Nc1ccccc1",
    "O=[N+]([O-])c1ccccc1",
    "c1ccc(O)cc1",
    "c1ccc(N)cc1",
    "c1ccncc1",
    "c1ccc2ccccc2c1",
    "C1CCCCC1",
    "C1CCCOC1",
    "CC1CCOC1",
    "O=C1CCCCC1",
    "CCOC=O",
    "CCC(=O)CC",
    "CC(C)=O",
    "CCCC=O",
    "CC=CCO",
    "OCCCO",
    "OCCCCO",
    "CC(C)(C)CO",
    "CC(O)CCO",
    "OC(CO)CO",
    "NCCO",
    "CC(N)=O",
    "NC(C)=O",
    "CCS",
    "CCSC",
    "CCSOCC",
    "CSC",
    "CCS",
    "CCOC(C)=O",
    "CCOC(=O)CC",
    "COC(=O)c1ccccc1",
    "O=C(O)CC(O)(CC(=O)O)C(=O)O",
    "O=C(O)C(O)C(O)C(O)C(O)CO",
    "CC(O)C(O)C(O)CO",
    "O=C(CO)CO",
    "CC(O)C=O",
    "CCOCOC",
    "COCOC",
    "C1COCCN1",
    "CN1CCCCC1",
    "CN(C)C",
    "CCN",
    "CCCNC",
    "O=C(NC)c1ccccc1",
    "CC(C)NC(C)C",
    "CCCOC",
    "CCCCOC",
    "CCC(O)CC",
    "CC(O)CCC",
    "CC(C)(C)OC",
    "CCOC(C)(C)C",
    "CCCC(O)CC",
    "CCCCC(O)CC",
    "CCOC(OCC)OCC",
    "c1ccoc1",
    "c1cnccn1",
    "c1ccc(O)c(O)c1",
    "COc1ccccc1",
    "COc1ccc(O)cc1",
    "Cc1ccc(O)cc1",
    "Cc1cccc(C)c1",
    "FC(F)(Oc1ccccc1)F",
    "CC(O)c1ccccc1",
    "O=C(C)c1ccccc1",
    "O=C(CC)c1ccccc1",
    "CC(=O)c1ccccc1",
    "CC(O)c1ccccc1",
    "c1ccc(CCO)cc1",
    "c1ccc(CCOCC)cc1",
    "OCc1ccc(O)cc1",
    "CC(O)c1ccc(O)cc1",
    "O=C(O)c1ccc(O)cc1",
    "O=C(O)c1ccc(N)cc1",
    "Nc1ccc(O)cc1",
    "c1ccc(CN)cc1",
    "c1ccc(CCO)nc1",
    "CCc1ccc(O)cc1",
    "CCCc1ccccc1",
    "CCCCc1ccccc1",
    "c1ccc2c(c1)OCO2",
    "c1ccc2scnc2c1",
    "O=S(=O)c1ccccc1",
    "CS(C)=O",
    "CCSO",
    "CCCCSCCCC",
    "CCCCCSC",
    "SCCCO",
    "OCCSCCO",
    "CC(C)OC(C)C",
    "CC(O)C(C)(C)C",
    "CCC(C)(C)C",
    "CCC(C)CC",
    "CCCC(C)C",
    "CC(C)CC(C)C",
]


def canonical_unique(smiles_list: list[str], limit: int) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            continue
        c = Chem.MolToSmiles(mol, isomericSmiles=True)
        if c not in seen:
            seen.add(c)
            out.append(c)
        if len(out) >= limit:
            break
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("-n", type=int, default=100, help="Number of compounds")
    p.add_argument("-o", type=Path, default=Path("batch_smiles_n100.csv"))
    args = p.parse_args()
    rows = canonical_unique(_RAW, args.n)
    if len(rows) < args.n:
        print(
            f"Only {len(rows)} unique valid SMILES available (need {args.n}).",
            file=sys.stderr,
        )
        sys.exit(1)
    args.o.parent.mkdir(parents=True, exist_ok=True)
    with args.o.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["smiles"])
        w.writeheader()
        for s in rows:
            w.writerow({"smiles": s})
    print(f"Wrote {len(rows)} rows to {args.o}")


if __name__ == "__main__":
    main()
