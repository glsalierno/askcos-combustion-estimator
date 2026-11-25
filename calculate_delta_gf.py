"""Read a prediction CSV and write its Gibbs free energy of formation lookup."""

import argparse
import logging
from ast import literal_eval
from pathlib import Path
from typing import List

import pandas as pd
from rdkit import Chem
from thermo import Chemical

logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        type=str.upper,
        default="INFO",
    )
    parser.add_argument("--log-file", type=Path, help="Also append logs to this file.")
    args = parser.parse_args()
    handlers: List[logging.Handler] = [logging.StreamHandler()]
    if args.log_file:
        try:
            handlers.append(logging.FileHandler(args.log_file, encoding="utf-8"))
        except OSError as error:
            parser.error(f"Could not open --log-file: {error}")
    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        handlers=handlers,
    )

    # Read the CSV file
    data = pd.read_csv("output_reactivity_oxygen_20250505_150317.csv")

    # Collect all unique SMILES (input and products)
    unique_smiles = set()

    # Add input SMILES
    for smiles in data["original_smiles"]:
        if pd.notna(smiles):
            unique_smiles.add(smiles)

    # Add product SMILES (parse the JSON-like string)
    for product_str in data["product_smiles"]:
        if pd.isna(product_str) or product_str == "":
            continue
        # Replace single quotes with double quotes for JSON parsing
        product_str = product_str.replace("'", '"')
        try:
            products = literal_eval(product_str)
            for product in products:
                product_smiles = product.get("outcome", "")
                if product_smiles:
                    unique_smiles.add(product_smiles)
        except Exception as e:
            logger.warning("Could not parse product SMILES (%s)", type(e).__name__)
            continue

    # Calculate Delta Gf for each unique SMILES
    delta_gf_data = []
    for smiles in unique_smiles:
        try:
            # Validate SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning("Invalid SMILES: %r", smiles)
                continue
            # Calculate Delta Gf using thermo
            chem = Chemical(smiles)
            delta_gf = chem.Gf  # Gibbs free energy of formation in J/mol
            if delta_gf is not None:
                delta_gf_kcal_mol = delta_gf / 4184  # Convert J/mol to kcal/mol
                delta_gf_data.append({"smiles": smiles, "delta_gf": delta_gf_kcal_mol})
            else:
                logger.warning("Delta Gf unavailable for SMILES %r", smiles)
        except Exception as e:
            logger.warning(
                "Delta Gf lookup failed for SMILES %r (%s)", smiles, type(e).__name__
            )

    # Save to CSV
    delta_gf_df = pd.DataFrame(delta_gf_data)
    delta_gf_df.to_csv("delta_gf_values.csv", index=False)
    logger.info("Saved Delta Gf values to delta_gf_values.csv")


if __name__ == "__main__":
    main()
