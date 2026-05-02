<p align="center">
  <img src="images/banner.svg" alt="ASKCOS Combustion Estimator — Predict fuel properties with the chemical kinetics API" width="100%" />
</p>

# ASKCOS Combustion Estimator

> ⚠️ **Requires a fully operational local ASKCOS instance**  
> This method calls the ASKCOS API. You must install and run ASKCOS (see [Installation](#installation)) before using the estimator.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- Run: pytest --cov=. --cov-config=.coveragerc --cov-report=term-missing once, then replace XX below with the total coverage percent. -->
![Coverage](https://img.shields.io/badge/coverage-59%25-yellowgreen)

Continuous integration runs **pytest** with coverage; optional local hooks are described in [CONTRIBUTING.md](CONTRIBUTING.md).

### Release v2.0.0

GitHub tag and notes: **[v2.0.0](https://github.com/glsalierno/askcos-combustion-estimator/releases/tag/v2.0.0)**.

Version **2.0.0** highlights **iterative** combustion modeling: recursive ASKCOS calls along oxidation pathways, **one row per pathway** with branch probabilities and **`cumulative_delta_g_kj_per_mol`**, tunable **`--max-depth`** / **`--prob-threshold`**, optional **`--combine-products`**, and pathway figures via **`--visualize`** or **`visualize_pathways.py`** (install with `pip install -e ".[viz]"`). The same CLI also supports the original single-shot mode (`askcos-iterative` / `process_reactivity_oxygen.py`). CI (GitHub Actions), `.coveragerc`, validation helpers under **`validation/`**, and **`scripts/monitor_run.sh`** (periodic run/output checks) round out the release. See [Usage](#usage) for flags and examples.

This repository provides a Python implementation that interfaces with the **ASKCOS API** (Automated System for Knowledge-based Continuous Organic Synthesis) to predict oxidation products for organic compounds and to estimate the **Gibbs free energy of combustion**. Combustion is treated as reaction with molecular oxygen (O=O) using ASKCOS forward-synthesis capabilities. The software is intended for **batch processing** of compounds supplied as a CSV list of SMILES strings.

---

## Scope and Capabilities

| Capability | Description |
|------------|-------------|
| **Product prediction** | The ASKCOS API is invoked to predict up to three candidate oxidation products per input SMILES string. |
| **Thermodynamic quantities** | Standard Gibbs free energies of formation (ΔG<sub>f</sub>) for reactants and products are computed using the `thermo` package; the user may then obtain ΔG<sub>combustion</sub> from the relationship given below. |
| **Batch input and output** | Input is a CSV file containing a column named `smiles`. Output is written to CSV with product assignments, rankings, scores, and molecular weights. |
| **Robustness** | API calls employ retry logic; SMILES strings are validated and normalized via RDKit. |

The Gibbs free energy of combustion is computed as:

**ΔG<sub>combustion</sub> = Σ(ΔG<sub>f</sub> products) − ΔG<sub>f</sub>(reactant) − ΔG<sub>f</sub>(O₂)**

where ΔG<sub>f</sub>(O₂) = 0 kcal/mol under standard conditions.

---

## System Requirements

- **Python:** version 3.8 or higher.
- **Dependencies:** Install with `pip install -r requirements.txt`. Required packages include `requests`, `pandas`, `rdkit`, and `thermo`, among others (see `requirements.txt`).
- **ASKCOS (required):** A **local** ASKCOS API server must be running so this package can call `http://<host>/api/forward/call-sync` (default base URL `http://0.0.0.0`). Deployment is typically via **Docker**. Plan for **at least ~8 GB RAM** for the ASKCOS stack on a workstation (more for large batches or concurrent jobs).
- **Upstream guides:** [**ASKCOS on GitHub**](https://github.com/ASKCOS/askcos) — official installation and deployment documentation; this repo also includes [ASKCOS_INSTALLATION.md](ASKCOS_INSTALLATION.md) with project-specific notes.

---

## Installation

**ASKCOS:** Install and start the ASKCOS server **before** running predictions. Follow the [**ASKCOS repository**](https://github.com/ASKCOS/askcos) instructions (Docker-based deployment is typical). See [ASKCOS_INSTALLATION.md](ASKCOS_INSTALLATION.md) for additional context.

**Python environment — from a clone:**

```bash
git clone https://github.com/glsalierno/askcos-combustion-estimator.git
cd askcos-combustion-estimator
pip install -r requirements.txt
```

**CLI package (includes `askcos-iterative` entry point):**

```bash
pip install -e ".[viz]"   # add [viz] for pathway PNGs (matplotlib + networkx)
# or from GitHub after release:
# pip install "askcos-iterative[viz] @ git+https://github.com/glsalierno/askcos-combustion-estimator.git"
askcos-iterative --help
```

**Development / tests:** `pip install -r requirements.txt -r requirements-dev.txt`

---

## Validation

Benchmarking scripts compare iterative pathway predictions to placeholder experimental ΔG°comb values and write a parity plot. See the **[validation/](validation/)** folder and [validation/README.md](validation/README.md) for scope, limitations, and how to run `python validation/benchmark.py`.

---

## Usage

### Original vs iterative prediction

| Aspect | **Original mode** (default) | **Iterative mode** (`--iterative`) |
|--------|-----------------------------|-------------------------------------|
| **Behavior** | **One** ASKCOS request per input SMILES; records **up to three** ranked products per compound from that single response. | **Many** requests per compound: recursively re-invokes ASKCOS on successive products until a **fully oxidized** terminal (e.g. CO₂, H₂O), **`--max-depth`**, or **`--prob-threshold`** stops expansion. |
| **Output** | Wide prediction CSV with `product_rank`, `product_smiles`, scores, raw API payload, etc. | Pathway CSV: **one row per pathway** with `pathway` (arrow-separated), branch **`probability`**, and **`cumulative_delta_g_kj_per_mol`** (sum of approximate step ΔG via `thermo`). |
| **Cost / speed** | ~**N** API calls for **N** compounds. | **Much higher** load and runtime—scales with branching and pathway length; avoid huge batches without tuning depth/threshold. |
| **Typical use** | Screening, quick oxidation suggestions, downstream **`calculate_delta_gf.py`**. | Multi-step combustion routes, pathway-level probabilities, and cumulative model ΔG along oxidation sequences. |

Prefer **original** mode unless you need full pathways and cumulative ΔG columns.

---

The workflow comprises two main scripts.

1. **Product prediction (original mode).** Execute:
   ```bash
   python process_reactivity_oxygen.py
   ```
   By default, the script reads from `DFP_smiles.csv` and writes to a timestamped file of the form `output_reactivity_oxygen_<timestamp>.csv`. The reagent is set to molecular oxygen (O=O). Progress may be saved at regular intervals by adjusting the `save_interval` parameter (default: every 10 compounds).

   **Command-line options (original mode):** `--input` / `-i`, `--output` / `-o`, `--reagent`, `--save-interval`, `--askcos-url` (default `http://0.0.0.0`).

### Iterative mode (`--iterative`)

Recursive ASKCOS expansion: each predicted product is fed back into ASKCOS until a **fully oxidized** terminal pattern is reached (e.g. CO₂, H₂O), **max depth** is hit, or branch probability falls below a **threshold**. Cycles along a pathway are skipped.

```bash
python process_reactivity_oxygen.py --iterative -i DFP_smiles.csv -o pathways_out \
  --max-depth 5 --prob-threshold 0.01 --askcos-url http://127.0.0.1
```

**Optional — co-products in one ASKCOS response:** `--combine-products` treats **all** returned products in a single response as one stoichiometric step (joint probability = product of probabilities; pathway segment shows `P1 + P2 + ...`). Use when the API lists co-products of a single reaction.

**Optional — pathway diagrams:** after an iterative run, add `--visualize` (requires `matplotlib` and `networkx`, e.g. `pip install -e ".[viz]"`):

```bash
python process_reactivity_oxygen.py --iterative -i DFP_smiles.csv -o pathways_out --visualize --visualize-dir ./plots
```

Or run `python visualize_pathways.py pathways_out_<timestamp>.csv -o ./plots`.

Output CSV columns:

| Column | Description |
|--------|-------------|
| `original_smiles` | Input fuel |
| `pathway` | Arrow-separated segments (e.g. `CCO -> CC=O -> O=C=O` or `CCO -> C=C + [OH2]` for a multi-product step) |
| `probability` | Product of ASKCOS branch probabilities along the pathway |
| `cumulative_delta_g_kj_per_mol` | Sum of approximate ΔG (kJ/mol) per balanced step using `thermo` |

### Balancing notes (single- and multi-product)

- **Single product:** `R + a O₂ + b H₂O → P` with carbon conserved (`C` in R equals `C` in P).
- **Multi-product:** `R + a O₂ + b H₂O → P₁ + P₂ + …` with sum of carbons in products equal to `C` in `R`. Coefficients are obtained from H and O atom balance; ΔG<sub>step</sub> = Σ(ΔG<sub>f</sub> products) − ΔG<sub>f</sub>(R) − *a*ΔG<sub>f</sub>(O₂) − *b*ΔG<sub>f</sub>(H₂O).

Core dependencies remain `requirements.txt`. Visualization and tests use extras / `requirements-dev.txt`.

2. **Computation of ΔG<sub>f</sub>.** After obtaining the prediction output, run:
   ```bash
   python calculate_delta_gf.py
   ```
   This script reads the prediction CSV, extracts unique SMILES (reactants and products), computes ΔG<sub>f</sub> in kcal/mol using the `thermo` package, and writes results to `delta_gf_values.csv`. Full ΔG<sub>combustion</sub> may be obtained by post-processing (e.g., reaction balancing and application of the formula above) or by extending the script.

Input CSV files must include a column named `smiles`. The user may modify the default input and output filenames in the scripts or extend the code to accept command-line arguments.

---

## Script Reference

| Script | Function |
|--------|----------|
| `process_reactivity_oxygen.py` | Calls the ASKCOS API with reagent O=O and writes predicted products and associated metadata to a CSV file. Supports **`--iterative`**, **`--combine-products`**, **`--visualize`**. |
| `visualize_pathways.py` | Reads iterative pathway CSV and writes **`{original}_tree_N.png`** (networkx + matplotlib). |
| `calculate_delta_gf.py` | Reads a prediction output CSV, computes ΔG<sub>f</sub> for all unique SMILES, and writes the results to CSV. |

---

## Example Output

Representative excerpt from the prediction output for ethanol (SMILES: `CCO`):

| original_smiles | normalized_smiles | product_rank | product_smiles | probability | score | molecular_weight |
|----------------|-------------------|--------------|----------------|-------------|-------|-------------------|
| CCO            | CCO               | 1            | CC=O           | 0.85        | 0.92  | 44.05             |

Representative ΔG<sub>f</sub> output:

| smiles | delta_gf |
|--------|----------|
| CCO    | -41.2    |
| CC=O   | -30.5    |

---

## Tests

Mocked ASKCOS (no live server) exercise **original** and **iterative** pipelines, **multi-product balancing**, and **visualization**:

```bash
pip install -r requirements.txt -r requirements-dev.txt
pip install -e ".[viz]"
pytest tests/ -v
```

**Coverage** (terminal summary + HTML report under `htmlcov/`):

```bash
pytest --cov=. --cov-config=.coveragerc --cov-report=term-missing --cov-report=html
```

After running once, update the coverage badge at the top of this README: replace `XX` with the reported total percentage.

---

## Known Limitations

- The ASKCOS predictions are based on machine-learning models and are not equivalent to detailed combustion or kinetic simulations. The results should be interpreted as **estimates** only.
- The default configuration assumes **complete oxidation** with O=O. For partial oxidation or other reaction conditions, the reagent and/or script logic must be modified accordingly.
- ΔG<sub>f</sub> values are obtained from the `thermo` package, which uses group-contribution methods. **Verification against experimental or high-level computational data is recommended** for applications requiring high accuracy.
- **Reaction balancing** in iterative mode is approximate (C/H/O with O₂ and H₂O only); exotic species, ions, or nitrogen/sulfur chemistry may need manual checks.
- The parsing of ASKCOS responses assumes a specific output structure. If the local ASKCOS deployment returns a different format, the parsing logic may require adaptation.
- The software is designed for use with a **local** ASKCOS instance; cloud or other remote deployments are not supported in the current version.

> **Extensibility**: Although designed for O₂ combustion, the framework can be adapted for other reagents by modifying the reagent SMILES in `process_reactivity_oxygen.py` and the balancing logic in `calculate_delta_gf.py`. See inline comments.

---

## Contributing

See **[CONTRIBUTING.md](CONTRIBUTING.md)** for environment setup, style, tests, coverage, and pull request guidelines.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Citation

Zenodo archives each GitHub release under its **own DOI**. For **v2.0.0** (iterative pathways release), cite the software record:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19964369.svg)](https://doi.org/10.5281/zenodo.19964369)

Landing page: [zenodo.org/records/19964369](https://zenodo.org/records/19964369).

---

## Acknowledgments

The authors acknowledge the **ASKCOS** project at MIT for the synthesis-prediction API; **RDKit** for SMILES validation and normalization; and the **thermo** package for thermodynamic property calculations.

*Document last updated: May 2026*
