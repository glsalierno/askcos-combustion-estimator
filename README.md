# ASKCOS-Based Combustion Product and Free Energy Estimator

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains a Python implementation that leverages the **[ASKCOS API](https://askcos.mit.edu/)** (Automated System for Knowledge-based Continuous Organic Synthesis) to estimate combustion products for organic compounds and calculate their Gibbs free energy of combustion. Combustion is modeled as oxidation with molecular oxygen (O=O) using ASKCOS's forward synthesis prediction capabilities.

The tool is designed for **batch processing**: Provide a CSV file containing a list of SMILES (Simplified Molecular Input Line Entry System) strings, and it will predict oxidation products and compute thermodynamic properties.

This is particularly useful for chemists, researchers, or engineers interested in predicting oxidation products (as proxies for combustion outcomes) and estimating thermodynamic properties like the free energy of combustion (ΔG_combustion), which is computed as:

**ΔG_combustion = Σ(ΔG_f products) - ΔG_f reactant - ΔG_f O₂**

*(Note: ΔG_f for O₂ is 0 kcal/mol at standard conditions.)*

## Key Features
- **Product Prediction**: Uses ASKCOS to predict up to 3 likely oxidation products per input SMILES.
- **Thermodynamic Calculations**: Computes Gibbs free energy of formation (ΔG_f) for reactants and products using the `thermo` package, then derives ΔG_combustion (extend as needed for full computation).
- **Batch Input**: Processes a CSV file with a 'smiles' column for multiple compounds.
- **Output**: Generates CSV files with detailed results, including normalized SMILES, product rankings, probabilities, scores, molecular weights, and raw API responses for debugging.
- **Retry Logic**: Built-in retries for API calls to handle transient errors.
- **Validation**: SMILES strings are validated and normalized using RDKit.

## Requirements
- **Python 3.8 or higher**
- **Required packages** (install via `pip install -r requirements.txt`):
  - `requests`
  - `json`
  - `pandas`
  - `typing`
  - `rdkit`
  - `thermo`
  - `datetime`
- A local installation of the ASKCOS API server (running on `http://0.0.0.0`). **Important**: This tool assumes ASKCOS is installed and running locally. For installation instructions, please refer to the official **[ASKCOS documentation](https://askcos.mit.edu/docs/)** or their **[GitHub repository](https://github.com/ASKCOS/askcos-core)**. ASKCOS typically runs via Docker for Linux environments—follow their setup guide for Linux-specific steps.

*(Add additional links here if needed, e.g., for advanced ASKCOS configuration: [ASKCOS Advanced Docs](https://example-link.com).)*

## Installation
1. Clone this repository:

   git clone https://github.com/glsalierno/askcos-combustion-estimator.git
cd askcos-combustion-estimator

2. Install dependencies:
   
   pip install -r requirements.txt

3. Set up ASKCOS locally:
- Follow the **[ASKCOS installation guide](https://askcos.mit.edu/docs/installation.html)** for Linux. This usually involves Docker Compose to spin up the API server.
- Ensure the ASKCOS server is running and accessible at `http://0.0.0.0/api/forward/call-sync`.

*(Add links here for Docker setup if relevant, e.g., [Docker Documentation](https://docs.docker.com/).)*

## Usage
The implementation consists of two main scripts:
- `process_reactivity_oxygen.py`: Predicts combustion products using ASKCOS.
- `calculate_delta_gf.py`: Calculates ΔG_f for unique SMILES from the prediction output.

### Processing a CSV File
Prepare a CSV file (e.g., `input_smiles.csv`) with a column named `smiles`, e.g.:

smiles
CCO
CC(=O)O
C1CCCCC1

Run the prediction script:

python process_reactivity_oxygen.py

- By default, it uses `input_csv = "DFP_smiles.csv"` and `output_csv = "output_reactivity_oxygen"`. Edit these in the script or add command-line arguments for flexibility (future enhancement).
- This processes all SMILES in the CSV, predicts products with reagent "O=O", and saves to `output_reactivity_oxygen_TIMESTAMP.csv`.
- Use the `save_interval` parameter in the script to save progress every N compounds (default: 10).

### Calculating Free Energy
After running the product prediction, compute ΔG_f:

python calculate_delta_gf.py

- By default, it reads from a hardcoded output CSV (e.g., `output_reactivity_oxygen_20250505_150317.csv`) and saves ΔG_f values to `delta_gf_values.csv`.
- This extracts unique SMILES (reactants + products), calculates ΔG_f in kcal/mol using `thermo`, and saves to CSV.
- To compute full ΔG_combustion, you can extend the script or post-process the CSV (e.g., via pandas) by balancing the reaction and applying the formula above.

*(Add links here for `thermo` package details, e.g., [thermo PyPI Page](https://pypi.org/project/thermo/) or [thermo Documentation](https://thermo.readthedocs.io/). For RDKit: [RDKit Docs](https://www.rdkit.org/docs/).)*

## Example Output
For an input CSV with SMILES "CCO" (ethanol), the prediction output CSV might include:

| original_smiles | normalized_smiles | product_rank | product_smiles | probability | score | molecular_weight | error | raw_response |
|-----------------|-------------------|--------------|----------------|-------------|-------|------------------|-------|--------------|
| CCO            | CCO              | 1            | CC=O           | 0.85       | 0.92 | 44.05           | None | {"result": [...]} |

ΔG_f CSV:

| smiles | delta_gf |
|--------|----------|
| CCO   | -41.2   |
| CC=O  | -30.5   |

## Limitations
- ASKCOS predictions are ML-based and may not always accurately model real combustion (which involves complex radical mechanisms). Use as an estimation tool only.
- The tool assumes complete oxidation; for partial oxidation or specific conditions, adjust the reagent in the script.
- ΔG_f calculations rely on the `thermo` package, which uses group contribution methods—verify with experimental data for critical applications.
- No automatic balancing of reactions in the current scripts; users may need to manually balance for accurate ΔG_combustion.
- Product parsing assumes a specific ASKCOS response structure (e.g., "smiles" or "outcome" keys)—test with your local setup and adjust if needed.
- Local ASKCOS setup required; cloud alternatives may exist but are not supported here.

## Contributing
Contributions are welcome! Feel free to open issues or pull requests for improvements, such as:
- Adding command-line arguments for input/output files and reagents.
- Implementing reaction balancing for automated ΔG_combustion calculation.
- Supporting more API endpoints or error handling.
- Enhancing with visualizations (e.g., via RDKit for molecule rendering).

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments
- **[ASKCOS](https://askcos.mit.edu/)** team at MIT for the powerful synthesis prediction API.
- **[RDKit](https://www.rdkit.org/)** for SMILES handling.
- **[thermo](https://pypi.org/project/thermo/)** for thermodynamic property calculations.

*(Add more acknowledgments or links here for related tools, e.g., [PubChem for SMILES Reference](https://pubchem.ncbi.nlm.nih.gov/).)*

*Last updated: February 2026*
