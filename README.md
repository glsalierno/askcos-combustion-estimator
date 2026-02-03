# askcos-combustion-estimator
ASKCOS API (Automated System for Knowledge-based Continuous Organic Synthesis) implemetation to predict combustion products for organic compounds and estimate a Gibbs free energy of combustion. The combustion is modeles as reaction with molecular oxygen using ASKCOS capabilities. The tool is designed for batch processing using a CSV list of SMILES.

This is particularly useful for chemists, researchers, or engineers interested in predicting oxidation products (as proxies for combustion outcomes) and estimating thermodynamic properties like the free energy of combustion (ΔG_combustion), which is computed as:

ΔG_combustion = Σ(ΔG_f products) - ΔG_f reactant - ΔG_f O₂

(Note: ΔG_f for O₂ is 0 kcal/mol at standard conditions.)

# Key Features

+ Product Prediction: Uses ASKCOS to predict up to 3 likely oxidation products per input SMILES.
+ Thermodynamic Calculations: Computes Gibbs free energy of formation (ΔG_f) for reactants and products using the thermo package, then derives ΔG_combustion (extend as needed for full computation).
+ Batch Input: Processes a CSV file with a 'smiles' column for multiple compounds.
+ Output: Generates CSV files with detailed results, including normalized SMILES, product rankings, probabilities, scores, molecular weights, and raw API responses for debugging.
+ Retry Logic: Built-in retries for API calls to handle transient errors.
+ Validation: SMILES strings are validated and normalized using RDKit.
