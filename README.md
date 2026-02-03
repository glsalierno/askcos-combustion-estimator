# askcos-combustion-estimator
ASKCOS API (Automated System for Knowledge-based Continuous Organic Synthesis) implemetation to predict combustion products for organic compounds and estimate a Gibbs free energy of combustion. The combustion is modeles as reaction with molecular oxygen using ASKCOS capabilities. The tool is designed for batch processing using a CSV list of SMILES.

This is particularly useful for chemists, researchers, or engineers interested in predicting oxidation products (as proxies for combustion outcomes) and estimating thermodynamic properties like the free energy of combustion (ΔG_combustion), which is computed as:
ΔG_combustion = Σ(ΔG_f products) - ΔG_f reactant - ΔG_f O₂
(Note: ΔG_f for O₂ is 0 kcal/mol at standard conditions.)
