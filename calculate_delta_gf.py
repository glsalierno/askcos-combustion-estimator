import pandas as pd
from thermo import Chemical
from rdkit import Chem

# Read the CSV file
data = pd.read_csv('output_reactivity_oxygen_20250505_150317.csv')

# Collect all unique SMILES (input and products)
unique_smiles = set()

# Add input SMILES
for smiles in data['original_smiles']:
    if pd.notna(smiles):
        unique_smiles.add(smiles)

# Add product SMILES (parse the JSON-like string)
for product_str in data['product_smiles']:
    if pd.isna(product_str) or product_str == '':
        continue
    # Replace single quotes with double quotes for JSON parsing
    product_str = product_str.replace("'", '"')
    try:
        products = eval(product_str)  # Using eval since jsondecode might not be available
        for product in products:
            product_smiles = product.get('outcome', '')
            if product_smiles:
                unique_smiles.add(product_smiles)
    except Exception as e:
        print(f"Error parsing product SMILES: {e}")
        continue

# Calculate Delta Gf for each unique SMILES
delta_gf_data = []
for smiles in unique_smiles:
    try:
        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES: {smiles}")
            continue
        # Calculate Delta Gf using thermo
        chem = Chemical(smiles)
        delta_gf = chem.Gf  # Gibbs free energy of formation in J/mol
        if delta_gf is not None:
            delta_gf_kcal_mol = delta_gf / 4184  # Convert J/mol to kcal/mol
            delta_gf_data.append({'smiles': smiles, 'delta_gf': delta_gf_kcal_mol})
        else:
            print(f"Failed to calculate Delta Gf for {smiles}")
    except Exception as e:
        print(f"Error calculating Delta Gf for {smiles}: {e}")

# Save to CSV
delta_gf_df = pd.DataFrame(delta_gf_data)
delta_gf_df.to_csv('delta_gf_values.csv', index=False)
print("Saved Delta Gf values to delta_gf_values.csv")