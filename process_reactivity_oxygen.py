import requests
import json
import pandas as pd
from typing import Dict, Any, List, Optional, Union
from rdkit import Chem
import time
from datetime import datetime

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

def extract_product_info(product: Union[Dict[str, Any], str]) -> Dict[str, Any]:
    """Extract product information handling both dictionary and string formats."""
    if isinstance(product, dict):
        return {
            "product_smiles": product.get("smiles", ""),
            "probability": product.get("probability", 0.0),
            "score": product.get("score", 0.0),
            "molecular_weight": product.get("molecular_weight", 0.0)
        }
    else:
        # If product is just a SMILES string
        return {
            "product_smiles": str(product),
            "probability": 0.0,
            "score": 0.0,
            "molecular_weight": 0.0
        }

def send_askcs_forward_synthesis(smiles: str, reagent: str = "O=O", max_retries: int = 3, retry_delay: int = 5) -> Dict[str, Any]:
    """Send request to ASKCOS API with retry logic and log full response."""
    base_url = "http://0.0.0.0"  # Using local ASKCOS server
    endpoint = "/api/forward/call-sync"
    
    headers = {
        "accept": "application/json",
        "Content-Type": "application/json"
    }
    
    # Validate and normalize SMILES first
    if not validate_smiles(smiles):
        return {
            "success": False,
            "error": "Invalid SMILES string",
            "original_smiles": smiles,
            "normalized_smiles": None,
            "products": [],
            "raw_response": None
        }
    
    normalized_smiles = normalize_smiles(smiles)
    if normalized_smiles is None:
        return {
            "success": False,
            "error": "Failed to normalize SMILES",
            "original_smiles": smiles,
            "normalized_smiles": None,
            "products": [],
            "raw_response": None
        }
    
    payload = {
        "smiles": [normalized_smiles],
        "reagents": reagent
    }
    
    # Retry logic
    for attempt in range(max_retries):
        try:
            response = requests.post(
                f"{base_url}{endpoint}",
                headers=headers,
                json=payload,
                timeout=30
            )
            
            raw_response_text = response.text
            if response.status_code == 200:
                try:
                    result = response.json()
                    # Try to extract products from different possible structures
                    products = []
                    if isinstance(result, dict):
                        if "result" in result:
                            products = result["result"]
                        elif "products" in result:
                            products = result["products"]
                        else:
                            print(f"[WARNING] Unexpected API response structure (no 'result' or 'products' key): {result}")
                    else:
                        print(f"[WARNING] API response is not a dict: {result}")
                    return {
                        "success": True,
                        "error": None,
                        "original_smiles": smiles,
                        "normalized_smiles": normalized_smiles,
                        "products": products,
                        "raw_response": raw_response_text
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
                        "raw_response": raw_response_text
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
                    "raw_response": raw_response_text
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
                "raw_response": None
            }
    
    return {
        "success": False,
        "error": "Max retries exceeded",
        "original_smiles": smiles,
        "normalized_smiles": normalized_smiles,
        "products": [],
        "raw_response": None
    }

def process_smiles_for_reactivity(input_csv: str, output_csv: str, reagent: str = "O=O", save_interval: int = 10) -> None:
    """Process SMILES strings and save results to CSV, including raw API response."""
    # Read CSV file
    df = pd.read_csv(input_csv)
    
    if 'smiles' not in df.columns:
        raise ValueError("CSV file must contain a 'smiles' column")
    
    # Initialize results list and counters
    results = []
    total = len(df)
    processed = 0
    successful = 0
    failed = 0
    
    # Create timestamp for this run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    print(f"\nStarting processing at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total SMILES to process: {total}\n")
    
    # Process each SMILES string
    for i, smiles in enumerate(df['smiles'], 1):
        print(f"Processing {i}/{total}: {smiles}")
        
        result = send_askcs_forward_synthesis(smiles, reagent)
        processed += 1
        
        if result["success"]:
            successful += 1
            # Process products
            products = result["products"]
            if products:
                for j, product in enumerate(products[:3], 1):  # Take up to 3 products
                    product_info = extract_product_info(product)
                    product_data = {
                        "original_smiles": result["original_smiles"],
                        "normalized_smiles": result["normalized_smiles"],
                        "product_rank": j,
                        **product_info,
                        "error": None,
                        "raw_response": result["raw_response"]
                    }
                    results.append(product_data)
            else:
                print(f"[WARNING] No products found for SMILES: {smiles}. Raw response: {result['raw_response']}")
                # No products but successful API call
                results.append({
                    "original_smiles": result["original_smiles"],
                    "normalized_smiles": result["normalized_smiles"],
                    "product_rank": 1,
                    "product_smiles": "",
                    "probability": 0.0,
                    "score": 0.0,
                    "molecular_weight": 0.0,
                    "error": "No products found",
                    "raw_response": result["raw_response"]
                })
        else:
            failed += 1
            print(f"[ERROR] Failed to process SMILES: {smiles}. Error: {result['error']}. Raw response: {result['raw_response']}")
            # Record the error
            results.append({
                "original_smiles": result["original_smiles"],
                "normalized_smiles": result["normalized_smiles"],
                "product_rank": 1,
                "product_smiles": "",
                "probability": 0.0,
                "score": 0.0,
                "molecular_weight": 0.0,
                "error": result["error"],
                "raw_response": result["raw_response"]
            })
        
        # Save results periodically
        if i % save_interval == 0 or i == total:
            df_results = pd.DataFrame(results)
            df_results.to_csv(f"{output_csv}_{timestamp}.csv", index=False)
            print(f"\nProgress update:")
            print(f"Processed: {processed}/{total} ({processed/total*100:.1f}%)")
            print(f"Successful: {successful}")
            print(f"Failed: {failed}")
            print(f"Results saved to: {output_csv}_{timestamp}.csv\n")
    
    # Final save and summary
    df_results = pd.DataFrame(results)
    df_results.to_csv(f"{output_csv}_{timestamp}.csv", index=False)
    
    print(f"\nProcessing complete at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Final statistics:")
    print(f"Total processed: {processed}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Success rate: {successful/total*100:.1f}%")
    print(f"Results saved to: {output_csv}_{timestamp}.csv")

if __name__ == "__main__":
    input_csv = "DFP_smiles.csv"
    output_csv = "output_reactivity_oxygen"
    process_smiles_for_reactivity(input_csv, output_csv)
