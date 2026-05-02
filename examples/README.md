# Example inputs and outputs

Files here are **small, illustrative** CSVs that match the shapes written by this repository. Numeric fields in `sample_output_*.csv` are **representative only**; live runs depend on your ASKCOS model and `thermo` data.

| File | Role |
|------|------|
| [`input_smiles.csv`](input_smiles.csv) | Minimal batch input: one column **`smiles`**. |
| [`sample_output_original.csv`](sample_output_original.csv) | Original (single-call) mode: ranked products per row. |
| [`sample_output_iterative.csv`](sample_output_iterative.csv) | Iterative mode: one row per **pathway**. |
| [`sample_delta_gf.csv`](sample_delta_gf.csv) | Shape of **`calculate_delta_gf.py`** output (`smiles`, `delta_gf` in kcal/mol). |

Try the CLI against the sample input (requires a running ASKCOS server):

```bash
python process_reactivity_oxygen.py -i examples/input_smiles.csv -o examples/run_demo
python process_reactivity_oxygen.py --iterative -i examples/input_smiles.csv -o examples/run_demo_pathways --askcos-url http://127.0.0.1
```

`calculate_delta_gf.py` ships with a **hard-coded input filename**; copy a real prediction CSV to that name or edit the script before running.
