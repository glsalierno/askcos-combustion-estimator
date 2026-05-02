# Validation benchmarks

This folder contains scripts to compare **model outputs** from the iterative ASKCOS pathway workflow against **experimental** standard Gibbs energies of combustion, ΔG°comb (kJ mol⁻¹).

## Experimental data

The numbers in `benchmark.py` (`BENCHMARK_MOLECULES`) are **rounded placeholder values** for demonstration. They are **not** intended as a definitive experimental database. For publication or regulatory use, replace them with values taken from primary literature or critically evaluated thermochemical tables, with explicit temperature, phase, and standard-state conventions.

## How to run

1. Start your local ASKCOS server (see the main [README](../README.md) and [ASKCOS_INSTALLATION.md](../ASKCOS_INSTALLATION.md)).
2. From the **repository root**:

   ```bash
   pip install -r requirements.txt
   python validation/benchmark.py --askcos-url http://127.0.0.1
   ```

   Optional flags: `--limit 3` (smoke test), `--max-depth`, `--prob-threshold`, `--combine-products`. Outputs default to `validation/benchmark_results.csv` and `validation/parity_plot.png` (override with `--output-table` / `--parity-plot`).

## Interpreting results

- **`delta_g_comb_pred_kj_per_mol`**: cumulative ΔG along the **selected** iterative oxidation pathway (sum of approximate step Gibbs energies using `thermo`), not necessarily a fully balanced global combustion reaction.
- **`selection`**: whether the row used a pathway that **heuristically** contains both CO₂ (`O=C=O`) and water (`[OH2]`); otherwise the highest-probability pathway was used (model may not reach full combustion).
- **MAE / RMSE**: mean absolute error and root mean square error of (predicted − experimental) in kJ mol⁻¹. Large errors often indicate incomplete pathways, thermo limitations, or mismatch between the pathway cumulative ΔG and tabulated ΔG°comb conventions.

The **parity plot** places experimental ΔG°comb on the *x*-axis and the model cumulative ΔG on the *y*-axis. Agreement along the diagonal indicates similar magnitude and sign under comparable conventions—interpret cautiously because the model quantity is pathway-dependent.
