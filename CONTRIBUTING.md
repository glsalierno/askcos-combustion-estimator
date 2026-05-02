# Contributing

Thank you for helping improve this project.

## Development setup

```bash
git clone https://github.com/glsalierno/askcos-combustion-estimator.git
cd askcos-combustion-estimator
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt -r requirements-dev.txt
pip install -e ".[viz]"
```

A **local ASKCOS instance** is required to exercise live API paths (see the main README). Tests use mocks and do not need ASKCOS running.

## Code style

- Follow [PEP 8](https://peps.python.org/pep-0008/).
- Formatting and linting use **Ruff** (`ruff format`, `ruff check`), matching CI and `.pre-commit-config.yaml`.

## Pre-commit (optional)

Install hooks once after cloning:

```bash
pip install -r requirements.txt -r requirements-dev.txt
pre-commit install
```

Hooks run the checks defined in `.pre-commit-config.yaml` (including **`pytest tests/`** with **`--cov`**). Run all hooks manually:

```bash
pre-commit run --all-files
```

## Testing and coverage

```bash
pytest tests/ -v
pytest --cov=. --cov-config=.coveragerc --cov-report=term-missing --cov-report=html
```

Open `htmlcov/index.html` for the HTML report. CI exposes **`coverage.xml`** as a workflow artifact for downstream tooling.

## Acknowledgments

This project relies on the **ASKCOS** software suite for reaction prediction and free energy estimation.
If you use this tool in your research, please also cite ASKCOS appropriately:

> **ASKCOS: An open source chemical synthesis planning software**
> https://askcos.mit.edu/
> Coley, C. W., et al. (2017). "A graph-convolutional neural network model for the prediction of chemical reactivity." *Chemical Science*, 8(4), 3190–3203.

Respect the ASKCOS license terms (MIT) and note that this estimator does **not** alter ASKCOS source code; it only calls its API.

## Issues and pull requests

- **Issues**: Describe expected vs. actual behavior, Python version, ASKCOS URL/version if relevant, and a minimal CSV or SMILES example.
- **Pull requests**: Keep changes focused; reference the issue number when applicable; ensure `pytest tests/` passes and, if you use pre-commit, `pre-commit run --all-files` is clean. Note limitations in the README or validation docs when behavior changes.

## Validation set

When adding molecules to `validation/benchmark.py`, include a short **citation or source comment** for experimental ΔG°comb (or mark clearly as a synthetic placeholder). Prefer SMILES that RDKit normalizes consistently with your ASKCOS workflow.
