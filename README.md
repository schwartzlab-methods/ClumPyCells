# ClumPyCells

ClumPyCells is a Python toolbox for spatial mark cross-correlation analysis in spatial omics data. It supports cell-size correction, long-running progress feedback, a Streamlit GUI, terminal/HPC execution, and downstream AUC, group comparison, decision-tree, survival, and permutation workflows.

Preprint: <https://www.biorxiv.org/content/10.64898/2026.03.26.714529v1>

## Installation

From this repository:

```bash
python -m venv CLPC
source CLPC/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

After publication to PyPI:

```bash
python -m pip install clumpycells
```

## GUI

Run the packaged Streamlit app:

```bash
clumpycells-gui
```

Or from a source checkout:

```bash
streamlit run streamlit_app.py
```

The GUI lets users map arbitrary CSV column names, preview images, run KMM/markcorr, generate HPC commands, and perform downstream analysis. The Cohort analysis tab uses the full CSV for the selected mapping, while Data preview remains intentionally sampled for responsiveness. The Image preview tab draws cell areas in the same coordinate scale as x/y and can filter visible cell types.

## Python API

```python
from ClumPyCells import MarkcorrResult, runSpatial

runSpatial(
    csv_path="cells.csv",
    savefolder="results/",
    xrange=[0, 1056],
    yrange=[0, 642],
    x_col="centroid_x",
    y_col="centroid_y",
    image_col="sample_id",
    area_col="cell_area",
    mark_columns=["cell_type", "CD3", "MPO"],
    sizeCorrection=True,
    pp_criterion=lambda df: df["Area"] > 100,
    max_workers=4,
)

result = MarkcorrResult(
    groups={"A": ["1", "2"], "B": ["3", "4"]},
    resultFolder="results/",
    axisName={"cell_type_on": "on", "cell_type_off": "off"},
)
auc, heatmaps = result.getAUC()
```

## Terminal and HPC

The installable package provides a terminal command for offline and cluster runs:

```bash
clumpycells inspect-csv --csv cells.csv

clumpycells run-markcorr \
    --csv cells.csv \
    --out results_hpc \
    --x-col centroid_x \
    --y-col centroid_y \
    --image-col sample_id \
    --area-col cell_area \
    --mark cell_type \
    --xrange 0 1056 \
    --yrange 0 642 \
    --max-workers 4

clumpycells run-permutation \
    --run-config results_hpc \
    --out results_hpc \
    --seed 42 \
    --n-perm 100 \
    --max-workers 4
```

`run-markcorr` writes `clumpycells_run_config.json` under the result folder; `run-permutation --run-config` reuses those KMM settings. For scheduler arrays, pass image chunks with `--image-id`. The GUI also generates generic Slurm templates for markcorr and permutation jobs.

## Documentation

Detailed usage, API examples, GUI notes, HPC templates, demo data descriptions, and testing commands are in [Documentation/UserGuide.md](Documentation/UserGuide.md).

The older function reference is in [Documentation/Documentation.md](Documentation/Documentation.md).

## Repository Layout

```text
ClumPyCells/                # importable Python package
streamlit_app.py            # Streamlit GUI entry point
clumpycells_cli.py          # terminal/HPC command implementation
Documentation/              # detailed docs
Data/Demo/                  # small mock data for UI and API testing
tests/                      # pytest coverage for core workflows
```

## Development Checks

```bash
python -m py_compile streamlit_app.py clumpycells_cli.py
pytest tests/ -v
```

Build metadata is defined in [pyproject.toml](pyproject.toml). For a local release check, install the development extras and run:

```bash
python -m pip install -e ".[dev]"
python -m build
python -m twine check dist/*
```
