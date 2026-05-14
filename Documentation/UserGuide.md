# ClumPyCells User Guide

This guide covers the Streamlit GUI, terminal/HPC workflows, and Python imports exposed by the installable package.

## Installation

From a local checkout:

```bash
python -m venv CLPC
source CLPC/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

After publication to PyPI, users can install the same package with:

```bash
python -m pip install clumpycells
```

For development and release checks:

```bash
python -m pip install -e ".[dev]"
python -m build
python -m twine check dist/*
```

## GUI

Start the GUI from an installed package:

```bash
clumpycells-gui
```

Or from a source checkout:

```bash
streamlit run streamlit_app.py
```

The GUI has two task areas:

| Area                | Purpose                                                                                                                                                    |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Run KMM / markcorr  | Load a cell CSV, map columns, inspect the cohort and individual images, choose size-correction settings, run markcorr, and generate terminal/HPC commands. |
| Downstream analysis | Open a result folder for KMM curves, AUC heatmaps, group comparisons, permutation runs, decision trees, survival analysis, and artifact browsing.          |

The Data preview tab only displays a small sample so that large CSVs remain responsive. The Cohort analysis tab reads the full CSV for the currently mapped image, area, and mark columns, so counts and distributions reflect the entire file. The Image preview tab uses equal x/y scaling and draws cell areas as coordinate-scale circles, with cell-type color and visibility filters.

## Input CSV

Each row is one cell. The GUI, CLI, and `runSpatial` can map arbitrary input names into the internal roles below.

| Role         | Required                 | Notes                                                                       |
| ------------ | ------------------------ | --------------------------------------------------------------------------- |
| x coordinate | Yes                      | Numeric cell centroid x coordinate.                                         |
| y coordinate | Yes                      | Numeric cell centroid y coordinate.                                         |
| image ID     | No                       | If omitted, the full CSV is treated as one image.                           |
| area         | Only for size correction | Used to derive cell diameters and large-cell occluder rules.                |
| mark columns | Yes                      | Categorical cell types or numeric features used for mark cross-correlation. |

## Python API

Top-level imports are available for common workflows:

```python
from ClumPyCells import MarkcorrResult, prepare_cell_table, runSpatial

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
auc, heatmaps = result.getAUC(norm="min_mid_max", plot=True)
```

For lower-level markcorr primitives:

```python
from ClumPyCells.Markcorr import markcorr, pointPattern, window
```

Decision-tree and survival helpers remain available from their modules:

```python
from ClumPyCells.Analysis.decisionTree import decision_tree_from_feature_table
from ClumPyCells.Analysis.survivalAnalysis import run_user_survival_analysis
```

## Terminal and HPC

The package installs a terminal command named `clumpycells`.

Inspect a CSV:

```bash
clumpycells inspect-csv --csv cells.csv
```

Run markcorr:

```bash
clumpycells run-markcorr \
    --csv cells.csv \
    --out results_hpc \
    --x-col centroid_x \
    --y-col centroid_y \
    --image-col sample_id \
    --area-col cell_area \
    --mark cell_type \
    --mark CD3 \
    --xrange 0 1056 \
    --yrange 0 642 \
    --size-correction \
    --pp-area-threshold 100 \
    --max-workers "${SLURM_CPUS_PER_TASK:-1}" \
    --threads-per-worker 1
```

Run within-image mark-label permutations:

```bash
clumpycells run-permutation \
    --run-config results_hpc \
    --out results_hpc \
    --seed 42 \
    --n-perm 100 \
    --max-workers "${SLURM_CPUS_PER_TASK:-1}" \
    --threads-per-worker 1
```

`run-markcorr` writes `clumpycells_run_config.json` under the result folder. Passing the result folder or that JSON file with `--run-config` makes permutation runs reuse the original CSV mapping, mark columns, x/y window, and size-correction settings. You can still override individual settings on the command line when needed.

For scheduler arrays, split work by image ID:

```bash
clumpycells run-markcorr \
    --csv /path/to/cells.csv \
    --out results_hpc \
    --x-col centroid_x \
    --y-col centroid_y \
    --image-col sample_id \
    --mark cell_type \
    --xrange 0 1056 \
    --yrange 0 642 \
    --max-workers 4 \
    --image-id "$IMAGE_IDS"
```

The GUI generates generic Slurm-array templates for both markcorr and permutation runs. Replace project paths, module loads, queues, memory, and time limits with values that match the local cluster.

## Outputs

`runSpatial` and `run-markcorr` write one folder per image:

```text
results/
├── image_1/
│   ├── r.csv
│   └── iso.csv
├── clumpycells_run_config.json
├── image_2/
└── ...
```

Permutation runs write under:

```text
results/perm_runs/perm_<seed>/image_<id>/
```

## Demo Data

Small mock files live in `Data/Demo/`:

| File                         | Use                                           |
| ---------------------------- | --------------------------------------------- |
| `demo_cells.csv`             | Cell-level CSV for GUI and CLI markcorr runs. |
| `demo_clinical.csv`          | Mock clinical table for survival UI testing.  |
| `demo_image_to_clinical.csv` | Image-to-clinical ID mapping.                 |
| `demo_feature_table.csv`     | Standalone decision-tree feature table.       |
| `demo_streamlit_config.json` | Importable GUI configuration.                 |

## Tests

```bash
pytest tests/ -v
```

The tests cover the markcorr engine, the CSV-to-result pipeline, size-correction invariants, and downstream AUC/group comparison behavior.
