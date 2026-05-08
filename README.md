# ClumPyCells: Spatial Correlation Analysis Toolbox with Cell Size Correction for Spatial Omics Data

ClumPyCells is an open-source toolbox for spatial correlation analysis based on point-process methods. By accounting for cell-size differences, ClumPyCells avoids biases that arise from tissue heterogeneity. The package provides:

- a Python apckeage (`Markcorr`) with cell-size correction, and
- a post-hoc analysis layer (`Analysis`) that summarises spatial features across images and groups, with publication-ready visualisations.

Preprint: (https://www.biorxiv.org/content/10.64898/2026.03.26.714529v1)

---

## 1. Installation

```bash
git clone git@github.com:schwartzlab-methods/ClumPyCells.git
cd ClumPyCells

# (recommended) create a virtual environment
python3.12 -m venv CLPC
source CLPC/bin/activate

pip install -r requirement.txt
```

The repository ships with a `config.json` at the root that the bundled analysis scripts (`AML_analysis.py`, `Melanoma_analysis.py`, `runPermutation.py`) consult to locate dataset folders. Update its `HOMEDIR` field to the absolute path of your ClumPyCells checkout (the trailing slash is required):

```json
{ "HOMEDIR": "/absolute/path/to/ClumPyCells/" }
```

The library itself can find `config.json` automatically; you can also override the search by setting the `CLUMPYCELLS_CONFIG` environment variable to an explicit JSON file.

---

## 2. Repository layout

```
ClumPyCells/
├── ClumPyCells/                # importable Python package
│   ├── ClumPyCells.py          # high-level driver (analyzeImage, runSpatial)
│   ├── Markcorr/               # mark cross-correlation engine
│   │   ├── window.py           # rectangular observation window
│   │   ├── pointPattern.py     # point pattern container with marks/diameters
│   │   ├── breakpts.py         # radius (r) discretisation
│   │   ├── closepairs.py       # close-pair enumeration with size correction
│   │   ├── unnormdensity.py    # weighted KDE used by sewsmod
│   │   └── markcorr.py         # main markcorr() routine
│   └── Analysis/               # post-hoc analysis & plotting
│       ├── markcorrResult.py   # AUC heatmaps, find_diff, box plots
│       ├── metadata.py         # dataset metadata (axisName, group splits)
│       ├── decisionTree.py     # downstream classifier helpers
│       ├── survivalAnalysis.py # Kaplan–Meier overlays
│       └── interactivePlot.py  # altair-based image viewer
├── Data/                       # bundled datasets and test fixtures
├── Result/                     # default output tree (created on first run)
├── tests/                      # pytest suite (see §6)
├── config.json                 # HOMEDIR for bundled analysis scripts
├── AML_analysis.py             # example end-to-end study (AML)
├── Melanoma_analysis.py        # example end-to-end study (Melanoma)
└── runPermutation.py           # group-label permutation null distribution
```

---

## 3. Input format

A single CSV file describing one or more images. The core package uses canonical column names internally, but the Streamlit app, CLI, and `runSpatial` can map your file's column names into that format:

| Internal role | Required?                      | Description                                                                                                                                           |
| ------------- | ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `x`, `y`      | Yes                            | Cell centroid coordinates (must lie inside the supplied window). Your CSV can call these columns something else, such as `centroid_x` / `centroid_y`. |
| `ImageNum`    | No                             | Image identifier. If no image column is selected, the whole CSV is treated as one image and written to `image_1/`.                                    |
| `Area`        | Only for `sizeCorrection=True` | Cell area in the same units as `x`/`y` squared. Your CSV can use a different name, such as `cell_area`.                                               |
| mark columns  | Yes                            | Continuous intensities, categorical cell types, or any other selected features used for mark cross-correlation.                                       |

Only the columns selected as marks become the rows/columns of the resulting mark cross-correlation matrix. Categorical (non-integer) columns are automatically one-hot encoded.

---

## 4. End-to-end usage

### 4.1 Run mark correlation on every image

```python
from ClumPyCells.ClumPyCells import runSpatial

runSpatial(
    csv_path="cells.csv",
    savefolder="results/",
    xrange=[0, 1056],
    yrange=[0, 642],
    sizeCorrection=True,                          # True ⇒ derive cell diameters from Area
    pp_criterion=lambda df: df["Area"] > 100,     # optional: cells used as occluders
    max_workers=4,                                # parallel images; use 1 for serial/HPC memory limits
    x_col="centroid_x",                           # optional input-column mapping
    y_col="centroid_y",
    image_col="sample_id",                       # use None if the CSV is one image
    area_col="cell_area",
    mark_columns=["cell_type", "CD3", "MPO"],    # optional; defaults to all non-metadata columns
)
```

Outputs (one folder per image):

```
results/
├── image_1/
│   ├── r.csv      # the radius vector used for kmm(r)
│   └── iso.csv    # isotropic-corrected kmm(r) for every (mark_i, mark_j) pair
├── image_2/
└── …
```

`runSpatial` is a multithreaded wrapper around `analyzeImage`, which can also be used for a single image when you already have a DataFrame in memory:

```python
from ClumPyCells.ClumPyCells import analyzeImage

analyzeImage(
    imageNum=1,
    imageData=df,
    savefolder="results/",
    xrange=[0, 1056],
    yrange=[0, 642],
    sizeCorrection=True,
)
```

### 4.2 Aggregate, compare, visualise

```python
from ClumPyCells.Analysis.markcorrResult import MarkcorrResult

groups   = {"AML": list(range(36)), "NBM": list(range(36, 51))}
axisName = {
    "Intensity_CD34": "CD34",
    "Intensity_CD3":  "CD3",
    "Intensity_MPO":  "MPO",
    # …
}

mc = MarkcorrResult(groups=groups, resultFolder="results/", axisName=axisName)

# AUC under the kmm(r) curve, per (image, mark-pair); returns dict per group + altair heatmaps
auc, heatmaps = mc.getAUC(norm="min_mid_max")

# Group comparison: Mann–Whitney U with Benjamini–Hochberg FDR (or "perm" for Fisher–Pitman)
diff_plot = MarkcorrResult.find_diff(
    auc, "AML", "NBM",
    method="MW",
    axisName=axisName,
    saveCsv="results/AML_vs_NBM_diffchart.csv",  # optional, no longer hardcoded
)

# Box plot of AUC per pair, coloured by group
boxplot = mc.getBoxPlot(auc, ["AML", "NBM"], axisName=axisName)
```

### 4.3 Working with the bundled datasets

Two ready-made subclasses, `AMLResult` and the Melanoma equivalent, preset `groups`, `axisName`, and `resultFolder` for the studies described in the manuscript. See `AML_analysis.py` and `Melanoma_analysis.py` for the full reproductions.

---

## 5. How the algorithm works

1. **Point pattern.** Every image is wrapped in a `pointPattern(x, y, diameter, marks, window)`. Points outside the rectangular `window` are dropped; missing diameters default to `0`.
2. **Radii.** `r` is built by `handle_r_b` using Ripley’s rule (`min(0.25·minEdge, sqrt(1000 / (π·λ)))`) and discretised into 513 evenly spaced values.
3. **Close pairs.** `closepairs` enumerates all `(i, j)` with distance ≤ `rmax`. With `sizeCorrection=True`, the centre-to-centre distance is reduced by `(d_i + d_j) / 2`; with `pp_criterion`, the segment is further shortened by any large "occluder" cell that intersects the line between the two points.
4. **Mark correlation.** For every mark pair `(coli, colj)`,

    $$kmm(r) = \frac{\widehat{\mathrm{KDE}}_{\,wt\cdot f}(r)}{E_f \cdot \widehat{\mathrm{KDE}}_{\,wt}(r)}, \qquad f = m_i \cdot m_j, \qquad E_f = \overline{m_i}\,\overline{m_j}$$

    with edge-correction weights `wt` from Ripley's isotropic correction (`edgecorrection`) or the translation correction (`edgetrans`).
5. **Aggregation.** `MarkcorrResult.getAUC` integrates `log(kmm(r))` (or `min_mid_max`-normalised values) over the radius range to obtain a single AUC per image and pair. `find_diff` then runs Mann–Whitney U or a Fisher–Pitman permutation test across the groups and applies BH-FDR.

### Numerical invariant

> Running `analyzeImage` (or `runSpatial`) with `sizeCorrection=True` on data whose `Area` column is identically zero **must** produce the exact same `iso.csv` as running it with `sizeCorrection=False`.

This invariant is enforced by `tests/test_pipeline.py::test_runSpatial_size_correction_with_zero_area_matches_no_size_correction` and `tests/test_markcorr.py::test_markcorr_size_correction_zero_equals_no_size_correction`, which compare the two code paths element-wise with `atol=1e-12`. The exact algebra of the size-correction step (radius subtraction, occluder chord, segment-membership check, clamp-at-zero) is pinned down by the eight unit tests in [tests/test_size_correction.py](tests/test_size_correction.py).

### Performance

The close-pair enumeration is fully vectorised (NumPy), so a 5000-cell image now runs in seconds rather than the minutes the original double-Python-loop took. Both `runSpatial` (across images) and `markcorr` (across mark pairs) print a tqdm progress bar so you can see how a long run is progressing.

---

## 6. Streamlit UI and Terminal CLI

A point-and-click front-end for the whole pipeline lives in [streamlit_app.py](streamlit_app.py):

```bash
source CLPC/bin/activate
streamlit run streamlit_app.py
```

The UI has three pages:

| Page                       | What it does                                                                                                                                                                                                                                                     |
| -------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **1. Run markcorr**        | Upload a CSV, choose which columns represent x/y/image/area, select mark columns, set parallel image workers, toggle size correction (with an optional area-threshold occluder rule), and launch `runSpatial`. A live progress bar reports each completed image. |
| **2. Downstream analysis** | Point at any markcorr output folder, define groups by selecting their image numbers, edit display labels, and explore the AUC heatmaps, group-vs-group `find_diff` (Mann–Whitney + BH FDR or Fisher–Pitman permutation), and per-pair box plots.                 |
| **3. Permutation test**    | Within-image label shuffling for any subset of images and any number of permutations; outputs land under `<savefolder>/perm_runs/perm_<seed>/image_<n>/iso.csv`.                                                                                                 |

For HPC clusters or offline environments where a browser/server UI is inconvenient, use the terminal entry point [clumpycells_cli.py](clumpycells_cli.py). It needs only the installed Python environment; no internet access is required at run time.

Inspect a CSV before choosing columns:

```bash
source CLPC/bin/activate
python clumpycells_cli.py inspect-csv --csv cells.csv
```

Run markcorr from a terminal or batch script:

```bash
python clumpycells_cli.py run-markcorr \
    --csv cells.csv \
    --out results_hpc \
    --x-col centroid_x \
    --y-col centroid_y \
    --image-col sample_id \
    --area-col cell_area \
    --mark cell_type \
    --mark CD3 \
    --mark MPO \
    --xrange 0 1056 \
    --yrange 0 642 \
    --size-correction \
    --pp-area-threshold 100 \
    --max-workers "${SLURM_CPUS_PER_TASK:-1}"
```

If the CSV represents one image, pass an empty image column:

```bash
python clumpycells_cli.py run-markcorr \
    --csv cells.csv \
    --out results_single_image \
    --image-col "" \
    --mark cell_type \
    --max-workers 1
```

Set `--max-workers` to the number of CPU cores allocated by the scheduler. Use `--max-workers 1` for serial execution or memory-constrained jobs.

---

## 7. Tests

The pytest suite covers the markcorr engine, the full CSV → AUC pipeline, and the analysis layer:

```bash
source CLPC/bin/activate
pytest tests/ -v
```

What the suite exercises:

| File                               | What it tests                                                                                                                                                                                                                      |
| ---------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `tests/test_markcorr.py`           | `window`, `pointPattern` filtering, frozen baseline regression of the BetaCells curves, the size-correction equivalence invariant, behaviour on a Poisson pattern, pickle-cache cleanup.                                           |
| `tests/test_pipeline.py`           | `runSpatial` produces the expected per-image artefacts; identical inputs yield identical outputs; `sizeCorrection=True` with zero areas matches `sizeCorrection=False`; `MarkcorrResult.getAUC`, `find_diff`, `getCombinedResult`. |
| `tests/test_legacy_smoke.py`       | The legacy top-level `Test.py` is importable without crashing.                                                                                                                                                                     |
| `tests/baseline/betacells_iso.npz` | Frozen reference iso curves used by the regression test. Regenerate only after a deliberate numerical change.                                                                                                                      |

All tests stage their outputs under pytest’s `tmp_path`; they never touch `Result/`.

### Re-generating the regression baseline

If you intentionally change the numerics (e.g. a faster KDE) and need to refresh the captured curves:

```python
import json, os, shutil, numpy as np, pandas as pd
from ClumPyCells.Markcorr.markcorr import markcorr
from ClumPyCells.Markcorr.pointPattern import pointPattern
from ClumPyCells.Markcorr.window import window

with open("Data/TestData/BetaCells.json") as f:
    data = json.load(f)
W = window([28.08, 778.08], [16.2, 1007.02])
m = pd.DataFrame({"type": data["type"]}); m["type"] = m["type"].astype("category")
shutil.rmtree("tests/baseline/_tmp", ignore_errors=True)
os.makedirs("tests/baseline/_tmp")
p = pointPattern(list(data["x"]), list(data["y"]), None, W, m)
r, funs = markcorr(p, savefolder="tests/baseline/_tmp/", remove_zeros=False,
                   correction=["isotropic"], saveImage=False, saveCache=False)
np.savez_compressed("tests/baseline/betacells_iso.npz",
                    r=np.array(r, dtype=float),
                    keys=np.array(list(funs.keys())),
                    values=np.array([funs[k][0] for k in funs], dtype=float))
shutil.rmtree("tests/baseline/_tmp")
```

---

## 8. Bug fixes shipped with this revision

While instrumenting the test suite the following bugs were identified and fixed:

- **`closePpairs`**: the cell-size correction had three issues: (i) it was a pure-Python `O(n^2 · |pp|)` double loop, (ii) the chord through an occluder was only subtracted when `chord < d2` rather than clamping `d2` to zero, and (iii) it had no segment-membership check, so an occluder lying on the *infinite* line through `(p_i, p_j)` but outside the segment still removed distance. The function is now fully vectorised, clamps to zero, and checks the foot-of-perpendicular falls within `[0, ||p_2 - p_1||]`. ([ClumPyCells/Markcorr/closepairs.py](ClumPyCells/Markcorr/closepairs.py))
- **`pointPattern`**: the in-window filter used `pop` while iterating and reassigned `marks = marks.drop(..., inplace=True)` (which returns `None`). It now uses a vectorised boolean mask. ([ClumPyCells/Markcorr/pointPattern.py](ClumPyCells/Markcorr/pointPattern.py))
- **`analyzeImage`**: when `sizeCorrection=True` and `pp_criterion=None`, `pp` was never bound, raising `UnboundLocalError`. The helper now initialises `pp = None` upfront and avoids constructing the point pattern twice. ([ClumPyCells/ClumPyCells.py](ClumPyCells/ClumPyCells.py))
- **`MarkcorrResult.find_diff`**: the result CSV was always written to `HOMEDIR + "/Result/Test/diffchart.csv"`. The path is now controlled by an optional `saveCsv=` argument. ([ClumPyCells/Analysis/markcorrResult.py](ClumPyCells/Analysis/markcorrResult.py))
- **`Analysis.metadata`**: `config.json` was opened with a relative path at import time, breaking any `from ClumPyCells.Analysis... import` issued from a different cwd. The loader now searches `$CLUMPYCELLS_CONFIG`, the cwd, and the package root, and the `altairThemes` import no longer relies on a typo'd absolute path. ([ClumPyCells/Analysis/metadata.py](ClumPyCells/Analysis/metadata.py))
- **`Test.py`**: a side-effecting `modify_csv_index(...)` call at module top level crashed on import. It now lives behind an `if __name__ == "__main__":` guard. ([Test.py](Test.py))
- **`unnormdensity.density`**: the noisy `WARNING: sum(weights) != 1` was emitted on every pair; the warning is no longer raised because the caller-driven KDE intentionally does not normalise weights.

---

## 9. Publicly available datasets and example studies

- AML IMC dataset: <https://zenodo.org/records/14711407>
- Reproducible figures: <https://github.com/schwartzlab-methods/ClumPyCells_paper_figure>

For more detailed function references, see [Documentation/Documentation.md](Documentation/Documentation.md).
