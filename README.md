# ClumPyCells: Spatial Correlation Analysis Toolbox for Spatial Omics Data 

ClumPyCells is an open-source toolbox for spatial correlation analysis based on point processing methods. By addressing cell-size differences, ClumPyCells helps to avoid biases that arise from tissue heterogeneity. The tool package provides calculations for cell aggregations using mark correlation functions, as well as post-hoc analysis tools that summarize spatial features and provides proper visualizations. 

Preprint: link to the paper

## 2. Download

ClumPyCells is publicly available at [GitHub](https://github.com/schwartzlab-methods/ClumPyCells).
```bash
git clone git@github.com:schwartzlab-methods/ClumPyCells.git
```


## 3. Required Packages

The following dependencies are required to run ClumPyCells. Install them using the `requirements.txt` file:

```bash
pip install -r requirements.txt
```

## 4. General Workflow

### Import File Format

The file need to be in `.CSV` format with the following required columns:

- `imageNum` – Indicates the image number that the cell belongs to.
- `x` – The x-coordinate of the cell.
- `y` – The y-coordinate of the cell.
- `Area` – Required if the Size Correction Algorithm is to be used.

All other columns will be treated as marks (as explained in part of the paper related to mark correlation function).

### Analyzing Single Image Using Mark Correlation Function

The `analyzeImage` function processes spatial data for a single image and computes mark cross-correlation.

**Example Usage:**

```python
from ClumPyCells import analyzeImage

analyzeImage(
    imageNum=1,
    imageData=csv_data,
    savefolder="./results",
    xrange=(0, 1000),
    yrange=(0, 1000),
    sizeCorrection=True,
    pp_criterion=lambda df: df["Area"] > 100,
    dropArea=True
)
```

---

### Running Batch Analysis for Multiple Images

The `runSpatial` function processes multiple images from a CSV file.

**Example Usage:**

```python
from ClumPyCells import runSpatial

runSpatial(
    csv_path="input_data.csv",
    savefolder="./results",
    xrange=(0, 1000),
    yrange=(0, 1000),
    sizeCorrection=True,
    pp_criterion=lambda df: df["Area"] > 100
)
```

**Output:** Each image's results are saved in the specified folder, including mark correlation analysis in `iso.csv` files.

---

## 5. Example Post-Analysis

The result saved in the `.iso` files can be used by `MarkcorrResult` for post analysis, providing summarization of different classes and corresponding visualizations. 

**Example Usage:**

```python
from ClumPyCell.Analysis.markcorrResult import MarkcorrResult

groups = {"Group1": [1, 2, 3], "Group2": [4, 5, 6]}
result_folder = "./results/"
axis_name = {"Type1": "Marker1", "Type2": "Marker2"}

markcorr = MarkcorrResult(groups, result_folder, axis_name)
combined_result = markcorr.getCombinedResult()
auc, plots = markcorr.getAUC()
boxplot = markcorr.getBoxPlot(auc, ["Group1", "Group2"])
diff_plot = markcorr.find_diff(auc, "Group1", "Group2")
```

**Output:**
- Combined result dataframe with mark correlation values.
- AUC plots visualizing spatial correlation.
- Box plots for group comparisons.
- Statistical difference plots between groups.

---

## 6. Publicly Available Datasets and Sample Codes

AML IMC dataset: [AML Data](https://zenodo.org/records/14711407)
Analysis Code: [Github](https://github.com/schwartzlab-methods/ClumPyCells_paper_figure)

For more in-depth explanations of functions, please refer to the dedicated markdown files in the [documentation](https://github.com/schwartzlab-methods/ClumPyCells/tree/main/Documentation).

---



