# ClumPyCells Function Documentation

## 1. ClumPyCells.py

### `analyzeImage`

**Description:**
Processes spatial data for a single image and computes mark cross-correlation, saving results in the specified directory.

**Parameters:**
- `imageNum` (int): Unique identifier for the image.
- `imageData` (DataFrame): A Pandas DataFrame containing spatial data with required columns (`imageNum`, `x`, `y`, and optional `Area`).
- `savefolder` (str): Directory path to save the output results.
- `xrange` (tuple): Tuple of min and max x-coordinates, e.g., `(0, 1000)`.
- `yrange` (tuple): Tuple of min and max y-coordinates, e.g., `(0, 1000)`.
- `sizeCorrection` (bool, optional): If `True`, applies size correction based on cell area; default is `False`.
- `pp_criterion` (function, optional): A function to filter specific data based on user criteria, default is `None`.
- `dropArea` (bool, optional): If `True`, removes the 'Area' column from data; default is `False`.

**Returns:**
- None (Results are saved in the specified directory).

**Example Usage:**
```python
analyzeImage(
    imageNum=1,
    imageData=df,
    savefolder="./results/",
    xrange=(0, 1000),
    yrange=(0, 1000),
    sizeCorrection=True,
    pp_criterion=lambda df: df["Area"] > 100,
    dropArea=True
)
```

---

### `runSpatial`

**Description:**
Processes multiple images from a CSV file and performs spatial analysis, generating correlation results.

**Parameters:**
- `csv_path` (str): File path to the input CSV containing spatial data.
- `savefolder` (str): Directory path to store output files.
- `xrange` (tuple): Min and max x-coordinates.
- `yrange` (tuple): Min and max y-coordinates.
- `sizeCorrection` (bool, optional): If `True`, applies size correction; default is `False`.
- `pp_criterion` (function, optional): Function to filter data; default is `None`.

**Returns:**
- None (Results are saved to the specified folder).

**Example Usage:**
```python
runSpatial(
    csv_path="input_data.csv",
    savefolder="./results/",
    xrange=(0, 1000),
    yrange=(0, 1000),
    sizeCorrection=True,
    pp_criterion=lambda df: df["Area"] > 100
)
```

---

## 2. MarkCorrResult.py

### `MarkcorrResult`

**Description:**
Handles the analysis of spatial correlation results and provides aggregation insights.

**Parameters:**
- `groups` (dict): Dictionary mapping group names to image indices.
- `resultFolder` (str): Path to the directory with analysis results.
- `axisName` (dict): Mapping of marker types to their readable names.

**Methods:**
- `getCombinedResult() -> DataFrame`: Merges results from multiple images into a single DataFrame.
- `getAUC(norm="min_mid_max", takeMean=True, plot=True, r_range=None) -> tuple`: Computes AUC of spatial correlation.
  - `norm` (str): Options: `"min_mid_max"` (default), `"log"`.
  - `takeMean` (bool): Compute mean across groups; default is `True`.
  - `plot` (bool): Whether to generate plots; default is `True`.
  - `r_range` (tuple, optional): Range of `r` values; default is `None`.
- `plotCurve(imageNum, type1, type2) -> alt.Chart`: Generates a visualization of the correlation curve.
- `getBoxPlot(auc, cols, axisName={}) -> alt.Chart`: Produces a box plot visualization.
- `find_diff(auc, col1, col2, takeMean=True, method="MW", axisName=None) -> alt.Chart`: Identifies statistical differences between groups.
- `displayImages(singlePlots: dict, order) -> alt.Chart`: Displays visualized images in a grid layout.

**Example Usage:**
```python
groups = {"Group1": [1, 2, 3], "Group2": [4, 5, 6]}
result_folder = "./results/"
axis_name = {"Type1": "Marker1", "Type2": "Marker2"}

markcorr = MarkcorrResult(groups, result_folder, axis_name)
combined_result = markcorr.getCombinedResult()
auc, plots = markcorr.getAUC()
image_grid = markcorr.displayImages(plots, [["Group1", "Group2"]])
```

---

## 3. DecisionTree.py

### Decision Tree Functions

**Description:**
Functions to perform classification using decision trees on AML and NBM datasets.

**Functions:**
- `fit_decision_tree(X, y, feature_names, bo=False, saveFig=True, saveFolder="./") -> DecisionTreeClassifier`: Trains a decision tree classifier.
  - `X` (array-like): Feature matrix.
  - `y` (array-like): Target labels.
  - `feature_names` (list): Column names of input features.
  - `bo` (bool): Enable Bayesian optimization; default is `False`.
  - `saveFig` (bool): Save decision tree visualization; default is `True`.
  - `saveFolder` (str): Directory to store the visualization.

**Example Usage:**
```python
X, y, feature_names = get_decision_tree_data(result, clinical=False, blastPercentage=False)
clf = fit_decision_tree(X, y, feature_names, bo=False, saveFig=True, saveFolder="./")
```

**Additional Functions:**
- `cross_validation(model, _X, _y, _cv=3) -> dict`: Perform cross-validation with accuracy, precision, recall, and F1 score.
- `view_decision_tree(X, y, clf, saveFolder) -> None`: Generates and saves a visual representation of the decision tree.
- `view_node_stats(X, clf, saveFolder) -> None`: Generates visual statistics of leaf nodes.
- `decision_tree(intensity=True, saveFolder="./") -> None`: Runs the entire decision tree pipeline.

**Example Usage:**
```python
decision_tree(intensity=True, saveFolder="./")
```

---

