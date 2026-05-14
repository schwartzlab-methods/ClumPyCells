import os
import shutil

import altair as alt
import dtreeviz
import numpy as np
from bayes_opt import BayesianOptimization
from sklearn.impute import SimpleImputer
from sklearn.metrics import accuracy_score, classification_report
from sklearn.model_selection import cross_validate
from sklearn.preprocessing import LabelEncoder
from sklearn.tree import DecisionTreeClassifier

from . import altairThemes as altthm

from .markcorrResult import *
from .metadata import *


def cross_validation(model, _X, _y, _cv=3):
    _scoring = ["accuracy", "precision", "recall", "f1"]
    results = cross_validate(
        estimator=model, X=_X, y=_y, cv=_cv, scoring=_scoring, return_train_score=True
    )
    return {
        "Training Accuracy scores": results["train_accuracy"],
        "Mean Training Accuracy": results["train_accuracy"].mean() * 100,
        "Training Precision scores": results["train_precision"],
        "Mean Training Precision": results["train_precision"].mean(),
        "Training Recall scores": results["train_recall"],
        "Mean Training Recall": results["train_recall"].mean(),
        "Training F1 scores": results["train_f1"],
        "Mean Training F1 Score": results["train_f1"].mean(),
        "Validation Accuracy scores": results["test_accuracy"],
        "Mean Validation Accuracy": results["test_accuracy"].mean() * 100,
        "Validation Precision scores": results["test_precision"],
        "Mean Validation Precision": results["test_precision"].mean(),
        "Validation Recall scores": results["test_recall"],
        "Mean Validation Recall": results["test_recall"].mean(),
        "Validation F1 scores": results["test_f1"],
        "Mean Validation F1 Score": results["test_f1"].mean(),
    }


def write_node(tree, node, left, leafNum):
    outputStr = ""
    if tree.children_left[node] == tree.children_right[node]:
        outputStr += '[{"_item": ['
        for i in range(int(leafNum.loc[node, "aml"])):
            outputStr += '{"_barcode": { "unCell": "aml"},"_cellRow": {"unRow": 5}},'
        for i in range(int(leafNum.loc[node, "nbm"])):
            outputStr += '{"_barcode": { "unCell": "nbm"},"_cellRow": {"unRow": 5}},'
        outputStr = outputStr.rstrip(",")
        outputStr += '],"_significance": null,"_distance": 5},[]]'
        if left:
            outputStr += ","
        return outputStr
    else:
        outputStr += '[{"_item": null,"_significance": null,"_distance": 5},['
        if tree.children_left[node] != -1:
            outputStr += write_node(tree, tree.children_left[node], True, leafNum)
        if tree.children_right[node] != -1:
            outputStr += write_node(tree, tree.children_right[node], False, leafNum)
        outputStr += "]]"
        if left:
            outputStr += ","
        return outputStr


def combineBlastPercentage(data):
    blastData = pd.read_csv(os.path.join(HOMEDIR, "Data/AML_blast_stroma_database.csv"))
    # blast_percentage = pd.cut(blastData["Blast_percentatge"], bins=[0, 75, 100], labels= ["low", "high"])
    blast_percentage = blastData["Blast_percentage"]
    # blast_percentage = np.where(blast_percentage > 75, 1, 0)
    blast_percentage = pd.DataFrame(blast_percentage, columns=["Blast_percentage"])
    blast_percentage["Biopsy_number"] = blastData["Biopsy_code"].apply(
        lambda x: x.split("-")[1]
    )
    image2biop = imageNum_to_biopNum()
    data["Biopsy_number"] = data["imageNum"].apply(
        lambda x: image2biop[x] if 0 <= x < len(image2biop) else None
    )
    data = data.merge(
        blast_percentage, left_on="Biopsy_number", right_on="Biopsy_number"
    )
    return data


def view_decision_tree(X, y, clf, saveFolder):
    f = open(
        os.path.join(HOMEDIR, os.path.join(saveFolder, "cluster_list.json")),
        "w",
    )
    f.write(
        '[[{"_barcode":{"unCell":"node1"},"_cellRow":{"unRow":5}},[{"unCluster":2},{"unCluster":1},{"unCluster":0}]],[{"_barcode":{"unCell":"node2"},"_cellRow":{"unRow":1406}},[{"unCluster":3},{"unCluster":1},{"unCluster":0}]],[{"_barcode":{"unCell":"node1"},"_cellRow":{"unRow":1406}},[{"unCluster":5},{"unCluster":4},{"unCluster":0}]],[{"_barcode":{"unCell":"node2"},"_cellRow":{"unRow":1406}},[{"unCluster":6},{"unCluster":4},{"unCluster":0}]]]'
    )
    f.close()

    f = open(
        os.path.join(HOMEDIR, os.path.join(saveFolder, "cluster_tree.json")),
        "w",
    )
    leafNum = get_leaf_number(X, y, clf)
    init_str = write_node(clf.tree_, 0, False, leafNum)
    f.write(init_str)
    f.close()


def get_leaf_number(X, y, clf):
    node_stat = pd.DataFrame()
    node_stat["leaf"] = clf.apply(X)
    node_stat["y"] = y
    aml = node_stat.loc[node_stat["y"] == 1]
    nbm = node_stat.loc[node_stat["y"] == 0]
    aml = aml.groupby("leaf").count().rename({"y": "aml"}, axis=1)
    nbm = nbm.groupby("leaf").count().rename({"y": "nbm"}, axis=1)
    tot_stat = pd.concat([aml, nbm], axis=1).sort_values(by="leaf").fillna(0)
    return tot_stat


def fit_decision_tree(
    X,
    y,
    feature_names,
    bo=False,
    saveFig=True,
    saveFolder="./",
    impute=True,
    max_depth=4,
    max_features=0.6701,
    bo_init_points=100,
    bo_n_iter=30,
    class_names=None,
    target_name="group",
    random_state=123,
):
    params = {"max_depth": int(max_depth), "max_features": float(max_features)}
    if impute:
        X = SimpleImputer(strategy="constant", fill_value=0).fit_transform(X)
    else:
        combined = np.hstack((X, y.reshape(-1, 1)))
        mask = ~np.isnan(combined).any(axis=1)
        X = X[mask]
        y = y[mask]
    if bo:

        def cart_bo(max_depth, max_features):
            params_cart = {}
            params_cart["max_depth"] = round(max_depth)
            params_cart["max_features"] = max_features
            score = cross_validate(
                estimator=DecisionTreeClassifier(
                    random_state=random_state, **params_cart, criterion="entropy"
                ),
                X=X,
                y=y,
                scoring=["accuracy"],
                cv=3,
                return_train_score=False,
            )["test_accuracy"].mean()
            return score

        # Run Bayesian Optimization
        params_cart = {"max_depth": (3, 10), "max_features": (0.6, 1)}
        bo_result = BayesianOptimization(cart_bo, params_cart, random_state=111)
        bo_result.maximize(init_points=int(bo_init_points), n_iter=int(bo_n_iter))
        max_result = bo_result.max
        if max_result is None:
            raise RuntimeError("Bayesian optimization did not return parameters")
        params_tuned = max_result["params"]
        params_tuned["max_depth"] = round(params_tuned["max_depth"])
        params = params_tuned

    decision_tree_model = DecisionTreeClassifier(
        criterion="entropy",
        max_depth=params["max_depth"],
        max_features=params["max_features"],
        random_state=random_state,
    )

    clf = decision_tree_model.fit(X, y)

    if saveFig:
        viz = dtreeviz.model(
            clf,
            X,
            y,
            target_name=target_name,
            feature_names=feature_names,
            class_names=class_names or [str(value) for value in np.unique(y)],
        )
        viz.view(orientation="LR", scale=0.8, fancy=True).save(
            os.path.join(saveFolder, "tree.svg")
        )

        leaves = np.asarray(clf.apply(X))
        leaf_counts = pd.Series(leaves.tolist()).value_counts()
        most_common_leaf = leaf_counts.idxmax()
        sample_idx = np.where(leaves == most_common_leaf)[0][0]
        print(f"Most common leaf: {most_common_leaf}, sample index: {sample_idx}")
        # Visualize the path for the largest leaf

        viz_path = os.path.join(saveFolder, "largest_leaf_path.svg")
        viz.view(x=X[sample_idx], show_just_path=True).save(viz_path)
        print(f"dtreeviz path for largest leaf saved to: {viz_path}")

    return clf


def decision_tree_from_feature_table(
    feature_table,
    target_col,
    saveFolder="./",
    feature_columns=None,
    bo=False,
    impute=True,
    save_fig=True,
    max_depth=4,
    max_features=0.6701,
    bo_init_points=25,
    bo_n_iter=10,
    random_state=123,
):
    os.makedirs(saveFolder, exist_ok=True)
    data = feature_table.copy()
    if target_col not in data.columns:
        raise ValueError(f"target column '{target_col}' was not found")

    if feature_columns is None:
        feature_columns = [column for column in data.columns if column != target_col]
    feature_columns = [column for column in feature_columns if column in data.columns]
    if not feature_columns:
        raise ValueError("Select at least one feature column for decision tree")

    X_df = data[feature_columns].apply(pd.to_numeric, errors="coerce")
    if impute:
        valid_mask = pd.Series(True, index=X_df.index)
    else:
        valid_mask = ~X_df.isna().any(axis=1)
        data = data.loc[valid_mask].copy()
        X_df = X_df.loc[valid_mask].copy()
    y_raw = data[target_col].astype(str)
    label_encoder = LabelEncoder()
    y = label_encoder.fit_transform(y_raw)
    class_names = label_encoder.classes_.tolist()

    clf = fit_decision_tree(
        X_df.to_numpy(),
        y,
        list(X_df.columns),
        bo=bool(bo),
        saveFig=bool(save_fig),
        saveFolder=saveFolder,
        impute=bool(impute),
        max_depth=int(max_depth),
        max_features=float(max_features),
        bo_init_points=int(bo_init_points),
        bo_n_iter=int(bo_n_iter),
        class_names=class_names,
        target_name=target_col,
        random_state=int(random_state),
    )

    X_fit = SimpleImputer(strategy="constant", fill_value=0).fit_transform(X_df)
    predictions = clf.predict(X_fit)
    feature_importance = pd.Series(clf.feature_importances_, index=X_df.columns)
    feature_importance = feature_importance.sort_values(ascending=False)
    feature_importance.to_csv(os.path.join(saveFolder, "feature_importance.csv"))

    data.to_csv(os.path.join(saveFolder, "training_samples.csv"), index=False)
    metrics = pd.DataFrame(
        {
            "metric": ["training_accuracy", "n_samples", "n_features"],
            "value": [accuracy_score(y, predictions), len(data), len(feature_columns)],
        }
    )
    metrics.to_csv(os.path.join(saveFolder, "model_metrics.csv"), index=False)
    report = classification_report(
        y,
        predictions,
        target_names=class_names,
        output_dict=True,
        zero_division=0,
    )
    pd.DataFrame(report).transpose().to_csv(
        os.path.join(saveFolder, "classification_report.csv")
    )
    return clf, feature_importance, metrics


def decision_tree_from_markcorr_groups(
    groups,
    resultFolder,
    axisName=None,
    norm="min_mid_max",
    saveFolder="./",
    **kwargs,
):
    os.makedirs(saveFolder, exist_ok=True)
    axisName = axisName or {}
    result = MarkcorrResult(groups=groups, resultFolder=resultFolder, axisName=axisName)
    auc_by_group, _ = result.getAUC(norm=norm, plot=False)
    rows = []
    for group_name, auc_table in auc_by_group.items():
        table = (
            auc_table.transpose().reset_index().rename(columns={"index": "image_id"})
        )
        table["image_id"] = (
            table["image_id"].astype(str).str.replace("image_", "", regex=False)
        )
        table["target"] = group_name
        rows.append(table)
    if not rows:
        raise ValueError("No group/image data available for decision tree")
    feature_table = pd.concat(rows, ignore_index=True)
    feature_table.to_csv(
        os.path.join(saveFolder, "decision_tree_input_table.csv"), index=False
    )
    feature_columns = [
        column
        for column in feature_table.columns
        if column not in {"image_id", "target"}
    ]
    return decision_tree_from_feature_table(
        feature_table,
        target_col="target",
        saveFolder=saveFolder,
        feature_columns=feature_columns,
        **kwargs,
    )


def get_decision_tree_data(
    result,
    clinical=False,
    blastPercentage=False,
    removeAdipocytes=True,
    simplifyNames=True,
):
    data = result.getCombinedResult()
    if blastPercentage:
        data = combineBlastPercentage(data=data)
        data = data.drop(["Biopsy_number"], axis=1)
    if removeAdipocytes:
        data = data.loc[:, ~data.columns.str.contains("Intensity_Adipocytes")]

    data["y"] = data["imageNum"].apply(lambda x: 1 if x < 36 else 0)
    y = data["y"].values
    data = data.drop(["imageNum", "y"], axis=1)

    if simplifyNames:
        data.columns = data.columns.str.replace("Intensity_", "", regex=True)

    X = data.values
    feature_names = data.columns
    return X, y, feature_names


def view_node_stats(X, clf, saveFolder):
    # Split the dataset into chunks and collect leaf indices
    num_images = len(X) // 513
    indices = np.array_split(X, num_images)
    leaf_indices = [clf.apply(chunk) for chunk in indices]
    node_by_image_full = pd.DataFrame(leaf_indices).transpose()

    # Reshape data for plotting
    node_by_image = node_by_image_full.melt(var_name="image_num", value_name="leaf_num")

    # Create color schemes using get_colour_scheme
    num_images = node_by_image["image_num"].nunique()
    num_leaves = node_by_image["leaf_num"].nunique()

    image_colors = altthm.get_colour_scheme("tab20", num_images)
    leaf_colors = altthm.get_colour_scheme("tab20", num_leaves)

    # First bar chart (leaf number by count, colored by image number)
    p1 = (
        alt.Chart(node_by_image)
        .mark_bar()
        .encode(
            x="leaf_num:N",
            y="count()",
            color=alt.Color(
                "image_num:Q",
                scale=alt.Scale(domain=list(range(num_images)), range=image_colors),
                legend=alt.Legend(gradientLength=200, title="Image Number"),
            ),
        )
        .properties(title="Leaf Node Distribution per Image", width=370)
    )

    # Second bar chart (image number by count, colored by leaf number)
    p2 = (
        alt.Chart(node_by_image)
        .mark_bar()
        .encode(
            y="image_num:N",
            x="count()",
            color=alt.Color(
                "leaf_num:Q",
                scale=alt.Scale(domain=list(range(num_leaves)), range=leaf_colors),
                legend=alt.Legend(gradientLength=200, title="Leaf Number"),
            ),
        )
        .properties(title="Image Distribution per Leaf Node")
    )

    # Combine the two charts horizontally
    plot = (
        alt.hconcat(p1, p2)
        .resolve_legend(color="independent")
        .resolve_scale(color="independent")
    )

    # Save the plot as an HTML file for interactivity
    plot_path = os.path.join(saveFolder, "nodeStats.html")
    plot.save(plot_path)

    print(f"Combined plot saved to {plot_path}")


def decision_tree(
    intensity=True,
    saveFolder="./",
    bo=True,
    impute=True,
    save_fig=True,
    max_depth=4,
    max_features=0.6701,
    bo_init_points=100,
    bo_n_iter=30,
):
    if intensity:
        result = AMLResult(sizeCorrection=True, intensity=True)
        folder = os.path.join(saveFolder, "intensity/")
        folder = saveFolder + "intensity/"
    else:
        result = AMLResult(sizeCorrection=True, intensity=False)
        folder = os.path.join(saveFolder, "cellType/")

    X, y, feature_names = get_decision_tree_data(
        result,
        clinical=False,
        blastPercentage=False,
        removeAdipocytes=True,
        simplifyNames=True,
    )

    clf = fit_decision_tree(
        X,
        y,
        feature_names,
        bo=bool(bo),
        saveFig=bool(save_fig),
        saveFolder=folder,
        impute=bool(impute),
        max_depth=int(max_depth),
        max_features=float(max_features),
        bo_init_points=int(bo_init_points),
        bo_n_iter=int(bo_n_iter),
        class_names=["NBM", "AML"],
        target_name="type",
    )
    feature_importance = pd.Series(clf.feature_importances_, index=feature_names)
    feature_importance = feature_importance.sort_values(ascending=False)
    feature_importance.to_csv(os.path.join(folder, "feature_importance.csv"))
    view_decision_tree(X, y, clf, folder)
    view_node_stats(X, clf, folder)
    # get_leaf_number(X, y, clf)
