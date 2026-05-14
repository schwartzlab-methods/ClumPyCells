import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter, statistics
from sksurv.preprocessing import OneHotEncoder
from statsmodels.stats.multitest import multipletests

from .markcorrResult import *
from .metadata import *


def _list_result_images(result_folder):
    images = []
    for child in sorted(os.listdir(result_folder)):
        if child.startswith("image_") and os.path.isdir(
            os.path.join(result_folder, child)
        ):
            images.append(child.split("_", 1)[1])
    return images


def _coerce_event(series, positive_value=None):
    if positive_value is not None and str(positive_value).strip() != "":
        return series.astype(str) == str(positive_value)
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().sum() == len(series):
        return numeric > 0
    positives = {"1", "true", "yes", "y", "dead", "deceased", "event", "relapse"}
    return series.astype(str).str.strip().str.lower().isin(positives)


def read_markcorr_auc_features(
    result_folder,
    image_numbers=None,
    norm="min_mid_max",
    r_range=None,
):
    result_folder = result_folder.rstrip("/") + "/"
    if image_numbers is None:
        image_numbers = _list_result_images(result_folder)
    if not image_numbers:
        raise ValueError("No image_<id>/iso.csv folders were found")

    iso_tables = {}
    combined_values = []
    for image_id in image_numbers:
        image_key = str(image_id).replace("image_", "")
        iso_path = os.path.join(result_folder, f"image_{image_key}", "iso.csv")
        if not os.path.exists(iso_path):
            continue
        iso = pd.read_csv(iso_path).drop(["Unnamed: 0"], axis=1, errors="ignore")
        iso_min = iso.where(iso > 0).min().min()
        if pd.isna(iso_min):
            iso_min = 1
        iso = iso.where(iso > 0, iso_min)
        if r_range is not None:
            r_path = os.path.join(result_folder, f"image_{image_key}", "r.csv")
            r = pd.read_csv(r_path).drop(["Unnamed: 0"], axis=1, errors="ignore")
            iso = iso.loc[(r["r"] >= r_range[0]) & (r["r"] <= r_range[1])]
        iso_tables[image_key] = iso
        combined_values.extend(iso.to_numpy().flatten().tolist())

    if not iso_tables:
        raise ValueError("No readable iso.csv files were found")

    max_data = min_data = None
    if norm == "min_mid_max":
        values = np.asarray(
            [value for value in combined_values if value != 1], dtype=float
        )
        if values.size:
            q1 = np.nanpercentile(values, 25)
            q3 = np.nanpercentile(values, 75)
            iqr = q3 - q1
            filtered = values[(values >= q1 - 5 * iqr) & (values <= q3 + 5 * iqr)]
            max_data = np.nanmax(filtered)
            min_data = np.nanmin(filtered)

    rows = []
    for image_key, iso in iso_tables.items():
        if norm == "log":
            values = np.log2(iso).sum()
        elif norm == "min_mid_max" and max_data is not None and min_data is not None:
            scaled = iso.map(
                lambda x: (
                    -np.abs((1 - x)) / np.abs(1 - min_data)
                    if x <= 1
                    else np.abs(x - 1) / np.abs(max_data - 1)
                )
            )
            scaled.mask(scaled > 1, 1, inplace=True)
            scaled.mask(scaled < -1, -1, inplace=True)
            values = scaled.sum()
        else:
            values = iso.sum()
        row = values.to_dict()
        row["image_id"] = image_key
        rows.append(row)
    return pd.DataFrame(rows)


def prepare_user_survival_table(
    spatial_features,
    clinical_data,
    clinical_id_col,
    duration_col,
    event_col,
    event_positive_value=None,
    covariate_cols=None,
    image_to_clinical=None,
    mapping_image_col=None,
    mapping_clinical_col=None,
):
    spatial = spatial_features.copy()
    spatial["image_id"] = spatial["image_id"].astype(str)
    clinical = clinical_data.copy()
    clinical[clinical_id_col] = clinical[clinical_id_col].astype(str)

    if image_to_clinical is not None:
        mapping = image_to_clinical[[mapping_image_col, mapping_clinical_col]].copy()
        mapping[mapping_image_col] = mapping[mapping_image_col].astype(str)
        mapping[mapping_clinical_col] = mapping[mapping_clinical_col].astype(str)
        spatial = spatial.merge(
            mapping,
            left_on="image_id",
            right_on=mapping_image_col,
            how="left",
        )
        spatial["clinical_id"] = spatial[mapping_clinical_col]
    else:
        spatial["clinical_id"] = spatial["image_id"]

    clinical = clinical.rename(
        columns={
            clinical_id_col: "clinical_id",
            duration_col: "duration",
            event_col: "event_raw",
        }
    )
    clinical["duration"] = pd.to_numeric(clinical["duration"], errors="coerce")
    clinical["event"] = _coerce_event(clinical["event_raw"], event_positive_value)

    covariate_cols = covariate_cols or []
    covariates = clinical[["clinical_id", "duration", "event"]].copy()
    for column in covariate_cols:
        if column in clinical_data.columns:
            covariates[column] = clinical_data[column]
    covariates = pd.get_dummies(
        covariates,
        columns=[
            column
            for column in covariate_cols
            if column in covariates.columns
            and not pd.api.types.is_numeric_dtype(covariates[column])
        ],
        drop_first=True,
    )

    merged = spatial.merge(covariates, on="clinical_id", how="inner")
    merged = merged.dropna(subset=["duration", "event"])
    merged["event"] = merged["event"].astype(bool)
    return merged.fillna(0)


def get_clinical_data():
    clinical = pd.read_csv(
        os.path.join(HOMEDIR, "Data/Clinical_data.csv"),
        usecols=[
            "Age",
            "Biopsy Number",
            "MSC density",
            "RELAPSE",
            "OVERALL.SURVIVAL",
            "TTR",
            "status.at.last.follow.up",
        ],
    ).rename({"OVERALL.SURVIVAL": "OST", "MSC density": "MSC_density"}, axis=1)
    biop = [num.split("-")[1] for num in clinical["Biopsy Number"]]
    clinical["Biopsy Number"] = biop

    # split the age group into young middle and old
    clinical["Age_group"] = pd.cut(
        clinical["Age"], bins=[0, 50, 70, 100], labels=["young", "middle", "old"]
    )

    # convert columns to categorical columns
    clinical["MSC_density"] = clinical["MSC_density"].astype(
        "category",
    )
    clinical["RELAPSE"] = clinical["RELAPSE"].astype("category")
    clinical["status.at.last.follow.up"] = clinical["status.at.last.follow.up"].astype(
        "category"
    )

    cat_data = clinical[["Age_group"]]
    cat_one_hot = OneHotEncoder().fit_transform(cat_data)
    clinical = pd.concat(
        [
            clinical[
                [
                    "status.at.last.follow.up",
                    "Biopsy Number",
                    "OST",
                    "TTR",
                    "Age",
                ]
            ],
            cat_one_hot,
        ],
        axis=1,
    )
    clinical["survival_status"] = np.where(
        clinical["status.at.last.follow.up"] == "dead", True, False
    )
    clinical = clinical.drop(["status.at.last.follow.up"], axis=1)
    clinical = clinical.drop(["TTR"], axis=1).dropna()
    return clinical


def get_spatial_clinical_combined(clinical, auc_t):
    biopsy_num = imageNum_to_biopNum()
    auc_t["Biopsy Number"] = biopsy_num[:36]
    survival_data = auc_t.merge(
        clinical, left_on="Biopsy Number", right_on="Biopsy Number"
    )
    survival_data["survival_status"] = survival_data["survival_status"].astype(bool)
    survival_data = survival_data.fillna(0)
    return survival_data


def select_surv_feature(
    clinical,
    savefolder="./",
    save_csv="sig_features_bh.csv",
    duration_col="OST",
    event_col="survival_status",
    candidate_cols=None,
    penalizer=0.0,
    min_variance=0.0,
):
    ost = clinical[[event_col, duration_col]].rename(
        columns={event_col: "survival_status", duration_col: "OST"}
    )
    if candidate_cols is None:
        candidate_cols = [
            column
            for column in clinical.columns
            if column not in {duration_col, event_col}
        ]
    clinical_x = (
        clinical[candidate_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    )

    sigPair = {}
    hazardRank = {}
    p_values = []
    all_p = {}
    all_coefs = {}
    all_exp_coefs = {}

    for col in clinical_x.columns:
        if clinical_x[col].var() <= float(min_variance):
            continue
        selected_x = clinical_x[[col]]
        cox = CoxPHFitter(penalizer=float(penalizer))
        df = pd.concat([ost, selected_x], axis=1)
        try:
            cox.fit(
                df,
                duration_col="OST",
                event_col="survival_status",
            )
            schoenfeld = statistics.proportional_hazard_test(
                cox, df, time_transform="rank"
            )
        except Exception as exc:
            print(f"Skipping {col}: {exc}")
            continue
        if schoenfeld.p_value < 0.05:
            print(f"{col} is not proportional hazard")
            continue
        coef = cox.summary["coef"].iloc[0]
        exp_coef = np.exp(coef)
        p = cox.summary["p"].iloc[0]
        all_coefs[col] = coef
        all_exp_coefs[col] = exp_coef
        if p < 0.05:
            sigPair[col] = (p, coef)
            p_values.append(p)
        hazardRank[col] = abs(coef)
        all_p[col] = p

    all_p = pd.Series(all_p, name="p_value")
    save_csv = os.path.join(savefolder, save_csv)
    all_p.to_csv(save_csv, index=True)

    # Perform Benjamini–Hochberg correction on all features
    if len(all_p) > 0:
        pval_array = pd.Series(all_p, name="p_value")
        reject, pvals_corrected, _, _ = multipletests(pval_array, method="fdr_bh")
        coefs = {col: all_coefs[col] for col in pval_array.index}
        exp_coefs = {col: all_exp_coefs[col] for col in pval_array.index}
        corrected_results = pd.DataFrame(
            {
                "feature": pval_array.index,
                "raw_p": pval_array.values,
                "coef": [coefs[c] for c in pval_array.index],
                "exp(coef)": [exp_coefs[c] for c in pval_array.index],
                "bh_corrected_p": pvals_corrected,
                "significant": reject,
            }
        )
        corrected_results.to_csv(save_csv, index=False)
    else:
        corrected_results = pd.DataFrame(
            columns=[
                "feature",
                "raw_p",
                "coef",
                "exp(coef)",
                "bh_corrected_p",
                "significant",
            ]
        )
        corrected_results.to_csv(save_csv, index=False)

    hazardRank = sorted(hazardRank.items(), key=lambda x: x[1])
    hazardRank = hazardRank[-9:-1]
    most_hazard = [pair[0] for pair in hazardRank]

    return most_hazard, sigPair, corrected_results


def KM_age(data, saveFolder="."):
    kmf = KaplanMeierFitter()

    # Define the age groups
    age_groups = {
        "Age_group=middle": "Age = middle",
        "Age_group=old": "Age = old",
        "Age_group=young": "Age = young",
    }

    for age_group, label in age_groups.items():
        if age_group == "Age_group=young":
            age_group_data = data[
                (data["Age_group=middle"] == 0) & (data["Age_group=old"] == 0)
            ]
        else:
            age_group_data = data[data[age_group] == 1]

        # Fit the Kaplan-Meier model
        timeline = np.linspace(data["OST"].min(), data["OST"].max(), 1000)
        kmf.fit(
            age_group_data["OST"],
            event_observed=age_group_data["survival_status"],
            timeline=timeline,
        )

        # Plot the survival function
        kmf.plot_survival_function(ci_show=False, label=label)

    plt.ylim(0, 1)
    plt.ylabel(r"est. probability of survival $\hat{S}(t)$")
    plt.xlabel(r"time $t$")
    plt.title("Survival Curves by Age Group")
    plt.legend(loc="best")
    plt.savefig(saveFolder + "KM_by_Age.png")
    plt.close()


def KM_median(
    data,
    col,
    plotCurve=True,
    saveFolder="./",
    duration_col="OST",
    event_col="survival_status",
) -> tuple[float, float, float, float]:
    # Calculate the median of the specified column
    kmm_median = data[col].median()
    # Split the data based on whether they are greater than or equal to the median
    gtm = data[data[col] > kmm_median]
    stm = data[data[col] < kmm_median]
    # Find rows exactly equal to the median
    equal_to_median = data[data[col] == kmm_median]

    # Assign median-equal rows to the group with fewer elements to balance the split
    lab = {}
    if len(gtm) < len(stm):
        gtm = pd.concat([gtm, equal_to_median])
        lab[0] = "greater than or equal to median"
        lab[1] = "smaller than median"
    else:
        stm = pd.concat([stm, equal_to_median])
        lab[1] = "smaller than or equal to median"
        lab[0] = "greater than median"

    # Initialize KaplanMeierFitter
    kmf = KaplanMeierFitter()
    timeline = np.linspace(data[duration_col].min(), data[duration_col].max(), 1000)
    # Fit data for greater than or equal to median
    kmf.fit(gtm[duration_col], event_observed=gtm[event_col], timeline=timeline)
    half_gtm = kmf.percentile(0.5)
    res = statistics.logrank_test(
        gtm[duration_col], stm[duration_col], gtm[event_col], stm[event_col]
    )

    if plotCurve:
        plt.figure(figsize=(6, 4))
        kmf.plot_survival_function(ci_show=False, label=lab[0])

    # Fit data for less than median
    kmf.fit(stm[duration_col], event_observed=stm[event_col], timeline=timeline)
    half_stm = kmf.percentile(0.5)

    if plotCurve:
        x_position = plt.xlim()[1] * 0.7
        y_position = 0.5
        kmf.plot_survival_function(ci_show=False, label=lab[1])
        plt.text(x_position, y_position, f"log rank test: {res.p_value:.2f}")
        plt.ylim(0, 1)
        plt.ylabel(r"est. probability of survival $\hat{S}(t)$")
        plt.xlabel("time $t$")
        plt.title(f"Survival curve based on {col}")
        plt.legend()
        safe_col = str(col).replace("/", "_").replace(os.sep, "_")
        plt.savefig(f"{saveFolder}{safe_col}.svg", format="svg")
        plt.close()

    return abs(half_gtm - half_stm), half_gtm, half_stm, res.p_value


def KM_median_patients(
    data,
    col,
    plot_curve=True,
    saveFolder="./",
    patient_id_col="Biopsy Number",
    duration_col="OST",
    event_col="survival_status",
) -> tuple[float, float, float, float]:
    # group by patient
    data = data.groupby(by=[patient_id_col]).mean(numeric_only=True)
    data[event_col] = data[event_col].apply(bool)
    data.to_csv(saveFolder + "clinical_with_spatial_patient.csv", index=False)
    return KM_median(
        data=data,
        col=col,
        plotCurve=plot_curve,
        saveFolder=saveFolder,
        duration_col=duration_col,
        event_col=event_col,
    )


def find_most_splitted_curves(
    data,
    cols,
    by_patient=False,
    num=10,
    saveFolder="./",
    save_logrank=True,
    logrank_csv="logrank_results.csv",
    patient_id_col="Biopsy Number",
    duration_col="OST",
    event_col="survival_status",
):
    half_surv_time = {}
    logrank_results = []

    for col in cols:
        try:
            if by_patient:
                surv_time_diff, half_gtm, half_stm, p_value = KM_median_patients(
                    data,
                    col,
                    False,
                    patient_id_col=patient_id_col,
                    duration_col=duration_col,
                    event_col=event_col,
                )
            else:
                surv_time_diff, half_gtm, half_stm, p_value = KM_median(
                    data,
                    col,
                    False,
                    duration_col=duration_col,
                    event_col=event_col,
                )
        except Exception as exc:
            logrank_results.append(
                {
                    "feature": col,
                    "survival_time_diff": np.nan,
                    "median_gtm": np.nan,
                    "median_stm": np.nan,
                    "logrank_p": np.nan,
                    "error": str(exc),
                }
            )
            continue
        half_surv_time[col] = surv_time_diff
        # Collect results for logging
        logrank_results.append(
            {
                "feature": col,
                "survival_time_diff": surv_time_diff,
                "median_gtm": half_gtm,
                "median_stm": half_stm,
                "logrank_p": p_value,
            }
        )

    # Save log-rank test results to CSV
    if save_logrank and logrank_results:
        df_logrank = pd.DataFrame(logrank_results)
        # Apply BH correction to the p-values
        reject, bh_pvals, _, _ = multipletests(df_logrank["logrank_p"], method="fdr_bh")
        df_logrank["logrank_p_bh"] = bh_pvals
        df_logrank["significant"] = reject
        df_logrank = df_logrank.sort_values("logrank_p_bh")
        df_logrank.to_csv(os.path.join(saveFolder, logrank_csv), index=False)

    # Sort and plot top curves
    sorted_surv_diff = sorted(half_surv_time.items(), key=lambda x: x[1], reverse=True)
    plot_count = len(sorted_surv_diff) if (num is None or num <= 0) else num
    for sur_col_time in sorted_surv_diff[:plot_count]:
        try:
            if by_patient:
                KM_median_patients(
                    data,
                    sur_col_time[0],
                    saveFolder=saveFolder,
                    patient_id_col=patient_id_col,
                    duration_col=duration_col,
                    event_col=event_col,
                )
            else:
                KM_median(
                    data,
                    sur_col_time[0],
                    saveFolder=saveFolder,
                    duration_col=duration_col,
                    event_col=event_col,
                )
        except Exception as exc:
            print(f"Skipping KM plot for {sur_col_time[0]}: {exc}")


def run_user_survival_analysis(
    result_folder,
    clinical_csv,
    clinical_id_col,
    duration_col,
    event_col,
    saveFolder="./",
    event_positive_value=None,
    covariate_cols=None,
    image_numbers=None,
    image_to_clinical_csv=None,
    mapping_image_col=None,
    mapping_clinical_col=None,
    norm="min_mid_max",
    km_top_n=10,
    run_km_by_patient=True,
    run_km_by_roi=True,
    penalizer=0.0,
    min_feature_variance=0.0,
):
    os.makedirs(saveFolder, exist_ok=True)
    spatial = read_markcorr_auc_features(
        result_folder=result_folder,
        image_numbers=image_numbers,
        norm=norm,
    )
    clinical = pd.read_csv(clinical_csv)
    mapping = pd.read_csv(image_to_clinical_csv) if image_to_clinical_csv else None
    survival_data = prepare_user_survival_table(
        spatial,
        clinical,
        clinical_id_col=clinical_id_col,
        duration_col=duration_col,
        event_col=event_col,
        event_positive_value=event_positive_value,
        covariate_cols=covariate_cols,
        image_to_clinical=mapping,
        mapping_image_col=mapping_image_col,
        mapping_clinical_col=mapping_clinical_col,
    )
    survival_data.to_csv(
        os.path.join(saveFolder, "clinical_with_spatial.csv"), index=False
    )

    feature_cols = [
        column
        for column in survival_data.columns
        if column not in {"image_id", "clinical_id", "duration", "event"}
    ]
    cols_by_hazard, _, bh_results_df = select_surv_feature(
        survival_data,
        savefolder=saveFolder,
        save_csv="sig_features_bh.csv",
        duration_col="duration",
        event_col="event",
        candidate_cols=feature_cols,
        penalizer=penalizer,
        min_variance=min_feature_variance,
    )

    bh_sig_features = bh_results_df.query("significant")["feature"].tolist()
    for selected_cols, filename in (
        (bh_sig_features, "cox_by_significance.csv"),
        (cols_by_hazard, "cox_by_hazard.csv"),
    ):
        selected_cols = [
            column for column in selected_cols if column in survival_data.columns
        ]
        if not selected_cols:
            continue
        cox_cols = selected_cols + ["duration", "event"]
        try:
            cox = CoxPHFitter(penalizer=float(penalizer))
            cox.fit(survival_data[cox_cols], duration_col="duration", event_col="event")
            cox.summary.to_csv(os.path.join(saveFolder, filename))
        except Exception as exc:
            pd.DataFrame({"error": [str(exc)]}).to_csv(
                os.path.join(saveFolder, filename), index=False
            )

    top_n = int(km_top_n)
    km_features = cols_by_hazard if cols_by_hazard else feature_cols[:top_n]
    if run_km_by_patient:
        patient_folder = os.path.join(saveFolder, "KM_by_patient") + "/"
        os.makedirs(patient_folder, exist_ok=True)
        find_most_splitted_curves(
            survival_data,
            km_features,
            by_patient=True,
            num=top_n,
            saveFolder=patient_folder,
            patient_id_col="clinical_id",
            duration_col="duration",
            event_col="event",
        )
    if run_km_by_roi:
        roi_folder = os.path.join(saveFolder, "KM_by_ROI") + "/"
        os.makedirs(roi_folder, exist_ok=True)
        find_most_splitted_curves(
            survival_data,
            km_features,
            by_patient=False,
            num=top_n,
            saveFolder=roi_folder,
            duration_col="duration",
            event_col="event",
        )
    return survival_data, bh_results_df


def run_survival_analysis(
    intensity=True,
    saveFolder="./",
    km_top_n=10,
    run_km_age=True,
    run_km_by_patient=True,
    run_km_by_roi=True,
    run_example_km_pairs=True,
):
    if intensity:
        result = AMLResult(sizeCorrection=True, intensity=True)
        folder = saveFolder + "intensity/"
    else:
        result = AMLResult(sizeCorrection=True, intensity=False)
        folder = saveFolder + "cellType/"

    marks = list(result.axisName.keys())
    idd_cols = result.get_idd_columns(marks)
    intensity_auc, _ = result.getAUC(norm="min_mid_max", plot=False)
    clinical_data = get_clinical_data()
    clinical_with_spatial = get_spatial_clinical_combined(
        clinical=clinical_data, auc_t=intensity_auc["AML"].transpose()[idd_cols]
    )

    clinical_with_spatial.to_csv(folder + "clinical_with_spatial.csv", index=False)
    cols_by_hazard, sigPair, bh_results_df = select_surv_feature(
        clinical_with_spatial, savefolder=folder, save_csv="sig_features_bh.csv"
    )

    # Use Benjamini–Hochberg significant features
    bh_sig_features = bh_results_df.query("significant")["feature"].tolist()
    selected_cols = bh_sig_features.copy()  # Safe copy for manipulation
    if selected_cols:
        selected_cols.extend(["OST", "survival_status"])
        cox = CoxPHFitter()
        cox.fit(
            clinical_with_spatial[selected_cols],
            duration_col="OST",
            event_col="survival_status",
        )
        results = statistics.proportional_hazard_test(
            cox, clinical_with_spatial[selected_cols], time_transform="rank"
        )
        results.print_summary()
        cox.summary.to_csv(folder + "cox_by_significance.csv")
    else:
        print("No significant features after Benjamini–Hochberg correction.")

    # Use top features by hazard
    selected_cols = cols_by_hazard.copy()
    if selected_cols:
        selected_cols.extend(["OST", "survival_status"])
        cox = CoxPHFitter()
        cox.fit(
            clinical_with_spatial[selected_cols],
            duration_col="OST",
            event_col="survival_status",
        )
        cox.summary.to_csv(folder + "cox_by_hazard.csv")
    else:
        print("No features selected by hazard ranking.")

    # KM curves
    if bool(run_km_age):
        KM_age(clinical_with_spatial, saveFolder=folder)

    top_n = int(km_top_n)
    if bool(run_km_by_patient):
        create_folder(folder + "KM_by_patient/")
        find_most_splitted_curves(
            clinical_with_spatial,
            idd_cols,
            num=top_n,
            by_patient=True,
            saveFolder=folder + "KM_by_patient/",
        )
    if bool(run_km_by_roi):
        create_folder(folder + "KM_by_ROI/")
        find_most_splitted_curves(
            clinical_with_spatial,
            idd_cols,
            num=top_n,
            by_patient=False,
            saveFolder=folder + "KM_by_ROI/",
        )

    if bool(run_example_km_pairs):
        KM_median_patients(
            clinical_with_spatial,
            "Intensity_Erythroids vs. Intensity_CD163",
            plot_curve=True,
            saveFolder=folder,
        )
        KM_median_patients(
            clinical_with_spatial,
            "Intensity_MPO vs. Intensity_MPO",
            plot_curve=True,
            saveFolder=folder,
        )

        KM_median_patients(
            clinical_with_spatial,
            "Intensity_CD31 vs. Intensity_CD163",
            plot_curve=True,
            saveFolder=folder,
        )
