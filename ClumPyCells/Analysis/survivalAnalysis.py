import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter, statistics
from sksurv.preprocessing import OneHotEncoder
from statsmodels.stats.multitest import multipletests

from .markcorrResult import *
from .metadata import *


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


def select_surv_feature(clinical, savefolder="./", save_csv="sig_features_bh.csv"):
    ost = clinical[["survival_status", "OST"]]
    clinical_x = clinical.drop(["OST", "survival_status"], axis=1)

    sigPair = {}
    hazardRank = {}
    p_values = []
    all_p = {}
    all_coefs = {}
    all_exp_coefs = {}

    for col in clinical_x.columns:
        selected_x = clinical_x[[col]]
        cox = CoxPHFitter()
        df = pd.concat([ost, selected_x], axis=1)
        cox.fit(
            df,
            duration_col="OST",
            event_col="survival_status",
        )
        schoenfeld = statistics.proportional_hazard_test(cox, df, time_transform="rank")
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


def KM_median(data, col, plotCurve=True, saveFolder="./"):
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
    timeline = np.linspace(data["OST"].min(), data["OST"].max(), 1000)
    # Fit data for greater than or equal to median
    kmf.fit(gtm["OST"], event_observed=gtm["survival_status"], timeline=timeline)
    half_gtm = kmf.percentile(0.5)
    res = statistics.logrank_test(
        gtm["OST"], stm["OST"], gtm["survival_status"], stm["survival_status"]
    )

    if plotCurve:
        plt.figure(figsize=(6, 4))
        kmf.plot_survival_function(ci_show=False, label=lab[0])

    # Fit data for less than median
    kmf.fit(stm["OST"], event_observed=stm["survival_status"], timeline=timeline)
    half_stm = kmf.percentile(0.5)

    if plotCurve:
        x_position = plt.xlim()[1] * 0.7
        y_position = 0.5
        kmf.plot_survival_function(ci_show=False, label=lab[1])
        plt.text(x_position, y_position, f"log rank test: {res.p_value}")
        plt.ylim(0, 1)
        plt.ylabel(r"est. probability of survival $\hat{S}(t)$")
        plt.xlabel("time $t$")
        plt.title(f"Survival curve based on {col}")
        plt.legend()
        plt.savefig(f"{saveFolder}{col}.svg", format="svg")
        plt.close()

    return abs(half_gtm - half_stm), half_gtm, half_stm, res.p_value


def KM_median_patients(data, col, plot_curve=True, saveFolder="./"):
    # group by patient
    data = data.groupby(by=["Biopsy Number"]).mean()
    data["survival_status"] = data["survival_status"].apply(bool)
    return KM_median(data=data, col=col, plotCurve=plot_curve, saveFolder=saveFolder)


def find_most_splitted_curves(
    data,
    cols,
    by_patient=False,
    num=10,
    saveFolder="./",
    save_logrank=True,
    logrank_csv="logrank_results.csv",
):
    func = KM_median_patients if by_patient else KM_median
    half_surv_time = {}
    logrank_results = []

    for col in cols:
        surv_time_diff, half_gtm, half_stm, p_value = func(data, col, False)
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
        func(data, sur_col_time[0], saveFolder=saveFolder)


def run_survival_analysis(intensity=True, saveFolder="./"):
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
    KM_age(clinical_with_spatial, saveFolder=folder)
    create_folder(folder + "KM_by_patient/")
    create_folder(folder + "KM_by_ROI/")
    find_most_splitted_curves(
        clinical_with_spatial,
        idd_cols,
        num=10,
        by_patient=True,
        saveFolder=folder + "KM_by_patient/",
    )
    find_most_splitted_curves(
        clinical_with_spatial,
        idd_cols,
        num=10,
        by_patient=False,
        saveFolder=folder + "KM_by_ROI/",
    )
