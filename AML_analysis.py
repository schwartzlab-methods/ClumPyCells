import altair as alt
from ClumPyCell.Analysis.markcorrResult import *
from ClumPyCell.Analysis.metadata import *


def plot_images(AMLvsNBM=True, MinMax=True, BlastPercentage=True):
    aml = AMLResult(sizeCorrection=True, intensity=True)
    aml_metadata = AML_metadata()
    if AMLvsNBM:
        auc, plts = aml.getAUC()
        auc_plts = AMLResult.displayImages(plts, [["NBM", "AML"]])
        auc_plts.save(HOMEDIR + "/Result/images/NBMvsAML.html")
        diff_plts = AMLResult.find_diff(auc, "AML", "NBM", aml.axisName)
        diff_plts.save(HOMEDIR + "/Result/images/AML_NBM_diff.html")
    if MinMax:
        id = auc["AML"].stack().idxmin()
        print(id)
        imageNumber = 1
        if imageNumber > 36:
            dataFile = aml_metadata.nbm_file
        else:
            dataFile = aml_metadata.aml_file
        dataFile = pd.read_csv(dataFile)
        plotImage(
            dataFile,
            aml_metadata.colInfo,
            imageName=imageNumber,
            area_colName="Area",
            selected_cellTypes=["Tcells", "Adipocytes"],
            saveName=HOMEDIR + "/Result/images/AML_min.svg",
        )
        imageNumber = 34
        if imageNumber > 36:
            dataFile = aml_metadata.nbm_file
        else:
            dataFile = aml_metadata.aml_file
        dataFile = pd.read_csv(dataFile)
        plotImage(
            dataFile,
            aml_metadata.colInfo,
            imageName=imageNumber,
            area_colName="Area",
            selected_cellTypes=["Bcells", "Adipocytes"],
            saveName=HOMEDIR + "/Result/images/AML_max.svg",
        )
    if BlastPercentage:
        BlastPercentage_groups = aml_metadata.get_blast_percentage_split()
        blast_group_result = AMLResult(
            sizeCorrection=True, groups=BlastPercentage_groups, intensity=True
        )
        auc, plots = blast_group_result.getAUC()
        auc_plts = blast_group_result.displayImages(plots, [["aml_high", "aml_low"]])
        auc_plts.save(HOMEDIR + "/Result/images/blastpercentage.html")
        diff_plts = AMLResult.find_diff(
            auc=auc,
            col1="aml_high",
            col2="aml_low",
            method="perm",
            axisName=blast_group_result.axisName,
        )
        diff_plts.save(HOMEDIR + "/Result/images/blastPercentage_diff.html")


def permutation_result():
    aml = AMLResult(sizeCorrection=True, intensity=True)
    aml_permutation = pd.DataFrame()
    nbm_permutation = pd.DataFrame()
    for i in range(100):
        permutation = MarkcorrResult(
            {"AML": range(36), "NBM": range(36, 51)},
            HOMEDIR + f"/Result/AML/Permutation/perm_{i}/",
        )
        auc, _ = permutation.getAUC(plot=False, axisName=aml.axisName)
        aml_permutation[f"perm_{i}"] = auc["AML"].mean(axis=1)
        print(aml_permutation)
        nbm_permutation[f"perm_{i}"] = auc["NBM"].mean(axis=1)
    aml_permutation.to_csv(
        "/cluster/home/t114231uhn/AML_Public/Result/AML/Permutation/perm_aml.csv"
    )
    nbm_permutation.to_csv(
        "/cluster/home/t114231uhn/AML_Public/Result/AML/Permutation/perm_nbm.csv"
    )


def permutation_p():
    aml = AMLResult(sizeCorrection=True, intensity=False)
    auc_obs, _ = aml.getAUC(plot=False)
    auc_obs_mean = auc_obs["AML"].mean(axis=1)
    auc_obs_mean = auc_obs_mean.transpose()
    perm_aml = pd.read_csv(
        "/cluster/home/t114231uhn/AML_Public/Result/AML/Permutation/perm_aml.csv"
    ).drop(["Unnamed: 0"], axis=1)
    perm_nbm = pd.read_csv(
        "/cluster/home/t114231uhn/AML_Public/Result/AML/Permutation/perm_nbm.csv"
    ).drop(["Unnamed: 0"], axis=1)
    auc_obs_mean_reordered = auc_obs_mean[perm_aml.columns]
    count_condition_met = pd.Series(index=perm_aml.columns, dtype=int)

    for column in perm_aml.columns:
        if auc_obs_mean_reordered[column] > 0:
            count_condition_met[column] = (
                perm_aml[column] > auc_obs_mean_reordered[column]
            ).sum()
        else:
            count_condition_met[column] = (
                perm_aml[column] < auc_obs_mean_reordered[column]
            ).sum()
    p_aml = count_condition_met / 100
    heatmapData = pd.DataFrame()
    heatmapData["auc"] = auc_obs_mean
    heatmapData["p"] = p_aml
    x, y = aml.get_XY(heatmapData, aml.axisName)
    heatmapData["from"] = x
    heatmapData["to"] = y

    all_combinations = pd.MultiIndex.from_product(
        [list(aml.axisName.values()), list(aml.axisName.values())],
        names=["from", "to"],
    ).to_frame(index=False)
    heatmapData = pd.merge(
        all_combinations, heatmapData, on=["from", "to"], how="outer"
    )
    heatmapData["p"] = heatmapData["p"].where(pd.notnull(heatmapData["p"]), 0)

    print(heatmapData)

    inSig = heatmapData.loc[heatmapData["p"] < 0.05]

    plot = (
        alt.Chart(heatmapData, title="AML")
        .mark_rect()
        .encode(
            x=alt.X("from").sort(list(aml.axisName.values())),
            y=alt.Y("to").sort(list(aml.axisName.values())),
            color=alt.Color(
                "auc",
                scale=alt.Scale(
                    scheme="redblue", domainMid=0, reverse=True, type="symlog"
                ),
            ),
        )
        .properties(height=500, width=500)
    )

    sig = (
        alt.Chart(inSig)
        .mark_point(size=500, filled=False, shape="triangle-up")
        .encode(
            x=alt.X("from").sort(list(aml.axisName.values())),
            y=alt.Y("to").sort(list(aml.axisName.values())),
            color=alt.value("black"),
        )
    )

    text = (
        alt.Chart(heatmapData)
        .mark_text()
        .encode(
            x=alt.X("from").sort(list(aml.axisName.values())),
            y=alt.Y("to").sort(list(aml.axisName.values())),
            color=alt.value("black"),
            text="p",
        )
    )
    aml_plot = plot + sig
    aml_plot.save(HOMEDIR + f"/Result/images/permutationAML.html")
    aml_text = plot + text
    aml_text.save(HOMEDIR + f"/Result/images/permText_aml.html")

    auc_obs_mean = auc_obs["NBM"].mean(axis=1)
    auc_obs_mean = auc_obs_mean.transpose()
    auc_obs_mean_reordered_nbm = auc_obs_mean[perm_nbm.columns]
    count_condition_met_nbm = pd.Series(index=perm_nbm.columns, dtype=int)

    for column in perm_nbm.columns:
        if auc_obs_mean_reordered_nbm[column] > 0:
            count_condition_met_nbm[column] = (
                perm_nbm[column] > auc_obs_mean_reordered_nbm[column]
            ).sum()
        else:
            count_condition_met_nbm[column] = (
                perm_nbm[column] < auc_obs_mean_reordered_nbm[column]
            ).sum()

    # Calculate p-value for perm_nbm
    p_nbm = count_condition_met_nbm / 100

    heatmapData = pd.DataFrame()
    heatmapData["auc"] = auc_obs_mean
    heatmapData["p"] = p_nbm
    # heatmapData["p"] = multi.fdrcorrection(heatmapData["p"])[1]
    x, y = aml.get_XY(heatmapData, aml.axisName)
    heatmapData["from"] = x
    heatmapData["to"] = y
    inSig = heatmapData.loc[heatmapData["p"] < 0.05]
    plot = (
        alt.Chart(heatmapData, title="NBM")
        .mark_rect()
        .encode(
            x=alt.X("from").sort(list(aml.axisName.values())),
            y=alt.Y("to").sort(list(aml.axisName.values())),
            color=alt.Color(
                "auc",
                scale=alt.Scale(
                    scheme="redblue", domainMid=0, reverse=True, type="symlog"
                ),
            ),
        )
        .properties(height=500, width=500)
    )
    sig = (
        alt.Chart(inSig)
        .mark_point(size=500, filled=False, shape="triangle-up")
        .encode(
            x=alt.X("from").sort(list(aml.axisName.values())),
            y=alt.Y("to").sort(list(aml.axisName.values())),
            color=alt.value("black"),
        )
    )

    text = (
        alt.Chart(heatmapData)
        .mark_text()
        .encode(
            x=alt.X("from").sort(list(aml.axisName.values())),
            y=alt.Y("to").sort(list(aml.axisName.values())),
            color=alt.value("black"),
            text="p",
        )
    )
    print(p_aml)
    nbm_plot = plot + sig
    nbm_plot.save(HOMEDIR + f"/Result/images/permutationNBM.html")
    nbm_text = plot + text
    nbm_text.save(HOMEDIR + f"/Result/images/permText_nbm.html")

    (aml_plot | nbm_plot).save(HOMEDIR + f"/Result/images/permutation_combined.html")
    (aml_text | nbm_text).save(HOMEDIR + f"/Result/images/perm_text_combined.html")


plot_images(False, False, True)
