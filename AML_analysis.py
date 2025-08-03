import sys

import altair as alt

from ClumPyCells.Analysis.decisionTree import *
from ClumPyCells.Analysis.markcorrResult import *
from ClumPyCells.Analysis.metadata import *
from ClumPyCells.Analysis.survivalAnalysis import *

IMAGEFOLDER = HOMEDIR + "Result/images/AML/"
sys.path.append(HOMEDIR + "altairThemes.py")

if True:  # In order to bypass isort when saving
    import altairThemes

# register the custom theme under a chosen name
alt.themes.register("publishTheme", altairThemes.publishTheme)

# enable the newly registered theme
alt.themes.enable("publishTheme")


def plot_images_cellType():
    aml = AMLResult(sizeCorrection=True, intensity=False)
    aml_metadata = AML_metadata()
    auc, plts = aml.getAUC()
    auc_plts = AMLResult.displayImages(plts, [["NBM", "AML"]])
    auc_plts.save(IMAGEFOLDER + "NBMvsAML_cellType.html")
    diff_plts = AMLResult.find_diff(
        auc=auc, col1="AML", col2="NBM", axisName=aml.axisName
    )
    diff_plts.save(IMAGEFOLDER + "AML_NBM_diff_cellType.html")


def plot_images(AMLvsNBM=True, MinMax=True, BlastPercentage=True):
    aml = AMLResult(sizeCorrection=True, intensity=True)
    aml_metadata = AML_metadata()
    auc, plts = aml.getAUC()
    if AMLvsNBM:
        auc_plts = AMLResult.displayImages(plts, [["NBM", "AML"]])
        auc_plts.save(IMAGEFOLDER + "NBMvsAML.html")
        diff_plts = AMLResult.find_diff(
            auc=auc, col1="AML", col2="NBM", axisName=aml.axisName
        )
        diff_plts.save(IMAGEFOLDER + "AML_NBM_diff.html")
    if MinMax:
        # Find the image with the highest and lowest AUC
        id_min = auc["AML"].stack().idxmin()
        print(f"min: {id_min}")
        id_max = auc["AML"].stack().idxmax()
        print(f"max: {id_max}")

        min_imageNum = int(id_min[1].split("_")[1])
        max_imageNum = int(id_max[1].split("_")[1])
        min_types = id_min[0].split(" vs. ")
        max_types = id_max[0].split(" vs. ")

        print(min_types, max_types, min_imageNum, max_imageNum)
        if min_imageNum > 36:
            dataFile = aml_metadata.nbm_file
        else:
            dataFile = aml_metadata.aml_file

        dataFile = pd.read_csv(dataFile)
        plotImage(
            dataFile,
            aml_metadata.colInfo,
            imageName=min_imageNum,
            area_colName="Area",
            selected_cellTypes=min_types,
            saveName=IMAGEFOLDER + "AML_min.svg",
        )
        if max_imageNum > 36:
            dataFile = aml_metadata.nbm_file
        else:
            dataFile = aml_metadata.aml_file
        dataFile = pd.read_csv(dataFile)
        plotImage(
            dataFile,
            aml_metadata.colInfo,
            imageName=max_imageNum,
            area_colName="Area",
            selected_cellTypes=["Bcells", "Bcells"],
            saveName=IMAGEFOLDER + "AML_max.svg",
        )
    if BlastPercentage:
        BlastPercentage_groups = aml_metadata.get_blast_percentage_split()
        blast_group_result = AMLResult(
            sizeCorrection=True, groups=BlastPercentage_groups, intensity=True
        )
        auc, plots = blast_group_result.getAUC()
        auc_plts = blast_group_result.displayImages(plots, [["aml_high", "aml_low"]])
        auc_plts.save(IMAGEFOLDER + "blastpercentage.html")
        diff_plts = AMLResult.find_diff(
            auc=auc,
            col1="aml_high",
            col2="aml_low",
            method="perm",
            axisName=blast_group_result.axisName,
        )
        diff_plts.save(IMAGEFOLDER + "blastPercentage_diff.html")


def permutation_result(permFolder, perm_num, groups, axisName):
    permResultMaster = {}
    for group in groups:
        permResultMaster[group] = pd.DataFrame()
    for i in range(perm_num):
        # Get result from each permutation
        permResFolder = os.path.join(permFolder + f"/perm_{i}/")
        perm_result = MarkcorrResult(
            groups=groups, resultFolder=permResFolder, axisName=axisName
        )
        auc, _ = perm_result.getAUC(plot=False)
        for group in groups:
            permResultMaster[group][f"perm_{i}"] = auc[group].mean(axis=1)
    for group in groups:
        permResultMaster[group].to_csv(f"{permFolder}{group}.csv")


def get_permuation_p(kmmResFolder, groups, axisName, permFolder, permNum, plot=True):
    kmmResult = MarkcorrResult(
        groups=groups, resultFolder=kmmResFolder, axisName=axisName
    )
    auc_obs, obs_plots = kmmResult.getAUC()
    permP = {}
    plots = {}
    for group in groups:

        auc_obs_mean = auc_obs[group].mean(axis=1)
        perm = pd.read_csv(f"{permFolder}{group}.csv", index_col=0)
        auc_obs_mean = auc_obs_mean.transpose()
        perm = perm.transpose()
        auc_obs_mean_reordered = auc_obs_mean[perm.columns]
        count_condition_met = pd.Series(index=perm.columns, dtype=int)
        for column in perm.columns:
            if auc_obs_mean_reordered[column] > 0:
                count_condition_met[column] = (
                    perm[column] > auc_obs_mean_reordered[column]
                ).sum()
            else:
                count_condition_met[column] = (
                    perm[column] < auc_obs_mean_reordered[column]
                ).sum()
        permP[group] = count_condition_met / permNum

        if plot:
            heatmapData = pd.DataFrame()
            heatmapData["p"] = permP[group]
            x, y = kmmResult.get_XY(heatmapData, kmmResult.axisName)
            heatmapData["from"] = x
            heatmapData["to"] = y

            all_combinations = pd.MultiIndex.from_product(
                [list(kmmResult.axisName.values()), list(kmmResult.axisName.values())],
                names=["from", "to"],
            ).to_frame(index=False)
            heatmapData = pd.merge(
                all_combinations, heatmapData, on=["from", "to"], how="outer"
            )
            heatmapData["p"] = heatmapData["p"].where(pd.notnull(heatmapData["p"]), 0)

            inSig = heatmapData.loc[heatmapData["p"] < 0.05]

            plot = obs_plots[group]

            sig = (
                alt.Chart(inSig)
                .mark_text(filled=True, text="*", dy=3)
                .encode(
                    x=alt.X(
                        "from",
                        axis=alt.Axis(
                            labelAngle=-45, labelFontSize=11, titleFontSize=11
                        ),
                    ).sort(list(kmmResult.axisName.values())),
                    y=alt.Y(
                        "to", axis=alt.Axis(labelFontSize=11, titleFontSize=11)
                    ).sort(list(kmmResult.axisName.values())),
                    color=alt.value("black"),
                )
            )

            text = (
                alt.Chart(heatmapData)
                .mark_text()
                .encode(
                    x=alt.X(
                        "from",
                        axis=alt.Axis(
                            labelAngle=-45, labelFontSize=11, titleFontSize=11
                        ),
                    ).sort(list(kmmResult.axisName.values())),
                    y=alt.Y(
                        "to", axis=alt.Axis(labelFontSize=11, titleFontSize=11)
                    ).sort(list(kmmResult.axisName.values())),
                    color=alt.value("black"),
                    text="p",
                )
            )
            perm_plot = plot + sig
            perm_plot.save(f"{IMAGEFOLDER}perm_{group}.html")
            perm_text = plot + text
            perm_text.save(f"{IMAGEFOLDER}perm_{group}_text.html")

            plots[group] = perm_plot
    if plot:
        # Combine all group plots into a single HTML file
        plots = {g: plots[g].properties(height=170, width=170) for g in groups}
        combined = alt.hconcat(*[plots[g] for g in groups])
        combined = alt.hconcat(plots["NBM"], plots["AML"])
        combined.save(f"{IMAGEFOLDER}perm_combined.html")
    return permP


def perm_aml():
    metadata = AML_metadata()
    aml_axisName = metadata.axisName
    permFolder = HOMEDIR + "Result/AML/Permutation/"
    kmmResFolder = HOMEDIR + "Result/AML/intensity_withSize/"
    # permutation_result(
    #     permFolder=permFolder,
    #     perm_num=50,
    #     groups=ICIgroups,
    #     axisName=mel_axisName,
    # )
    get_permuation_p(
        kmmResFolder=kmmResFolder,
        groups={"AML": range(36), "NBM": range(36, 51)},
        axisName=aml_axisName,
        permFolder=permFolder,
        permNum=100,
    )


def plot_AML(imageNum):
    aml_meta = AML_metadata()
    col_info = aml_meta.colInfo
    i = imageNum
    if imageNum > 36:
        dataFile = pd.read_csv(aml_meta.nbm_file)
        i -= 36  # Adjust image number for NBM images
    else:
        dataFile = pd.read_csv(aml_meta.aml_file)
    plotImage(
        dataFile=dataFile,
        dataFile_colNames=col_info,
        imageName=i,
        saveName=IMAGEFOLDER + f"cell_images/image_{imageNum - 1}_image.svg",
    )


def surv_result():
    run_survival_analysis(intensity=True, saveFolder=HOMEDIR + "/Result/AML/survival/")


def plot_decision_tree():
    decision_tree(intensity=True, saveFolder=IMAGEFOLDER)


def plot_each_auc():
    imageGroups = {f"image_{i}": [i] for i in range(51)}
    imageResult = AMLResult(groups=imageGroups, intensity=True, sizeCorrection=True)
    auc, plots = imageResult.getAUC(plot=True)
    labels = [f"image_{i}" for i in range(51)]
    # First 3 rows: groups of 12
    rows_1_3 = [labels[i : i + 12] for i in range(0, 36, 12)]
    # Last 3 rows: groups of 5
    rows_4_6 = [labels[i : i + 5] for i in range(36, 51, 5)]
    # Combine all rows
    image_array = rows_1_3 + rows_4_6
    imageResult.displayImages(plots, image_array).save(
        IMAGEFOLDER + "AML_auc_each_plot_fixed.html"
    )


surv_result()
