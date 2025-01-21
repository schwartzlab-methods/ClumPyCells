import os

import altair as alt

from ClumPyCells.Analysis.markcorrResult import *
from ClumPyCells.Analysis.metadata import *

IMGF = HOMEDIR + "/Result/images/Melanoma/"


def generateImages():
    metadata = Melanoma_metadata()
    ICI_response = metadata.ICI_response_group
    nonICI_source = metadata.nonICI_source
    ICI_source = metadata.ICI_source
    mel_dataFile = pd.read_csv(metadata.cellTypeFile, sep="\t")
    biopCode = metadata.image2biop
    mel_colInfo = metadata.colInfo

    mel_cellType = MelanomaResult(intensity=False, groups=ICI_response)
    auc, plots = mel_cellType.getAUC()
    plot = MelanomaResult.displayImages(
        plots, [["nonICI", "ICI_nonResponse", "ICI_response"]]
    )
    plot.save(IMGF + "ICIres_nonRes_nonICI.html")

    plot = MelanomaResult.find_diff(
        auc,
        col1="ICI_response",
        col2="ICI_nonResponse",
        method="perm",
        axisName=mel_cellType.axisName,
    )
    plot.save(IMGF + "ICIres_nonRes_diff.html")

    plot = MelanomaResult.find_diff(
        auc, col1="ICI", col2="nonICI", method="perm", axisName=mel_cellType.axisName
    )
    plot.save(IMGF + "ICI_nonICI_diff.html")

    mel_cellType = MelanomaResult(intensity=False, groups=ICI_response)
    auc, plots = mel_cellType.getAUC(takeMean=False)
    plot = MelanomaResult.displayImages(
        plots, [["ICI", "nonICI"], ["ICI_response", "ICI_nonResponse"]]
    )
    plot.save(IMGF + "ICI_nonICI_four.html")

    # mel_cellType = MelanomaResult(intensity=False, groups=nonICI_source)
    # _, plots = mel_cellType.getAUC()
    # plot1 = MelanomaResult.displayImages(plots, [list(nonICI_source.keys())])
    # mel_cellType = MelanomaResult(intensity=False, groups=ICI_source)
    # _, plots = mel_cellType.getAUC()
    # plot2 = MelanomaResult.displayImages(plots, [list(ICI_source.keys())])
    # plot = alt.vconcat(plot1, plot2)
    # plot.save(IMGF + "nonICI_source.html")

    # mel_cellType = MelanomaResult(intensity=False, groups={"all": range(72)})
    # auc, _ = mel_cellType.getAUC()
    # id = auc["all"].stack().idxmax()
    # print(f"max: {id}")
    # id = auc["all"].stack().idxmin()
    # print(f"min: {id}")

    # plotImage(
    #     dataFile=mel_dataFile,
    #     dataFile_colNames=mel_colInfo,
    #     imageName=biopCode[4],
    #     selected_cellTypes=["melano"],
    #     saveName=HOMEDIR + "/Result/images/mel_max.svg",
    # )
    # plotImage(
    #     dataFile=mel_dataFile,
    #     dataFile_colNames=mel_colInfo,
    #     imageName=biopCode[45],
    #     selected_cellTypes=["melano", "B"],
    #     saveName=HOMEDIR + "/Result/images/mel_min.svg",
    # )


def ICIvsNonICI():
    metadata = Melanoma_metadata()
    ICI_response = metadata.ICI_response_group
    mel_cellType = MelanomaResult(intensity=False, groups=ICI_response)

    auc, plots = mel_cellType.getAUC()
    plot = MelanomaResult.displayImages(
        plots, [["nonICI", "ICI_nonResponse", "ICI_response", "ICI"]]
    )
    plot.save(IMGF + "ICIres_nonRes_nonICI.html")
    plot = MelanomaResult.find_diff(
        auc,
        col1="ICI_response",
        col2="ICI_nonResponse",
        method="perm",
        axisName=mel_cellType.axisName,
        takeMean=True,
    )
    plot.save(IMGF + "ICIres_nonRes_diff.html")

    auc, plots = mel_cellType.getAUC(r_range=[120, 170])
    plot = MelanomaResult.find_diff(
        auc,
        col1="ICI_response",
        col2="ICI_nonResponse",
        method="perm",
        axisName=mel_cellType.axisName,
        takeMean=True,
    )
    plot.save(IMGF + "ICIres_nonRes_diff_Rrange.html")

    # plot = MelanomaResult.find_diff(
    #     auc,
    #     col1="ICI",
    #     col2="nonICI",
    #     method="perm",
    #     axisName=mel_cellType.axisName,
    #     takeMean=True,
    # )
    # plot.save(IMGF + "ICI_nonICI_diff.html")

    # mel_cellType = MelanomaResult(intensity=False, groups=ICI_response)
    # auc, plots = mel_cellType.getAUC(takeMean=False)
    # plot = MelanomaResult.displayImages(
    #     plots, [["ICI", "nonICI"], ["ICI_response", "ICI_nonResponse"]]
    # )
    # plot.save(IMGF + "ICI_nonICI_four.html")


def plotMedianAggregation():
    metadata = Melanoma_metadata()
    ICI_response = metadata.ICI_response_group
    mel_cellType = MelanomaResult(intensity=False, groups=ICI_response)
    auc, _ = mel_cellType.getAUC()


def heatmap_all():
    metadata = Melanoma_metadata()
    mel_cellType = MelanomaResult(intensity=False, groups={"all": range(0, 72)})
    auc, plots = mel_cellType.getAUC(r_range=[120, 170])
    plot = MelanomaResult.displayImages(plots, [["all"]])
    plot.save(IMGF + "heatmap_all.html")


def boxplots():
    metadata = Melanoma_metadata()
    ICI_response = metadata.ICI_response_group
    mel_cellType = MelanomaResult(intensity=False, groups=ICI_response)
    auc, plots = mel_cellType.getAUC()

    box_plots = mel_cellType.getBoxPlot(auc, ["ICI_response", "ICI_nonResponse"])
    box_plots.save(IMGF + "ICI_nonICI_boxplot.html")


def tcae_vs_melano():
    metadata = Melanoma_metadata()
    ICI_response = metadata.ICI_response_group
    mel_cellType = MelanomaResult(intensity=False, groups=ICI_response)
    auc, plots = mel_cellType.getAUC(r_range=[220, 225])

    res = auc["ICI_response"].transpose()["Cluster_Tc.ae vs. Cluster_melano"]
    non_res = auc["ICI_nonResponse"].transpose()["Cluster_Tc.ae vs. Cluster_melano"]

    res.dropna(inplace=True)
    non_res.dropna(inplace=True)
    # res.reset_index(drop=True, inplace=True)
    # non_res.reset_index(drop=True, inplace=True)

    res = pd.DataFrame(res)
    non_res = pd.DataFrame(non_res)
    res["group"] = "Responders"
    non_res["group"] = "Non-responders"

    boxData = pd.concat([res, non_res], axis=0)
    boxData.rename(columns={"Cluster_Tc.ae vs. Cluster_melano": "value"}, inplace=True)
    print(boxData)
    plot = (
        alt.Chart(boxData)
        .mark_boxplot()
        .encode(
            x=alt.X("group", title="", axis=alt.Axis(labelAngle=-45)),
            y=alt.Y("value", title="AUC"),
            color=alt.Color("group"),
        )
    )
    plot.save(IMGF + "tcae_vs_melano_boxplot.html")


def plotPointDistribution():
    metadata = Melanoma_metadata()
    mel_dataFile = pd.read_csv(metadata.cellTypeFile, sep="\t")
    biopCode = metadata.image2biop
    mel_colInfo = metadata.colInfo
    mel_cellType = MelanomaResult(intensity=False, groups={"all": range(72)})
    auc, _ = mel_cellType.getAUC()
    plotImage(
        dataFile=mel_dataFile,
        dataFile_colNames=mel_colInfo,
        imageName=biopCode[2],
        selected_cellTypes=["melano", "Tc.ae"],
        saveName=HOMEDIR + "/Result/images/melanoVsTcae.svg",
    )


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
    auc_obs, _ = kmmResult.getAUC(plot=False)
    permP = {}
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
            heatmapData["auc"] = auc_obs_mean
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

            plot = (
                alt.Chart(heatmapData, title=group)
                .mark_rect()
                .encode(
                    x=alt.X("from", axis=alt.Axis(labelAngle=-45)).sort(
                        list(kmmResult.axisName.values())
                    ),
                    y=alt.Y("to").sort(list(kmmResult.axisName.values())),
                    color=alt.Color(
                        "auc",
                        scale=alt.Scale(
                            scheme="redblue", domainMid=0, reverse=True, type="symlog"
                        ),
                    ),
                )
            )

            sig = (
                alt.Chart(inSig)
                .mark_text(filled=True, text="*")
                .encode(
                    x=alt.X("from", axis=alt.Axis(labelAngle=-45)).sort(
                        list(kmmResult.axisName.values())
                    ),
                    y=alt.Y("to").sort(list(kmmResult.axisName.values())),
                    color=alt.value("black"),
                )
            )

            text = (
                alt.Chart(heatmapData)
                .mark_text()
                .encode(
                    x=alt.X("from", axis=alt.Axis(labelAngle=-45)).sort(
                        list(kmmResult.axisName.values())
                    ),
                    y=alt.Y("to").sort(list(kmmResult.axisName.values())),
                    color=alt.value("black"),
                    text="p",
                )
            )
            perm_plot = plot + sig
            perm_plot.save(f"{permFolder}perm_{group}.html")
            perm_text = plot + text
            perm_text.save(f"{permFolder}perm_{group}_text.html")

    return permP


def mel_permuation():
    metadata = Melanoma_metadata()
    mel_axisName = metadata.axisName
    ICIgroups = metadata.ICI_response_group
    permFolder = HOMEDIR + "Result/Melanoma/Permutation/"
    kmmResFolder = HOMEDIR + "Result/Melanoma/Melanoma_cellType_more/"
    # permutation_result(
    #     permFolder=permFolder,
    #     perm_num=50,
    #     groups=ICIgroups,
    #     axisName=mel_axisName,
    # )
    get_permuation_p(
        kmmResFolder=kmmResFolder,
        groups=ICIgroups,
        axisName=mel_axisName,
        permFolder=permFolder,
        permNum=51,
    )


tcae_vs_melano()
plotPointDistribution()
