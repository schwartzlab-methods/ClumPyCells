import altair as alt
from ClumPyCell.Analysis.markcorrResult import *
from ClumPyCell.Analysis.metadata import *
import os


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
    plot.save(HOMEDIR + "/Result/images/ICIres_nonRes_nonICI.html")

    plot = MelanomaResult.find_diff(
        auc,
        col1="ICI_response",
        col2="ICI_nonResponse",
        method="perm",
        axisName=mel_cellType.axisName,
    )
    plot.save(HOMEDIR + "/Result/images/ICIres_nonRes_diff.html")

    plot = MelanomaResult.find_diff(
        auc, col1="ICI", col2="nonICI", method="perm", axisName=mel_cellType.axisName
    )
    plot.save(HOMEDIR + "/Result/images/ICI_nonICI_diff.html")

    mel_cellType = MelanomaResult(intensity=False, groups=ICI_response)
    auc, plots = mel_cellType.getAUC(takeMean=False)
    plot = MelanomaResult.displayImages(
        plots, [["ICI", "nonICI"], ["ICI_response", "ICI_nonResponse"]]
    )
    plot.save(HOMEDIR + "/Result/images/ICI_nonICI_four.html")

    # mel_cellType = MelanomaResult(intensity=False, groups=nonICI_source)
    # _, plots = mel_cellType.getAUC()
    # plot1 = MelanomaResult.displayImages(plots, [list(nonICI_source.keys())])
    # mel_cellType = MelanomaResult(intensity=False, groups=ICI_source)
    # _, plots = mel_cellType.getAUC()
    # plot2 = MelanomaResult.displayImages(plots, [list(ICI_source.keys())])
    # plot = alt.vconcat(plot1, plot2)
    # plot.save(HOMEDIR + "/Result/images/nonICI_source.html")

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
    auc, plots = mel_cellType.getAUC(r_range=[120, 170])
    plot = MelanomaResult.displayImages(
        plots, [["nonICI", "ICI_nonResponse", "ICI_response"]]
    )
    plot.save(HOMEDIR + "/Result/images/ICIres_nonRes_nonICI.html")

    plot = MelanomaResult.find_diff(
        auc,
        col1="ICI_response",
        col2="ICI_nonResponse",
        method="perm",
        axisName=mel_cellType.axisName,
        takeMean=True,
    )
    plot.save(HOMEDIR + "/Result/images/ICIres_nonRes_diff.html")

    # plot = MelanomaResult.find_diff(
    #     auc,
    #     col1="ICI",
    #     col2="nonICI",
    #     method="perm",
    #     axisName=mel_cellType.axisName,
    #     takeMean=True,
    # )
    # plot.save(HOMEDIR + "/Result/images/ICI_nonICI_diff.html")

    # mel_cellType = MelanomaResult(intensity=False, groups=ICI_response)
    # auc, plots = mel_cellType.getAUC(takeMean=False)
    # plot = MelanomaResult.displayImages(
    #     plots, [["ICI", "nonICI"], ["ICI_response", "ICI_nonResponse"]]
    # )
    # plot.save(HOMEDIR + "/Result/images/ICI_nonICI_four.html")


def boxplots():
    metadata = Melanoma_metadata()
    ICI_response = metadata.ICI_response_group
    mel_cellType = MelanomaResult(intensity=False, groups=ICI_response)
    auc, plots = mel_cellType.getAUC()

    box_plots = mel_cellType.getBoxPlot(auc, ["ICI_response", "ICI_nonResponse"])
    box_plots.save(HOMEDIR + "/Result/images/ICI_nonICI_boxplot.html")


def plotimage():
    metadata = Melanoma_metadata()
    mel_dataFile = pd.read_csv(metadata.cellTypeFile, sep="\t")
    biopCode = metadata.image2biop
    mel_colInfo = metadata.colInfo
    mel_cellType = MelanomaResult(intensity=False, groups={"all": range(72)})
    auc, _ = mel_cellType.getAUC()
    plotImage(
        dataFile=mel_dataFile,
        dataFile_colNames=mel_colInfo,
        imageName=biopCode[0],
        selected_cellTypes=["melano", "Tc.ae"],
        saveName=HOMEDIR + "/Result/images/temp.svg",
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
        print(perm.columns)
        count_condition_met = pd.Series(index=perm.columns, dtype=int)
        for column in perm.columns:
            print(column)
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
                alt.Chart(heatmapData, title="AML")
                .mark_rect()
                .encode(
                    x=alt.X("from").sort(list(kmmResult.axisName.values())),
                    y=alt.Y("to").sort(list(kmmResult.axisName.values())),
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
                    x=alt.X("from").sort(list(kmmResult.axisName.values())),
                    y=alt.Y("to").sort(list(kmmResult.axisName.values())),
                    color=alt.value("black"),
                )
            )

            text = (
                alt.Chart(heatmapData)
                .mark_text()
                .encode(
                    x=alt.X("from").sort(list(kmmResult.axisName.values())),
                    y=alt.Y("to").sort(list(kmmResult.axisName.values())),
                    color=alt.value("black"),
                    text="p",
                )
            )
            aml_plot = plot + sig
            aml_plot.save(f"{permFolder}perm_{group}.html")
            aml_text = plot + text
            aml_text.save(f"{permFolder}perm_{group}_text.html")
            break

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
    get_permuation_p(kmmResFolder, ICIgroups, mel_axisName, permFolder, 1)


mel_permuation()
