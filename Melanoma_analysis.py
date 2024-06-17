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


ICIvsNonICI()
