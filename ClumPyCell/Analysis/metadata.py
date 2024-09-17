import json
import logging
import os
import shutil
import sys

import altair as alt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

with open("config.json", "r") as file:
    config = json.load(file)
HOMEDIR = config["HOMEDIR"]

sys.path.append(HOMEDIR + "ClumPyCell/Analysis/altairThemes.py")

if True:  # In order to bypass isort when saving
    import altairThemes

# register the custom theme under a chosen name
alt.themes.register("publishTheme", altairThemes.publishTheme)

# enable the newly registered theme
alt.themes.enable("publishTheme")


class AML_metadata:
    def __init__(self) -> None:
        self.nbm_file = f"{HOMEDIR}Data/output/AML.csv"
        self.aml_file = f"{HOMEDIR}Data/output/Normal.csv"
        self.imageSize_x = 1056
        self.imageSize_y = 642
        self.cellTypes_antibody = {
            "Adipocytes": "Perilipin",
            "B_cells": "CD20",
            "HSC": "CD34",
            "Erythroids": "ECAD",
            "Ki67": "Ki67",
            "MSC": "NGFR",
            "Macrophages": "CD163",
            "Megakaryocytes": "CD31",
            "Monocytes": "CD68",
            "Myeloids": "MPO",
            "T_cells": "CD3",
        }

        self.colInfo = [
            "ImageNumber",
            "x",
            "y",
            "CellType",
            self.imageSize_x,
            self.imageSize_y,
        ]

    def get_blast_percentage_split(self, splitPercentage=75):
        blastData = pd.read_csv(
            os.path.join(HOMEDIR, "Data/AML_blast_stroma_database.csv")
        )
        blastData["Blast_percentage_group"] = pd.cut(
            blastData["Blast_percentage"],
            bins=[0, splitPercentage, 100],
            labels=["low", "high"],
        )
        blastData["Biopsy_number"] = blastData["Biopsy_code"].apply(
            lambda x: x.split("-")[1]
        )
        image2biop = imageNum_to_biopNum()
        aml_high, aml_low = [], []
        print(len(image2biop))
        for i, biopCode in enumerate(image2biop):
            biop = blastData.loc[blastData["Biopsy_number"] == biopCode]
            try:
                group = biop["Blast_percentage_group"].values[0]
            except:
                logging.warning("biopsy code data is not available and will be skipped")
                continue
            if group == "high":
                aml_high.append(i)
            else:
                aml_low.append(i)
        return {
            "aml_high": aml_high,
            "aml_low": aml_low,
        }


class Melanoma_metadata:
    def __init__(self) -> None:
        # Melanoma cell location file based on Cell type
        self.cellTypeFile = HOMEDIR + "Data/ST4_cell_data.txt"
        self.colInfo = ["sample.id", "coord.x", "coord.y", "Cluster", 1200, 1200]
        self.axisName = {
            "Cluster_Tc.ae": "Tc.ae",
            "Cluster_Tc.naive": "Tc.naive",
            "Cluster_Th.naive": "Th.naive",
            "Cluster_B": "B",
            "Cluster_Th.ae": "Th.ae",
            "Cluster_Treg": "Treg",
            "Cluster_CD31": "CD31",
            "Cluster_melano": "melano",
            "Cluster_macro.mono": "macro.mono",
            "Cluster_others": "others",
        }
        mel_cellType = pd.read_csv(self.cellTypeFile, sep="\t")

        # image number to biopsy code reference dict, use index() to get the reverse
        self.image2biop = list(mel_cellType["sample.id"].unique())
        imageInfo = pd.read_csv(
            f"{HOMEDIR}Data/ST1_sample_info.txt",
            sep="\t",
            usecols=["Sample_ID", "Cohort", "Response", "Tissue_Source"],
        )

        # frist group
        ICI_response_group = {}
        ICI = imageInfo[imageInfo["Cohort"] == "ICI"]
        ICI_response_group["ICI"] = sorted(
            [self.image2biop.index(i) for i in ICI["Sample_ID"]]
        )

        nonICI = imageInfo[imageInfo["Cohort"] == "non-ICI"]
        ICI_response_group["nonICI"] = sorted(
            [self.image2biop.index(i) for i in nonICI["Sample_ID"]]
        )

        ICI_response = imageInfo[imageInfo["Response"] == "Yes"]
        ICI_response_group["ICI_response"] = sorted(
            [self.image2biop.index(i) for i in ICI_response["Sample_ID"]]
        )

        ICI_nonResponse = imageInfo[imageInfo["Response"] == "No"]
        ICI_response_group["ICI_nonResponse"] = sorted(
            [self.image2biop.index(i) for i in ICI_nonResponse["Sample_ID"]]
        )
        self.ICI_response_group = ICI_response_group

        nonICI_source = {}
        unique_sources = sorted(nonICI["Tissue_Source"].unique())

        # Iterate over each unique source and populate the dictionary
        for source in unique_sources:
            sample_ids = nonICI[nonICI["Tissue_Source"] == source]["Sample_ID"]
            nonICI_source[source] = sorted(
                [self.image2biop.index(i) for i in sample_ids]
            )
        self.nonICI_source = nonICI_source

        ICI_source = {}
        unique_sources = sorted(ICI["Tissue_Source"].unique())

        # Iterate over each unique source and populate the dictionary
        for source in unique_sources:
            sample_ids = ICI[ICI["Tissue_Source"] == source]["Sample_ID"]
            ICI_source[source] = sorted([self.image2biop.index(i) for i in sample_ids])
        self.ICI_source = ICI_source


def create_folder(folder_path: str):
    if os.path.exists(folder_path):
        # Clear the folder's contents
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")
        print(f"Folder cleared at: {folder_path}")
    else:
        # Create the folder if it does not exist
        os.makedirs(folder_path)
        print(f"Folder created at: {folder_path}")


def imageNum_to_biopNum():
    image_to_imageNum = pd.read_csv(os.path.join(HOMEDIR, "Data/AML_image.csv"))[
        ["FileName_Original"]
    ]
    aml_imageNum = []
    for i in range(len(image_to_imageNum)):
        aml_imageNum.append(image_to_imageNum.iloc[i, 0].split("_")[0])

    image_to_imageNum = pd.read_csv(os.path.join(HOMEDIR, "Data/NBM_image.csv"))[
        ["FileName_Original"]
    ]
    nbm_imageNum = []
    for i in range(len(image_to_imageNum)):
        nbm_imageNum.append(image_to_imageNum.iloc[i, 0].split("_")[0])

    print(len(aml_imageNum), len(nbm_imageNum))
    image2biop = aml_imageNum + nbm_imageNum

    return image2biop


def plotImage(
    dataFile: pd.DataFrame,
    dataFile_colNames: list,
    imageName,
    area_colName=None,
    selected_cellTypes: list = None,
    saveName=None,
):
    imageName_colName = dataFile_colNames[0]
    x_colName = dataFile_colNames[1]
    y_colName = dataFile_colNames[2]
    cellType_colName = dataFile_colNames[3]
    x_range = dataFile_colNames[4]
    y_range = dataFile_colNames[5]
    cellTypes = dataFile[cellType_colName].unique()

    fig, ax = plt.subplots()
    plt.xticks(range(0, x_range, 100))
    plt.yticks(range(0, y_range, 100))
    plt.xlim((-5, x_range))
    plt.ylim((-5, y_range))
    image = dataFile.loc[dataFile[imageName_colName] == imageName]

    colors = np.array(
        [
            "#c82423",
            "tab:green",
            "tab:orange",
            "tab:purple",
            "tab:brown",
            "tab:pink",
            "tab:gray",
            "tab:olive",
            "#3B5284",
            "tab:cyan",
            "k",
        ]
    )
    if selected_cellTypes:
        selectedTypes = np.in1d(image[cellType_colName], selected_cellTypes)
    else:
        selectedTypes = np.in1d(image[cellType_colName], cellTypes)

    image = image.iloc[selectedTypes]
    x = image[x_colName].to_list()
    y = image[y_colName].to_list()
    if area_colName:
        radius = np.sqrt(image[area_colName] / 3.14)
    for j in range(len(image)):
        index = np.where(cellTypes == image[cellType_colName].iloc[j])
        color = colors[index][0]
        r = radius.iloc[j] if area_colName else 5
        c = plt.Circle(xy=(x[j], y[j]), radius=r, color=color)
        ax.add_artist(c)
    plt.title(imageName)
    plt.gca().set_aspect("equal")
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    customeLines = []
    for j in range(len(cellTypes)):
        customeLines.append(
            Line2D([0], [0], marker="o", color=colors[j], label=cellTypes[j])
        )

    ax.legend(handles=customeLines, loc="center left", bbox_to_anchor=(1, 0.5))
    if saveName:
        plt.savefig(saveName)
    else:
        plt.show()
