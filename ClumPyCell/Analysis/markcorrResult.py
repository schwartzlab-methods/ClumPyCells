import itertools
import logging
import sys

import altair as alt
import numpy as np
import pandas as pd
import scipy.stats as ss
import statsmodels.stats.multitest as multi
from permute.core import two_sample

from .metadata import *

# sys.path.append(HOMEDIR + "ClumPyCell/Analysis/altairThemes.py")

# if True:  # In order to bypass isort when saving
#     import altairThemes

# # register the custom theme under a chosen name
# alt.themes.register("publishTheme", altairThemes.publishTheme)

# # enable the newly registered theme
# alt.themes.enable("publishTheme")


class MarkcorrResult:
    def __init__(self, groups: dict, resultFolder: str, axisName: dict) -> None:
        """
        groups is a dictionary with the keys as the experiment names and the value as the number of images in that group
        """
        self.resultFolder = resultFolder
        self.groups = groups
        self.axisName = axisName

    def getCombinedResult(self):
        combinedResult = pd.DataFrame()
        for groupName in self.groups:
            for i in self.groups[groupName]:
                isoResult = pd.read_csv(f"{self.resultFolder}image_{i}/iso.csv").drop(
                    ["Unnamed: 0"], axis=1
                )
                isoMin = isoResult.where(isoResult > 0).min().min()
                isoResult = isoResult.where(isoResult > 0, isoMin)
                isoResult["imageNum"] = i
                combinedResult = pd.concat(
                    [combinedResult, isoResult], ignore_index=True
                )
        return combinedResult

    def getAUC(
        self,
        norm: str = "min_mid_max",
        takeMean=True,
        plot=True,
        r_range=None,
    ):
        axisName = self.axisName
        tot_auc = {}
        # find max and min data based on the whole dataset
        if norm == "min_mid_max":
            combinedResult = pd.DataFrame()
            for groupName in self.groups:
                for i in self.groups[groupName]:
                    isoResult = pd.read_csv(
                        f"{self.resultFolder}image_{i}/iso.csv"
                    ).drop(["Unnamed: 0"], axis=1)
                    isoMin = isoResult.where(isoResult > 0).min().min()
                    isoResult = isoResult.where(isoResult > 0, isoMin)
                    combinedResult = pd.concat([combinedResult, isoResult])
            combinedResult = combinedResult.to_numpy().flatten()
            combinedResult = combinedResult[combinedResult != 1]
            max_data, min_data = self.find_max_min(combinedResult)

        # Generate all possible combinations of axis values
        axis_keys = list(axisName.keys())
        combinations = list(itertools.product(axis_keys, repeat=2))
        formatted_combinations = [f"{comb[0]} vs. {comb[1]}" for comb in combinations]

        for groupName in self.groups:
            auc = pd.DataFrame(index=formatted_combinations)
            r = pd.read_csv(f"{self.resultFolder}image_{i}/r.csv").drop(
                ["Unnamed: 0"], axis=1
            )
            for i in self.groups[groupName]:
                isoResult = pd.read_csv(f"{self.resultFolder}image_{i}/iso.csv").drop(
                    ["Unnamed: 0"], axis=1
                )
                isoMin = isoResult.where(isoResult > 0).min().min()
                isoResult = isoResult.where(isoResult > 0, isoMin)
                if r_range:
                    isoResult = isoResult.loc[r["r"] <= r_range[1]].loc[
                        r["r"] >= r_range[0]
                    ]

                if norm == "log":
                    isoResult = np.log2(isoResult)
                elif norm == "min_mid_max":
                    isoResult = self.min_mid_max_normalization(
                        isoResult, max_data, 1, min_data
                    )
                else:
                    logging.error(
                        "Only min_mid_max and log are supported for the norm entry"
                    )
                auc["image_%d" % i] = isoResult.sum()

            auc = auc[auc != 0]
            tot_auc[groupName] = auc

        # plot the heat map with vega altair
        tot_plot = {}
        if plot:
            for groupName in self.groups:
                if takeMean:
                    auc_representation = tot_auc[groupName].mean(axis=1)
                else:
                    auc_representation = tot_auc[groupName].median(axis=1)
                heatmap_x, heatmap_y = MarkcorrResult.get_XY(
                    auc_representation, axisName
                )

                auc = pd.DataFrame(
                    {"from": heatmap_x, "to": heatmap_y, "auc": auc_representation}
                )

                heatmap = (
                    alt.Chart(auc, title=groupName)
                    .mark_rect()
                    .encode(
                        x=alt.X("from", axis=alt.Axis(labelAngle=-45)).sort(
                            list(axisName.values())
                        ),
                        y=alt.Y("to").sort(list(axisName.values())),
                        color=alt.Color(
                            "auc",
                            scale=alt.Scale(
                                scheme="redblue",
                                domainMid=0,
                                reverse=True,
                                type="symlog",
                            ),
                        ),
                    )
                )

                tot_plot[groupName] = heatmap
        return tot_auc, tot_plot

    def plotCurve(self, imageNum, type1, type2):
        col = f"{type1} vs. {type2}"
        isoResult = pd.read_csv(f"{self.resultFolder}image_{imageNum}/iso.csv").drop(
            ["Unnamed: 0"], axis=1
        )
        r = pd.read_csv(f"{self.resultFolder}image_{imageNum}/r.csv").drop(
            ["Unnamed: 0"], axis=1
        )
        colResult = isoResult[col]
        lineData = pd.concat([colResult, r], axis=1)
        lineData["ref"] = [0] * 513
        lineData["log"] = np.log2(lineData[col])
        lineData.rename(columns={col: "kmm"}, inplace=True)
        maxData = colResult.max() + 0.1
        minData = colResult.min() - 0.1

        line = (
            alt.Chart(lineData)
            .mark_area()
            .encode(
                x=alt.X("r", scale=alt.Scale(domain=[0, 160])),
                # y=alt.Y("kmm", title="kmm", scale=alt.Scale(domain=[minData, maxData])),
                y=alt.Y("log", title="log", scale=alt.Scale(domain=[-0.7, 0.3])),
            )
        )
        ref = (
            alt.Chart(lineData)
            .mark_line(strokeDash=[8, 8], color="red")
            .encode(x="r", y="ref")
        )
        return (line + ref).properties(width=480, height=320).configure_axis(grid=False)

    @staticmethod
    def get_XY(auc: pd.DataFrame, axisName=None):
        x = []
        y = []
        if not axisName:
            for i in range(len(auc)):
                x.append(auc.index[i].split(" vs. ")[0])
                y.append(auc.index[i].split(" vs. ")[1])
        else:
            for i in range(len(auc)):
                x.append(axisName[auc.index[i].split(" vs. ")[0]])
                y.append(axisName[auc.index[i].split(" vs. ")[1]])
        return x, y

    @staticmethod
    def find_max_min(data):
        q1 = np.nanpercentile(data, 25)
        q3 = np.nanpercentile(data, 75)
        iqr = q3 - q1
        lower_bound = q1 - 5 * iqr
        upper_bound = q3 + 5 * iqr
        filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
        max_data = np.max(filtered_data)
        min_data = np.min(filtered_data)
        return max_data, min_data

    @staticmethod
    def get_idd_columns(marks):
        idd_cols = []
        for i in range(len(marks)):
            for j in range(i, len(marks)):
                idd_cols.append(f"{marks[i]} vs. {marks[j]}")
        return idd_cols

    @staticmethod
    def min_mid_max_normalization(data, max, mid, min):
        data = data.map(
            lambda x: (
                -np.abs((mid - x)) / np.abs(mid - min)
                if x <= mid
                else np.abs(x - mid) / np.abs(max - mid)
            )
        )
        data.mask(data > 1, 1, inplace=True)
        data.mask(data < -1, -1, inplace=True)
        return data

    def getBoxPlot(self, auc: dict, cols, axisName: dict = {}):
        combined_auc = pd.DataFrame()
        for group in cols:
            group_auc = auc.get(group)
            auc_t = group_auc.transpose()
            auc_t = auc_t.melt(var_name="pair", value_name="auc")
            auc_t["group"] = group
            combined_auc = pd.concat([combined_auc, auc_t], ignore_index=True)

        selection = alt.selection_interval(bind="scales")

        # Boxplot
        boxplot = alt.Chart(combined_auc).mark_boxplot().encode(x="group:N", y="auc:Q")

        # Points for mean
        points = (
            alt.Chart(combined_auc)
            .mark_point(color="red")
            .encode(x="group:N", y="mean(auc):Q", tooltip="mean(auc):Q")
        )

        # Text for mean and std

        std_text = (
            alt.Chart(combined_auc)
            .mark_text(align="left", dx=5, dy=5)
            .encode(
                x="group:N",
                y="mean(auc):Q",
                text=alt.Text("stdev(auc):Q", format=".2f"),
            )
        )

        # Layering the plots
        layeredPlot = alt.layer(boxplot, points, std_text).facet(
            "pair:N", columns=len(axisName)
        )

        return layeredPlot

    @staticmethod
    def displayImages(singlePlots: dict, order):
        rowPlots = []
        for row in order:
            singleRowImage = alt.hconcat(*[singlePlots[i] for i in row])
            rowPlots.append(singleRowImage)
        combinedPlot = alt.vconcat(*rowPlots)
        return combinedPlot

    @staticmethod
    def find_diff(auc: dict, col1, col2, takeMean=True, method="MW", axisName=None):
        auc1_t = auc[col1].transpose()
        auc2_t = auc[col2].transpose()
        MWresult = {}

        if method == "MW":
            # perform mannwhitney U test
            for col in auc1_t:
                _, pValue = ss.mannwhitneyu(
                    auc1_t.get(col, pd.Series([0])),
                    auc2_t.get(col, pd.Series([0])),
                    nan_policy="omit",
                )
                MWresult[col] = pValue

        elif method == "perm":
            # perform fisher-pitman permuation test
            for col in auc1_t:
                pValue, _ = two_sample(
                    auc1_t.get(col, pd.Series([0])).dropna(),
                    auc2_t.get(col, pd.Series([0])).dropna(),
                    reps=10000,
                    alternative="two-sided",
                )
                MWresult[col] = pValue
        # Generate all possible combinations of axis values
        axis_keys = list(axisName.keys())
        combinations = list(itertools.product(axis_keys, repeat=2))
        formatted_combinations = [f"{comb[0]} vs. {comb[1]}" for comb in combinations]
        diffChart = pd.DataFrame(index=formatted_combinations)
        diffChart["p"] = MWresult

        if method == "MW":
            # apply fdr correction with benjamini hochberg correction
            diffChart["adjusted"] = multi.fdrcorrection(diffChart["p"])[1]
        else:
            diffChart["adjusted"] = diffChart["p"]

        heatmap_x, heatmap_y = MarkcorrResult.get_XY(diffChart, axisName)
        diffChart["from"] = heatmap_x
        diffChart["to"] = heatmap_y

        if takeMean:
            auc_rep1 = auc[col1].mean(axis=1)
            auc_rep2 = auc[col2].mean(axis=1)
        else:
            auc_rep1 = auc[col1].median(axis=1)
            auc_rep2 = auc[col2].median(axis=1)

        diff = auc_rep1 - auc_rep2
        diffChart["diff"] = diff
        diffChart["GT0"] = diffChart["diff"] > 0
        sig = diffChart.loc[diffChart["adjusted"] < 0.05]
        diffChart.to_csv(HOMEDIR + "/Result/Test/diffchart.csv")

        heatmap_shape = (
            alt.Chart(diffChart)
            .mark_point(size=200, filled=True)
            .encode(
                x=alt.X("from", axis=alt.Axis(labelAngle=-45)).sort(
                    list(axisName.values())
                ),
                y=alt.Y("to").sort(list(axisName.values())),
                color=alt.Color(
                    "diff",
                    scale=alt.Scale(
                        scheme="redblue", domainMid=0, reverse=True, type="symlog"
                    ),
                ),
                shape="GT0",
            )
        )
        triangle = (
            alt.Chart(sig)
            .mark_point(size=100, filled=False, shape="triangle-up", strokeWidth=0.756)
            .encode(
                x=alt.X("from").sort(list(axisName.values())),
                y=alt.Y("to").sort(list(axisName.values())),
                color=alt.value("black"),
            )
        )

        plot = alt.layer(heatmap_shape, triangle)

        return plot


class AMLResult(MarkcorrResult):
    def __init__(
        self,
        sizeCorrection: bool = False,
        intensity: bool = True,
        groups={"AML": range(36), "NBM": range(36, 51)},
    ) -> None:
        self.intensity = intensity
        if sizeCorrection:
            if intensity:
                resultFolder = HOMEDIR + "Result/AML/intensity_withSize/"
            else:
                resultFolder = HOMEDIR + "Result/AML/cellType_withSize/"
        else:
            if intensity:
                resultFolder = HOMEDIR + "Result/AML/intensity_withoutSize/"
            else:
                logging.error(
                    "cell type result without size correction has not been run yet"
                )

        if self.intensity:
            axisName = {
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
        else:
            axisName = {
                "Adipocytes": "Adipocytes",
                "B_cells": "B_cells",
                "HSC": "HSC",
                "Erythroids": "Erythroids",
                "Ki67": "Ki67",
                "MSC": "MSC",
                "Macrophages": "Macrophages",
                "Megakaryocytes": "Megakaryocytes",
                "Monocytes": "Monocytes",
                "Myeloids": "Myeloids",
                "T_cells": "T_cells",
            }
        super().__init__(groups, resultFolder=resultFolder, axisName=axisName)

    def getAUC(self, norm: str = "min_mid_max", plot=True, r_range=None, takeMean=True):
        axisName = self.axisName
        auc, plot = super().getAUC(
            norm=norm, takeMean=takeMean, plot=plot, r_range=r_range, axisName=axisName
        )
        return auc, plot


class MelanomaResult(MarkcorrResult):
    def __init__(self, intensity: bool = False, groups={"All": range(72)}) -> None:
        if intensity:
            resultFolder = HOMEDIR + "Result/Melanoma/Melanoma_intensity/"
        else:
            resultFolder = HOMEDIR + "Result/Melanoma/Melanoma_cellType_more/"
            axisName = {
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
        super().__init__(groups=groups, resultFolder=resultFolder, axisName=axisName)

    def getAUC(self, norm: str = "min_mid_max", plot=True, r_range=None, takeMean=True):
        axisName = self.axisName
        return super().getAUC(norm=norm, takeMean=takeMean, plot=plot, r_range=r_range)

    def getBoxPlot(self, auc: dict, cols):
        axisName = self.axisName
        return super().getBoxPlot(auc, cols, axisName)
