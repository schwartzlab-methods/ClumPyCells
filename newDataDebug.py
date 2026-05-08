import pickle
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


def get_CD34_count():

    AML_file = pd.read_csv(HOMEDIR + "Data/output/AML.csv")
    NBM_file = pd.read_csv(HOMEDIR + "Data/output/Normal.csv")
    AML_file_old = pd.read_csv(HOMEDIR + "Data/output/AML_old.csv")
    NBM_file_old = pd.read_csv(HOMEDIR + "Data/output/Normal_old.csv")

    AML_counts = (
        AML_file[AML_file["CellType"] == "CD34"]
        .groupby(["ImageNumber", "CellType"])
        .size()
        .reset_index(name="CD34_Count")
    )

    NBM_counts = (
        NBM_file[NBM_file["CellType"] == "CD34"]
        .groupby(["ImageNumber", "CellType"])
        .size()
        .reset_index(name="CD34_Count")
    )

    AML_counts_old = (
        AML_file_old[AML_file_old["CellType"] == "CD34"]
        .groupby(["ImageNumber", "CellType"])
        .size()
        .reset_index(name="CD34_Count")
    )
    NBM_counts_old = (
        NBM_file_old[NBM_file_old["CellType"] == "CD34"]
        .groupby(["ImageNumber", "CellType"])
        .size()
        .reset_index(name="CD34_Count")
    )

    print("AML CD34 counts:")
    print(AML_counts)

    print("\nAML CD34 counts (old):")
    print(AML_counts_old)

    print("\nNBM CD34 counts:")
    print(NBM_counts)
    print("\nNBM CD34 counts (old):")
    print(NBM_counts_old)


def combine_pickle_to_csv():
    base_dir = "/Users/leo/ClumPyCells/Result/AML/intensity_withSize2"

    for i in range(51):
        folder_name = os.path.join(base_dir, f"image_{i}")
        iso_file_path = os.path.join(folder_name, "iso.csv")
        data = {}
        if os.path.exists(folder_name):
            for filename in os.listdir(folder_name):
                if filename.endswith(".pkl"):
                    file_path = os.path.join(folder_name, filename)

                    # Load the pickle data
                    with open(file_path, "rb") as file:
                        try:
                            file_data = pickle.load(file)
                        except Exception as e:
                            print(f"Error loading {file_path}: {e}")
                            exit(1)

                    data[filename[:-4]] = file_data
            if len(data) == 121:
                iso = {key: val[0] for key, val in data.items()}
                df = pd.DataFrame(iso)

                if len(df.columns) not in [36, 49, 64, 81, 100, 121]:
                    logging.error(
                        f"image_{i} contains unfinished data (wrong column count: {len(df.columns)})"
                    )
                else:
                    df.to_csv(os.path.join(folder_name, "iso.csv"))
                    print(f"saved at {os.path.join(folder_name, 'iso.csv')}")

        # print(i, end=" ")


def plot_heatmap_individual_image(image_index):
    imageGroup = {f"image_{image_index}": [image_index]}
    a = AMLResult(
        groups=imageGroup,
        sizeCorrection=True,
        intensity=True,
        resultFolder="/Users/leo/ClumPyCells/Result/AML/intensity_withSize2/",
    )
    _, plots = a.getAUC(plot=True)
    image_auc_plot = a.displayImages(plots, order=[[f"image_{image_index}"]])
    image_auc_plot.save(f"{IMAGEFOLDER}image_{image_index}_auc.html")
    print(f"Heatmap for image_{image_index} saved.")


plot_heatmap_individual_image(37)
