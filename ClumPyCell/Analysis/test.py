import os

from markcorrResult import *


def test_dataset():
    a = AMLResult(sizeCorrection=True, intensity=True)
    _, singlePlots = a.getAUC(plot=True, r_max=150)
    plot = a.displayImages(singlePlots=singlePlots, order=[["AML", "NBM"]])
    plot.save("/Users/leo/Schwartzlab/AML_project/rewriteImages/auc1.html")

    _, singlePlots = a.getAUC(plot=True)
    plot = a.displayImages(singlePlots=singlePlots, order=[["AML", "NBM"]])
    plot.save("/Users/leo/Schwartzlab/AML_project/rewriteImages/auc.html")


def test_mel_dataset():
    a = MelanomaResult(intensity=False)
    _, singlePlots = a.getAUC(plot=True, r_max=150)
    plot = a.displayImages(singlePlots=singlePlots, order=[["ICI", "non_ICI"]])
    plot.save("/Users/leo/Schwartzlab/AML_project/rewriteImages/auc1.html")

    _, singlePlots = a.getAUC(plot=True)
    plot = a.displayImages(singlePlots=singlePlots, order=[["ICI", "non_ICI"]])
    plot.save("/Users/leo/Schwartzlab/AML_project/rewriteImages/auc.html")


def test_p_value():
    a = AMLResult(sizeCorrection=True, intensity=True)
    auc, singlePlots = a.getAUC(plot=True)
    plot = a.find_diff(auc, "AML", "NBM")
    plot.save("/Users/leo/Schwartzlab/AML_project/rewriteImages/p-value.html")


def combine_cache():
    base_dir = (
        "/cluster/home/t114231uhn/AML_Public/Result/Melanoma/Melanoma_cellType_more"
    )
    for i in range(51, 72):
        folder_name = os.path.join(base_dir, f"image_{i}")
        # List to store data from each pickle file
        data = {}
        iso_file_path = os.path.join(folder_name, "iso.csv")
        # Iterate over each file in the directory
        if os.path.exists(iso_file_path) or os.path.exists(folder_name):
            for filename in os.listdir(folder_name):
                if filename.endswith(".pkl"):
                    # Construct the full file path
                    file_path = os.path.join(folder_name, filename)

                    # Load the list from the pickle file
                    with open(file_path, "rb") as file:
                        try:
                            file_data = pickle.load(file)
                        except:
                            print(file_path)
                            exit(1)
                    data[filename[:-4]] = file_data
                    iso = {}
            for j in data:
                iso[j] = data[j][0]
            iso = pd.DataFrame(iso)
            length = len(iso.columns)
            if length not in [36, 49, 64, 81, 100, 121]:
                errorMsg = f"image_{i} contains unfinished data detected by wrong number of output file"
                logging.error(errorMsg)
            if length > 36:
                iso.to_csv(folder_name + "/iso.csv")
        print(i, end="")


def convert_result():
    antibody_cell = {
        "Intensity_Adipocytes": "Adipocytes",
        "Intensity_CD163": "Macrophages",
        "Intensity_CD20": "B_cells",
        "Intensity_CD3": "T_cells",
        "Intensity_CD31": "Megakaryocytes",
        "Intensity_CD34": "HSC",
        "Intensity_CD68": "Monocytes",
        "Intensity_Erythroids": "Erythroids",
        "Intensity_Ki67": "Ki67",
        "Intensity_MPO": "Myeloids",
        "Intensity_NGFR": "MSC",
        "Intensity__Adipocytes": "Adipocytes",
        "Intensity__Bcells": "B_cells",
        "Intensity__CD34": "HSC",
        "Intensity__Erythroids": "Erythroids",
        "Intensity__Ki67": "Ki67",
        "Intensity__MSC": "MSC",
        "Intensity__Macrophages": "Macrophages",
        "Intensity__Megs": "Megakaryocytes",
        "Intensity__Monocytes": "Monocytes",
        "Intensity__Myeloid": "Myeloids",
        "Intensity__Tcells": "T_cells",
    }
    isoResult = pd.read_csv(
        "/Users/leo/Schwartzlab/AML_project/result1/Normal/image_0/iso.csv"
    ).drop("Unnamed: 0", axis=1)
    nameMapper = {}
    for cols in isoResult.columns:
        firstMark = antibody_cell[cols.split(" vs. ")[0]]
        secondMark = antibody_cell[cols.split(" vs. ")[1]]
        nameMapper[cols] = f"{firstMark} vs. {secondMark}"
    for i in range(15):
        isoResult = pd.read_csv(
            "/Users/leo/Schwartzlab/AML_project/result1/Normal/image_%d/iso.csv" % i
        ).drop("Unnamed: 0", axis=1)
        isoResult.rename(columns=nameMapper, inplace=True)
        os.mkdir(
            "/Users/leo/Schwartzlab/AML_project/result/AML/cellType_withSize/NBM/image_%d"
            % i
        )
        isoResult.to_csv(
            "/Users/leo/Schwartzlab/AML_project/result/AML/cellType_withSize/NBM/image_%d/iso.csv"
            % i
        )


def generate_level(level, max_level=4):
    # Base case: If the current level is the max level, return an empty structure
    if level > max_level:
        return ""

    # For simplicity, assume a pattern based on your provided structure to replicate
    node_template = """{
        "_item": [
            {
                "_barcode": {
                    "unCell": "node{node_id}"
                },
                "_cellRow": {
                    "unRow": {row}
                }
            }
        ],
        "_significance": null,
        "_distance": 5.972810594217861
    }"""

    # Generate nodes for the current level
    node1_str = node_template.format(
        node_id=1 + (level - 2) * 2, row=5 + (level - 2) * 1000
    )
    node2_str = node_template.format(
        node_id=2 + (level - 2) * 2, row=1406 + (level - 2) * 1000
    )

    # Generate the next level string
    next_level_str = generate_level(level + 1, max_level)

    # Combine current nodes and next level into the structure
    if next_level_str:
        combined_str = (
            f"[{node1_str}, [{next_level_str}], {node2_str}, [{next_level_str}]]"
        )
    else:
        combined_str = f"[{node1_str}, [], {node2_str}, []]"

    return combined_str


# Start from level 2 (as per the given structure) and generate up to level 6
level_2_structure = generate_level(2)

# The top-level structure, wrapping the generated structure for level 2 to 6
top_level_structure = f"[{level_2_structure}]"

print(top_level_structure)
