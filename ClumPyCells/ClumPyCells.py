import concurrent.futures
import json
import os
import shutil

from .Analysis.metadata import *
from .Markcorr.markcorr import *


def analyzeImage(
    imageNum,
    imageData,
    savefolder,
    xrange,
    yrange,
    sizeCorrection=False,
    pp_criterion=None,
    dropArea=True,
):
    # Create folder for the image
    image_folder = os.path.join(savefolder, f"image_{imageNum}/")
    os.makedirs(image_folder, exist_ok=True)

    # Create point pattern based on filtered image Data
    imageData = imageData[imageData["ImageNum"] == imageNum]
    x = imageData["x"].tolist()
    y = imageData["y"].tolist()
    W = window(xrange=xrange, yrange=yrange)
    mark = imageData.drop(["x", "y", "ImageNum"], axis=1)

    # Create point pattern for large cells if size correction is to be used
    if sizeCorrection:
        if not "Area" in imageData.columns:
            logging.error(
                "Size correction cannot be applied due to lack of Area column"
            )
            return

        area = imageData["Area"]
        d = (np.sqrt(area / np.pi) * 2).tolist()
        points = pointPattern(x, y, d, W, mark)
        if pp_criterion:
            pp = imageData[pp_criterion(imageData)]
            if len(pp) == 0:
                pp = None
            else:
                pp_x = pp["x"].tolist()
                pp_y = pp["y"].tolist()
                pp_d = (np.sqrt(pp["Area"] / np.pi) * 2).tolist()
                pp = pointPattern(pp_x, pp_y, pp_d, W)
    else:
        d = None
        pp = None

    points = pointPattern(x, y, d, W, mark)

    if dropArea:
        mark.drop(["Area"], axis=1, inplace=True)
    print("here")
    # Run mark cross correlation function
    # return
    _, funs = markcorr(
        X=points,
        savefolder=image_folder,
        saveImage=False,
        pp=pp,
        correction=["isotropic"],
        remove_zeros=False,
        saveCache=True,
    )
    iso = {}
    for i in funs:
        iso[i] = funs[i][0]
    iso = pd.DataFrame(iso)
    iso.to_csv(os.path.join(image_folder, "iso.csv"))


def runSpatial(
    csv_path, savefolder, xrange, yrange, sizeCorrection=False, pp_criterion=None
):
    # Read CSV file in chunks to handle large files efficiently
    chunks = pd.read_csv(csv_path, chunksize=10000)
    csv_data = pd.concat(chunks)

    # Drop unnecessary columns
    if "Unnamed: 0" in csv_data.columns:
        csv_data = csv_data.drop(columns=["Unnamed: 0"])

    # Detect columns that are not all integers and set them as categorical
    for column in csv_data.columns:
        if column not in [
            "ImageNum",
            "x",
            "y",
            "Area",
        ] and not pd.api.types.is_integer_dtype(csv_data[column]):
            csv_data[column] = csv_data[column].astype("category")

    # Get unique image numbers
    image_numbers = csv_data["ImageNum"].unique()

    # Use multithreading to process each image
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(
                analyzeImage,
                image_num,
                csv_data,
                savefolder,
                xrange,
                yrange,
                sizeCorrection,
                pp_criterion,
            )
            for image_num in image_numbers
        ]
    results = [f.result() for f in futures]
