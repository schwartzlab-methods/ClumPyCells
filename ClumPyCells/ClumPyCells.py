import concurrent.futures
import json
import logging
import os
import shutil

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from .Analysis.metadata import *
from .Markcorr.markcorr import *

CANONICAL_COLUMNS = {"ImageNum", "x", "y", "Area"}


def _empty_to_none(value):
    if value is None:
        return None
    if isinstance(value, str) and value.strip() == "":
        return None
    return value


def _unique_output_name(name, used):
    if name not in CANONICAL_COLUMNS and name not in used:
        return name
    base = f"mark_{name}" if name in CANONICAL_COLUMNS else name
    candidate = base
    counter = 2
    while candidate in used or candidate in CANONICAL_COLUMNS:
        candidate = f"{base}_{counter}"
        counter += 1
    return candidate


def prepare_cell_table(
    cell_table,
    x_col="x",
    y_col="y",
    image_col="ImageNum",
    area_col="Area",
    mark_columns=None,
):
    """Return a canonical ClumPyCells input table.

    The markcorr engine expects ``x``, ``y`` and ``ImageNum`` column names.
    User-facing entry points can pass arbitrary source column names here. If
    ``image_col`` is blank/``None``, every row is treated as belonging to image
    ``1``.
    """
    x_col = _empty_to_none(x_col)
    y_col = _empty_to_none(y_col)
    image_col = _empty_to_none(image_col)
    area_col = _empty_to_none(area_col)

    if x_col not in cell_table.columns:
        raise ValueError(f"x column '{x_col}' was not found in the CSV")
    if y_col not in cell_table.columns:
        raise ValueError(f"y column '{y_col}' was not found in the CSV")
    if image_col is not None and image_col not in cell_table.columns:
        raise ValueError(f"image column '{image_col}' was not found in the CSV")
    if area_col is not None and area_col not in cell_table.columns:
        area_col = None

    out = pd.DataFrame(
        {
            "ImageNum": cell_table[image_col].to_numpy() if image_col else 1,
            "x": pd.to_numeric(cell_table[x_col], errors="raise"),
            "y": pd.to_numeric(cell_table[y_col], errors="raise"),
        }
    )

    excluded = {x_col, y_col, image_col, area_col, "Unnamed: 0"}
    excluded.discard(None)
    if area_col is not None:
        out["Area"] = pd.to_numeric(cell_table[area_col], errors="raise")

    if mark_columns is None:
        mark_columns = [c for c in cell_table.columns if c not in excluded]
    else:
        mark_columns = [c for c in mark_columns if c not in excluded]

    if not mark_columns:
        raise ValueError("Select at least one mark column for markcorr")

    used = set(out.columns)
    for column in mark_columns:
        if column not in cell_table.columns:
            raise ValueError(f"mark column '{column}' was not found in the CSV")
        out_name = _unique_output_name(str(column), used)
        out[out_name] = cell_table[column].to_numpy()
        used.add(out_name)

    return out


def read_cell_table(
    csv_path,
    x_col="x",
    y_col="y",
    image_col="ImageNum",
    area_col="Area",
    mark_columns=None,
    chunksize=10000,
):
    chunks = pd.read_csv(csv_path, chunksize=chunksize)
    return pd.concat(
        [
            prepare_cell_table(
                chunk,
                x_col=x_col,
                y_col=y_col,
                image_col=image_col,
                area_col=area_col,
                mark_columns=mark_columns,
            )
            for chunk in chunks
        ],
        ignore_index=True,
    )


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

    # Optional cell-size correction: derive a per-cell diameter from Area.
    pp = None
    if sizeCorrection:
        if "Area" not in imageData.columns:
            logging.error(
                "Size correction cannot be applied due to lack of Area column"
            )
            return

        area = imageData["Area"]
        d = (np.sqrt(area / np.pi) * 2).tolist()
        if pp_criterion:
            pp_df = imageData[pp_criterion(imageData)]
            if len(pp_df) > 0:
                pp_x = pp_df["x"].tolist()
                pp_y = pp_df["y"].tolist()
                pp_d = (np.sqrt(pp_df["Area"] / np.pi) * 2).tolist()
                pp = pointPattern(pp_x, pp_y, pp_d, W)
    else:
        d = None

    points = pointPattern(x, y, d, W, mark)

    if dropArea and "Area" in mark.columns:
        mark.drop(["Area"], axis=1, inplace=True)

    # Run mark cross correlation function
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
    csv_path,
    savefolder,
    xrange,
    yrange,
    sizeCorrection=False,
    pp_criterion=None,
    max_workers=None,
    progress_callback=None,
    show_progress=True,
    x_col="x",
    y_col="y",
    image_col="ImageNum",
    area_col="Area",
    mark_columns=None,
    chunksize=10000,
):
    """Run mark cross-correlation on every image in ``csv_path``.

    Parameters
    ----------
    progress_callback : callable, optional
        Invoked as ``progress_callback(done, total, image_num)`` after each
        image completes. Useful for driving an external progress bar (e.g.
        Streamlit's ``st.progress``).
    show_progress : bool
        Display a tqdm progress bar over images. Disable when running in a
        non-interactive context such as a notebook export or test suite.
    """
    csv_data = read_cell_table(
        csv_path,
        x_col=x_col,
        y_col=y_col,
        image_col=image_col,
        area_col=area_col,
        mark_columns=mark_columns,
        chunksize=chunksize,
    )
    if sizeCorrection and "Area" not in csv_data.columns:
        raise ValueError("sizeCorrection=True requires selecting an area column")

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
    image_numbers = list(csv_data["ImageNum"].unique())
    total = len(image_numbers)

    # Use multithreading to process each image
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_image = {
            executor.submit(
                analyzeImage,
                image_num,
                csv_data,
                savefolder,
                xrange,
                yrange,
                sizeCorrection,
                pp_criterion,
            ): image_num
            for image_num in image_numbers
        }
        iterator = concurrent.futures.as_completed(future_to_image)
        if show_progress:
            iterator = tqdm(iterator, total=total, desc="markcorr images", unit="img")
        results = []
        for done_count, fut in enumerate(iterator, 1):
            image_num = future_to_image[fut]
            results.append(fut.result())
            if progress_callback is not None:
                try:
                    progress_callback(done_count, total, image_num)
                except Exception:
                    pass
    return results
