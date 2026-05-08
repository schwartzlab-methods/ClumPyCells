"""End-to-end tests that exercise the full ``runSpatial`` -> ``MarkcorrResult``
pipeline that downstream notebooks depend on.

The tests stage every artefact under ``tmp_path`` so they can run in parallel
without touching the production ``Result/`` tree.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ClumPyCells.ClumPyCells import runSpatial
from ClumPyCells.Analysis.markcorrResult import MarkcorrResult

XRANGE = [28.08, 778.08]
YRANGE = [16.2, 1007.02]
AXIS = {"type_ on ": "on", "type_ off ": "off"}


def _run_pipeline(csv_path: Path, out_dir: Path, *, sizeCorrection: bool = False):
    runSpatial(
        csv_path=str(csv_path),
        savefolder=str(out_dir) + "/",
        xrange=XRANGE,
        yrange=YRANGE,
        sizeCorrection=sizeCorrection,
    )


def test_runSpatial_produces_expected_outputs(betacells_csv: Path, out_dir: Path):
    _run_pipeline(betacells_csv, out_dir)
    image_dirs = sorted(
        p for p in out_dir.iterdir() if p.is_dir() and p.name.startswith("image_")
    )
    assert len(image_dirs) == 4
    for d in image_dirs:
        assert (d / "iso.csv").exists()
        assert (d / "r.csv").exists()
        iso = pd.read_csv(d / "iso.csv")
        # 513 r values, plus the unnamed index column
        assert len(iso) == 513
        # at least the 4 type-vs-type pairs are present
        cols = [c for c in iso.columns if "vs." in c]
        assert len(cols) >= 4


def test_runSpatial_identical_inputs_produce_identical_outputs(
    betacells_csv: Path, out_dir: Path
):
    """Each of the 4 images is a copy of the same point pattern, so all per-image
    iso.csv files should be numerically identical."""
    _run_pipeline(betacells_csv, out_dir)
    isos = []
    for i in range(1, 5):
        df = pd.read_csv(out_dir / f"image_{i}" / "iso.csv").drop(
            columns=["Unnamed: 0"], errors="ignore"
        )
        isos.append(df.sort_index(axis=1))
    for i in range(1, 4):
        pd.testing.assert_frame_equal(isos[0], isos[i], check_exact=False, atol=1e-12)


def test_runSpatial_size_correction_with_zero_area_matches_no_size_correction(
    betacells_csv: Path, tmp_path: Path
):
    """The size-correction code path must reduce to the no-size-correction path
    when the cell areas are zero (so the implied diameters are also zero)."""
    df = pd.read_csv(betacells_csv)
    zero_csv = tmp_path / "betacells_zero_area.csv"
    df_zero = df.copy()
    df_zero["Area"] = 0.0
    df_zero.to_csv(zero_csv, index=False)

    a = tmp_path / "no_size"
    a.mkdir()
    b = tmp_path / "size_zero"
    b.mkdir()

    _run_pipeline(betacells_csv, a, sizeCorrection=False)
    _run_pipeline(zero_csv, b, sizeCorrection=True)

    for i in range(1, 5):
        df_a = pd.read_csv(a / f"image_{i}" / "iso.csv").drop(
            columns=["Unnamed: 0"], errors="ignore"
        )
        df_b = pd.read_csv(b / f"image_{i}" / "iso.csv").drop(
            columns=["Unnamed: 0"], errors="ignore"
        )
        # After size correction the dataframe also drops Area; both paths
        # should agree on the type-vs-type curves.
        common = sorted(set(df_a.columns) & set(df_b.columns))
        assert common, "expected at least one shared mark-pair column"
        np.testing.assert_allclose(
            df_a[common].to_numpy(dtype=float),
            df_b[common].to_numpy(dtype=float),
            rtol=0,
            atol=1e-12,
            err_msg=f"size correction with zero area drifts for image_{i}",
        )


def test_runSpatial_accepts_custom_column_mapping_and_selected_marks(
    betacells_csv: Path, tmp_path: Path
):
    df = pd.read_csv(betacells_csv).rename(
        columns={
            "ImageNum": "sample_id",
            "x": "centroid_x",
            "y": "centroid_y",
            "Area": "cell_area",
            "type": "phenotype",
        }
    )
    df["unused_feature"] = df["phenotype"]
    csv_path = tmp_path / "renamed_columns.csv"
    df.to_csv(csv_path, index=False)

    out = tmp_path / "custom_mapping"
    out.mkdir()
    runSpatial(
        csv_path=str(csv_path),
        savefolder=str(out) + "/",
        xrange=XRANGE,
        yrange=YRANGE,
        sizeCorrection=False,
        show_progress=False,
        x_col="centroid_x",
        y_col="centroid_y",
        image_col="sample_id",
        area_col="cell_area",
        mark_columns=["phenotype"],
    )

    iso = pd.read_csv(out / "image_1" / "iso.csv")
    pair_columns = [column for column in iso.columns if "vs." in column]
    assert pair_columns
    assert all("phenotype" in column for column in pair_columns)
    assert all("unused_feature" not in column for column in pair_columns)


def test_runSpatial_without_image_column_treats_csv_as_one_image(
    betacells_csv: Path, tmp_path: Path
):
    df = pd.read_csv(betacells_csv)
    df = df[df["ImageNum"] == 1].drop(columns=["ImageNum"])
    csv_path = tmp_path / "single_image.csv"
    df.to_csv(csv_path, index=False)

    out = tmp_path / "single_image_out"
    out.mkdir()
    runSpatial(
        csv_path=str(csv_path),
        savefolder=str(out) + "/",
        xrange=XRANGE,
        yrange=YRANGE,
        sizeCorrection=False,
        show_progress=False,
        image_col=None,
        mark_columns=["type"],
    )

    image_dirs = sorted(
        path.name
        for path in out.iterdir()
        if path.is_dir() and path.name.startswith("image_")
    )
    assert image_dirs == ["image_1"]
    assert (out / "image_1" / "iso.csv").exists()


# ---------------------------------------------------------------------------
# Analysis layer
# ---------------------------------------------------------------------------


def test_markcorrResult_getAUC_and_findDiff(betacells_csv: Path, out_dir: Path):
    _run_pipeline(betacells_csv, out_dir)
    result = MarkcorrResult(
        groups={"g1": [1, 2], "g2": [3, 4]},
        resultFolder=str(out_dir) + "/",
        axisName=AXIS,
    )
    auc, plots = result.getAUC(plot=True)
    assert set(auc.keys()) == {"g1", "g2"}
    # 2 marks -> 4 ordered pairs in the AUC table
    assert auc["g1"].shape == (4, 2)

    # All 4 images are identical so the within-group spread is zero.
    for group in ("g1", "g2"):
        # Each row must have either equal columns or be all-NaN.
        for _, row in auc[group].iterrows():
            vals = row.dropna().to_numpy(dtype=float)
            if vals.size > 1:
                assert np.allclose(
                    vals, vals[0]
                ), "identical inputs should yield identical AUC"

    # And because g1 == g2 image-for-image, find_diff must return diff == 0
    # everywhere and the saveCsv option must materialise a file on disk.
    diff_csv = out_dir / "diffchart.csv"
    plot = MarkcorrResult.find_diff(
        auc, "g1", "g2", axisName=AXIS, saveCsv=str(diff_csv)
    )
    assert plot is not None
    assert diff_csv.exists()
    diff_df = pd.read_csv(diff_csv, index_col=0)
    assert "p" in diff_df.columns
    np.testing.assert_allclose(diff_df["diff"].fillna(0).to_numpy(), 0.0, atol=1e-12)


def test_markcorrResult_combinedResult_replaces_non_positive(
    betacells_csv: Path, out_dir: Path
):
    _run_pipeline(betacells_csv, out_dir)
    result = MarkcorrResult(
        groups={"all": [1, 2, 3, 4]},
        resultFolder=str(out_dir) + "/",
        axisName=AXIS,
    )
    combined = result.getCombinedResult()
    assert "imageNum" in combined.columns
    # Non-positive entries in any iso curve should have been replaced with the
    # smallest positive value, so the data are strictly positive everywhere.
    numeric = combined.drop(columns=["imageNum"]).select_dtypes(include="number")
    assert (numeric.to_numpy() > 0).all()
