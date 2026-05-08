"""Tests for the low-level Markcorr engine.

These tests cover the behaviours that downstream analyses rely on:
* the ``window`` / ``pointPattern`` building blocks (including the previously
  buggy in-window filter),
* a deterministic numerical regression for ``markcorr`` on the BetaCells
  fixture, and
* the equivalence guarantee: running with size correction enabled but with
  zero diameters must produce the exact same numbers as running without size
  correction at all.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ClumPyCells.Markcorr.markcorr import markcorr
from ClumPyCells.Markcorr.pointPattern import pointPattern
from ClumPyCells.Markcorr.window import window

BASELINE = Path(__file__).parent / "baseline" / "betacells_iso.npz"


# ---------------------------------------------------------------------------
# window
# ---------------------------------------------------------------------------


def test_window_inWindow_and_geometry():
    W = window([0, 10], [0, 5])
    assert W.inWindow(5, 2.5)
    assert W.inWindow(0, 0)
    assert W.inWindow(10, 5)
    assert not W.inWindow(-0.1, 2)
    assert not W.inWindow(2, 5.1)
    assert W.getSideLength() == [10, 5]
    assert W.minEdge() == 5


# ---------------------------------------------------------------------------
# pointPattern
# ---------------------------------------------------------------------------


def test_pointPattern_filters_points_outside_window():
    """Regression test for the pop-while-iterating bug in ``pointPattern``."""
    W = window([0, 10], [0, 10])
    x = [1.0, 5.0, -1.0, 11.0, 7.0]
    y = [1.0, 5.0, 5.0, 5.0, 7.0]
    marks = pd.DataFrame({"a": [1, 2, 3, 4, 5]})
    pp = pointPattern(list(x), list(y), [0] * 5, W, marks)
    assert pp.n == 3
    assert list(pp.getX()) == [1.0, 5.0, 7.0]
    assert list(pp.getY()) == [1.0, 5.0, 7.0]
    assert list(pp.getMarks()["a"]) == [1, 2, 5]


def test_pointPattern_default_diameter_is_zero():
    W = window([0, 10], [0, 10])
    pp = pointPattern([1.0, 2.0], [1.0, 2.0], None, W)
    assert pp.n == 2
    assert np.all(pp.getD() == 0)


# ---------------------------------------------------------------------------
# markcorr core regression
# ---------------------------------------------------------------------------


@pytest.fixture()
def betacells_pp(betacells_data):
    W = window([28.08, 778.08], [16.2, 1007.02])
    marks = pd.DataFrame({"type": betacells_data["type"]})
    marks["type"] = marks["type"].astype("category")
    return pointPattern(
        list(betacells_data["x"]),
        list(betacells_data["y"]),
        None,
        W,
        marks,
    )


def test_markcorr_baseline_regression(out_dir, betacells_pp):
    """The iso curves must remain bit-identical to the captured baseline."""
    r, funs = markcorr(
        betacells_pp,
        savefolder=str(out_dir) + "/",
        correction=["isotropic"],
        remove_zeros=False,
        saveImage=False,
        saveCache=False,
    )
    assert len(r) == 513

    baseline = np.load(BASELINE, allow_pickle=False)
    expected_keys = list(baseline["keys"])
    expected_values = baseline["values"]
    assert set(funs.keys()) >= set(expected_keys)
    for idx, key in enumerate(expected_keys):
        np.testing.assert_allclose(
            np.asarray(funs[key][0], dtype=float),
            expected_values[idx],
            rtol=0,
            atol=1e-12,
            err_msg=f"iso curve drift for {key!r}",
        )


def test_markcorr_size_correction_zero_equals_no_size_correction(
    out_dir: Path, betacells_data
):
    """sizeCorrection=True with zero diameters == sizeCorrection=False."""
    W = window([28.08, 778.08], [16.2, 1007.02])
    marks = pd.DataFrame({"type": betacells_data["type"]})
    marks["type"] = marks["type"].astype("category")

    a_dir = out_dir / "a"
    a_dir.mkdir()
    b_dir = out_dir / "b"
    b_dir.mkdir()

    pp_no_size = pointPattern(
        list(betacells_data["x"]), list(betacells_data["y"]), None, W, marks.copy()
    )
    pp_zero_size = pointPattern(
        list(betacells_data["x"]),
        list(betacells_data["y"]),
        [0.0] * len(betacells_data["x"]),
        W,
        marks.copy(),
    )

    _, fa = markcorr(
        pp_no_size,
        savefolder=str(a_dir) + "/",
        correction=["isotropic"],
        remove_zeros=False,
        saveImage=False,
        saveCache=False,
    )
    _, fb = markcorr(
        pp_zero_size,
        savefolder=str(b_dir) + "/",
        correction=["isotropic"],
        remove_zeros=False,
        saveImage=False,
        saveCache=False,
    )
    assert set(fa.keys()) == set(fb.keys())
    for key in fa:
        np.testing.assert_array_equal(
            np.asarray(fa[key][0]),
            np.asarray(fb[key][0]),
            err_msg=f"size correction with d=0 differs for {key!r}",
        )


def test_markcorr_isotropic_curve_around_unity_for_uncorrelated(out_dir: Path):
    """For a Poisson-like pattern with constant marks, kmm(r) ~ 1."""
    rng = np.random.default_rng(0)
    n = 200
    x = list(rng.uniform(0, 100, n))
    y = list(rng.uniform(0, 100, n))
    W = window([0, 100], [0, 100])
    marks = pd.DataFrame({"m": np.ones(n)})
    pp = pointPattern(x, y, None, W, marks)
    r, funs = markcorr(
        pp,
        savefolder=str(out_dir) + "/",
        correction=["isotropic"],
        remove_zeros=False,
        saveImage=False,
        saveCache=False,
    )
    iso = np.asarray(funs["m vs. m"][0], dtype=float)
    # For constant marks, the mark correlation function is identically 1.
    np.testing.assert_allclose(iso, 1.0, atol=5e-2)


def test_markcorr_pickle_cache_cleans_up(out_dir: Path, betacells_pp):
    """When saveCache=True, the per-pair pickles are removed at the end."""
    r, funs = markcorr(
        betacells_pp,
        savefolder=str(out_dir) + "/",
        correction=["isotropic"],
        remove_zeros=False,
        saveImage=False,
        saveCache=True,
    )
    leftover_pkls = list(out_dir.glob("*.pkl"))
    assert leftover_pkls == [], f"unexpected pickle leftovers: {leftover_pkls}"
    assert len(funs) > 0
