"""Correctness tests for the size-correction distance algorithm.

The cell-size correction in ClumPyCells reduces every pairwise distance by

    d_corrected = max(0, ||p_i - p_j|| - r_i - r_j - sum_k chord_k)

where ``r_i = dia_i / 2`` is the radius of point cell *i* and ``chord_k`` is
the length of the chord swept by an *occluder* cell ``pp_k`` whose centre's
foot-of-perpendicular lies on the segment between ``p_i`` and ``p_j``.

These tests pin down the mathematical contract of ``closePpairs`` so that
future optimisations cannot quietly drift away from the published algorithm.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ClumPyCells.Markcorr.closepairs import closePpairs
from ClumPyCells.Markcorr.pointPattern import pointPattern
from ClumPyCells.Markcorr.window import window

W = window([0.0, 100.0], [0.0, 100.0])


def _run(xx, yy, dia, rmax, pp=None):
    """Thin wrapper around ``closePpairs`` returning a dict keyed by (i, j)."""
    i, j, d, _ = closePpairs(xx, yy, dia, rmax, nguess=4096, pp=pp)
    # closePpairs returns 1-based indices.
    return {(int(a), int(b)): float(c) for a, b, c in zip(i, j, d)}


def test_closePpairs_no_correction_is_plain_distance():
    """With zero diameters and no occluders the algorithm reproduces the
    Euclidean distance exactly (this is the invariant that the pipeline-level
    test ``test_runSpatial_size_correction_with_zero_area_matches_no_size_correction``
    relies on)."""
    xx = np.array([0.0, 3.0, 7.0])
    yy = np.array([0.0, 4.0, 0.0])
    dia = np.zeros(3)
    out = _run(xx, yy, dia, rmax=10.0)
    assert out[(1, 2)] == pytest.approx(5.0)
    assert out[(1, 3)] == pytest.approx(7.0)
    assert out[(2, 3)] == pytest.approx(math.hypot(4.0, 4.0))


def test_closePpairs_subtracts_small_cell_radii():
    """Each pair distance must be reduced by ``(d_i + d_j) / 2``."""
    xx = np.array([0.0, 10.0])
    yy = np.array([0.0, 0.0])
    dia = np.array([2.0, 4.0])  # radii 1 and 2 -> 3 units of overlap removed
    out = _run(xx, yy, dia, rmax=20.0)
    assert out[(1, 2)] == pytest.approx(10.0 - 1.0 - 2.0)


def test_closePpairs_clamps_corrected_distance_to_zero():
    """If the small cells together are larger than the centre-to-centre
    distance, the corrected distance saturates at zero, never goes negative."""
    xx = np.array([0.0, 5.0])
    yy = np.array([0.0, 0.0])
    dia = np.array([6.0, 8.0])  # radii 3 + 4 = 7 > 5
    out = _run(xx, yy, dia, rmax=20.0)
    assert out[(1, 2)] == pytest.approx(0.0)


def test_closePpairs_subtracts_chord_through_occluder():
    """A large occluder centred *on* the segment removes its full diameter
    from the corrected distance."""
    xx = np.array([0.0, 20.0])
    yy = np.array([0.0, 0.0])
    dia = np.zeros(2)
    # Occluder of diameter 4 centred at (10, 0) -> chord through diameter = 4.
    pp = pointPattern([10.0], [0.0], [4.0], W)
    out = _run(xx, yy, dia, rmax=30.0, pp=pp)
    assert out[(1, 2)] == pytest.approx(20.0 - 4.0)


def test_closePpairs_chord_uses_perpendicular_distance():
    """A laterally-offset occluder removes only the chord, not its diameter."""
    xx = np.array([0.0, 20.0])
    yy = np.array([0.0, 0.0])
    dia = np.zeros(2)
    # Occluder radius 5, perpendicular offset 3 -> chord = 2*sqrt(25-9) = 8.
    pp = pointPattern([10.0], [3.0], [10.0], W)
    out = _run(xx, yy, dia, rmax=30.0, pp=pp)
    assert out[(1, 2)] == pytest.approx(20.0 - 8.0)


def test_closePpairs_offsegment_occluder_does_not_count():
    """An occluder whose foot-of-perpendicular lies *outside* the segment must
    not reduce the distance, even if the infinite line through the segment
    would intersect it.  The previous loop-based implementation got this wrong;
    the vectorised replacement now enforces a proper segment-membership check.
    """
    xx = np.array([0.0, 5.0])
    yy = np.array([0.0, 0.0])
    dia = np.zeros(2)
    # Occluder centred at x=20, well past the segment [0, 5].
    pp = pointPattern([20.0], [0.0], [4.0], W)
    out = _run(xx, yy, dia, rmax=10.0, pp=pp)
    assert out[(1, 2)] == pytest.approx(5.0)


def test_closePpairs_multiple_occluders_cumulative_chord():
    """Two non-overlapping occluders along the segment subtract the sum of
    their chords (cumulative size correction)."""
    xx = np.array([0.0, 30.0])
    yy = np.array([0.0, 0.0])
    dia = np.zeros(2)
    # Two occluders of diameter 4 centred on the segment; both contribute
    # their full diameter (4) to the chord sum.
    pp = pointPattern([10.0, 20.0], [0.0, 0.0], [4.0, 4.0], W)
    out = _run(xx, yy, dia, rmax=40.0, pp=pp)
    assert out[(1, 2)] == pytest.approx(30.0 - 4.0 - 4.0)


def test_closePpairs_chord_clamped_to_zero():
    """If the occluder chords would over-subtract, the corrected distance
    clamps to zero rather than becoming negative."""
    xx = np.array([0.0, 10.0])
    yy = np.array([0.0, 0.0])
    dia = np.zeros(2)
    # Single huge occluder spanning the whole segment.
    pp = pointPattern([5.0], [0.0], [50.0], W)
    out = _run(xx, yy, dia, rmax=20.0, pp=pp)
    assert out[(1, 2)] == pytest.approx(0.0)
