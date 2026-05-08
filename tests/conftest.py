"""Shared pytest fixtures for the ClumPyCells test suite.

The fixtures intentionally avoid touching the real ``Result/`` tree used by the
production pipelines; everything is staged in ``tmp_path`` so the tests can be
run repeatedly and in parallel.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

# Ensure the repo root is on ``sys.path`` so ``import ClumPyCells`` works when
# pytest is invoked from another directory.
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


# Geometry of the bundled BetaCells dataset.
BETACELLS_XRANGE = [28.08, 778.08]
BETACELLS_YRANGE = [16.2, 1007.02]


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return ROOT


@pytest.fixture(scope="session")
def betacells_data(repo_root: Path) -> dict:
    """Load the BetaCells JSON test fixture once per session."""
    path = repo_root / "Data" / "TestData" / "BetaCells.json"
    with path.open("r") as fh:
        return json.load(fh)


@pytest.fixture()
def betacells_csv(tmp_path: Path, betacells_data: dict) -> Path:
    """Build a 4-image CSV (2 per group) replicating the BetaCells fixture."""
    rows = []
    for img_num in range(1, 5):
        for i in range(len(betacells_data["x"])):
            rows.append(
                {
                    "ImageNum": img_num,
                    "x": betacells_data["x"][i],
                    "y": betacells_data["y"][i],
                    "Area": betacells_data["area"][i],
                    "type": betacells_data["type"][i],
                }
            )
    csv_path = tmp_path / "betacells.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    return csv_path


@pytest.fixture()
def out_dir(tmp_path: Path) -> Path:
    out = tmp_path / "out"
    out.mkdir()
    return out


@pytest.fixture(autouse=True)
def _chdir_repo_root(monkeypatch, repo_root: Path):
    """A few legacy modules read ``config.json`` from the cwd at import time;
    keep tests deterministic by chdir-ing to the repo root for every test."""
    monkeypatch.chdir(repo_root)
    yield
