"""Smoke tests verifying that ``Test.py`` can be imported without side effects
and that its public helper functions still work."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path


def test_legacy_Test_module_imports_cleanly(repo_root: Path):
    sys.path.insert(0, str(repo_root))
    if "Test" in sys.modules:
        del sys.modules["Test"]
    try:
        mod = importlib.import_module("Test")
    finally:
        sys.path.remove(str(repo_root))
    assert hasattr(mod, "test_markCorr_BetaCells")
    assert hasattr(mod, "test_pipeline")
    assert hasattr(mod, "modify_csv_index")
