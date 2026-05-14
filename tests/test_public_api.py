from __future__ import annotations

import importlib.util


def test_top_level_public_api_exports():
    import ClumPyCells as cpc

    assert cpc.__version__
    assert callable(cpc.runSpatial)
    assert callable(cpc.prepare_cell_table)
    assert callable(cpc.read_cell_table)
    assert callable(cpc.analyzeImage)
    assert cpc.MarkcorrResult.__name__ == "MarkcorrResult"


def test_packaged_streamlit_app_is_discoverable():
    assert importlib.util.find_spec("streamlit_app") is not None
