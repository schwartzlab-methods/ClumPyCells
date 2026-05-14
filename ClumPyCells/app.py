"""Console launcher for the ClumPyCells Streamlit application."""

from __future__ import annotations

import importlib.util
import sys


def main() -> int:
    """Run the packaged Streamlit UI."""
    spec = importlib.util.find_spec("streamlit_app")
    if spec is None or spec.origin is None:
        raise RuntimeError("Could not locate the packaged streamlit_app module")

    from streamlit.web import cli as streamlit_cli

    sys.argv = ["streamlit", "run", spec.origin, *sys.argv[1:]]
    return int(streamlit_cli.main() or 0)


__all__ = ["main"]
