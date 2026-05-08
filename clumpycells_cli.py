"""Command-line entry point for running ClumPyCells without Streamlit.

Examples
--------
Inspect a CSV before choosing columns::

    python clumpycells_cli.py inspect-csv --csv cells.csv

Run markcorr serially::

    python clumpycells_cli.py run-markcorr --csv cells.csv --out results \
        --x-col centroid_x --y-col centroid_y --image-col sample_id \
        --area-col cell_area --mark cell_type --mark CD3

Run one CSV as a single image on 8 allocated HPC cores::

    python clumpycells_cli.py run-markcorr --csv cells.csv --out results \
        --image-col "" --max-workers 8 --mark cell_type
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ClumPyCells.ClumPyCells import runSpatial  # noqa: E402


def _none_if_blank(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip()
    if value == "" or value.lower() in {"none", "null", "na"}:
        return None
    return value


def _parse_mark_columns(values: list[str] | None) -> list[str] | None:
    if not values:
        return None
    columns: list[str] = []
    for value in values:
        for item in value.split(","):
            item = item.strip()
            if item:
                columns.append(item)
    return columns or None


def _default_workers() -> int:
    for env_name in ("SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS"):
        try:
            value = int(os.environ.get(env_name, ""))
        except ValueError:
            continue
        if value > 0:
            return value
    return 1


def _infer_ranges(csv_path: str, x_col: str, y_col: str):
    coords = pd.read_csv(csv_path, usecols=[x_col, y_col])
    x_values = pd.to_numeric(coords[x_col], errors="raise")
    y_values = pd.to_numeric(coords[y_col], errors="raise")
    return [float(x_values.min()), float(x_values.max())], [
        float(y_values.min()),
        float(y_values.max()),
    ]


def inspect_csv(args: argparse.Namespace) -> int:
    preview = pd.read_csv(args.csv, nrows=args.rows)
    print("Columns:")
    for index, column in enumerate(preview.columns, 1):
        print(f"  {index:>3}. {column}")
    print()
    print(f"Preview (first {min(args.rows, len(preview))} rows):")
    print(preview.to_string(index=False, max_cols=None))
    return 0


def run_markcorr(args: argparse.Namespace) -> int:
    image_col = _none_if_blank(args.image_col)
    area_col = _none_if_blank(args.area_col)
    mark_columns = _parse_mark_columns(args.mark)

    if args.max_workers < 1:
        raise SystemExit("--max-workers must be >= 1")

    if args.xrange is None or args.yrange is None:
        inferred_xrange, inferred_yrange = _infer_ranges(
            args.csv, args.x_col, args.y_col
        )
        xrange = args.xrange if args.xrange is not None else inferred_xrange
        yrange = args.yrange if args.yrange is not None else inferred_yrange
    else:
        xrange = args.xrange
        yrange = args.yrange

    pp_criterion = None
    if args.pp_area_threshold is not None:
        pp_criterion = (
            lambda frame, threshold=args.pp_area_threshold: frame["Area"] > threshold
        )

    Path(args.out).mkdir(parents=True, exist_ok=True)
    print("ClumPyCells run-markcorr")
    print(f"  csv:          {args.csv}")
    print(f"  out:          {args.out}")
    print(f"  x/y:          {args.x_col}, {args.y_col}")
    print(f"  image column: {image_col or '<single image>'}")
    print(f"  area column:  {area_col or '<none>'}")
    print(
        f"  marks:        {', '.join(mark_columns) if mark_columns else '<all non-metadata columns>'}"
    )
    print(f"  xrange:       {xrange}")
    print(f"  yrange:       {yrange}")
    print(f"  workers:      {args.max_workers}")

    started = time.time()
    runSpatial(
        csv_path=args.csv,
        savefolder=str(args.out).rstrip("/") + "/",
        xrange=xrange,
        yrange=yrange,
        sizeCorrection=args.size_correction,
        pp_criterion=pp_criterion,
        max_workers=args.max_workers,
        show_progress=not args.no_progress,
        x_col=args.x_col,
        y_col=args.y_col,
        image_col=image_col,
        area_col=area_col,
        mark_columns=mark_columns,
        chunksize=args.chunksize,
    )
    elapsed = time.time() - started
    print(f"Done in {elapsed:.1f}s")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run ClumPyCells from a terminal or HPC batch job."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    inspect_parser = subparsers.add_parser(
        "inspect-csv", help="Print columns and a small preview of an input CSV."
    )
    inspect_parser.add_argument("--csv", required=True, help="Input cell-table CSV.")
    inspect_parser.add_argument(
        "--rows", type=int, default=5, help="Preview row count."
    )
    inspect_parser.set_defaults(func=inspect_csv)

    run_parser = subparsers.add_parser(
        "run-markcorr", help="Run markcorr on every image in a CSV."
    )
    run_parser.add_argument("--csv", required=True, help="Input cell-table CSV.")
    run_parser.add_argument("--out", required=True, help="Output result folder.")
    run_parser.add_argument(
        "--x-col", default="x", help="Column containing x coordinates."
    )
    run_parser.add_argument(
        "--y-col", default="y", help="Column containing y coordinates."
    )
    run_parser.add_argument(
        "--image-col",
        default="ImageNum",
        help="Image ID column. Use an empty string to treat the whole CSV as one image.",
    )
    run_parser.add_argument(
        "--area-col",
        default="Area",
        help="Area column for size correction. Use an empty string if unavailable.",
    )
    run_parser.add_argument(
        "--mark",
        action="append",
        help="Mark column to include. Repeat or provide comma-separated names. If omitted, all non-metadata columns are used.",
    )
    run_parser.add_argument("--xrange", nargs=2, type=float, metavar=("MIN", "MAX"))
    run_parser.add_argument("--yrange", nargs=2, type=float, metavar=("MIN", "MAX"))
    run_parser.add_argument("--size-correction", action="store_true")
    run_parser.add_argument(
        "--pp-area-threshold",
        type=float,
        help="Treat cells with Area greater than this value as occluders.",
    )
    run_parser.add_argument(
        "--max-workers",
        type=int,
        default=_default_workers(),
        help="Parallel image workers. Defaults to scheduler CPU allocation if available, else 1.",
    )
    run_parser.add_argument("--chunksize", type=int, default=10000)
    run_parser.add_argument("--no-progress", action="store_true")
    run_parser.set_defaults(func=run_markcorr)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
