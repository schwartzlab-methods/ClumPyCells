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
import concurrent.futures
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ClumPyCells.ClumPyCells import (
    analyzeImage,
    prepare_cell_table,
    runSpatial,
)  # noqa: E402

RUN_CONFIG_FILENAME = "clumpycells_run_config.json"


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


def _parse_image_ids(values: list[str] | None) -> list[str] | None:
    if not values:
        return None
    image_ids: list[str] = []
    for value in values:
        for item in value.split(","):
            item = item.strip()
            if item:
                image_ids.append(item)
    return image_ids or None


def _default_workers() -> int:
    for env_name in ("SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS"):
        try:
            value = int(os.environ.get(env_name, ""))
        except ValueError:
            continue
        if value > 0:
            return value
    return 1


def _recommended_workers(
    cpus_per_job: int, concurrent_jobs: int, threads_per_worker: int
) -> int:
    denom = max(1, int(concurrent_jobs) * int(threads_per_worker))
    return max(1, int(cpus_per_job) // denom)


def _set_math_threads(threads_per_worker: int) -> None:
    value = str(int(max(1, threads_per_worker)))
    for env_name in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[env_name] = value


def _jsonable(value):
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    return str(value)


def _write_run_config(out_folder: str | Path, config: dict) -> Path:
    path = Path(out_folder) / RUN_CONFIG_FILENAME
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(config), indent=2), encoding="utf-8")
    return path


def _load_run_config(path: str | None) -> dict:
    if not path:
        return {}
    config_path = Path(path)
    if config_path.is_dir():
        config_path = config_path / RUN_CONFIG_FILENAME
    if not config_path.exists():
        raise SystemExit(f"Run config not found: {config_path}")
    try:
        return json.loads(config_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SystemExit(f"Could not parse run config {config_path}: {exc}") from exc


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
    image_numbers = _parse_image_ids(args.image_id)

    cpus_per_job = int(args.cpus_per_job)
    concurrent_jobs = int(args.concurrent_jobs)
    threads_per_worker = int(args.threads_per_worker)
    if cpus_per_job < 1:
        raise SystemExit("--cpus-per-job must be >= 1")
    if concurrent_jobs < 1:
        raise SystemExit("--concurrent-jobs must be >= 1")
    if threads_per_worker < 1:
        raise SystemExit("--threads-per-worker must be >= 1")
    if args.max_workers is not None and args.max_workers < 1:
        raise SystemExit("--max-workers must be >= 1")

    max_workers = (
        int(args.max_workers)
        if args.max_workers is not None
        else _recommended_workers(cpus_per_job, concurrent_jobs, threads_per_worker)
    )

    _set_math_threads(threads_per_worker)

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
    print(
        f"  image IDs:    {', '.join(image_numbers) if image_numbers else '<all images>'}"
    )
    print(f"  xrange:       {xrange}")
    print(f"  yrange:       {yrange}")
    print(f"  cpus/job:     {cpus_per_job}")
    print(f"  jobs shared:  {concurrent_jobs}")
    print(f"  threads/worker: {threads_per_worker}")
    print(f"  workers:      {max_workers}")

    started = time.time()
    runSpatial(
        csv_path=args.csv,
        savefolder=str(args.out).rstrip("/") + "/",
        xrange=xrange,
        yrange=yrange,
        sizeCorrection=args.size_correction,
        pp_criterion=pp_criterion,
        max_workers=max_workers,
        show_progress=not args.no_progress,
        x_col=args.x_col,
        y_col=args.y_col,
        image_col=image_col,
        area_col=area_col,
        mark_columns=mark_columns,
        chunksize=args.chunksize,
        image_numbers=image_numbers,
    )
    config_path = _write_run_config(
        args.out,
        {
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "source": "clumpycells_cli run-markcorr",
            "csv_path": args.csv,
            "result_folder": args.out,
            "x_col": args.x_col,
            "y_col": args.y_col,
            "image_col": image_col,
            "area_col": area_col,
            "mark_columns": mark_columns,
            "image_numbers": image_numbers,
            "xrange": [float(xrange[0]), float(xrange[1])],
            "yrange": [float(yrange[0]), float(yrange[1])],
            "size_correction": bool(args.size_correction),
            "pp_area_threshold": args.pp_area_threshold,
            "max_workers": max_workers,
            "threads_per_worker": threads_per_worker,
            "cpus_per_job": cpus_per_job,
            "concurrent_jobs": concurrent_jobs,
        },
    )
    elapsed = time.time() - started
    print(f"Run config: {config_path}")
    print(f"Done in {elapsed:.1f}s")
    return 0


def run_permutation(args: argparse.Namespace) -> int:
    run_config = _load_run_config(args.run_config)
    csv_path = args.csv or run_config.get("csv_path")
    if not csv_path:
        raise SystemExit("--csv is required unless --run-config provides csv_path")

    x_col = args.x_col or run_config.get("x_col") or "x"
    y_col = args.y_col or run_config.get("y_col") or "y"
    image_col = _none_if_blank(
        args.image_col if args.image_col is not None else run_config.get("image_col")
    )
    area_col = _none_if_blank(
        args.area_col if args.area_col is not None else run_config.get("area_col")
    )
    mark_columns = _parse_mark_columns(args.mark) or run_config.get("mark_columns")
    requested_images = _parse_image_ids(args.image_id)
    size_correction = bool(
        args.size_correction or run_config.get("size_correction", False)
    )
    pp_area_threshold = (
        args.pp_area_threshold
        if args.pp_area_threshold is not None
        else run_config.get("pp_area_threshold")
    )

    cpus_per_job = int(args.cpus_per_job)
    concurrent_jobs = int(args.concurrent_jobs)
    threads_per_worker = int(args.threads_per_worker)
    if cpus_per_job < 1:
        raise SystemExit("--cpus-per-job must be >= 1")
    if concurrent_jobs < 1:
        raise SystemExit("--concurrent-jobs must be >= 1")
    if threads_per_worker < 1:
        raise SystemExit("--threads-per-worker must be >= 1")
    if args.max_workers is not None and args.max_workers < 1:
        raise SystemExit("--max-workers must be >= 1")

    max_workers = (
        int(args.max_workers)
        if args.max_workers is not None
        else _recommended_workers(cpus_per_job, concurrent_jobs, threads_per_worker)
    )
    _set_math_threads(threads_per_worker)

    configured_xrange = run_config.get("xrange")
    configured_yrange = run_config.get("yrange")
    if args.xrange is None and configured_xrange is not None:
        xrange = configured_xrange
    elif args.xrange is None:
        xrange, _ = _infer_ranges(csv_path, x_col, y_col)
    else:
        xrange = args.xrange

    if args.yrange is None and configured_yrange is not None:
        yrange = configured_yrange
    elif args.yrange is None:
        _, yrange = _infer_ranges(csv_path, x_col, y_col)
    else:
        yrange = args.yrange

    df_full = prepare_cell_table(
        pd.read_csv(csv_path),
        x_col=x_col,
        y_col=y_col,
        image_col=image_col,
        area_col=area_col,
        mark_columns=mark_columns,
    )
    if size_correction and "Area" not in df_full.columns:
        raise SystemExit("--size-correction requires selecting an area column")
    df_full["ImageNum"] = df_full["ImageNum"].astype(str)
    for column in df_full.columns:
        if column not in {
            "ImageNum",
            "x",
            "y",
            "Area",
        } and not pd.api.types.is_integer_dtype(df_full[column]):
            df_full[column] = df_full[column].astype("category")

    available_images = list(dict.fromkeys(df_full["ImageNum"].astype(str)))
    if requested_images is None:
        image_numbers = available_images
    else:
        requested = {str(image_id) for image_id in requested_images}
        image_numbers = [
            image_id for image_id in available_images if image_id in requested
        ]
    if not image_numbers:
        raise SystemExit("No requested image IDs were found in the input CSV")

    shuffle_cols = [
        column
        for column in df_full.columns
        if column not in {"ImageNum", "x", "y", "Area"}
    ]
    if not shuffle_cols:
        raise SystemExit("No mark columns are available to permute")

    pp_criterion = None
    if size_correction and pp_area_threshold is not None:
        pp_criterion = (
            lambda frame, threshold=pp_area_threshold: frame["Area"] > threshold
        )

    perm_root = Path(args.out) / "perm_runs"
    perm_root.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(int(args.seed))
    seeds = [int(rng.integers(0, 2**31 - 1)) for _ in range(int(args.n_perm))]

    print("ClumPyCells run-permutation")
    print(f"  csv:          {csv_path}")
    if args.run_config:
        print(f"  run config:   {args.run_config}")
    print(f"  out:          {perm_root}")
    print(f"  x/y:          {x_col}, {y_col}")
    print(f"  image column: {image_col or '<single image>'}")
    print(f"  area column:  {area_col or '<none>'}")
    print(f"  image IDs:    {', '.join(image_numbers)}")
    print(f"  permutations: {len(seeds)}")
    print(f"  workers:      {max_workers}")

    started = time.time()
    completed = 0
    total = len(seeds) * len(image_numbers)

    def _run_one(seed_value: int, image_id: str, perm_folder: Path):
        sub = df_full[df_full["ImageNum"] == str(image_id)].copy()
        permuted = sub[shuffle_cols].sample(
            frac=1, random_state=int(seed_value), ignore_index=True
        )
        sub.loc[:, shuffle_cols] = permuted.values
        analyzeImage(
            imageNum=str(image_id),
            imageData=sub,
            savefolder=str(perm_folder).rstrip("/") + "/",
            xrange=xrange,
            yrange=yrange,
            sizeCorrection=size_correction,
            pp_criterion=pp_criterion,
        )

    for seed_value in seeds:
        perm_folder = perm_root / f"perm_{seed_value}"
        perm_folder.mkdir(parents=True, exist_ok=True)
        if max_workers > 1 and len(image_numbers) > 1:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=max_workers
            ) as executor:
                futures = [
                    executor.submit(_run_one, seed_value, image_id, perm_folder)
                    for image_id in image_numbers
                ]
                for future in concurrent.futures.as_completed(futures):
                    future.result()
                    completed += 1
                    if not args.no_progress:
                        print(f"  completed {completed}/{total}", flush=True)
        else:
            for image_id in image_numbers:
                _run_one(seed_value, image_id, perm_folder)
                completed += 1
                if not args.no_progress:
                    print(f"  completed {completed}/{total}", flush=True)

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
    run_parser.add_argument(
        "--image-id",
        action="append",
        help="Image ID to process. Repeat or provide comma-separated IDs. If omitted, all images are processed.",
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
        default=None,
        help="Parallel image workers. If omitted, computed from --cpus-per-job, --concurrent-jobs, and --threads-per-worker.",
    )
    run_parser.add_argument(
        "--cpus-per-job",
        type=int,
        default=_default_workers(),
        help="CPUs allocated to this job/task (scheduler value).",
    )
    run_parser.add_argument(
        "--concurrent-jobs",
        type=int,
        default=1,
        help="Number of jobs sharing CPUs on the same node/allocation.",
    )
    run_parser.add_argument(
        "--threads-per-worker",
        type=int,
        default=1,
        help="Math-library threads per worker (OMP/MKL/OpenBLAS/NumExpr).",
    )
    run_parser.add_argument("--chunksize", type=int, default=10000)
    run_parser.add_argument("--no-progress", action="store_true")
    run_parser.set_defaults(func=run_markcorr)

    perm_parser = subparsers.add_parser(
        "run-permutation", help="Shuffle mark labels within images and run markcorr."
    )
    perm_parser.add_argument("--csv", help="Input cell-table CSV.")
    perm_parser.add_argument(
        "--run-config",
        help=f"Path to {RUN_CONFIG_FILENAME}, or to a result folder containing it. Values from this config are used when matching an original KMM run.",
    )
    perm_parser.add_argument(
        "--out",
        required=True,
        help="Result folder where perm_runs/perm_<seed>/image_<id>/ outputs are written.",
    )
    perm_parser.add_argument(
        "--x-col", default=None, help="Column containing x coordinates."
    )
    perm_parser.add_argument(
        "--y-col", default=None, help="Column containing y coordinates."
    )
    perm_parser.add_argument(
        "--image-col",
        default=None,
        help="Image ID column. Use an empty string to treat the whole CSV as one image.",
    )
    perm_parser.add_argument(
        "--area-col",
        default=None,
        help="Area column for optional size correction. Use an empty string if unavailable.",
    )
    perm_parser.add_argument(
        "--mark",
        action="append",
        help="Mark column to permute. Repeat or provide comma-separated names. If omitted, all non-metadata columns are used.",
    )
    perm_parser.add_argument(
        "--image-id",
        action="append",
        help="Image ID to process. Repeat or provide comma-separated IDs. If omitted, all images are processed.",
    )
    perm_parser.add_argument("--xrange", nargs=2, type=float, metavar=("MIN", "MAX"))
    perm_parser.add_argument("--yrange", nargs=2, type=float, metavar=("MIN", "MAX"))
    perm_parser.add_argument(
        "--seed", type=int, default=42, help="Seed used to generate permutation seeds."
    )
    perm_parser.add_argument(
        "--n-perm", type=int, default=10, help="Number of permutations to run."
    )
    perm_parser.add_argument("--size-correction", action="store_true")
    perm_parser.add_argument(
        "--pp-area-threshold",
        type=float,
        help="Treat cells with Area greater than this value as occluders during permutation runs.",
    )
    perm_parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Parallel image workers. If omitted, computed from --cpus-per-job, --concurrent-jobs, and --threads-per-worker.",
    )
    perm_parser.add_argument(
        "--cpus-per-job",
        type=int,
        default=_default_workers(),
        help="CPUs allocated to this job/task (scheduler value).",
    )
    perm_parser.add_argument(
        "--concurrent-jobs",
        type=int,
        default=1,
        help="Number of jobs sharing CPUs on the same node/allocation.",
    )
    perm_parser.add_argument(
        "--threads-per-worker",
        type=int,
        default=1,
        help="Math-library threads per worker (OMP/MKL/OpenBLAS/NumExpr).",
    )
    perm_parser.add_argument("--no-progress", action="store_true")
    perm_parser.set_defaults(func=run_permutation)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
