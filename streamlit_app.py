"""Streamlit UI for running ClumPyCells interactively or from HPC-friendly commands.

The app is organized by task rather than by a forced workflow. Users can either
run a new KMM/markcorr calculation or open downstream analysis tools for any
existing result folder.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

import altair as alt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ClumPyCells.Analysis.markcorrResult import MarkcorrResult  # noqa: E402
from ClumPyCells.ClumPyCells import prepare_cell_table, runSpatial  # noqa: E402

st.set_page_config(page_title="ClumPyCells", layout="wide")

RUN_CONFIG_FILENAME = "clumpycells_run_config.json"
DEMO_DIR = REPO_ROOT / "Data" / "Demo"
DEMO_CONFIG_PATH = DEMO_DIR / "demo_streamlit_config.json"
DEMO_CELLS_CSV = DEMO_DIR / "demo_cells.csv"
DEMO_CLINICAL_CSV = DEMO_DIR / "demo_clinical.csv"
DEMO_FEATURE_TABLE_CSV = DEMO_DIR / "demo_feature_table.csv"
DEMO_IMAGE_MAPPING_CSV = DEMO_DIR / "demo_image_to_clinical.csv"
DEMO_OUTPUT_FOLDER = REPO_ROOT / "Result" / "demo_streamlit_run"


def _ensure_state():
    defaults = {
        "uploaded_csv_path": None,
        "last_savefolder": None,
        "last_decision_tree_folder": None,
        "last_survival_folder": None,
    }
    for key, value in defaults.items():
        st.session_state.setdefault(key, value)


_ensure_state()


CONFIG_EXPORT_CATEGORIES = {
    "Input data and column mapping": (
        "active_task",
        "cell_csv_source",
        "cell_csv_path",
        "uploaded_csv_path",
        "run_",
        "perm_",
    ),
    "Run settings and size correction": (
        "xrange_",
        "yrange_",
        "size_correction",
        "pp_area_threshold",
        "cpus_per_job",
        "threads_per_worker",
        "max_workers",
        "save_folder",
        "show_terminal_progress",
    ),
    "HPC and terminal commands": ("hpc_",),
    "Downstream groups and heatmaps": (
        "downstream_",
        "axis_",
        "curve_",
        "auc_",
        "diff_",
        "boxplot_",
    ),
    "Decision tree and survival": (
        "dt_",
        "tree_",
        "sv_",
        "surv_",
        "mapping_",
        "last_decision_tree_folder",
        "last_survival_folder",
    ),
    "Result browsing and remembered folders": (
        "generic_",
        "last_savefolder",
        "downstream_result_folder",
        "survival_view_folder",
    ),
}


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


def _config_key_category(key) -> str:
    text = str(key)
    for category, prefixes in CONFIG_EXPORT_CATEGORIES.items():
        if any(text == prefix or text.startswith(prefix) for prefix in prefixes):
            return category
    return "Other settings"


def _export_config(categories=None) -> dict:
    selected_categories = (
        set(CONFIG_EXPORT_CATEGORIES.keys()) if categories is None else set(categories)
    )
    return {
        key: _jsonable(value)
        for key, value in st.session_state.items()
        if _is_exportable_config_key(key)
        and _config_key_category(key) in selected_categories
    }


def _is_exportable_config_key(key) -> bool:
    text = str(key)
    return not (
        text.startswith("_")
        or text.startswith("btn_")
        or text.startswith("train_")
        or text
        in {
            "load_config",
            "download_config",
            "config_import_upload",
            "dt_feature_table",
            "surv_clinical_upload",
            "surv_mapping_upload",
        }
    )


def _load_demo_config() -> dict:
    if DEMO_CONFIG_PATH.exists():
        config = json.loads(DEMO_CONFIG_PATH.read_text(encoding="utf-8"))
    else:
        config = {}
    config.update(
        {
            "active_task": "Run KMM / markcorr",
            "cell_csv_source": "Use demo mock data",
            "cell_csv_path": str(DEMO_CELLS_CSV),
            "uploaded_csv_path": str(DEMO_CELLS_CSV),
            "run_x_col": "x",
            "run_y_col": "y",
            "run_image_col": "ImageNum",
            "run_area_col": "Area",
            "run_mark_columns": ["type", "insulin", "glucagon"],
            "preview_cell_type_col": "type",
            "preview_cell_type_filter": ["off", "on"],
            "preview_show_area": True,
            "xrange_lo": 0.0,
            "xrange_hi": 780.0,
            "yrange_lo": 0.0,
            "yrange_hi": 1020.0,
            "size_correction": True,
            "pp_area_threshold": 360.0,
            "cpus_per_job": 2,
            "threads_per_worker": 1,
            "max_workers": 2,
            "save_folder": str(DEMO_OUTPUT_FOLDER),
            "perm_source_csv_path": str(DEMO_CELLS_CSV),
            "hpc_stable_csv_path": str(DEMO_CELLS_CSV),
            "hpc_output_folder": str(REPO_ROOT / "Result" / "demo_hpc_run"),
            "downstream_result_folder": str(DEMO_OUTPUT_FOLDER),
            "surv_clinical_path": str(DEMO_CLINICAL_CSV),
            "mapping_image_col": "ImageNum",
            "mapping_clinical_col": "sample_id",
        }
    )
    return config


def _demo_config_token() -> str:
    if not DEMO_CONFIG_PATH.exists():
        return "missing-demo-config"
    stat = DEMO_CONFIG_PATH.stat()
    return f"{stat.st_mtime_ns}:{stat.st_size}"


def _apply_demo_defaults_for_run():
    token = _demo_config_token()
    if st.session_state.get("_demo_defaults_token") == token:
        return
    for key, value in _load_demo_config().items():
        if key in {"active_task", "cell_csv_source"}:
            continue
        if _is_exportable_config_key(key):
            st.session_state[key] = value
    st.session_state["_demo_defaults_token"] = token


def _apply_pending_config():
    pending_config = st.session_state.pop("_pending_config", None)
    if not pending_config:
        return
    for key, value in pending_config.items():
        if _is_exportable_config_key(key):
            st.session_state[key] = value


_apply_pending_config()


def _section_title(title: str, info: str | None = None):
    st.markdown(f"#### {title}", help=info)


st.sidebar.title("ClumPyCells")
task = st.sidebar.radio(
    "Choose task",
    ["Run KMM / markcorr", "Downstream analysis"],
    key="active_task",
)

with st.sidebar.expander("App documentation", expanded=False):
    st.markdown("""
### What can I do here?
- **Run KMM / markcorr** starts a new spatial mark cross-correlation job.
- **Downstream analysis** opens existing KMM results for heatmaps, KMM curves,
  group tests, permutation tests, decision trees, survival analysis, and result
  browsing.

### Expected cell CSV
- One row per cell.
- Numeric x/y coordinate columns.
- Optional image ID column. If omitted, the whole CSV is one image.
- Optional area column. Required for size correction and large-cell cutoff.

### HPC / Slurm note
The app generates scheduler-neutral worker commands plus a Slurm array example.
Different clusters use different headers, queues, modules, and storage paths, so
the generated commands are designed to be copied into your local template rather
than assuming one universal batch-script format.
""")

with st.sidebar.expander("Configuration", expanded=False):
    config_upload = st.file_uploader(
        "Import app config JSON",
        type=["json"],
        key="config_import_upload",
        help=(
            "Import restores widget values saved in a previous config file. "
            "Export lets you choose broad groups such as input mappings, HPC settings, "
            "or downstream analysis options."
        ),
    )
    if config_upload is not None and st.button("Load config", key="load_config"):
        try:
            loaded_config = json.loads(config_upload.getvalue().decode("utf-8"))
            st.session_state["_pending_config"] = loaded_config
            st.rerun()
        except Exception as exc:
            st.error(f"Could not load config: {exc}")
    if hasattr(st, "popover"):
        export_container = st.popover("Export configuration")
    else:
        export_container = st.expander("Export configuration", expanded=False)
    with export_container:
        export_categories = st.multiselect(
            "Categories to include",
            list(CONFIG_EXPORT_CATEGORIES.keys()),
            default=list(CONFIG_EXPORT_CATEGORIES.keys()),
            key="config_export_categories",
            help="Choose broad areas of the app state to include in the exported JSON.",
        )
        selected_config = _export_config(export_categories)
        st.caption(f"{len(selected_config)} setting(s) will be exported.")
        if export_categories:
            st.download_button(
                "Download selected config",
                data=json.dumps(selected_config, indent=2),
                file_name="clumpycells_streamlit_config.json",
                mime="application/json",
                key="download_config",
            )
        else:
            st.warning("Choose at least one category to export.")


NO_IMAGE_OPTION = "No image column: treat the whole CSV as one image"
NO_AREA_OPTION = "No area column"


def _image_numeric_value(value):
    text = str(value).strip()
    if re.fullmatch(r"[+-]?\d+(?:\.0+)?", text):
        return int(float(text))
    return None


def _image_sort_key(value):
    text = str(value).strip()
    numeric_value = _image_numeric_value(text)
    if numeric_value is not None:
        return ((0, numeric_value),)
    return tuple(
        (0, int(part)) if part.isdigit() else (1, part.lower())
        for part in re.split(r"(\d+)", text)
        if part != ""
    )


def _sort_image_ids(image_ids) -> list[str]:
    return sorted(
        dict.fromkeys(str(image_id) for image_id in image_ids), key=_image_sort_key
    )


def _sort_label_values(values) -> list[str]:
    return sorted(
        dict.fromkeys(str(value) for value in values if pd.notna(value)),
        key=_image_sort_key,
    )


def _parse_image_selection_text(text: str, image_numbers: list[str]) -> list[str]:
    image_numbers = _sort_image_ids(image_numbers)
    available = {str(image_id): str(image_id) for image_id in image_numbers}
    numeric_lookup = {}
    for image_id in image_numbers:
        numeric_value = _image_numeric_value(image_id)
        if numeric_value is not None:
            numeric_lookup[numeric_value] = str(image_id)

    selected = []
    for token in re.split(r"[,;\s]+", text.strip()):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            start, end = [piece.strip() for piece in token.split("-", 1)]
            try:
                lo = int(start)
                hi = int(end)
            except ValueError:
                if token in available:
                    selected.append(available[token])
                continue
            if lo > hi:
                lo, hi = hi, lo
            selected.extend(
                numeric_lookup[number]
                for number in range(lo, hi + 1)
                if number in numeric_lookup
            )
        elif token in available:
            selected.append(available[token])
        else:
            numeric_value = _image_numeric_value(token)
            if numeric_value in numeric_lookup:
                selected.append(numeric_lookup[numeric_value])
    return _sort_image_ids(selected)


def _image_selector(
    label: str,
    image_numbers: list[str],
    key_prefix: str,
    default=None,
    allow_empty=False,
) -> list[str]:
    image_numbers = _sort_image_ids(image_numbers)
    default = _sort_image_ids(default if default is not None else image_numbers)
    default = [image_id for image_id in default if image_id in image_numbers]
    methods = ["All images", "Range/list", "Manual multi-select"]
    default_method_index = 0 if default == image_numbers else 2
    method = st.selectbox(
        f"{label} selection method",
        methods,
        index=default_method_index,
        key=f"{key_prefix}_selection_method",
        help="Use ranges like 1-25, 30, 42-60 or fall back to manual multi-select.",
    )
    if method == "All images":
        selected = image_numbers
    elif method == "Range/list":
        example = "1-25, 30, 42-60"
        if image_numbers:
            first = image_numbers[0]
            last = image_numbers[-1]
            example = f"{first}-{last}" if first != last else first
        text = st.text_input(
            f"{label} image IDs or ranges",
            value="",
            placeholder=example,
            key=f"{key_prefix}_range_text",
            help="Separate IDs/ranges with commas, spaces, semicolons, or new lines. Blank uses the default selection.",
        )
        selected = (
            _parse_image_selection_text(text, image_numbers)
            if text.strip()
            else default
        )
    else:
        selected = st.multiselect(
            label,
            image_numbers,
            default=default,
            key=f"{key_prefix}_manual_multiselect",
        )
        selected = _sort_image_ids(selected)

    if selected:
        st.caption(f"Selected {len(selected):,} of {len(image_numbers):,} image(s).")
    elif not allow_empty:
        st.warning("Select at least one image.")
    return selected


def _list_image_dirs(folder: str | Path) -> list[str]:
    path = Path(folder)
    if not path.is_dir():
        return []
    image_ids = []
    for child in path.iterdir():
        if child.is_dir() and child.name.startswith("image_"):
            image_ids.append(child.name.split("_", 1)[1])
    return _sort_image_ids(image_ids)


def _read_image_ids(csv_path: str | Path, image_col: str | None) -> list[str]:
    if image_col is None:
        return ["1"]
    values = pd.read_csv(csv_path, usecols=[image_col])[image_col].dropna().unique()
    return _sort_image_ids(values)


def _file_signature(path: str | Path) -> tuple[str, int, int]:
    resolved = Path(path).resolve()
    stat = resolved.stat()
    return str(resolved), int(stat.st_mtime_ns), int(stat.st_size)


def _run_config_path(result_folder: str | Path) -> Path:
    return Path(result_folder) / RUN_CONFIG_FILENAME


def _write_run_config(result_folder: str | Path, config: dict) -> Path:
    path = _run_config_path(result_folder)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(config), indent=2), encoding="utf-8")
    return path


def _load_run_config(result_folder: str | Path) -> dict:
    path = _run_config_path(result_folder)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _seed_session_defaults(seed_name: str, defaults: dict):
    token = json.dumps(_jsonable(defaults), sort_keys=True)
    token_key = f"_{seed_name}_defaults_token"
    if st.session_state.get(token_key) == token:
        return
    for key, value in defaults.items():
        if value is not None:
            st.session_state[key] = value
    st.session_state[token_key] = token


def _seed_mapping_state_from_run_config(
    key_prefix: str, columns: list[str], run_config: dict
):
    if not run_config:
        return
    token = json.dumps(
        _jsonable(
            {
                "columns": columns,
                "x_col": run_config.get("x_col"),
                "y_col": run_config.get("y_col"),
                "image_col": run_config.get("image_col"),
                "area_col": run_config.get("area_col"),
                "mark_columns": run_config.get("mark_columns"),
            }
        ),
        sort_keys=True,
    )
    token_key = f"_{key_prefix}_mapping_seed_token"
    if st.session_state.get(token_key) == token:
        return

    if run_config.get("x_col") in columns:
        st.session_state[f"{key_prefix}_x_col"] = run_config["x_col"]
    if run_config.get("y_col") in columns:
        st.session_state[f"{key_prefix}_y_col"] = run_config["y_col"]

    image_col = run_config.get("image_col")
    st.session_state[f"{key_prefix}_image_col"] = (
        image_col if image_col in columns else NO_IMAGE_OPTION
    )
    area_col = run_config.get("area_col")
    st.session_state[f"{key_prefix}_area_col"] = (
        area_col if area_col in columns else NO_AREA_OPTION
    )
    mark_columns = [
        column for column in run_config.get("mark_columns") or [] if column in columns
    ]
    if mark_columns:
        st.session_state[f"{key_prefix}_mark_columns"] = mark_columns
    st.session_state[token_key] = token


@st.cache_data(show_spinner=False)
def _load_cohort_table(
    csv_path: str,
    file_signature: tuple[str, int, int],
    columns: tuple[str, ...],
) -> pd.DataFrame:
    del file_signature
    return pd.read_csv(csv_path, usecols=list(columns))


def _peek_axis_keys(folder: str | Path) -> list[str]:
    path = Path(folder)
    for child in sorted(
        path.glob("image_*/iso.csv"),
        key=lambda item: _image_sort_key(item.parent.name.split("_", 1)[1]),
    ):
        try:
            df = pd.read_csv(child, nrows=1)
        except Exception:
            continue
        cols = [column for column in df.columns if " vs. " in column]
        return sorted(
            {column.split(" vs. ")[0] for column in cols}
            | {column.split(" vs. ")[1] for column in cols}
        )
    return []


def _list_pair_columns(folder: str | Path) -> list[str]:
    path = Path(folder)
    for child in sorted(
        path.glob("image_*/iso.csv"),
        key=lambda item: _image_sort_key(item.parent.name.split("_", 1)[1]),
    ):
        try:
            df = pd.read_csv(child, nrows=1)
        except Exception:
            continue
        return [column for column in df.columns if " vs. " in column]
    return []


def _default_index(options, candidates, fallback=0):
    lowered = {str(option).lower(): index for index, option in enumerate(options)}
    for candidate in candidates:
        index = lowered.get(str(candidate).lower())
        if index is not None:
            return index
    return fallback


def _scheduler_int(env_names, fallback):
    for env_name in env_names:
        try:
            value = int(os.environ.get(env_name, ""))
        except ValueError:
            continue
        if value > 0:
            return value
    return fallback


def _recommended_workers(
    cpus_per_job: int, concurrent_jobs: int, threads_per_worker: int
):
    denom = max(1, int(concurrent_jobs) * int(threads_per_worker))
    return max(1, int(cpus_per_job) // denom)


def _set_math_threads(threads_per_worker: int):
    value = str(int(max(1, threads_per_worker)))
    for env_name in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[env_name] = value


def _run_demo_markcorr(progress_callback=None) -> tuple[Path, Path]:
    config = _load_demo_config()
    csv_path = Path(config["cell_csv_path"])
    if not csv_path.exists():
        raise FileNotFoundError(f"Demo cell CSV not found: {csv_path}")
    output_folder = Path(config["save_folder"])
    output_folder.mkdir(parents=True, exist_ok=True)

    threshold = config.get("pp_area_threshold")
    pp_criterion = None
    if config.get("size_correction") and threshold is not None:
        pp_criterion = (
            lambda frame, _threshold=float(threshold): frame["Area"] > _threshold
        )

    _set_math_threads(int(config.get("threads_per_worker", 1)))
    runSpatial(
        csv_path=str(csv_path),
        savefolder=str(output_folder).rstrip("/") + "/",
        xrange=[float(config["xrange_lo"]), float(config["xrange_hi"])],
        yrange=[float(config["yrange_lo"]), float(config["yrange_hi"])],
        sizeCorrection=bool(config.get("size_correction", False)),
        pp_criterion=pp_criterion,
        max_workers=int(config.get("max_workers", 1)),
        progress_callback=progress_callback,
        show_progress=False,
        x_col=str(config["run_x_col"]),
        y_col=str(config["run_y_col"]),
        image_col=str(config.get("run_image_col") or ""),
        area_col=str(config.get("run_area_col") or ""),
        mark_columns=list(config.get("run_mark_columns") or []),
    )
    run_config_path = _write_run_config(
        output_folder,
        {
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "source": "streamlit_demo",
            "csv_path": str(csv_path),
            "result_folder": str(output_folder),
            "x_col": config["run_x_col"],
            "y_col": config["run_y_col"],
            "image_col": config.get("run_image_col"),
            "area_col": config.get("run_area_col"),
            "mark_columns": list(config.get("run_mark_columns") or []),
            "xrange": [float(config["xrange_lo"]), float(config["xrange_hi"])],
            "yrange": [float(config["yrange_lo"]), float(config["yrange_hi"])],
            "size_correction": bool(config.get("size_correction", False)),
            "pp_area_threshold": (float(threshold) if threshold is not None else None),
            "max_workers": int(config.get("max_workers", 1)),
            "threads_per_worker": int(config.get("threads_per_worker", 1)),
            "cpus_per_job": int(config.get("cpus_per_job", 1)),
        },
    )
    return output_folder, run_config_path


def _chunks(values: list[str], chunk_size: int) -> list[list[str]]:
    values = _sort_image_ids(values)
    chunk_size = max(1, int(chunk_size))
    return [
        values[index : index + chunk_size]
        for index in range(0, len(values), chunk_size)
    ]


def _area_distribution_chart(area_values: pd.Series, threshold: float):
    values = pd.to_numeric(area_values, errors="coerce").dropna().to_numpy(dtype=float)
    if values.size < 2:
        return None
    counts, edges = np.histogram(
        values, bins=min(60, max(12, values.size // 10)), density=True
    )
    centers = (edges[:-1] + edges[1:]) / 2.0
    smooth = pd.Series(counts).rolling(window=5, min_periods=1, center=True).mean()
    chart_df = pd.DataFrame({"Area": centers, "Density": smooth})
    rule_df = pd.DataFrame({"cutoff": [float(threshold)]})
    curve = (
        alt.Chart(chart_df)
        .mark_line(strokeWidth=2)
        .encode(
            x=alt.X("Area:Q", title="Area"),
            y=alt.Y("Density:Q", title="Density"),
            tooltip=[
                alt.Tooltip("Area:Q", format=".2f"),
                alt.Tooltip("Density:Q", format=".4f"),
            ],
        )
    )
    cutoff = (
        alt.Chart(rule_df)
        .mark_rule(color="#d62728", strokeDash=[5, 5], strokeWidth=2)
        .encode(x="cutoff:Q")
    )
    return (curve + cutoff).properties(height=260)


def _cohort_overview(df_cohort: pd.DataFrame, image_col, area_col, mark_columns):
    total_cells = int(len(df_cohort))
    total_images = 1 if image_col is None else int(df_cohort[image_col].nunique())
    total_marks = int(len(mark_columns))
    metric_cols = st.columns(3)
    metric_cols[0].metric("Cells in file", f"{total_cells:,}")
    metric_cols[1].metric("Images in file", f"{total_images:,}")
    metric_cols[2].metric("Selected mark columns", f"{total_marks:,}")

    if image_col is not None:
        image_counts = (
            df_cohort.groupby(image_col, dropna=False)
            .size()
            .rename("cell_count")
            .reset_index()
        )
        image_counts["image_id"] = image_counts[image_col].astype(str)
        image_counts = image_counts.sort_values(
            "image_id", key=lambda values: values.map(_image_sort_key)
        )
        display_counts = image_counts.head(40).copy()
        image_sort_order = display_counts["image_id"].tolist()
        st.altair_chart(
            alt.Chart(display_counts)
            .mark_bar()
            .encode(
                x=alt.X("image_id:N", title="Image ID", sort=image_sort_order),
                y=alt.Y("cell_count:Q", title="Cell count"),
                tooltip=["image_id:N", "cell_count:Q"],
            )
            .properties(height=220),
            use_container_width=True,
        )
        c1, c2, c3 = st.columns(3)
        c1.metric("Min cells / image", int(image_counts["cell_count"].min()))
        c2.metric("Median cells / image", int(image_counts["cell_count"].median()))
        c3.metric("Max cells / image", int(image_counts["cell_count"].max()))

    if mark_columns:
        if st.session_state.get("cohort_mark_summary") not in mark_columns:
            st.session_state.pop("cohort_mark_summary", None)
        mark_column = st.selectbox(
            "Mark column to summarize",
            mark_columns,
            key="cohort_mark_summary",
            help="This list comes from the currently selected mark columns in Data preview.",
        )
        mark_values = df_cohort[mark_column].dropna()
        if pd.api.types.is_numeric_dtype(mark_values):
            mark_frame = pd.DataFrame({mark_column: mark_values})
            st.altair_chart(
                alt.Chart(mark_frame)
                .mark_bar(opacity=0.8)
                .encode(
                    x=alt.X(
                        f"{mark_column}:Q", bin=alt.Bin(maxbins=40), title=mark_column
                    ),
                    y=alt.Y("count():Q", title="Cell count"),
                    tooltip=[alt.Tooltip("count():Q", title="Cells")],
                )
                .properties(height=220),
                use_container_width=True,
            )
        else:
            mark_counts = (
                mark_values.astype(str)
                .str.strip()
                .value_counts()
                .rename_axis(mark_column)
                .reset_index(name="cell_count")
            )
            st.altair_chart(
                alt.Chart(mark_counts)
                .mark_bar(opacity=0.8)
                .encode(
                    x=alt.X(f"{mark_column}:N", title=mark_column, sort="-y"),
                    y=alt.Y("cell_count:Q", title="Cell count"),
                    tooltip=[mark_column, "cell_count"],
                )
                .properties(height=220),
                use_container_width=True,
            )

    if area_col is not None:
        area_values = pd.to_numeric(df_cohort[area_col], errors="coerce").dropna()
        if not area_values.empty:
            hist_df = pd.DataFrame({"Area": area_values.to_numpy(dtype=float)})
            st.altair_chart(
                alt.Chart(hist_df)
                .mark_bar(opacity=0.8)
                .encode(
                    x=alt.X("Area:Q", bin=alt.Bin(maxbins=40), title="Area"),
                    y=alt.Y("count():Q", title="Cell count"),
                    tooltip=[alt.Tooltip("count():Q", title="Cells")],
                )
                .properties(height=220),
                use_container_width=True,
            )


def _find_artifacts(base_folder: str | Path, patterns: list[str]) -> list[Path]:
    base = Path(base_folder)
    if not base.exists():
        return []
    files = []
    for pattern in patterns:
        files.extend(base.rglob(pattern))
    return sorted(
        [path for path in files if path.is_file()], key=lambda path: str(path)
    )


def _render_file_preview(path: Path):
    suffix = path.suffix.lower()
    st.caption(f"Preview: `{path.name}`")
    if suffix == ".csv":
        st.dataframe(_format_display_table(pd.read_csv(path)), use_container_width=True)
    elif suffix in {".png", ".jpg", ".jpeg", ".gif", ".webp"}:
        st.image(str(path), use_container_width=True)
    elif suffix == ".svg":
        components.html(
            path.read_text(encoding="utf-8", errors="ignore"),
            height=700,
            scrolling=True,
        )
    elif suffix == ".html":
        components.html(
            path.read_text(encoding="utf-8", errors="ignore"),
            height=700,
            scrolling=True,
        )
    elif suffix in {".txt", ".log", ".json"}:
        st.code(path.read_text(encoding="utf-8", errors="ignore")[:20000])
    else:
        st.info("No inline preview for this file type. Use the download button.")


def _format_display_table(df: pd.DataFrame) -> pd.DataFrame:
    df = df.drop(
        columns=[column for column in df.columns if str(column).startswith("Unnamed")]
    )
    rename_map = {
        "p": "P Value",
        "p_value": "P Value",
        "raw_p": "P Value",
        "adjusted": "Adjusted P Value",
        "bh_corrected_p": "Adjusted P Value",
        "logrank_p": "Log-rank P Value",
        "logrank_p_bh": "Adjusted Log-rank P Value",
        "coef": "Cox Coefficient",
        "exp(coef)": "Hazard Ratio",
        "diff": "AUC Difference",
        "from": "From Mark",
        "to": "To Mark",
        "GT0": "Group A Greater Than Group B",
        "feature": "Feature",
        "significant": "Significant After Correction",
        "survival_time_diff": "Median Survival Difference",
        "median_gtm": "Median Survival Above Cutoff",
        "median_stm": "Median Survival Below Cutoff",
    }
    return df.rename(
        columns={column: rename_map.get(column, column) for column in df.columns}
    )


def _safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "chart"


def _cache_chart(chart, cache_folder: str | Path, filename: str) -> Path | None:
    try:
        cache_path = Path(cache_folder)
        cache_path.mkdir(parents=True, exist_ok=True)
        output_path = cache_path / f"{_safe_filename(filename)}.html"
        chart.save(str(output_path))
        return output_path
    except Exception:
        return None


def _column_mapping_controls(df_preview: pd.DataFrame, key_prefix: str):
    _section_title(
        "Column mapping",
        "Map your CSV column names into the roles ClumPyCells needs. If your CSV "
        "does not have an image ID column, choose the no-image option and the file "
        "will be treated as one image.",
    )
    columns = df_preview.columns.tolist()
    if not columns:
        st.error("The CSV has no columns.")
        st.stop()
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        x_col = st.selectbox(
            "x coordinate column",
            columns,
            index=_default_index(columns, ["x", "X", "centroid_x", "CenterX"]),
            key=f"{key_prefix}_x_col",
            help="Numeric cell centroid x coordinate.",
        )
    with c2:
        y_col = st.selectbox(
            "y coordinate column",
            columns,
            index=_default_index(columns, ["y", "Y", "centroid_y", "CenterY"]),
            key=f"{key_prefix}_y_col",
            help="Numeric cell centroid y coordinate.",
        )
    with c3:
        image_options = [NO_IMAGE_OPTION] + columns
        image_choice = st.selectbox(
            "image ID column",
            image_options,
            index=_default_index(
                image_options,
                ["ImageNum", "image", "image_id", "ImageID", "sample_id"],
                fallback=0,
            ),
            key=f"{key_prefix}_image_col",
            help="Choose the column that identifies each image, ROI, or sample. Select the no-image option for a single-image CSV.",
        )
        image_col = None if image_choice == NO_IMAGE_OPTION else image_choice
    with c4:
        area_options = [NO_AREA_OPTION] + columns
        area_choice = st.selectbox(
            "area column",
            area_options,
            index=_default_index(area_options, ["Area", "area", "cell_area"]),
            key=f"{key_prefix}_area_col",
            help="Area is optional unless you enable cell-size correction or area-scaled previews.",
        )
        area_col = None if area_choice == NO_AREA_OPTION else area_choice

    reserved = {x_col, y_col, image_col, area_col, "Unnamed: 0"}
    reserved.discard(None)
    default_marks = [column for column in columns if column not in reserved]
    selected_marks = st.multiselect(
        "Columns to use as marks",
        columns,
        default=default_marks,
        key=f"{key_prefix}_mark_columns",
        help="Coordinate, image ID, and area columns are excluded automatically.",
    )
    mark_columns = [column for column in selected_marks if column not in reserved]
    if len(mark_columns) != len(selected_marks):
        st.info("Coordinate, image ID, and area columns are not used as marks.")
    if not mark_columns:
        st.warning("Select at least one mark column before running markcorr.")
    return x_col, y_col, image_col, area_col, mark_columns


def _current_column_mapping(df_preview: pd.DataFrame, key_prefix: str):
    columns = df_preview.columns.tolist()
    x_col = st.session_state.get(f"{key_prefix}_x_col")
    y_col = st.session_state.get(f"{key_prefix}_y_col")
    image_choice = st.session_state.get(f"{key_prefix}_image_col")
    area_choice = st.session_state.get(f"{key_prefix}_area_col")

    if x_col not in columns:
        x_col = columns[_default_index(columns, ["x", "X", "centroid_x", "CenterX"])]
    if y_col not in columns:
        y_col = columns[_default_index(columns, ["y", "Y", "centroid_y", "CenterY"])]

    image_col = None
    if image_choice and image_choice != NO_IMAGE_OPTION and image_choice in columns:
        image_col = image_choice

    area_col = None
    if area_choice and area_choice != NO_AREA_OPTION and area_choice in columns:
        area_col = area_choice

    reserved = {x_col, y_col, image_col, area_col, "Unnamed: 0"}
    reserved.discard(None)
    default_marks = [column for column in columns if column not in reserved]
    selected_marks = st.session_state.get(f"{key_prefix}_mark_columns", default_marks)
    mark_columns = [
        column
        for column in selected_marks
        if column in columns and column not in reserved
    ]
    return x_col, y_col, image_col, area_col, mark_columns


def _default_cell_type_column(
    df_preview: pd.DataFrame, column_options: list[str]
) -> int:
    if not column_options:
        return 0
    lowered = {
        str(column).lower(): index for index, column in enumerate(column_options)
    }
    for candidate in ["type", "cell_type", "celltype", "phenotype", "class"]:
        if candidate.lower() in lowered:
            return int(lowered[candidate.lower()])
    for index, column in enumerate(column_options):
        if column in df_preview.columns and not pd.api.types.is_numeric_dtype(
            df_preview[column]
        ):
            return index
    return 0


def _load_image_preview_data(
    csv_path: str | Path,
    x_col: str,
    y_col: str,
    image_col: str | None,
    area_col: str | None,
    cell_type_col: str | None,
    image_id: str | None,
):
    preview_columns = [x_col, y_col]
    if image_col is not None:
        preview_columns.append(image_col)
    if area_col is not None:
        preview_columns.append(area_col)
    if cell_type_col is not None:
        preview_columns.append(cell_type_col)
    preview_columns = list(dict.fromkeys(preview_columns))

    frame = pd.read_csv(csv_path, usecols=preview_columns)
    if image_col is not None and image_id is not None:
        frame = frame[frame[image_col].astype(str) == str(image_id)].copy()
    total_points = len(frame)
    frame["_x"] = pd.to_numeric(frame[x_col], errors="coerce")
    frame["_y"] = pd.to_numeric(frame[y_col], errors="coerce")
    if area_col is not None:
        frame["_area"] = pd.to_numeric(frame[area_col], errors="coerce")
    if cell_type_col is not None:
        frame["_cell_type"] = frame[cell_type_col].astype(str).str.strip()
    else:
        frame["_cell_type"] = "Cells"
    return frame.dropna(subset=["_x", "_y"]), total_points


def _nice_tick_step(span: float, target_ticks: int = 8) -> float:
    if span <= 0:
        return 1.0
    raw_step = span / max(1, target_ticks)
    exponent = np.floor(np.log10(raw_step))
    fraction = raw_step / 10**exponent
    if fraction <= 1:
        nice = 1
    elif fraction <= 2:
        nice = 2
    elif fraction <= 5:
        nice = 5
    else:
        nice = 10
    return float(nice * 10**exponent)


def _image_preview_figure(
    frame: pd.DataFrame,
    x_col: str,
    y_col: str,
    area_col: str | None,
    show_area: bool,
    x_range: list[float],
    y_range: list[float],
):
    x_span = max(float(x_range[1]) - float(x_range[0]), 1.0)
    y_span = max(float(y_range[1]) - float(y_range[0]), 1.0)
    longest_axis_inches = 7.5
    longest_span = max(x_span, y_span)
    figure_width = max(3.0, longest_axis_inches * x_span / longest_span)
    figure_height = max(3.0, longest_axis_inches * y_span / longest_span)
    fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=130)

    cell_types = _sort_label_values(frame["_cell_type"].unique())
    cmap = plt.get_cmap("tab10")
    color_map = {
        cell_type: cmap(index % cmap.N) for index, cell_type in enumerate(cell_types)
    }
    for cell_type in cell_types:
        type_frame = frame[frame["_cell_type"] == cell_type].copy()
        color = color_map[cell_type]
        ax.scatter(
            type_frame["_x"],
            type_frame["_y"],
            s=6,
            color=color,
            alpha=0.95,
            label=str(cell_type),
            linewidths=0,
        )
        if show_area and area_col is not None and "_area" in type_frame.columns:
            for _, row in type_frame.iterrows():
                if pd.isna(row["_area"]) or float(row["_area"]) < 0:
                    continue
                radius = float(np.sqrt(float(row["_area"]) / np.pi))
                ax.add_patch(
                    Circle(
                        (float(row["_x"]), float(row["_y"])),
                        radius=radius,
                        facecolor=color,
                        edgecolor=color,
                        linewidth=0.45,
                        alpha=0.32,
                    )
                )

    tick_step = _nice_tick_step(max(x_span, y_span))
    ax.set_xlim(float(x_range[0]), float(x_range[1]))
    ax.set_ylim(float(y_range[0]), float(y_range[1]))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.xaxis.set_major_locator(MultipleLocator(tick_step))
    ax.yaxis.set_major_locator(MultipleLocator(tick_step))
    ax.grid(True, color="#d9d9d9", linewidth=0.55, alpha=0.8)
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=str(cell_type),
            markerfacecolor=color_map[cell_type],
            markersize=6,
        )
        for cell_type in cell_types
    ]
    if legend_handles:
        ax.legend(
            handles=legend_handles,
            title="Cell type",
            loc="upper center",
            bbox_to_anchor=(0.5, -0.11),
            ncol=min(4, len(legend_handles)),
            frameon=False,
        )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return fig


def _cli_command(
    csv_path,
    save_folder,
    x_col,
    y_col,
    image_col,
    area_col,
    mark_columns,
    xrange,
    yrange,
    size_correction,
    pp_threshold,
    max_workers,
    threads_per_worker,
    cpus_per_job,
    concurrent_jobs,
    image_ids=None,
):
    parts = [
        "CLPC/bin/python",
        "clumpycells_cli.py",
        "run-markcorr",
        "--csv",
        str(csv_path),
        "--out",
        str(save_folder),
        "--x-col",
        str(x_col),
        "--y-col",
        str(y_col),
        "--image-col",
        "" if image_col is None else str(image_col),
        "--area-col",
        "" if area_col is None else str(area_col),
        "--xrange",
        str(xrange[0]),
        str(xrange[1]),
        "--yrange",
        str(yrange[0]),
        str(yrange[1]),
        "--max-workers",
        str(int(max_workers)),
        "--threads-per-worker",
        str(int(threads_per_worker)),
        "--cpus-per-job",
        str(int(cpus_per_job)),
        "--concurrent-jobs",
        str(int(concurrent_jobs)),
    ]
    for column in mark_columns:
        parts.extend(["--mark", str(column)])
    if image_ids:
        parts.extend(["--image-id", ",".join(map(str, image_ids))])
    if size_correction:
        parts.append("--size-correction")
        if pp_threshold is not None:
            parts.extend(["--pp-area-threshold", str(float(pp_threshold))])
    command = " ".join(shlex.quote(part) for part in parts)
    exports = (
        f"export OMP_NUM_THREADS={int(threads_per_worker)}\n"
        f"export MKL_NUM_THREADS={int(threads_per_worker)}\n"
        f"export OPENBLAS_NUM_THREADS={int(threads_per_worker)}\n"
        f"export NUMEXPR_NUM_THREADS={int(threads_per_worker)}"
    )
    return f"{exports}\n{command}"


def _permutation_cli_command(
    csv_path,
    out_folder,
    x_col,
    y_col,
    image_col,
    area_col,
    mark_columns,
    xrange,
    yrange,
    seed,
    n_perm,
    size_correction,
    pp_threshold,
    max_workers,
    threads_per_worker,
    cpus_per_job,
    concurrent_jobs,
    image_ids=None,
):
    parts = [
        "CLPC/bin/python",
        "clumpycells_cli.py",
        "run-permutation",
        "--csv",
        str(csv_path),
        "--out",
        str(out_folder),
        "--x-col",
        str(x_col),
        "--y-col",
        str(y_col),
        "--image-col",
        "" if image_col is None else str(image_col),
        "--area-col",
        "" if area_col is None else str(area_col),
        "--xrange",
        str(xrange[0]),
        str(xrange[1]),
        "--yrange",
        str(yrange[0]),
        str(yrange[1]),
        "--seed",
        str(int(seed)),
        "--n-perm",
        str(int(n_perm)),
        "--max-workers",
        str(int(max_workers)),
        "--threads-per-worker",
        str(int(threads_per_worker)),
        "--cpus-per-job",
        str(int(cpus_per_job)),
        "--concurrent-jobs",
        str(int(concurrent_jobs)),
    ]
    for column in mark_columns:
        parts.extend(["--mark", str(column)])
    if image_ids:
        parts.extend(["--image-id", ",".join(map(str, image_ids))])
    if size_correction:
        parts.append("--size-correction")
        if pp_threshold is not None:
            parts.extend(["--pp-area-threshold", str(float(pp_threshold))])
    command = " ".join(shlex.quote(part) for part in parts)
    exports = (
        f"export OMP_NUM_THREADS={int(threads_per_worker)}\n"
        f"export MKL_NUM_THREADS={int(threads_per_worker)}\n"
        f"export OPENBLAS_NUM_THREADS={int(threads_per_worker)}\n"
        f"export NUMEXPR_NUM_THREADS={int(threads_per_worker)}"
    )
    return f"{exports}\n{command}"


def _auc_x_marker_count(auc_by_group, groups, min_valid_percent) -> int:
    total = 0
    for group_name, auc_table in (auc_by_group or {}).items():
        if "count" not in auc_table.columns:
            continue
        threshold = len(groups.get(group_name, [])) * float(min_valid_percent) / 100.0
        counts = pd.to_numeric(auc_table["count"], errors="coerce")
        total += int((counts < threshold).sum())
    return total


def _slurm_array_script(
    base_command: str,
    image_chunks: list[list[str]],
    cpus_per_job: int,
    mem_gb: int,
    hours: int,
    project_dir: str = "/path/to/ClumPyCells",
):
    chunks = " ".join(shlex.quote(",".join(chunk)) for chunk in image_chunks)
    array_max = max(0, len(image_chunks) - 1)
    return f"""#!/usr/bin/env bash
#SBATCH --job-name=clumpycells
#SBATCH --cpus-per-task={int(cpus_per_job)}
#SBATCH --mem={int(mem_gb)}G
#SBATCH --time={int(hours):02d}:00:00
#SBATCH --array=0-{array_max}
#SBATCH --output=logs/clumpycells_%A_%a.out
#SBATCH --error=logs/clumpycells_%A_%a.err

set -euo pipefail
cd {shlex.quote(project_dir)}
mkdir -p logs

IMAGE_CHUNKS=({chunks})
IMAGE_CHUNK="${{IMAGE_CHUNKS[$SLURM_ARRAY_TASK_ID]}}"

{base_command} --image-id "$IMAGE_CHUNK"
"""


def _kmm_curve_chart(
    result_folder: str | Path, image_id: str, pair_col: str, transform: str
):
    image_folder = Path(result_folder) / f"image_{image_id}"
    iso = pd.read_csv(image_folder / "iso.csv").drop(
        ["Unnamed: 0"], axis=1, errors="ignore"
    )
    r = pd.read_csv(image_folder / "r.csv").drop(
        ["Unnamed: 0"], axis=1, errors="ignore"
    )
    if pair_col not in iso.columns:
        raise ValueError(f"{pair_col} was not found in image_{image_id}/iso.csv")
    line_data = pd.DataFrame(
        {"r": r["r"], "kmm": pd.to_numeric(iso[pair_col], errors="coerce")}
    )
    if transform == "log2":
        positive_min = line_data.loc[line_data["kmm"] > 0, "kmm"].min()
        line_data["value"] = np.log2(
            line_data["kmm"].where(line_data["kmm"] > 0, positive_min)
        )
        y_title = "log2(KMM)"
        ref_value = 0
    else:
        line_data["value"] = line_data["kmm"]
        y_title = "KMM"
        ref_value = 1
    line_data["reference"] = ref_value
    curve = (
        alt.Chart(line_data)
        .mark_line(strokeWidth=2)
        .encode(
            x=alt.X("r:Q", title="Radius"),
            y=alt.Y("value:Q", title=y_title),
            tooltip=[
                alt.Tooltip("r:Q", format=".2f"),
                alt.Tooltip("value:Q", format=".4f"),
            ],
        )
    )
    ref = (
        alt.Chart(line_data)
        .mark_rule(strokeDash=[6, 4], color="#d62728")
        .encode(y="reference:Q")
    )
    return (curve + ref).properties(height=360)


def _artifact_browser(
    base_folder: str | Path,
    label: str,
    key_prefix: str,
    patterns=None,
    preferred_names=None,
):
    patterns = patterns or ["*.csv", "*.svg", "*.png", "*.html", "*.json", "*.txt"]
    files = _find_artifacts(base_folder, patterns)
    if not files:
        st.info(f"No {label} artifacts found for the selected folder.")
        return
    base = Path(base_folder)
    preferred_names = preferred_names or []
    default_index = 0
    for preferred_name in preferred_names:
        for index, file_path in enumerate(files):
            if file_path.name == preferred_name:
                default_index = index
                break
        if files and files[default_index].name == preferred_name:
            break
    chosen = st.selectbox(
        f"{label} artifact",
        files,
        index=default_index,
        format_func=lambda path: str(path.relative_to(base)),
        key=f"{key_prefix}_artifact",
    )
    _render_file_preview(chosen)
    st.download_button(
        f"Download selected {label} artifact",
        data=chosen.read_bytes(),
        file_name=chosen.name,
        mime="application/octet-stream",
        key=f"{key_prefix}_download",
    )


def _image_groups_controls(image_numbers: list[str], key_prefix: str):
    image_numbers = _sort_image_ids(image_numbers)
    groups: dict[str, list[str]] = {}
    with st.expander(
        "Image groups for group tests and supervised models", expanded=False
    ):
        group_help = (
            "Groups define which images belong together for AUC heatmaps, statistical "
            "comparisons, and supervised models. Use manual selection for real study arms "
            "or quick splits for testing the workflow."
        )
        n_groups = st.number_input(
            "Number of groups",
            min_value=1,
            max_value=10,
            value=2,
            step=1,
            key=f"{key_prefix}_n_groups",
            help=group_help,
        )
        assignment_mode = st.selectbox(
            "Group assignment method",
            ["Consecutive ranges", "Alternating split", "Manual per group"],
            key=f"{key_prefix}_assignment_mode",
            help="Use quick splits for large cohorts, or switch to manual per-group selection for custom assignments.",
        )
        cols = st.columns(min(int(n_groups), 4))
        if assignment_mode == "Consecutive ranges":
            default_groups = [
                _sort_image_ids(chunk.tolist())
                for chunk in np.array_split(
                    np.array(image_numbers, dtype=object), int(n_groups)
                )
            ]
        else:
            default_groups = [
                image_numbers[index :: int(n_groups)] for index in range(int(n_groups))
            ]

        for group_index in range(int(n_groups)):
            with cols[group_index % len(cols)]:
                name = st.text_input(
                    f"Group {group_index + 1} name",
                    value=f"group_{group_index + 1}",
                    key=f"{key_prefix}_group_name_{group_index}",
                )
                if assignment_mode == "Manual per group":
                    picked = _image_selector(
                        f"Images in {name}",
                        image_numbers,
                        f"{key_prefix}_group_images_{group_index}",
                        default=default_groups[group_index],
                    )
                else:
                    picked = default_groups[group_index]
                    st.caption(
                        f"{len(picked):,} image(s): "
                        f"{', '.join(picked[:8])}{' ...' if len(picked) > 8 else ''}"
                    )
                if name and picked:
                    groups[name] = list(map(str, picked))
    return groups


def _render_permutation_panel(result_folder, image_numbers):
    run_config = _load_run_config(result_folder)
    if run_config:
        st.success(f"Loaded KMM run settings from `{_run_config_path(result_folder)}`.")
        run_defaults = {
            "perm_source_csv_path": run_config.get("csv_path"),
            "perm_xrange_lo": (run_config.get("xrange") or [None, None])[0],
            "perm_xrange_hi": (run_config.get("xrange") or [None, None])[1],
            "perm_yrange_lo": (run_config.get("yrange") or [None, None])[0],
            "perm_yrange_hi": (run_config.get("yrange") or [None, None])[1],
            "perm_size_correction": bool(run_config.get("size_correction", False)),
            "perm_pp_area_threshold": run_config.get("pp_area_threshold"),
            "perm_hpc_output_folder": str(result_folder),
        }
        _seed_session_defaults("perm_run_config", run_defaults)
    else:
        st.info(
            "No saved KMM run settings were found in this result folder. "
            "Use the same CSV, column mapping, window, and size-correction settings "
            "that were used for the original KMM run."
        )
    _section_title(
        "Permutation test",
        "Shuffle selected mark columns within each image, then rerun KMM to create "
        "a null set under result_folder/perm_runs. The same settings can be copied "
        "to a terminal or Slurm array job for HPC runs.",
    )
    perm_csv_source = st.radio(
        "Permutation source CSV",
        ["Use saved/source path", "Use demo cell CSV"],
        horizontal=True,
        key="perm_csv_source",
    )
    if perm_csv_source == "Use demo cell CSV":
        csv_path = str(DEMO_CELLS_CSV)
        st.session_state["perm_source_csv_path"] = csv_path
        st.success(f"Using bundled demo cell CSV: `{csv_path}`")
    else:
        csv_path = st.text_input(
            "Original cell-table CSV",
            value=st.session_state.get("uploaded_csv_path") or "",
            key="perm_source_csv_path",
            help="Use the same cell-level CSV that produced the selected KMM result folder.",
        )
    if not csv_path:
        st.info("Provide the original cell-table CSV to run permutations.")
        return
    if not Path(csv_path).exists():
        st.error("The source CSV path does not exist.")
        return

    df_perm_preview = pd.read_csv(csv_path, nrows=2000)
    _seed_mapping_state_from_run_config(
        "perm", df_perm_preview.columns.tolist(), run_config
    )
    x_col, y_col, image_col, area_col, mark_columns = _column_mapping_controls(
        df_perm_preview, "perm"
    )
    if not mark_columns:
        st.warning("Select at least one mark column before running permutations.")
        return

    try:
        x_values = pd.to_numeric(df_perm_preview[x_col], errors="raise")
        y_values = pd.to_numeric(df_perm_preview[y_col], errors="raise")
    except Exception as exc:
        st.error(f"Selected x/y columns must be numeric: {exc}")
        return

    default_xrange = [float(x_values.min()), float(x_values.max())]
    default_yrange = [float(y_values.min()), float(y_values.max())]
    c1, c2 = st.columns(2)
    with c1:
        perm_xrange_lo = st.number_input(
            "Permutation xrange min", value=default_xrange[0], key="perm_xrange_lo"
        )
        perm_xrange_hi = st.number_input(
            "Permutation xrange max", value=default_xrange[1], key="perm_xrange_hi"
        )
    with c2:
        perm_yrange_lo = st.number_input(
            "Permutation yrange min", value=default_yrange[0], key="perm_yrange_lo"
        )
        perm_yrange_hi = st.number_input(
            "Permutation yrange max", value=default_yrange[1], key="perm_yrange_hi"
        )

    option_cols = st.columns(4)
    with option_cols[0]:
        seed = st.number_input("Random seed", value=42, step=1, key="perm_seed")
    with option_cols[1]:
        n_perm = st.number_input(
            "Number of permutations",
            min_value=1,
            max_value=1000,
            value=10,
            step=1,
            key="perm_n_perm",
        )
    with option_cols[2]:
        perm_size_correction = st.checkbox(
            "Use size correction",
            value=False,
            disabled=area_col is None,
            key="perm_size_correction",
            help="Enable this if the null runs should match a size-corrected KMM analysis.",
        )
    perm_pp_threshold = None
    if perm_size_correction and area_col is not None:
        area_values = pd.to_numeric(df_perm_preview[area_col], errors="coerce").dropna()
        with option_cols[3]:
            perm_pp_threshold = st.number_input(
                "Large-cell cutoff",
                value=(
                    float(area_values.quantile(0.9)) if not area_values.empty else 0.0
                ),
                key="perm_pp_area_threshold",
                help="Cells above this area are treated as occluders during permutation runs.",
            )

    chosen_images = _image_selector(
        "Images to permute",
        image_numbers,
        "perm_images",
        default=image_numbers,
    )
    if not chosen_images:
        return

    local_tab, hpc_tab = st.tabs(["Run in app", "HPC / terminal"])
    with local_tab:
        if st.button("Run permutations", type="primary", key="btn_perm_run"):
            from ClumPyCells.ClumPyCells import analyzeImage

            if float(perm_xrange_hi) <= float(perm_xrange_lo) or float(
                perm_yrange_hi
            ) <= float(perm_yrange_lo):
                st.error("Invalid permutation window: max must be greater than min.")
                return
            df_full = prepare_cell_table(
                pd.read_csv(csv_path),
                x_col=x_col,
                y_col=y_col,
                image_col=image_col or "",
                area_col=area_col or "",
                mark_columns=mark_columns,
            )
            if perm_size_correction and "Area" not in df_full.columns:
                st.error("Size correction requires selecting an area column.")
                return
            df_full["ImageNum"] = df_full["ImageNum"].astype(str)
            for column in df_full.columns:
                if column not in {
                    "ImageNum",
                    "x",
                    "y",
                    "Area",
                } and not pd.api.types.is_integer_dtype(df_full[column]):
                    df_full[column] = df_full[column].astype("category")
            progress = st.progress(0.0, text="Permuting")
            total = int(n_perm) * len(chosen_images)
            done = 0
            rng = np.random.default_rng(int(seed))
            perm_root = Path(result_folder) / "perm_runs"
            perm_root.mkdir(exist_ok=True)
            pp_criterion = None
            if perm_size_correction and perm_pp_threshold is not None:
                threshold = float(perm_pp_threshold)
                pp_criterion = (
                    lambda frame, _threshold=threshold: frame["Area"] > _threshold
                )
            for _ in range(int(n_perm)):
                seed_p = int(rng.integers(0, 2**31 - 1))
                perm_folder = perm_root / f"perm_{seed_p}"
                perm_folder.mkdir(exist_ok=True)
                for image_id in chosen_images:
                    sub = df_full[df_full["ImageNum"] == str(image_id)].copy()
                    shuffle_cols = [
                        column
                        for column in sub.columns
                        if column not in {"ImageNum", "x", "y", "Area"}
                    ]
                    permuted = sub[shuffle_cols].sample(
                        frac=1, random_state=seed_p, ignore_index=True
                    )
                    sub.loc[:, shuffle_cols] = permuted.values
                    analyzeImage(
                        imageNum=str(image_id),
                        imageData=sub,
                        savefolder=str(perm_folder) + "/",
                        xrange=[perm_xrange_lo, perm_xrange_hi],
                        yrange=[perm_yrange_lo, perm_yrange_hi],
                        sizeCorrection=bool(perm_size_correction),
                        pp_criterion=pp_criterion,
                    )
                    done += 1
                    progress.progress(done / total, text=f"{done}/{total}")
            progress.progress(1.0, text="Done")
            st.success(f"Permutation outputs written under `{perm_root}`.")

    with hpc_tab:
        st.caption(
            "Use these commands for offline or scheduler-based permutation runs. "
            "Replace placeholder paths with paths visible from compute nodes."
        )
        stable_csv_path = st.text_input(
            "Stable CSV path visible from compute nodes",
            value="/path/to/cells.csv",
            key="perm_hpc_stable_csv_path",
        )
        hpc_output_folder = st.text_input(
            "Permutation output root on the cluster",
            value="results_clumpycells",
            key="perm_hpc_output_folder",
        )
        project_dir = st.text_input(
            "Project directory on the cluster",
            value="/path/to/ClumPyCells",
            key="perm_hpc_project_dir",
        )
        hpc_c1, hpc_c2, hpc_c3, hpc_c4 = st.columns(4)
        with hpc_c1:
            hpc_cpus = st.number_input(
                "CPUs per array task",
                min_value=1,
                value=max(1, int(_scheduler_int(("SLURM_CPUS_PER_TASK",), 8))),
                step=1,
                key="perm_hpc_cpus",
            )
        with hpc_c2:
            hpc_workers = st.number_input(
                "Processes per task",
                min_value=1,
                value=max(1, min(int(hpc_cpus), len(chosen_images))),
                step=1,
                key="perm_hpc_workers",
            )
        with hpc_c3:
            hpc_threads = st.number_input(
                "Threads per process",
                min_value=1,
                value=1,
                step=1,
                key="perm_hpc_threads",
            )
        with hpc_c4:
            images_per_job = st.number_input(
                "Images per array task",
                min_value=1,
                value=max(1, int(hpc_workers)),
                step=1,
                key="perm_hpc_images_per_job",
            )
        mem_gb = st.number_input(
            "Memory per task (GB)", min_value=1, value=32, step=1, key="perm_hpc_mem_gb"
        )
        hours = st.number_input(
            "Wall time per task (hours)",
            min_value=1,
            value=4,
            step=1,
            key="perm_hpc_hours",
        )

        single_cmd = _permutation_cli_command(
            stable_csv_path,
            hpc_output_folder,
            x_col,
            y_col,
            image_col,
            area_col,
            mark_columns,
            [perm_xrange_lo, perm_xrange_hi],
            [perm_yrange_lo, perm_yrange_hi],
            seed,
            n_perm,
            perm_size_correction,
            perm_pp_threshold,
            hpc_workers,
            hpc_threads,
            hpc_cpus,
            1,
            image_ids=chosen_images,
        )
        st.subheader("Single terminal command")
        st.code(single_cmd, language="bash")

        hpc_chunks = _chunks(chosen_images, int(images_per_job))
        st.metric("Array tasks needed", len(hpc_chunks))
        chunk_table = pd.DataFrame(
            {
                "task_id": list(range(len(hpc_chunks))),
                "image_ids": [",".join(chunk) for chunk in hpc_chunks],
                "command_suffix": [
                    f"--image-id {','.join(chunk)}" for chunk in hpc_chunks
                ],
            }
        )
        st.dataframe(chunk_table, use_container_width=True)
        slurm_base = _permutation_cli_command(
            stable_csv_path,
            hpc_output_folder,
            x_col,
            y_col,
            image_col,
            area_col,
            mark_columns,
            [perm_xrange_lo, perm_xrange_hi],
            [perm_yrange_lo, perm_yrange_hi],
            seed,
            n_perm,
            perm_size_correction,
            perm_pp_threshold,
            hpc_workers,
            hpc_threads,
            hpc_cpus,
            1,
        ).replace("\n", " && ")
        st.subheader("Slurm array template")
        st.code(
            _slurm_array_script(
                slurm_base,
                hpc_chunks,
                int(hpc_cpus),
                int(mem_gb),
                int(hours),
                project_dir=project_dir,
            ),
            language="bash",
        )


if task == "Run KMM / markcorr":
    st.header("Run KMM / markcorr")
    st.caption("Start a new mark cross-correlation run from a cell-level CSV.")

    source = st.radio(
        "Cell CSV source",
        ["Upload CSV", "Use path on this machine or HPC", "Use demo mock data"],
        horizontal=True,
        key="cell_csv_source",
    )
    csv_path = None
    if source == "Upload CSV":
        uploaded = st.file_uploader("Cell-table CSV", type=["csv"])
        if uploaded is not None:
            tmp_csv = Path(tempfile.gettempdir()) / f"cpc_input_{uploaded.name}"
            tmp_csv.write_bytes(uploaded.getvalue())
            csv_path = str(tmp_csv)
            st.session_state["uploaded_csv_path"] = csv_path
    elif source == "Use path on this machine or HPC":
        csv_path = st.text_input(
            "CSV path",
            value=st.session_state.get("uploaded_csv_path") or "",
            key="cell_csv_path",
        )
        if csv_path:
            st.session_state["uploaded_csv_path"] = csv_path
    else:
        _apply_demo_defaults_for_run()
        csv_path = str(DEMO_CELLS_CSV)
        st.session_state["uploaded_csv_path"] = csv_path
        st.success(f"Using bundled demo data: `{csv_path}`")

    if not csv_path:
        st.info("Upload a CSV or provide a path to begin.")
        st.stop()
    if not Path(csv_path).exists():
        st.error("The selected CSV path does not exist.")
        st.stop()

    df_preview = pd.read_csv(csv_path, nrows=2000)
    tab_data, tab_cohort, tab_image_preview, tab_settings, tab_hpc = st.tabs(
        [
            "Data preview",
            "Cohort analysis",
            "Image preview",
            "Run settings",
            "HPC / terminal",
        ]
    )

    with tab_data:
        _section_title(
            "Data preview",
            "Use this tab to verify the raw CSV and choose the columns that define "
            "coordinates, image IDs, areas, and marks. The cohort and image preview "
            "tabs use the same current selections.",
        )
        st.dataframe(df_preview.head(50), use_container_width=True)
        _column_mapping_controls(df_preview, "run")

    x_col, y_col, image_col, area_col, mark_columns = _current_column_mapping(
        df_preview, "run"
    )

    try:
        x_values = pd.to_numeric(df_preview[x_col], errors="raise")
        y_values = pd.to_numeric(df_preview[y_col], errors="raise")
    except Exception as exc:
        st.error(f"Selected x/y columns must be numeric: {exc}")
        st.stop()
    x_min, x_max = float(x_values.min()), float(x_values.max())
    y_min, y_max = float(y_values.min()), float(y_values.max())
    image_ids = _read_image_ids(csv_path, image_col)

    with tab_cohort:
        _section_title(
            "Cohort analysis",
            "These summaries update from the current column mapping in Data preview, "
            "including the selected image ID, area, and mark columns. The full CSV is "
            "used here, not just the displayed preview rows.",
        )
        st.caption(
            f"Using image column `{image_col or 'single image'}`, area column "
            f"`{area_col or 'none'}`, and {len(mark_columns)} selected mark column(s)."
        )
        cohort_columns = [x_col]
        if image_col is not None:
            cohort_columns.append(image_col)
        if area_col is not None:
            cohort_columns.append(area_col)
        cohort_columns.extend(mark_columns)
        cohort_columns = tuple(dict.fromkeys(cohort_columns))
        try:
            with st.spinner("Loading full CSV for cohort summaries"):
                cohort_df = _load_cohort_table(
                    str(csv_path), _file_signature(csv_path), cohort_columns
                )
            _cohort_overview(
                cohort_df,
                image_col=image_col,
                area_col=area_col,
                mark_columns=mark_columns,
            )
        except Exception as exc:
            st.error(f"Could not load the full CSV for cohort analysis: {exc}")
        st.caption(f"Detected {len(image_ids):,} image(s) in the selected CSV.")

    with tab_image_preview:
        _section_title(
            "Image preview",
            "Inspect one image at a time before running KMM. Dots are cells; you can "
            "color/filter by cell type and draw area circles in the same coordinate "
            "scale as the x/y axes.",
        )
        preview_cols = st.columns(4)
        with preview_cols[0]:
            preview_image = st.selectbox(
                "Image",
                image_ids,
                key="preview_image_id",
                help="Image IDs come from the currently selected image ID column.",
            )
        with preview_cols[1]:
            cell_type_options = ["No cell type column"] + df_preview.columns.tolist()
            if st.session_state.get("preview_cell_type_col") not in cell_type_options:
                st.session_state.pop("preview_cell_type_col", None)
            default_cell_type_index = 1 + _default_cell_type_column(
                df_preview, df_preview.columns.tolist()
            )
            cell_type_col = st.selectbox(
                "Cell type column",
                cell_type_options,
                index=default_cell_type_index,
                key="preview_cell_type_col",
                help="Cells are colored and filtered by values in this selected CSV column. This does not need to be one of the selected KMM mark columns.",
            )
            cell_type_col = (
                None if cell_type_col == "No cell type column" else cell_type_col
            )
        with preview_cols[2]:
            show_area = st.checkbox(
                "Draw true-scale area circles",
                value=area_col is not None,
                disabled=area_col is None,
                key="preview_show_area",
                help="When enabled, each cell is drawn as a circle with radius sqrt(area / pi) in x/y coordinate units.",
            )
        with preview_cols[3]:
            max_preview_points = st.number_input(
                "Max preview cells",
                min_value=100,
                max_value=50000,
                value=5000,
                step=100,
                key="preview_max_points",
                help="Large images are sampled for faster plotting; analysis still uses all cells.",
            )
        try:
            image_frame, total_points = _load_image_preview_data(
                csv_path,
                x_col,
                y_col,
                image_col,
                area_col,
                cell_type_col,
                preview_image,
            )
            available_cell_types = _sort_label_values(
                image_frame["_cell_type"].unique()
            )
            existing_cell_types = st.session_state.get(
                "preview_cell_type_filter", available_cell_types
            )
            default_cell_types = [
                value for value in existing_cell_types if value in available_cell_types
            ] or available_cell_types
            if existing_cell_types != default_cell_types:
                st.session_state["preview_cell_type_filter"] = default_cell_types
            selected_cell_types = st.multiselect(
                "Cell types to show",
                available_cell_types,
                default=default_cell_types,
                key="preview_cell_type_filter",
                help="Clear a cell type here to hide it from the image preview.",
            )
            image_frame = image_frame[
                image_frame["_cell_type"].isin(selected_cell_types)
            ].copy()
            filtered_points = len(image_frame)
            if filtered_points > int(max_preview_points):
                image_frame = image_frame.sample(
                    n=int(max_preview_points), random_state=17
                ).copy()
            preview_xrange = [
                float(st.session_state.get("xrange_lo", x_min)),
                float(st.session_state.get("xrange_hi", x_max)),
            ]
            preview_yrange = [
                float(st.session_state.get("yrange_lo", y_min)),
                float(st.session_state.get("yrange_hi", y_max)),
            ]
            if preview_xrange[1] <= preview_xrange[0]:
                preview_xrange = [
                    float(image_frame["_x"].min()),
                    float(image_frame["_x"].max()),
                ]
            if preview_yrange[1] <= preview_yrange[0]:
                preview_yrange = [
                    float(image_frame["_y"].min()),
                    float(image_frame["_y"].max()),
                ]
            st.caption(
                f"Displaying {len(image_frame):,} of {filtered_points:,} selected cell(s) "
                f"from {total_points:,} total cell(s) for image `{preview_image}`."
            )
            if image_frame.empty:
                st.warning("No plottable cells were found for the selected image.")
            else:
                figure = _image_preview_figure(
                    image_frame,
                    x_col,
                    y_col,
                    area_col,
                    bool(show_area),
                    preview_xrange,
                    preview_yrange,
                )
                st.pyplot(figure, clear_figure=True, use_container_width=False)
        except Exception as exc:
            st.exception(exc)

    with tab_settings:
        _section_title(
            "Run settings",
            "Set the rectangular analysis window, cell-size correction, and CPU layout. "
            "The x/y ranges should cover the coordinates in the selected CSV.",
        )
        col1, col2 = st.columns(2)
        with col1:
            xrange_lo = st.number_input(
                "xrange min", value=float(x_min), key="xrange_lo"
            )
            xrange_hi = st.number_input(
                "xrange max", value=float(x_max), key="xrange_hi"
            )
        with col2:
            yrange_lo = st.number_input(
                "yrange min", value=float(y_min), key="yrange_lo"
            )
            yrange_hi = st.number_input(
                "yrange max", value=float(y_max), key="yrange_hi"
            )

        size_correction = st.checkbox(
            "Enable size correction with large-cell cutoff",
            value=area_col is not None,
            disabled=area_col is None,
            help="When enabled, the area cutoff below defines which cells are large occluders for pp_criterion.",
            key="size_correction",
        )
        pp_threshold = None
        if size_correction:
            area_values = pd.to_numeric(df_preview[area_col], errors="coerce").dropna()
            if area_values.empty:
                st.error("The selected area column has no numeric values.")
                st.stop()
            area_min = float(area_values.min())
            area_max = float(area_values.max())
            q90 = float(area_values.quantile(0.9))
            pp_threshold = st.slider(
                "Large-cell area cutoff",
                min_value=area_min,
                max_value=area_max,
                value=q90,
                step=max((area_max - area_min) / 200.0, 1e-6),
                help="Cells with area above this value are treated as large cells in pp_criterion.",
                key="pp_area_threshold",
            )
            chart = _area_distribution_chart(area_values, float(pp_threshold))
            if chart is not None:
                st.altair_chart(chart, use_container_width=True)

        c1, c2, c3 = st.columns(3)
        with c1:
            cpus_per_job = st.number_input(
                "CPUs available",
                min_value=1,
                value=int(
                    _scheduler_int(
                        ("SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS"),
                        max(1, min(len(image_ids), 8)),
                    )
                ),
                step=1,
                key="cpus_per_job",
            )
        with c2:
            threads_per_worker = st.number_input(
                "Threads per worker",
                min_value=1,
                value=1,
                step=1,
                key="threads_per_worker",
            )
        with c3:
            recommended_workers = _recommended_workers(
                int(cpus_per_job), 1, int(threads_per_worker)
            )
            max_workers = st.number_input(
                "Worker processes",
                min_value=1,
                value=int(min(recommended_workers, max(1, len(image_ids)))),
                step=1,
                key="max_workers",
            )

        default_out = (
            REPO_ROOT
            / "Result"
            / f"streamlit_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )
        save_folder = st.text_input(
            "Output folder", value=str(default_out), key="save_folder"
        )
        show_terminal_progress = st.checkbox(
            "Also print terminal progress", value=False, key="show_terminal_progress"
        )

        if st.button("Run KMM / markcorr", type="primary"):
            if not mark_columns:
                st.error("Select at least one mark column before running markcorr.")
                st.stop()
            if float(xrange_hi) <= float(xrange_lo) or float(yrange_hi) <= float(
                yrange_lo
            ):
                st.error(
                    "Invalid window range: max must be greater than min for both axes."
                )
                st.stop()
            os.makedirs(save_folder, exist_ok=True)
            progress = st.progress(0.0, text="Starting")

            def _cb(done, total, image_num):
                progress.progress(
                    done / total, text=f"image_{image_num} done ({done}/{total})"
                )

            pp_criterion = None
            if size_correction and pp_threshold is not None:
                threshold = float(pp_threshold)
                pp_criterion = (
                    lambda frame, _threshold=threshold: frame["Area"] > _threshold
                )  # noqa: E731

            started = time.time()
            try:
                _set_math_threads(int(threads_per_worker))
                runSpatial(
                    csv_path=csv_path,
                    savefolder=str(save_folder).rstrip("/") + "/",
                    xrange=[xrange_lo, xrange_hi],
                    yrange=[yrange_lo, yrange_hi],
                    sizeCorrection=size_correction,
                    pp_criterion=pp_criterion,
                    max_workers=int(max_workers),
                    progress_callback=_cb,
                    show_progress=show_terminal_progress,
                    x_col=x_col,
                    y_col=y_col,
                    image_col=image_col or "",
                    area_col=area_col or "",
                    mark_columns=mark_columns,
                )
            except Exception as exc:
                st.exception(exc)
            else:
                elapsed = time.time() - started
                st.session_state["last_savefolder"] = str(save_folder)
                run_config_path = _write_run_config(
                    save_folder,
                    {
                        "created_at": datetime.now().isoformat(timespec="seconds"),
                        "source": "streamlit",
                        "csv_path": str(csv_path),
                        "result_folder": str(save_folder),
                        "x_col": x_col,
                        "y_col": y_col,
                        "image_col": image_col,
                        "area_col": area_col,
                        "mark_columns": mark_columns,
                        "xrange": [float(xrange_lo), float(xrange_hi)],
                        "yrange": [float(yrange_lo), float(yrange_hi)],
                        "size_correction": bool(size_correction),
                        "pp_area_threshold": (
                            float(pp_threshold) if pp_threshold is not None else None
                        ),
                        "max_workers": int(max_workers),
                        "threads_per_worker": int(threads_per_worker),
                        "cpus_per_job": int(cpus_per_job),
                        "show_terminal_progress": bool(show_terminal_progress),
                    },
                )
                progress.progress(1.0, text="Done")
                st.success(
                    f"Finished in {elapsed:.1f}s. Results saved to `{save_folder}`."
                )
                st.caption(f"Saved KMM run settings to `{run_config_path}`.")

    with tab_hpc:
        _section_title(
            "HPC / terminal",
            "These commands are for non-Streamlit use. Replace placeholder paths with "
            "paths visible from the compute node, then paste the command into your shell "
            "or scheduler template.",
        )
        st.caption(
            "Use these commands when Streamlit is not available or when splitting work across scheduler array jobs."
        )
        hpc_project_dir = st.text_input(
            "Project directory on the cluster",
            value="/path/to/ClumPyCells",
            key="hpc_project_dir",
        )
        stable_csv_path = st.text_input(
            "Stable CSV path visible from compute nodes",
            value="/path/to/cells.csv",
            key="hpc_stable_csv_path",
        )
        hpc_output_folder = st.text_input(
            "Output folder on the cluster",
            value="results_clumpycells",
            key="hpc_output_folder",
        )
        hpc_col1, hpc_col2, hpc_col3, hpc_col4 = st.columns(4)
        with hpc_col1:
            hpc_cpus = st.number_input(
                "CPUs per array task",
                min_value=1,
                value=max(1, int(_scheduler_int(("SLURM_CPUS_PER_TASK",), 8))),
                step=1,
                key="hpc_cpus",
            )
        with hpc_col2:
            hpc_workers = st.number_input(
                "Processes per task",
                min_value=1,
                value=max(1, min(int(hpc_cpus), len(image_ids))),
                step=1,
                key="hpc_workers",
            )
        with hpc_col3:
            hpc_threads = st.number_input(
                "Threads per process", min_value=1, value=1, step=1, key="hpc_threads"
            )
        with hpc_col4:
            images_per_job = st.number_input(
                "Images per array task",
                min_value=1,
                value=max(1, int(hpc_workers)),
                step=1,
                key="hpc_images_per_job",
            )
        hpc_mem = st.number_input(
            "Memory per task (GB)", min_value=1, value=32, step=1, key="hpc_mem_gb"
        )
        hpc_hours = st.number_input(
            "Wall time per task (hours)", min_value=1, value=4, step=1, key="hpc_hours"
        )

        hpc_chunks = _chunks(image_ids, int(images_per_job))
        st.metric("Array tasks needed", len(hpc_chunks))
        base_cmd = _cli_command(
            stable_csv_path,
            hpc_output_folder,
            x_col,
            y_col,
            image_col,
            area_col,
            mark_columns,
            [xrange_lo, xrange_hi],
            [yrange_lo, yrange_hi],
            size_correction,
            pp_threshold,
            hpc_workers,
            hpc_threads,
            hpc_cpus,
            1,
        )
        st.subheader("Single terminal command")
        st.code(base_cmd, language="bash")

        st.subheader("Scheduler-neutral image chunks")
        chunk_table = pd.DataFrame(
            {
                "task_id": list(range(len(hpc_chunks))),
                "image_ids": [",".join(chunk) for chunk in hpc_chunks],
                "command_suffix": [
                    f"--image-id {','.join(chunk)}" for chunk in hpc_chunks
                ],
            }
        )
        st.dataframe(chunk_table, use_container_width=True)

        st.subheader("Slurm array template")
        slurm_base = base_cmd.replace("\n", " && ")
        st.code(
            _slurm_array_script(
                slurm_base,
                hpc_chunks,
                int(hpc_cpus),
                int(hpc_mem),
                int(hpc_hours),
                project_dir=hpc_project_dir,
            ),
            language="bash",
        )

else:
    st.header("Downstream analysis")
    st.caption(
        "Open existing KMM/markcorr outputs for visualization, statistics, permutation tests, model fitting, and survival analysis."
    )

    result_folder = st.text_input(
        "KMM / markcorr result folder",
        value=st.session_state.get("last_savefolder") or "",
        key="downstream_result_folder",
        help="Folder containing image_<id>/iso.csv and r.csv outputs from KMM / markcorr.",
    )
    image_numbers = _list_image_dirs(result_folder) if result_folder else []
    folder_ready = bool(result_folder and image_numbers)
    if folder_ready:
        st.success(f"Detected {len(image_numbers)} image(s).")
    elif result_folder:
        st.warning("No image_<id>/iso.csv files were found in this folder.")
    else:
        st.info("Paste a result folder to enable result-based downstream tools.")

    axis_keys = _peek_axis_keys(result_folder) if folder_ready else []
    axis_name = {key: key for key in axis_keys}
    if axis_keys:
        with st.expander("Display labels for marks", expanded=False):
            for key in axis_keys:
                axis_name[key] = st.text_input(key, value=key, key=f"axis_{key}")
    groups = _image_groups_controls(image_numbers, "downstream") if folder_ready else {}

    tab_curves, tab_auc, tab_tree, tab_survival, tab_results = st.tabs(
        [
            "KMM curves",
            "AUC and groups",
            "Decision tree",
            "Survival",
            "Results",
        ]
    )

    with tab_curves:
        if not folder_ready:
            st.info("Provide a result folder first.")
        else:
            pair_columns = _list_pair_columns(result_folder)
            c1, c2, c3 = st.columns(3)
            with c1:
                curve_image = st.selectbox("Image", image_numbers, key="curve_image")
            with c2:
                curve_pair = st.selectbox("KMM pair", pair_columns, key="curve_pair")
            with c3:
                curve_transform = st.selectbox(
                    "Y transform", ["log2", "raw"], key="curve_transform"
                )
            if curve_pair:
                try:
                    curve_chart = _kmm_curve_chart(
                        result_folder, curve_image, curve_pair, curve_transform
                    )
                    st.altair_chart(curve_chart, use_container_width=True)
                    curve_cache = _cache_chart(
                        curve_chart,
                        Path(result_folder) / "streamlit_cache",
                        f"kmm_curve_{curve_image}_{curve_pair}_{curve_transform}",
                    )
                    if curve_cache is not None:
                        st.caption(f"Cached KMM curve: `{curve_cache}`")
                except Exception as exc:
                    st.exception(exc)

    with tab_auc:
        _section_title(
            "AUC and group analysis",
            "AUC condenses each KMM curve into one value per image and mark pair. "
            "Group comparison then tests whether those values differ between image groups.",
        )
        if not folder_ready:
            st.info("Provide a result folder first.")
        elif not groups:
            st.info("Define image groups above to run AUC and group comparison.")
        else:
            folder_norm = str(result_folder).rstrip("/") + "/"
            chart_cache_folder = Path(folder_norm) / "streamlit_cache"
            result = MarkcorrResult(
                groups=groups, resultFolder=folder_norm, axisName=axis_name
            )
            norm = st.selectbox(
                "AUC normalization", ["min_mid_max", "log"], index=0, key="auc_norm"
            )
            min_valid_percent = st.slider(
                "Minimum valid-image percentage for heatmap cells",
                min_value=0,
                max_value=100,
                value=50,
                step=5,
                key="auc_min_valid_percent",
                help="Cells with valid KMM/AUC values in fewer images than this threshold are marked with X.",
            )
            auc_context = json.dumps(
                {
                    "folder": folder_norm,
                    "groups": groups,
                    "axis_name": axis_name,
                    "norm": norm,
                    "min_valid_percent": int(min_valid_percent),
                },
                sort_keys=True,
            )
            has_current_auc = st.session_state.get("_auc_context") == auc_context
            sub_auc, sub_diff, sub_box, sub_perm = st.tabs(
                ["AUC heatmap", "Group comparison", "Box plot", "Permutation"]
            )
            with sub_auc:
                if st.button("Compute AUC heatmaps", key="btn_auc"):
                    with st.spinner("Computing AUC"):
                        auc, plots = result.getAUC(
                            norm=norm,
                            plot=True,
                            min_nanPercentile=float(min_valid_percent) / 100.0,
                        )
                    st.session_state["_auc"] = auc
                    st.session_state["_auc_plots"] = plots
                    st.session_state["_auc_context"] = auc_context
                    has_current_auc = True
                    auc_cache_files = {}
                    for group_name, chart in plots.items():
                        cache_path = _cache_chart(
                            chart,
                            chart_cache_folder,
                            f"auc_heatmap_{group_name}",
                        )
                        if cache_path is not None:
                            auc_cache_files[group_name] = str(cache_path)
                    st.session_state["_auc_cache_files"] = auc_cache_files
                if "_auc_plots" in st.session_state and has_current_auc:
                    for group_name, chart in st.session_state["_auc_plots"].items():
                        st.markdown(
                            f"#### {group_name} (n = {len(groups[group_name])})"
                        )
                        st.altair_chart(chart, use_container_width=False)
                    cache_files = st.session_state.get("_auc_cache_files") or {}
                    if cache_files:
                        st.caption(
                            "Cached heatmap HTML files: "
                            + ", ".join(path for path in cache_files.values() if path)
                        )
                    x_marker_count = _auc_x_marker_count(
                        st.session_state.get("_auc"), groups, min_valid_percent
                    )
                    if x_marker_count > 0:
                        st.info(
                            f"In these AUC heatmaps, `X` marks {x_marker_count} cell(s) "
                            f"where fewer than {min_valid_percent}% of images had valid "
                            "values for that mark pair. Adjust the threshold above if needed."
                        )
                elif "_auc_plots" in st.session_state:
                    st.info(
                        "AUC results are cached for a different folder or setting. "
                        "Compute AUC heatmaps for the current configuration."
                    )
            with sub_diff:
                if "_auc" not in st.session_state or not has_current_auc:
                    st.info("Compute AUC heatmaps first.")
                else:
                    group_names = list(groups.keys())
                    if len(group_names) < 2:
                        st.info("Need at least two groups to compare.")
                    else:
                        c1, c2, c3 = st.columns(3)
                        with c1:
                            g1 = st.selectbox("Group A", group_names, index=0)
                        with c2:
                            g2 = st.selectbox("Group B", group_names, index=1)
                        with c3:
                            method = st.selectbox(
                                "Test", ["MW", "perm"], index=0, key="diff_test_method"
                            )
                        diff_context = json.dumps(
                            {
                                "auc_context": auc_context,
                                "group_a": g1,
                                "group_b": g2,
                                "method": method,
                            },
                            sort_keys=True,
                        )
                        if g1 != g2 and st.button(
                            "Run group comparison", key="btn_diff"
                        ):
                            out_csv = Path(folder_norm) / f"diff_{g1}_vs_{g2}.csv"
                            chart = MarkcorrResult.find_diff(
                                st.session_state["_auc"],
                                g1,
                                g2,
                                method=method,
                                axisName=axis_name,
                                saveCsv=str(out_csv),
                            )
                            st.session_state["_diff_chart"] = chart
                            st.session_state["_diff_csv"] = str(out_csv)
                            st.session_state["_diff_context"] = diff_context
                            diff_cache_path = _cache_chart(
                                chart,
                                chart_cache_folder,
                                f"diff_heatmap_{g1}_vs_{g2}_{method}",
                            )
                            st.session_state["_diff_cache_file"] = (
                                str(diff_cache_path)
                                if diff_cache_path is not None
                                else None
                            )
                        if (
                            "_diff_chart" in st.session_state
                            and st.session_state.get("_diff_context") == diff_context
                        ):
                            st.altair_chart(
                                st.session_state["_diff_chart"],
                                use_container_width=False,
                            )
                            diff_csv = st.session_state.get("_diff_csv")
                            if diff_csv and Path(diff_csv).exists():
                                st.dataframe(
                                    _format_display_table(
                                        pd.read_csv(diff_csv, index_col=0)
                                    ),
                                    use_container_width=True,
                                )
                            diff_cache = st.session_state.get("_diff_cache_file")
                            if diff_cache:
                                st.caption(
                                    f"Cached group-comparison heatmap: `{diff_cache}`"
                                )
            with sub_box:
                if "_auc" not in st.session_state or not has_current_auc:
                    st.info("Compute AUC heatmaps first.")
                else:
                    chosen = st.multiselect(
                        "Groups for box plot",
                        list(groups.keys()),
                        default=list(groups.keys()),
                        key="boxplot_groups",
                    )
                    box_context = json.dumps(
                        {"auc_context": auc_context, "groups": chosen}, sort_keys=True
                    )
                    if chosen and st.button("Render box plot", key="btn_box"):
                        box_chart = result.getBoxPlot(
                            st.session_state["_auc"], chosen, axisName=axis_name
                        )
                        st.session_state["_box_chart"] = box_chart
                        st.session_state["_box_context"] = box_context
                        box_cache_path = _cache_chart(
                            box_chart,
                            chart_cache_folder,
                            f"boxplot_{'_vs_'.join(chosen)}",
                        )
                        st.session_state["_box_cache_file"] = (
                            str(box_cache_path) if box_cache_path is not None else None
                        )
                    if (
                        "_box_chart" in st.session_state
                        and st.session_state.get("_box_context") == box_context
                    ):
                        st.altair_chart(
                            st.session_state["_box_chart"], use_container_width=True
                        )
                        box_cache = st.session_state.get("_box_cache_file")
                        if box_cache:
                            st.caption(f"Cached box plot: `{box_cache}`")
            with sub_perm:
                _render_permutation_panel(result_folder, image_numbers)

    with tab_tree:
        _section_title(
            "Decision tree",
            "Train a classifier from AUC features or an uploaded feature table, then preview "
            "the dtreeviz SVG outputs saved by the analysis.",
        )
        tree_source = st.radio(
            "Training data source",
            [
                "Use image groups from this result folder",
                "Upload feature table",
                "Use demo feature table",
            ],
            horizontal=True,
            key="dt_training_data_source",
        )
        tree_out = st.text_input(
            "Decision tree output folder",
            value=str(REPO_ROOT / "Result" / "streamlit_models" / "decision_tree"),
        )
        c1, c2, c3 = st.columns(3)
        with c1:
            dt_use_bo = st.checkbox(
                "Use Bayesian optimization", value=False, key="dt_use_bo"
            )
            dt_impute = st.checkbox(
                "Impute missing values", value=True, key="dt_impute"
            )
        with c2:
            dt_max_depth = st.number_input(
                "Base max depth", min_value=1, max_value=30, value=4, step=1
            )
            dt_max_features = st.slider(
                "Base max features",
                min_value=0.05,
                max_value=1.0,
                value=0.67,
                step=0.01,
            )
        with c3:
            dt_bo_init = st.number_input(
                "BO initial points", min_value=1, max_value=500, value=25, step=1
            )
            dt_bo_iter = st.number_input(
                "BO iterations", min_value=1, max_value=500, value=10, step=1
            )
        if tree_source == "Use image groups from this result folder":
            dt_norm = st.selectbox(
                "Feature normalization", ["min_mid_max", "log"], index=0, key="dt_norm"
            )
            if st.button(
                "Train decision tree", type="primary", key="train_tree_groups"
            ):
                if not folder_ready or not groups:
                    st.error("Provide a result folder and define image groups first.")
                else:
                    from ClumPyCells.Analysis.decisionTree import (
                        decision_tree_from_markcorr_groups,
                    )

                    os.makedirs(tree_out, exist_ok=True)
                    with st.spinner("Training decision tree"):
                        decision_tree_from_markcorr_groups(
                            groups=groups,
                            resultFolder=str(result_folder).rstrip("/") + "/",
                            axisName=axis_name,
                            norm=dt_norm,
                            saveFolder=str(tree_out).rstrip("/") + "/",
                            bo=dt_use_bo,
                            impute=dt_impute,
                            save_fig=True,
                            max_depth=int(dt_max_depth),
                            max_features=float(dt_max_features),
                            bo_init_points=int(dt_bo_init),
                            bo_n_iter=int(dt_bo_iter),
                        )
                    st.session_state["last_decision_tree_folder"] = str(tree_out)
                    st.success(f"Decision tree outputs saved under `{tree_out}`.")
        elif tree_source == "Use demo feature table":
            feature_df = pd.read_csv(DEMO_FEATURE_TABLE_CSV)
            st.success(f"Using bundled demo feature table: `{DEMO_FEATURE_TABLE_CSV}`")
            st.dataframe(feature_df.head(20), use_container_width=True)
            target_col = st.selectbox(
                "Target/label column",
                feature_df.columns,
                index=_default_index(
                    feature_df.columns, ["phenotype", "target", "label", "group"]
                ),
                key="dt_demo_target_col",
            )
            feature_cols = st.multiselect(
                "Feature columns",
                [column for column in feature_df.columns if column != target_col],
                default=[
                    column for column in feature_df.columns if column != target_col
                ],
                key="dt_demo_feature_cols",
            )
            if st.button(
                "Train decision tree", type="primary", key="train_tree_demo_upload"
            ):
                from ClumPyCells.Analysis.decisionTree import (
                    decision_tree_from_feature_table,
                )

                os.makedirs(tree_out, exist_ok=True)
                decision_tree_from_feature_table(
                    feature_df,
                    target_col=target_col,
                    saveFolder=str(tree_out).rstrip("/") + "/",
                    feature_columns=feature_cols,
                    bo=dt_use_bo,
                    impute=dt_impute,
                    save_fig=True,
                    max_depth=int(dt_max_depth),
                    max_features=float(dt_max_features),
                    bo_init_points=int(dt_bo_init),
                    bo_n_iter=int(dt_bo_iter),
                )
                st.session_state["last_decision_tree_folder"] = str(tree_out)
                st.success(f"Decision tree outputs saved under `{tree_out}`.")
        else:
            feature_upload = st.file_uploader(
                "Feature table CSV", type=["csv"], key="dt_feature_table"
            )
            if feature_upload is not None:
                feature_path = (
                    Path(tempfile.gettempdir()) / f"dt_features_{feature_upload.name}"
                )
                feature_path.write_bytes(feature_upload.getvalue())
                feature_df = pd.read_csv(feature_path)
                st.dataframe(feature_df.head(20), use_container_width=True)
                target_col = st.selectbox(
                    "Target/label column",
                    feature_df.columns,
                    index=_default_index(
                        feature_df.columns, ["phenotype", "target", "label", "group"]
                    ),
                )
                feature_cols = st.multiselect(
                    "Feature columns",
                    [column for column in feature_df.columns if column != target_col],
                    default=[
                        column for column in feature_df.columns if column != target_col
                    ],
                )
                if st.button(
                    "Train decision tree", type="primary", key="train_tree_upload"
                ):
                    from ClumPyCells.Analysis.decisionTree import (
                        decision_tree_from_feature_table,
                    )

                    os.makedirs(tree_out, exist_ok=True)
                    decision_tree_from_feature_table(
                        feature_df,
                        target_col=target_col,
                        saveFolder=str(tree_out).rstrip("/") + "/",
                        feature_columns=feature_cols,
                        bo=dt_use_bo,
                        impute=dt_impute,
                        save_fig=True,
                        max_depth=int(dt_max_depth),
                        max_features=float(dt_max_features),
                        bo_init_points=int(dt_bo_init),
                        bo_n_iter=int(dt_bo_iter),
                    )
                    st.session_state["last_decision_tree_folder"] = str(tree_out)
                    st.success(f"Decision tree outputs saved under `{tree_out}`.")
        tree_view = st.text_input(
            "Decision tree results folder to display",
            value=st.session_state.get("last_decision_tree_folder") or tree_out,
            key="tree_view_folder",
        )
        if Path(tree_view, "tree.svg").exists():
            st.markdown("#### dtreeviz decision tree")
            _render_file_preview(Path(tree_view) / "tree.svg")
        if Path(tree_view, "largest_leaf_path.svg").exists():
            st.markdown("#### dtreeviz largest-leaf path")
            _render_file_preview(Path(tree_view) / "largest_leaf_path.svg")
        _artifact_browser(
            tree_view,
            "decision tree",
            "dt",
            ["*.svg", "*.csv", "*.html", "*.json", "*.txt"],
            preferred_names=["tree.svg", "largest_leaf_path.svg"],
        )

    with tab_survival:
        _section_title(
            "Survival analysis",
            "Join image-level spatial features with clinical data. The event/status column "
            "should identify whether the endpoint occurred for each sample.",
        )
        if not folder_ready:
            st.info(
                "Provide a result folder first so spatial features can be joined to clinical data."
            )
        clinical_source = st.radio(
            "Clinical data source",
            ["Upload clinical CSV", "Use clinical CSV path", "Use demo clinical CSV"],
            horizontal=True,
            key="surv_clinical_source",
        )
        clinical_path = ""
        if clinical_source == "Upload clinical CSV":
            clinical_upload = st.file_uploader(
                "Clinical CSV", type=["csv"], key="surv_clinical_upload"
            )
            if clinical_upload is not None:
                clinical_tmp = (
                    Path(tempfile.gettempdir()) / f"clinical_{clinical_upload.name}"
                )
                clinical_tmp.write_bytes(clinical_upload.getvalue())
                clinical_path = str(clinical_tmp)
        elif clinical_source == "Use clinical CSV path":
            clinical_path = st.text_input("Clinical CSV path", key="surv_clinical_path")
        else:
            clinical_path = str(DEMO_CLINICAL_CSV)
            st.session_state["surv_clinical_path"] = clinical_path
            st.success(f"Using bundled demo clinical CSV: `{clinical_path}`")

        if clinical_path and Path(clinical_path).exists():
            clinical_df = pd.read_csv(clinical_path, nrows=2000)
            st.dataframe(clinical_df.head(20), use_container_width=True)
            clinical_cols = clinical_df.columns.tolist()
            clinical_id_col = st.selectbox(
                "Clinical/sample ID column", clinical_cols, index=0
            )
            duration_col = st.selectbox(
                "Survival time column",
                clinical_cols,
                index=_default_index(
                    clinical_cols, ["OST", "OS", "survival", "time", "duration"]
                ),
            )
            event_col = st.selectbox(
                "Event/status column",
                clinical_cols,
                index=_default_index(
                    clinical_cols,
                    ["survival_status", "status", "event", "dead", "relapse"],
                ),
            )
            event_values = ["Auto infer"] + [
                str(value)
                for value in pd.Series(clinical_df[event_col].dropna().unique()).head(
                    20
                )
            ]
            event_positive = st.selectbox("Value meaning event occurred", event_values)
            event_positive_value = (
                None if event_positive == "Auto infer" else event_positive
            )
            covariate_options = [
                column
                for column in clinical_cols
                if column not in {clinical_id_col, duration_col, event_col}
            ]
            covariates = st.multiselect(
                "Clinical covariates to include", covariate_options, default=[]
            )

            with st.expander("Image to clinical ID mapping", expanded=False):
                st.caption(
                    "Leave this empty if image IDs in the result folder match the clinical/sample ID column."
                )
                mapping_source = st.radio(
                    "Mapping data source",
                    ["No mapping", "Upload mapping CSV", "Use demo mapping CSV"],
                    horizontal=True,
                    key="surv_mapping_source",
                )
                mapping_path = None
                mapping_image_col = None
                mapping_clinical_col = None
                if mapping_source == "Upload mapping CSV":
                    mapping_upload = st.file_uploader(
                        "Optional mapping CSV", type=["csv"], key="surv_mapping_upload"
                    )
                else:
                    mapping_upload = None
                if mapping_source == "Use demo mapping CSV":
                    mapping_path = str(DEMO_IMAGE_MAPPING_CSV)
                    mapping_df = pd.read_csv(mapping_path, nrows=2000)
                    st.success(f"Using bundled demo mapping CSV: `{mapping_path}`")
                    st.dataframe(mapping_df.head(20), use_container_width=True)
                    mapping_image_col = st.selectbox(
                        "Mapping image ID column",
                        mapping_df.columns,
                        index=_default_index(
                            mapping_df.columns, ["ImageNum", "image_id", "image"]
                        ),
                        key="mapping_image_col",
                    )
                    mapping_clinical_col = st.selectbox(
                        "Mapping clinical ID column",
                        mapping_df.columns,
                        index=_default_index(
                            mapping_df.columns,
                            ["sample_id", "patient_id", "clinical_id"],
                        ),
                        key="mapping_clinical_col",
                    )
                elif mapping_upload is not None:
                    mapping_tmp = (
                        Path(tempfile.gettempdir()) / f"mapping_{mapping_upload.name}"
                    )
                    mapping_tmp.write_bytes(mapping_upload.getvalue())
                    mapping_path = str(mapping_tmp)
                    mapping_df = pd.read_csv(mapping_path, nrows=2000)
                    st.dataframe(mapping_df.head(20), use_container_width=True)
                    mapping_image_col = st.selectbox(
                        "Mapping image ID column",
                        mapping_df.columns,
                        index=_default_index(
                            mapping_df.columns, ["ImageNum", "image_id", "image"]
                        ),
                        key="mapping_image_col",
                    )
                    mapping_clinical_col = st.selectbox(
                        "Mapping clinical ID column",
                        mapping_df.columns,
                        index=_default_index(
                            mapping_df.columns,
                            ["sample_id", "patient_id", "clinical_id"],
                        ),
                        key="mapping_clinical_col",
                    )

            surv_out = st.text_input(
                "Survival output folder",
                value=str(REPO_ROOT / "Result" / "streamlit_models" / "survival"),
            )
            sv_col1, sv_col2, sv_col3 = st.columns(3)
            with sv_col1:
                sv_norm = st.selectbox(
                    "Spatial feature normalization",
                    ["min_mid_max", "log"],
                    index=0,
                    key="sv_norm",
                )
                sv_top_n = st.number_input(
                    "Top N KM features", min_value=1, max_value=200, value=10, step=1
                )
            with sv_col2:
                sv_km_patient = st.checkbox("Generate KM by patient/sample", value=True)
                sv_km_roi = st.checkbox("Generate KM by image/ROI", value=True)
            with sv_col3:
                sv_penalizer = st.number_input(
                    "Cox penalizer", min_value=0.0, value=0.0, step=0.01
                )
                sv_min_var = st.number_input(
                    "Minimum feature variance",
                    min_value=0.0,
                    value=0.0,
                    step=0.0001,
                    format="%.4f",
                )
            if st.button("Run survival analysis", type="primary"):
                if not folder_ready:
                    st.error("Provide a KMM result folder first.")
                else:
                    from ClumPyCells.Analysis.survivalAnalysis import (
                        run_user_survival_analysis,
                    )

                    os.makedirs(surv_out, exist_ok=True)
                    with st.spinner("Running survival analysis"):
                        run_user_survival_analysis(
                            result_folder=str(result_folder).rstrip("/") + "/",
                            clinical_csv=clinical_path,
                            clinical_id_col=clinical_id_col,
                            duration_col=duration_col,
                            event_col=event_col,
                            saveFolder=str(surv_out).rstrip("/") + "/",
                            event_positive_value=event_positive_value,
                            covariate_cols=covariates,
                            image_numbers=image_numbers,
                            image_to_clinical_csv=mapping_path,
                            mapping_image_col=mapping_image_col,
                            mapping_clinical_col=mapping_clinical_col,
                            norm=sv_norm,
                            km_top_n=int(sv_top_n),
                            run_km_by_patient=bool(sv_km_patient),
                            run_km_by_roi=bool(sv_km_roi),
                            penalizer=float(sv_penalizer),
                            min_feature_variance=float(sv_min_var),
                        )
                    st.session_state["last_survival_folder"] = str(surv_out)
                    st.success(f"Survival outputs saved under `{surv_out}`.")
        elif clinical_path:
            st.error("The clinical CSV path does not exist.")
        else:
            st.info("Provide clinical data to enable survival analysis.")

        survival_view = st.text_input(
            "Survival results folder to display",
            value=st.session_state.get("last_survival_folder")
            or str(REPO_ROOT / "Result" / "streamlit_models" / "survival"),
            key="survival_view_folder",
        )
        _artifact_browser(survival_view, "survival", "survival")
        st.markdown("#### KM curve gallery")
        km_files = [
            path
            for path in _find_artifacts(survival_view, ["*.svg", "*.png"])
            if "km" in path.name.lower() or "km" in str(path.parent).lower()
        ]
        if km_files:
            selected_km = st.multiselect(
                "KM curves to display",
                km_files,
                default=km_files[: min(4, len(km_files))],
                format_func=lambda path: str(path.relative_to(Path(survival_view))),
            )
            for km_path in selected_km:
                st.markdown(f"**{km_path.name}**")
                _render_file_preview(km_path)
        else:
            st.info("No KM curve images found yet.")

    with tab_results:
        browse_folder = st.text_input(
            "Folder to browse",
            value=result_folder or str(REPO_ROOT / "Result"),
            key="generic_browse_folder",
        )
        _artifact_browser(browse_folder, "result", "generic")
