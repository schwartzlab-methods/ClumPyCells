"""Streamlit UI for the ClumPyCells workflow.

Run with::

    streamlit run streamlit_app.py

The app walks the user through the canonical ClumPyCells pipeline:

1. **Run markcorr** — upload a CSV, configure the observation window, optionally
   enable cell-size correction (with an occluder threshold), and launch
   ``runSpatial``. Per-image progress is reported via ``st.progress`` and the
   raw ``iso.csv`` artefacts are written to a chosen output folder.

2. **Downstream analysis** — point at any folder produced in step 1, define
   the groups to compare, and run ``MarkcorrResult.getAUC`` /
   ``MarkcorrResult.find_diff`` / ``getBoxPlot``. Both Mann-Whitney with BH
   FDR and the Fisher-Pitman permutation test are exposed.

3. **Permutation test** — a within-image label-shuffling null distribution
   for one or more images, mirroring ``runPermutation.py``.
"""

from __future__ import annotations

import io
import os
import shlex
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

# Make the package importable regardless of where streamlit is launched from.
REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ClumPyCells.Analysis.markcorrResult import MarkcorrResult  # noqa: E402
from ClumPyCells.ClumPyCells import prepare_cell_table, runSpatial  # noqa: E402

st.set_page_config(page_title="ClumPyCells", layout="wide")


# ---------------------------------------------------------------------------
# Session state defaults
# ---------------------------------------------------------------------------
def _ensure_state():
    defaults = {
        "uploaded_csv_path": None,
        "csv_columns": [],
        "image_numbers": [],
        "last_savefolder": None,
        "groups": {},
        "axis_name": {},
    }
    for k, v in defaults.items():
        st.session_state.setdefault(k, v)


_ensure_state()


# ---------------------------------------------------------------------------
# Sidebar navigation
# ---------------------------------------------------------------------------
st.sidebar.title("ClumPyCells")
page = st.sidebar.radio(
    "Workflow step",
    [
        "1. Run markcorr",
        "2. Downstream analysis",
        "3. Permutation test",
    ],
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _list_image_dirs(folder: str) -> list[int]:
    p = Path(folder)
    if not p.is_dir():
        return []
    nums = []
    for child in p.iterdir():
        if child.is_dir() and child.name.startswith("image_"):
            raw = child.name.split("_", 1)[1]
            try:
                nums.append(int(raw))
            except ValueError:
                nums.append(raw)
    return sorted(nums, key=lambda value: str(value))


def _peek_axis_keys(folder: str) -> list[str]:
    """Return the mark-pair column names (e.g. ``'CD3 vs. CD20'``) found in
    the first iso.csv inside ``folder``."""
    p = Path(folder)
    for child in sorted(p.glob("image_*/iso.csv")):
        try:
            df = pd.read_csv(child, nrows=1)
        except Exception:
            continue
        cols = [c for c in df.columns if "vs." in c]
        keys = sorted(
            {c.split(" vs. ")[0] for c in cols} | {c.split(" vs. ")[1] for c in cols}
        )
        return keys
    return []


NO_IMAGE_OPTION = "No image column: treat the whole CSV as one image"
NO_AREA_OPTION = "No area column"


def _default_index(options, candidates, fallback=0):
    lowered = {str(option).lower(): i for i, option in enumerate(options)}
    for candidate in candidates:
        idx = lowered.get(candidate.lower())
        if idx is not None:
            return idx
    return fallback


def _hpc_default_workers():
    for env_name in ("SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS"):
        try:
            value = int(os.environ.get(env_name, ""))
        except ValueError:
            continue
        if value > 0:
            return value
    return 1


def _column_mapping_controls(df_preview: pd.DataFrame, key_prefix: str):
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
        )
    with c2:
        y_col = st.selectbox(
            "y coordinate column",
            columns,
            index=_default_index(columns, ["y", "Y", "centroid_y", "CenterY"]),
            key=f"{key_prefix}_y_col",
        )
    with c3:
        image_options = [NO_IMAGE_OPTION] + columns
        image_col_choice = st.selectbox(
            "image ID column",
            image_options,
            index=_default_index(
                image_options,
                ["ImageNum", "image", "image_id", "ImageID", "sample_id"],
                fallback=0,
            ),
            key=f"{key_prefix}_image_col",
        )
        image_col = None if image_col_choice == NO_IMAGE_OPTION else image_col_choice
    with c4:
        area_options = [NO_AREA_OPTION] + columns
        area_col_choice = st.selectbox(
            "area column",
            area_options,
            index=_default_index(area_options, ["Area", "area", "cell_area"]),
            key=f"{key_prefix}_area_col",
        )
        area_col = None if area_col_choice == NO_AREA_OPTION else area_col_choice

    reserved = {x_col, y_col, image_col, area_col, "Unnamed: 0"}
    reserved.discard(None)
    default_marks = [column for column in columns if column not in reserved]
    selected_marks = st.multiselect(
        "Columns to use as marks",
        columns,
        default=default_marks,
        key=f"{key_prefix}_mark_columns",
        help="Only selected columns become markcorr features; coordinate, image, and area columns are excluded automatically.",
    )
    mark_columns = [column for column in selected_marks if column not in reserved]
    if len(mark_columns) != len(selected_marks):
        st.info("Coordinate, image ID, and area columns are not used as marks.")
    if not mark_columns:
        st.warning("Select at least one mark column before running markcorr.")

    return x_col, y_col, image_col, area_col, mark_columns


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
    ]
    for column in mark_columns:
        parts.extend(["--mark", str(column)])
    if size_correction:
        parts.append("--size-correction")
    if pp_threshold is not None:
        parts.extend(["--pp-area-threshold", str(float(pp_threshold))])
    return " ".join(shlex.quote(part) for part in parts)


# ===========================================================================
# Page 1 — Run markcorr
# ===========================================================================
if page.startswith("1"):
    st.header("Step 1 · Run mark cross-correlation")
    st.caption(
        "Upload a CSV (one row per cell), map its coordinate/image columns, "
        "choose which columns become marks, and set the amount of parallel "
        "work to use on this machine or HPC allocation."
    )

    uploaded = st.file_uploader("Cell-table CSV", type=["csv"])
    if uploaded is not None:
        # Persist to a temp file so runSpatial (which uses pandas chunks) can
        # re-read it without keeping the entire frame in memory twice.
        tmp_csv = Path(tempfile.gettempdir()) / f"cpc_input_{uploaded.name}"
        tmp_csv.write_bytes(uploaded.getvalue())
        st.session_state["uploaded_csv_path"] = str(tmp_csv)
        df_preview = pd.read_csv(tmp_csv, nrows=2000)
        st.session_state["csv_columns"] = df_preview.columns.tolist()
        st.write("Preview (first 2000 rows):")
        st.dataframe(df_preview.head(20), use_container_width=True)

        st.subheader("Input columns")
        x_col, y_col, image_col, area_col, mark_columns = _column_mapping_controls(
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

        if image_col is None:
            n_images_preview = 1
        else:
            n_images_preview = int(df_preview[image_col].nunique())
        st.caption(
            f"Preview contains {n_images_preview} image(s). "
            "If no image ID column is selected, the whole CSV is processed as image_1."
        )

        col1, col2 = st.columns(2)
        with col1:
            xrange_lo = st.number_input("xrange min", value=float(x_min))
            xrange_hi = st.number_input("xrange max", value=float(x_max))
        with col2:
            yrange_lo = st.number_input("yrange min", value=float(y_min))
            yrange_hi = st.number_input("yrange max", value=float(y_max))

        size_correction = st.checkbox(
            "Enable cell-size correction",
            value=area_col is not None,
            disabled=area_col is None,
            help=(
                "Subtract the radii of the two cells from every centre-to-"
                "centre distance, and additionally remove the chord of any "
                "large *occluder* cell crossing the segment between them. "
                "Requires selecting an area column."
            ),
        )

        pp_threshold = None
        if size_correction:
            use_pp = st.checkbox(
                "Treat large cells as occluders (pp_criterion)",
                value=False,
                help=(
                    "Cells whose selected area exceeds the threshold below "
                    "will additionally subtract their chord from every "
                    "intersecting pairwise segment."
                ),
            )
            if use_pp:
                pp_threshold = st.number_input(
                    "Occluder area threshold",
                    value=float(pd.to_numeric(df_preview[area_col]).quantile(0.9)),
                    min_value=0.0,
                )

        st.subheader("Compute resources")
        max_workers = st.number_input(
            "Parallel image workers",
            min_value=1,
            value=max(1, min(_hpc_default_workers(), max(1, n_images_preview))),
            step=1,
            help=(
                "Set this to the number of CPU cores allocated by your scheduler. "
                "Use 1 for serial execution or memory-constrained jobs."
            ),
        )
        show_terminal_progress = st.checkbox(
            "Also print tqdm progress in the terminal",
            value=False,
            help="Useful for HPC logs; the Streamlit progress bar remains visible either way.",
        )

        default_out = (
            REPO_ROOT
            / "Result"
            / f"streamlit_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )
        save_folder = st.text_input("Output folder", value=str(default_out))

        with st.expander("Terminal / HPC equivalent command", expanded=False):
            st.code(
                _cli_command(
                    st.session_state["uploaded_csv_path"],
                    save_folder,
                    x_col,
                    y_col,
                    image_col,
                    area_col,
                    mark_columns,
                    [xrange_lo, xrange_hi],
                    [yrange_lo, yrange_hi],
                    size_correction,
                    pp_threshold,
                    max_workers,
                ),
                language="bash",
            )
            st.caption(
                "On a cluster, replace the temporary upload path with the stable CSV path "
                "visible from the compute node. The command does not require Streamlit."
            )

        if st.button("Run markcorr", type="primary"):
            if not mark_columns:
                st.error("Select at least one mark column before running markcorr.")
                st.stop()
            os.makedirs(save_folder, exist_ok=True)
            progress = st.progress(0.0, text="Starting…")

            def _cb(done, total, image_num):
                progress.progress(
                    done / total,
                    text=f"image_{image_num} done ({done}/{total})",
                )

            pp_criterion = None
            if pp_threshold is not None:
                thr = float(pp_threshold)
                pp_criterion = lambda d, _t=thr: d["Area"] > _t  # noqa: E731

            t0 = time.time()
            try:
                runSpatial(
                    csv_path=st.session_state["uploaded_csv_path"],
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
                    image_col=image_col,
                    area_col=area_col,
                    mark_columns=mark_columns,
                )
            except Exception as e:  # surface user-relevant errors in the UI
                st.exception(e)
            else:
                progress.progress(1.0, text="Done")
                elapsed = time.time() - t0
                st.session_state["last_savefolder"] = str(save_folder)
                st.session_state["image_numbers"] = _list_image_dirs(save_folder)
                st.success(
                    f"markcorr finished in {elapsed:.1f}s — "
                    f"{len(st.session_state['image_numbers'])} images written to "
                    f"`{save_folder}`."
                )


# ===========================================================================
# Page 2 — Downstream analysis
# ===========================================================================
elif page.startswith("2"):
    st.header("Step 2 · Downstream analysis")
    folder = st.text_input(
        "Result folder (containing `image_<n>/iso.csv`)",
        value=st.session_state.get("last_savefolder") or "",
    )
    if not folder:
        st.info("Run step 1 first or paste a path here.")
        st.stop()
    image_numbers = _list_image_dirs(folder)
    if not image_numbers:
        st.error("No `image_*/iso.csv` files found in this folder.")
        st.stop()
    st.write(f"Detected **{len(image_numbers)}** images: {image_numbers}")

    # Group definitions
    st.subheader("Groups")
    st.caption("Assign images to one or more groups for downstream comparison.")
    n_groups = st.number_input(
        "Number of groups", min_value=1, max_value=10, value=2, step=1
    )
    groups: dict[str, list[int]] = {}
    cols = st.columns(min(int(n_groups), 4))
    for g in range(int(n_groups)):
        with cols[g % len(cols)]:
            name = st.text_input(
                f"Group {g + 1} name", value=f"group_{g + 1}", key=f"gn_{g}"
            )
            picked = st.multiselect(
                f"Images in `{name}`",
                options=image_numbers,
                default=image_numbers[g :: int(n_groups)],
                key=f"gp_{g}",
            )
            if name and picked:
                groups[name] = list(map(int, picked))
    if not groups:
        st.warning("Define at least one non-empty group to proceed.")
        st.stop()

    # Axis labels
    st.subheader("Mark labels (optional)")
    axis_keys = _peek_axis_keys(folder)
    axis_name = {}
    if axis_keys:
        with st.expander(
            f"Edit display labels for {len(axis_keys)} marks", expanded=False
        ):
            for k in axis_keys:
                axis_name[k] = st.text_input(k, value=k, key=f"ax_{k}")
    else:
        st.warning("No mark columns detected.")
        axis_name = {}

    # Build the result wrapper
    folder_norm = folder.rstrip("/") + "/"
    result = MarkcorrResult(groups=groups, resultFolder=folder_norm, axisName=axis_name)

    norm = st.selectbox("AUC normalisation", options=["min_mid_max", "log"], index=0)

    tab_auc, tab_diff, tab_box = st.tabs(
        ["AUC heatmap", "Group comparison", "Box plot"]
    )

    with tab_auc:
        if st.button("Compute AUC heatmaps", key="btn_auc"):
            with st.spinner("Computing AUC…"):
                auc, plots = result.getAUC(norm=norm, plot=True)
            st.session_state["_auc"] = auc
            st.session_state["_auc_plots"] = plots
        if "_auc_plots" in st.session_state:
            for gname, chart in st.session_state["_auc_plots"].items():
                st.markdown(f"#### `{gname}` (n = {len(groups[gname])})")
                st.altair_chart(chart, use_container_width=True)

    with tab_diff:
        if "_auc" not in st.session_state:
            st.info("Compute the AUC heatmaps first.")
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
                        "Test",
                        ["MW", "perm"],
                        index=0,
                        help="MW = Mann-Whitney + BH FDR; perm = Fisher-Pitman permutation",
                    )
                if g1 != g2 and st.button("Run find_diff", key="btn_diff"):
                    out_csv = Path(folder_norm) / f"diff_{g1}_vs_{g2}.csv"
                    with st.spinner("Running test…"):
                        chart = MarkcorrResult.find_diff(
                            st.session_state["_auc"],
                            g1,
                            g2,
                            method=method,
                            axisName=axis_name,
                            saveCsv=str(out_csv),
                        )
                    st.altair_chart(chart, use_container_width=True)
                    if out_csv.exists():
                        st.dataframe(pd.read_csv(out_csv, index_col=0))
                        st.caption(f"Saved CSV: `{out_csv}`")

    with tab_box:
        if "_auc" not in st.session_state:
            st.info("Compute the AUC heatmaps first.")
        else:
            chosen = st.multiselect(
                "Groups for the box plot",
                options=list(groups.keys()),
                default=list(groups.keys()),
            )
            if chosen and st.button("Render box plot", key="btn_box"):
                chart = result.getBoxPlot(
                    st.session_state["_auc"], chosen, axisName=axis_name
                )
                st.altair_chart(chart, use_container_width=True)


# ===========================================================================
# Page 3 — Permutation test
# ===========================================================================
else:
    st.header("Step 3 · Mark-label permutation test")
    st.caption(
        "Shuffle the categorical mark labels within each image and re-run "
        "markcorr to build a null distribution. The shuffled outputs land "
        "next to the originals in `<savefolder>/perm_<seed>/image_<n>/iso.csv`."
    )

    folder = st.text_input(
        "Result folder (the markcorr output)",
        value=st.session_state.get("last_savefolder") or "",
    )
    csv_path = st.text_input(
        "Original cell-table CSV",
        value=st.session_state.get("uploaded_csv_path") or "",
    )
    image_numbers = _list_image_dirs(folder) if folder else []
    if not image_numbers or not csv_path:
        st.info(
            "Provide both the markcorr result folder and the source CSV "
            "to enable the permutation step."
        )
        st.stop()
    if not Path(csv_path).exists():
        st.error("The source CSV path does not exist.")
        st.stop()

    df_perm_preview = pd.read_csv(csv_path, nrows=2000)
    st.subheader("Source CSV columns")
    x_col, y_col, image_col, area_col, mark_columns = _column_mapping_controls(
        df_perm_preview, "perm"
    )

    seed = st.number_input("Random seed", value=42, step=1)
    n_perm = st.number_input(
        "Number of permutations",
        min_value=1,
        max_value=1000,
        value=10,
        step=1,
    )
    chosen_images = st.multiselect(
        "Images to permute", options=image_numbers, default=image_numbers
    )

    if st.button("Run permutations", type="primary"):
        from ClumPyCells.ClumPyCells import analyzeImage  # local import for speed

        if not mark_columns:
            st.error("Select at least one mark column before running permutations.")
            st.stop()

        df_full = prepare_cell_table(
            pd.read_csv(csv_path),
            x_col=x_col,
            y_col=y_col,
            image_col=image_col,
            area_col=area_col,
            mark_columns=mark_columns,
        )
        for col in df_full.columns:
            if col not in {
                "ImageNum",
                "x",
                "y",
                "Area",
            } and not pd.api.types.is_integer_dtype(df_full[col]):
                df_full[col] = df_full[col].astype("category")
        x_min, x_max = float(df_full["x"].min()), float(df_full["x"].max())
        y_min, y_max = float(df_full["y"].min()), float(df_full["y"].max())

        progress = st.progress(0.0, text="Permuting…")
        total = int(n_perm) * len(chosen_images)
        done = 0
        rng = np.random.default_rng(int(seed))
        perm_root = Path(folder) / f"perm_runs"
        perm_root.mkdir(exist_ok=True)
        for p in range(int(n_perm)):
            seed_p = int(rng.integers(0, 2**31 - 1))
            perm_folder = perm_root / f"perm_{seed_p}"
            perm_folder.mkdir(exist_ok=True)
            for img in chosen_images:
                sub = df_full[df_full["ImageNum"] == img].copy()
                # Shuffle every non-coordinate / non-ID column.
                shuffle_cols = [
                    c for c in sub.columns if c not in {"ImageNum", "x", "y", "Area"}
                ]
                permuted = sub[shuffle_cols].sample(
                    frac=1, random_state=seed_p, ignore_index=True
                )
                sub.loc[:, shuffle_cols] = permuted.values
                analyzeImage(
                    imageNum=img,
                    imageData=sub,
                    savefolder=str(perm_folder) + "/",
                    xrange=[x_min, x_max],
                    yrange=[y_min, y_max],
                    sizeCorrection=False,
                )
                done += 1
                progress.progress(done / total, text=f"{done}/{total}")
        progress.progress(1.0, text="Done")
        st.success(f"Permutation outputs written under `{perm_root}`.")
