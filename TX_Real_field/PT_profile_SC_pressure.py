# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 12:11:52 2026

@author: jpauya1
"""


import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import CoolProp.CoolProp as CP
import numpy as np


def save_figure_png_and_eps(fig, filename, dpi=300):
    """Save a Matplotlib figure in both PNG and EPS formats."""
    png_path = Path(filename)
    eps_path = png_path.with_suffix(".eps")

    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(eps_path, format="eps", dpi=dpi, bbox_inches="tight")

    print(f"Saved: {png_path}")
    print(f"Saved: {eps_path}")
    
    
_TIME_RE = re.compile(r"^\s*Time\s*\[s\]\s*:\s*([+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)\s*$")


def _split_row(line: str) -> list[str]:
    """
    Split a data row that may be tab-separated or whitespace-aligned.
    """
    line = line.strip()
    if not line:
        return []
    # Prefer tabs if present; otherwise fall back to generic whitespace split.
    if "\t" in line:
        return [x for x in line.split("\t") if x != ""]
    return line.split()


def parse_ledaflow_profile_tsv(
    file_path: str | Path,
    mesh_col: int = 1,
    pressure_col: int = 25,
    temperature_col: int = 35,
    columns_are_1_based: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse a LedaFlow "Profile result file" exported as TSV/text blocks by time.

    Returns:
      pressure_df: index = mesh boundary (depth/position), columns = time [s]
      temperature_df: index = mesh boundary (depth/position), columns = time [s]

    Notes:
    - mesh_col is the "Mesh boundaries" column (your 1st column).
    - pressure_col is your 25th column.
    - temperature_col is your 35th column.
    - Times are read from lines like:  Time [s]:  0.0
    """

    file_path = Path(file_path)

    # Convert to 0-based indices for Python lists.
    if columns_are_1_based:
        mesh_i = mesh_col - 1
        p_i = pressure_col - 1
        t_i = temperature_col - 1
    else:
        mesh_i, p_i, t_i = mesh_col, pressure_col, temperature_col

    need_cols = max(mesh_i, p_i, t_i)

    pressure_series_by_time: dict[float, pd.Series] = {}
    temp_series_by_time: dict[float, pd.Series] = {}

    current_time: float | None = None
    current_mesh: list[float] = []
    current_p: list[float] = []
    current_T: list[float] = []

    def flush_block():
        nonlocal current_time, current_mesh, current_p, current_T
        if current_time is None or not current_mesh:
            return

        # Build Series indexed by mesh boundaries.
        idx = pd.Index(current_mesh, name="mesh_boundary")
        pressure_series_by_time[current_time] = pd.Series(current_p, index=idx, name=current_time)
        temp_series_by_time[current_time] = pd.Series(current_T, index=idx, name=current_time)

        # Reset for next block.
        current_mesh, current_p, current_T = [], [], []

    with file_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")

            # Detect start of a new time block.
            m = _TIME_RE.match(line)
            if m:
                flush_block()
                current_time = float(m.group(1))
                continue

            # If we haven't hit a time block yet, ignore everything.
            if current_time is None:
                continue

            # Try to parse a data row.
            parts = _split_row(line)
            if not parts:
                # Blank lines inside a block: ignore.
                continue

            # Skip obvious header/unit rows (non-numeric first token).
            try:
                float(parts[0])
            except ValueError:
                continue

            # Must have enough columns to reach your requested indices.
            if len(parts) <= need_cols:
                continue

            # Parse needed fields.
            try:
                mesh_val = float(parts[mesh_i])
                p_val = float(parts[p_i])
                T_val = float(parts[t_i])
            except ValueError:
                # If a row is malformed, skip it.
                continue

            current_mesh.append(mesh_val)
            current_p.append(p_val)
            current_T.append(T_val)

    # Flush last block at EOF.
    flush_block()

    if not pressure_series_by_time:
        raise ValueError(
            "No time blocks / data rows were parsed. "
            "Check that the file contains lines like 'Time [s]: 0.0' "
            "and that the column indices are correct."
        )

    # Combine series into DataFrames (columns are times).
    pressure_df = pd.concat(pressure_series_by_time.values(), axis=1)
    temperature_df = pd.concat(temp_series_by_time.values(), axis=1)

    # Sort columns by time (numeric).
    pressure_df = pressure_df.reindex(sorted(pressure_df.columns), axis=1)
    temperature_df = temperature_df.reindex(sorted(temperature_df.columns), axis=1)

    # Optional: ensure monotonic index ordering.
    pressure_df = pressure_df.sort_index()
    temperature_df = temperature_df.sort_index()
    
    # Keep mesh boundaries order; reverse only the P/T values down the well for each time
    pressure_df = pd.DataFrame(
        {t: pressure_df[t].to_numpy()[::-1] for t in pressure_df.columns},
        index=pressure_df.index,
    )
    temperature_df = pd.DataFrame(
        {t: temperature_df[t].to_numpy()[::-1] for t in temperature_df.columns},
        index=temperature_df.index,
    )
    pressure_df.columns.name = "time_s"
    temperature_df.columns.name = "time_s"

    # Name the columns meaningfully.
    pressure_df.columns.name = "time_s"
    temperature_df.columns.name = "time_s"

    return pressure_df, temperature_df


if __name__ == "__main__":
    file_path = "Pure_SC_CO2_injection_real_field_radial_TX_Ledaflow_pressure_for_Feb25_profile.tsv"

    P_df, T_df = parse_ledaflow_profile_tsv(
        file_path,
        mesh_col=1,
        pressure_col=25,
        temperature_col=35,
        columns_are_1_based=True,
    )

    print("Pressure DF shape:", P_df.shape)
    print("Temperature DF shape:", T_df.shape)

    # Save if you want
    P_df.to_csv("pressure_by_time.csv")       # index is mesh_boundary
    T_df.to_csv("temperature_by_time.csv")

    # -----------------------------
    # Phase contour (bar, °C)
    # -----------------------------
    # depth axis from index (mesh boundaries)
    depth_values = P_df.index.to_numpy(dtype=float)

    # time axis from columns (already floats: seconds)
    time_values = P_df.columns.to_numpy(dtype=float)

    # Make sure time is sorted and aligned in both dfs
    sort_idx = np.argsort(time_values)
    time_values = time_values[sort_idx]
    P_mat_bar = P_df.to_numpy()[:, sort_idx]
    T_mat_C = T_df.to_numpy()[:, sort_idx]

    n_depth, n_time = P_mat_bar.shape

    # Identify CO2 phase from the local pressure (Pa) and temperature (K).
    P_mat_Pa = P_mat_bar * 1e5
    # Preserve plot codes; CoolProp's liquid-like/gas-like regions use 1/2.
    phase_mapping = {
        "supercritical": 0,
        "liquid": 1,
        "supercritical_liquid": 1,
        "gas": 2,
        "supercritical_gas": 2,
        "twophase": np.nan,
        "critical_point": np.nan,
    }
    phase_numeric = np.full((n_depth, n_time), np.nan)
    for i, j in np.ndindex(phase_numeric.shape):
        if not (np.isfinite(P_mat_Pa[i, j]) and np.isfinite(T_mat_C[i, j])):
            continue
        phase = CP.PhaseSI("P", P_mat_Pa[i, j], "T", T_mat_C[i, j] + 273.15, "CO2")
        if phase not in phase_mapping:
            raise ValueError(f"CO2 phase unresolved at depth {depth_values[i]}, time {time_values[j]} s: {phase}")
        phase_numeric[i, j] = phase_mapping[phase]

    # Plot
    time_mesh, depth_mesh = np.meshgrid(time_values, depth_values)

    cmap = mcolors.ListedColormap(["green", "blue", "red"])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # perforations given in ft → convert to m
    perf_top_ft = 4479
    perf_bot_ft = 4977
    ft_to_m = 0.3048
    perf_top_m = perf_top_ft * ft_to_m
    perf_bot_m = perf_bot_ft * ft_to_m

    plt.figure(figsize=(12, 8))
    time_hours = time_mesh / 3600.0
    time_days = time_mesh / (24*3600)

    cf = plt.contourf(
        time_days, depth_mesh, phase_numeric,
        levels=[-0.5, 0.5, 1.5, 2.5],
        cmap=cmap, norm=norm
    )
    cbar = plt.colorbar(cf, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(["Supercritical", "Liquid", "Gas"])

    plt.xlabel("Time (days)")
    plt.ylabel("Depth (m)")
    #plt.title("CO₂ Phase Evolution in the Wellbore")
    plt.gca().invert_yaxis()  # if increasing depth goes downward
    plt.grid(True, alpha=0.3)
    
    # --- horizontal perforation lines ---
    plt.axhline(perf_top_m, color="brown", linestyle="--", linewidth=2)
    plt.axhline(perf_bot_m, color="brown", linestyle="--", linewidth=2)
    
    # --- label placement (a bit inside the left edge) ---
    y_offset_m = 35.0  # tweak this
    x_min, x_max = plt.xlim()
    label_x_pos = x_min + 0.05 * (x_max - x_min)
    
    plt.text(label_x_pos, perf_top_m - y_offset_m, f"Perf Top: {perf_top_m} m",
             color="brown", va="center", fontsize=10, fontweight="bold")
    plt.text(label_x_pos, perf_bot_m - y_offset_m, f"Perf Bottom: {perf_bot_m} m",
             color="brown", va="center", fontsize=10, fontweight="bold")
    
    # Save the contour plot
    save_figure_png_and_eps(plt.gcf(), "CO2_phase_evolution_contour.png", dpi=300)
    print("Contour plot saved as: CO2_phase_evolution_contour.png")
    
    plt.show()