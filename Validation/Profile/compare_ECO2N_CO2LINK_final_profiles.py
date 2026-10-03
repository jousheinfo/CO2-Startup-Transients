# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 12:23:26 2026

@author: jpauya1
"""


"""
Compare final wellbore profiles from:
  1) T2WELL / ECO2N .listing
  2) CO2LINK / LedaFlow profile .tsv

Paper Figure 9: one combined 2x3 figure showing temperature, common-reference
specific enthalpy change, specific heat capacity, Joule-Thomson coefficient,
density, and dynamic viscosity versus depth.

Thermodynamic treatment:
  - ECO2N: enthalpy and viscosity are interpolated directly from CO2TAB.
            cp is derived from the CO2TAB enthalpy surface as (dh/dT)_P.
            mu_JT is derived from the CO2TAB enthalpy surface:

                mu_JT = (dT/dP)_h
                      = -(dh/dP)_T / (dh/dT)_P

  - CO2LINK: enthalpy, viscosity, and mu_JT are calculated with
             CoolProp/HEOS::CO2 at each CO2LINK pressure-temperature point.
             CoolProp's CO2 implementation uses the Span-Wagner EOS.

Important enthalpy note:
  The enthalpy plot shows the native absolute specific enthalpy returned by
  each thermodynamic implementation. ECO2N/CO2TAB and CoolProp may use
  different enthalpy reference states, so an absolute offset between the two
  curves must NOT by itself be interpreted as a physical enthalpy discrepancy.

Required packages:
    pip install numpy scipy pandas matplotlib CoolProp

The TSV is NOT loaded into memory. The code seeks directly to the final
"Time [s]:" block, which is important for very large profile files.

For profile consistency with the trend comparison, the T2WELL/ECO2N profile
is truncated at cell a85 (inclusive), corresponding to the selected
1472.5 m reservoir-layer center. The CO2LINK profile is left unchanged.
"""

from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from t2listing import t2listing

try:
    import CoolProp.CoolProp as CP
except ImportError:
    CP = None


# =============================================================================
# USER SETTINGS
# =============================================================================

BASE_DIR = Path.cwd()

LISTING_FILE = BASE_DIR / "CO2_injection_no_caprock_ECO2N_final_test_no_salinity_pytough.listing"
TSV_FILE = BASE_DIR / "CO2_injection_no_caprock_no_salinity_Ledaflow_comparison_ECO2N_SW_SP_final_profile.tsv"
CO2TAB_FILE = BASE_DIR / "CO2TAB"
FFLOW_FILE = BASE_DIR / "FFlow"
MESH_FILE = BASE_DIR / "MESH"

OUTPUT_DIR = BASE_DIR / "final_profile_comparison"

# CO2LINK profile columns. These names match the uploaded TSV.
TSV_DEPTH_COLUMN = "Mesh centers"
TSV_PRESSURE_COLUMN = "Pressure"
TSV_TEMPERATURE_COLUMN = "Temperature - average"
TSV_GRAVITY_GRADIENT_COLUMN = "Pressure gradient - gravity"
TSV_FRICTION_GRADIENT_COLUMN = "Pressure gradient - friction"

# Units in the uploaded TSV:
# Pressure = bara, Temperature = degC, Mesh centers = m
#
# IMPORTANT FOR THIS CO2LINK EXPORT:
# The raw Mesh-centre coordinate increases from the BOTTOM of the well
# toward the WELLHEAD. Therefore it is NOT already "depth below wellhead".
# The profile plots use these mesh-center values directly; the coordinate is
# only transformed to depth below the wellhead and then sorted so that the
# reverse order in the TSV is handled without separating P/T/property values.
#
#     depth_below_wellhead = wellhead_coordinate - raw_mesh_center
#
# Set to "top_to_bottom" only if a future export already increases downward.
CO2LINK_MESH_DIRECTION = "bottom_to_top"

# Known wellhead coordinate for this 1500 m CO2LINK model. This is used only
# as the reference for converting the reverse Mesh-centre coordinate to
# physical depth; Mesh boundaries are not used for the profile coordinates.
CO2LINK_WELLHEAD_DEPTH_M = 1500.0
CO2LINK_DEPTH_OFFSET_M = 0.0

# T2WELL profile endpoint.
# The trend comparison uses the selected 1472.5 m reservoir-layer center,
# represented by cell a85 in this T2WELL model. Keep the T2WELL profile only
# (inclusive). CO2LINK is left unchanged.
T2WELL_LAST_PROFILE_CELL = "a85"

# Plot settings
SAVE_DPI = 300
SHOW_PLOTS = True

# Thermodynamic fluid used for CO2LINK.
# HEOS::CO2 is CoolProp's Helmholtz EOS implementation for CO2.
COOLPROP_FLUID = "HEOS::CO2"


# =============================================================================
# GENERAL HELPERS
# =============================================================================

def normalize_block_name(name: str) -> str:
    """Normalize a TOUGH block name for matching while preserving uniqueness."""
    return "".join(str(name).split()).lower()


def ensure_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found:\n  {path}")


def validate_numeric_array(name, arr):
    arr = np.asarray(arr, dtype=float)
    if not np.all(np.isfinite(arr)):
        bad = np.where(~np.isfinite(arr))[0]
        raise ValueError(f"{name} contains non-finite values at indices {bad[:10].tolist()}")
    return arr


def truncate_t2well_profile_at_cell(df, last_cell):
    """
    Keep the T2WELL profile from the wellhead through ``last_cell`` inclusive.

    The dataframe is already sorted by physical depth, so this truncates the
    profile at the physical depth of the requested T2WELL cell. All quantities
    (P, T, h, viscosity, JT, etc.) therefore use the same endpoint.
    """
    target = normalize_block_name(last_cell)
    normalized_blocks = df["block"].map(normalize_block_name)

    matches = df.index[normalized_blocks == target].tolist()
    if not matches:
        available = ", ".join(df["block"].astype(str).tail(10).tolist())
        raise ValueError(
            f"T2WELL endpoint cell '{last_cell}' was not found in the extracted "
            f"well profile. Last available cells are: {available}"
        )

    target_index = matches[0]
    target_depth = float(df.loc[target_index, "depth_m"])

    # Keep every well cell from the wellhead down to the selected a85 cell,
    # including a85 itself.
    truncated = df[df["depth_m"] <= target_depth].copy()
    truncated = truncated.sort_values("depth_m").reset_index(drop=True)

    return truncated, target_depth


# =============================================================================
# 1. T2WELL / ECO2N: READ FINAL WELL PROFILE FROM .listing
# =============================================================================

def read_t2well_final_profile(listing_path: Path):
    """
    Read the well-cell geometry and the LAST printed element-state block
    directly from a T2WELL/ECO2N listing file.

    This avoids hard-coding the number of well cells or their spacing.

    Returns
    -------
    df : pandas.DataFrame
        Columns:
            cell_id, block, depth_m, pressure_Pa, pressure_bar, temperature_C,
            density_kgm3_listing
    final_time_days : float
    """

    # Example geometry line in this listing:
    #   a90       90    0.14975E+04 E   1 ...
    #   a 2        2    0.30003E+02 I   2 ...
    geometry_pattern = re.compile(
        r"^\s*(.*?)\s+(\d+)\s+"
        r"([+-]?(?:\d+\.?\d*|\.\d+)[EeDd][+-]?\d+)\s+"
        r"([SIE])\s+"
    )

    time_pattern = re.compile(
        r"THE TIME IS\s*([+-]?[0-9.]+(?:[EeDd][+-]?\d+)?)\s*DAYS",
        re.IGNORECASE,
    )

    well_geometry = {}
    reading_geometry = False

    # This dictionary is reset every time a new OUTPUT DATA block begins.
    # At EOF it therefore contains the final printed timestep.
    current_profile = {}
    final_time_days = None
    output_blocks_seen = 0

    with open(listing_path, "r", encoding="latin1", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")

            # -------- Well geometry block --------
            if "Cell  ID   Distance" in line:
                reading_geometry = True
                continue

            if reading_geometry:
                if "Totally" in line:
                    reading_geometry = False
                    continue

                m = geometry_pattern.match(line)
                if m:
                    block_raw = m.group(1).strip()
                    cell_id = int(m.group(2))
                    depth_m = float(m.group(3).replace("D", "E").replace("d", "e"))
                    key = normalize_block_name(block_raw)
                    well_geometry[key] = {
                        "cell_id": cell_id,
                        "block": block_raw,
                        "depth_m": depth_m,
                    }
                continue

            # -------- Start of a new printed timestep --------
            if "OUTPUT DATA AFTER" in line:
                output_blocks_seen += 1
                current_profile = {}

                mt = time_pattern.search(line)
                if mt:
                    final_time_days = float(
                        mt.group(1).replace("D", "E").replace("d", "e")
                    )
                continue

            # -------- Element rows --------
            # T2WELL prints a 5-character block name beginning at column 2.
            # Example:
            #  *wa 1   1 0.91929E+07  35.00 ...
            #    a 2  51 0.92806E+07  35.20 ...
            if well_geometry and len(line) >= 7:
                block_raw = line[1:6]
                key = normalize_block_name(block_raw)

                if key not in well_geometry:
                    continue

                fields = line[6:].split()
                if len(fields) < 3:
                    continue

                try:
                    global_index = int(fields[0])
                    pressure_Pa = float(fields[1].replace("D", "E").replace("d", "e"))
                    temperature_C = float(fields[2].replace("D", "E").replace("d", "e"))

                    # T2WELL/ECO2N element-state output:
                    # fields[10] = DG = gas density [kg/m3].
                    density_kgm3 = float(
                        fields[10].replace("D", "E").replace("d", "e")
                    )
                except (ValueError, IndexError):
                    # This rejects connection/generator rows that may start
                    # with the same block name but are not element-state rows.
                    continue

                current_profile[key] = {
                    "global_index": global_index,
                    "pressure_Pa": pressure_Pa,
                    "temperature_C": temperature_C,
                    "density_kgm3": density_kgm3,
                }

    if not well_geometry:
        raise RuntimeError(
            "Could not find the T2WELL well geometry block in the listing file."
        )

    if output_blocks_seen == 0 or not current_profile:
        raise RuntimeError(
            "Could not find a final OUTPUT DATA element-state block in the listing file."
        )

    rows = []
    missing = []

    for key, geom in well_geometry.items():
        state = current_profile.get(key)
        if state is None:
            missing.append(geom["block"])
            continue

        rows.append(
            {
                "cell_id": geom["cell_id"],
                "block": geom["block"],
                "depth_m": geom["depth_m"],
                "global_index": state["global_index"],
                "pressure_Pa": state["pressure_Pa"],
                "pressure_bar": state["pressure_Pa"] / 1.0e5,
                "temperature_C": state["temperature_C"],
                "density_kgm3_listing": state["density_kgm3"],
            }
        )

    if missing:
        warnings.warn(
            f"{len(missing)} well cells were present in the geometry block but "
            f"not found in the final output block. First missing cells: {missing[:10]}"
        )

    df = pd.DataFrame(rows).sort_values("depth_m").reset_index(drop=True)

    if df.empty:
        raise RuntimeError("No wellbore cells were extracted from the final listing output.")

    return df, final_time_days



# =============================================================================
# 1A. T2WELL / MESH: READ WELLBORE CELL-CENTER DEPTHS
# =============================================================================

def read_t2well_mesh_wellbore_centers(mesh_path: Path):
    """
    Read the vertical wellbore cell-center depths from the MESH ELEME section.

    The current MESH stores the vertical coordinate as the last scientific-
    notation field in each ELEME record. The wellbore cells are read directly from the MESH ELEME section. The MESH coordinate is negative downward, so

        depth [m] = -z [m].

    No interpolation is performed.
    """
    ensure_file(mesh_path, "T2WELL MESH file")

    sci_pattern = re.compile(
        r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)[Ee][+-]\d{2}"
    )

    centers = {}
    in_eleme = False

    with open(mesh_path, "r", encoding="latin1", errors="replace") as f:
        for line in f:
            stripped = line.strip()

            if stripped == "ELEME":
                in_eleme = True
                continue

            if in_eleme and stripped == "CONNE":
                break

            if not in_eleme or len(line) < 6:
                continue

            block_raw = line[:5].strip()
            block_key = normalize_block_name(block_raw)

            # Wellbore cells are *wa 1 and a 2 ... a 91.
            if block_key == normalize_block_name("*wa1"):
                is_wellbore = True
            else:
                m = re.fullmatch(r"a(\d+)", block_key)
                is_wellbore = m is not None

            if not is_wellbore:
                continue

            numbers = sci_pattern.findall(line[6:])
            if len(numbers) < 6:
                continue

            z_m = float(numbers[-1].replace("D", "E").replace("d", "e"))
            centers[block_key] = -z_m

    if not centers:
        raise RuntimeError(
            "No T2WELL wellbore cell centers were found in the MESH ELEME section."
        )

    if normalize_block_name("*wa1") not in centers:
        raise RuntimeError("MESH wellbore cell *wa 1 was not found.")

    numeric_a_cells = [
        key for key in centers
        if re.fullmatch(r"a\d+", key)
    ]
    if len(numeric_a_cells) < 1:
        raise RuntimeError("No numeric T2WELL wellbore cells were found in MESH.")

    depths = np.asarray(
        [centers[key] for key in numeric_a_cells],
        dtype=float,
    )
    if not np.all(np.isfinite(depths)):
        raise RuntimeError("Non-finite T2WELL MESH wellbore depths were found.")

    return centers


# =============================================================================
# 1AA. T2WELL / ECO2N: READ FINAL WELLBORE FRICTION TERM
# =============================================================================

def read_t2well_final_friction(listing_path: Path, well_geometry_df: pd.DataFrame):
    """
    Read the FINAL T2WELL wellbore friction-gradient values from the
    momentum-equation table.

    Important details for the T2WELL listing format used here:
      * The final momentum table is split across two printed sections, each
        repeating the ELEM1/ELEM2/.../fric header. Therefore both sections
        must be accumulated for the same final timestep.
      * The first connection is printed as ``a 2*wa 1`` without the normal
        spacing between the two five-character element names. Therefore the
        element names are not parsed by fixed character positions.
      * The friction value itself is a fixed-column numeric field and can be
        read safely from the header-defined ``fric`` column.
      * The listing's own well-cell geometry is used for the connection
        midpoint. This is intentional: the supplied MESH file does not have
        the same wellbore geometry as this listing.

    The friction value belongs to a connection between adjacent wellbore
    cells. It is plotted at the midpoint of the two corresponding cell
    centers. No interpolation is performed.
    """

    # Build the physical well-cell ordering from the listing geometry.
    geometry = well_geometry_df.copy()
    geometry["block_key"] = geometry["block"].map(normalize_block_name)
    geometry = geometry.sort_values("depth_m").reset_index(drop=True)

    if len(geometry) < 2:
        raise RuntimeError(
            "The T2WELL listing geometry contains fewer than two wellbore cells."
        )

    # The momentum table in this listing contains one connection per adjacent
    # pair of wellbore cells. We use the listing geometry itself to establish
    # the connection depths, rather than relying on the element-name formatting
    # in the momentum table.
    n_expected_connections = len(geometry) - 1

    # The final momentum-equation section begins after the final ELEMENT SOURCE
    # table and ends before the following volume/mass-balance section.
    with open(listing_path, "r", encoding="latin1", errors="replace") as f:
        text = f.read()

    source_matches = list(re.finditer(r"^\s*ELEMENT SOURCE INDEX\b", text, re.MULTILINE))
    if not source_matches:
        raise RuntimeError(
            "Could not locate the final T2WELL ELEMENT SOURCE section."
        )

    section_start = source_matches[-1].start()

    balance_match = re.search(
        r"^\s*\*{10,}\s*VOLUME- AND MASS-BALANCES",
        text[section_start:],
        re.MULTILINE,
    )
    if balance_match:
        section_end = section_start + balance_match.start()
    else:
        section_end = len(text)

    final_section = text[section_start:section_end]

    # Locate every momentum-table header inside the final section. The final
    # table is split into multiple printed pages, so all headers in this
    # section belong to the same final timestep.
    header_pattern = re.compile(
        r"^[ \t]*ELEM1\s+ELEM2\s+INDEX.*\bfric\b.*\bMomDL\b",
        re.MULTILINE,
    )
    headers = list(header_pattern.finditer(final_section))

    if not headers:
        raise RuntimeError(
            "Could not find the momentum-equation table header containing "
            "the T2WELL 'fric' column in the final timestep."
        )

    # The header gives the exact fixed-column position of the fric field.
    # Use the same column positions for all repeated headers.
    header_line = headers[-1].group(0)
    fric_start = header_line.index("fric")
    momdl_start = header_line.index("MomDL")

    rows_by_index = {}

    # Each repeated momentum header is followed by data rows until the next
    # header or the end of the section. Parse only rows with an integer INDEX.
    for hnum, header_match in enumerate(headers):
        data_start = header_match.end()
        data_end = (
            headers[hnum + 1].start()
            if hnum + 1 < len(headers)
            else len(final_section)
        )

        block = final_section[data_start:data_end]

        for line in block.splitlines():
            if not line.strip():
                continue

            # The element-name area is not safely tokenizable because the
            # first connection is printed as "a 2*wa 1", and the character
            # position of INDEX changes when the index has more digits.
            # Instead, identify INDEX immediately before the first scientific-
            # notation value (Vmix). The remaining numerical fields then are:
            #
            # INDEX, Vmix, Vdrift, Vsugas, Vsuliq, MomDT, Efpg, Rey, fric, MomDL
            #
            index_vmix_pattern = re.search(
                r"(\d+)\s+"
                r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)[Ee][+-]?\d+)",
                line,
            )
            if index_vmix_pattern is None:
                continue

            try:
                index = int(index_vmix_pattern.group(1))
                remaining = line[index_vmix_pattern.end():].split()

                # After Vmix, fric is the seventh remaining numerical value:
                # Vdrift, Vsugas, Vsuliq, MomDT, Efpg, Rey, fric.
                if len(remaining) < 7:
                    continue

                fric = float(
                    remaining[6].replace("D", "E").replace("d", "e")
                )
            except (ValueError, IndexError):
                continue

            # Keep only the wellbore connections expected from this listing.
            if not (1 <= index <= n_expected_connections):
                continue

            if not np.isfinite(fric):
                continue

            # The final occurrence of an index is retained defensively.
            rows_by_index[index] = fric

    expected_indices = set(range(1, n_expected_connections + 1))
    found_indices = set(rows_by_index)

    missing = sorted(expected_indices - found_indices)
    if missing:
        raise RuntimeError(
            "The final T2WELL momentum table is incomplete. Missing friction "
            f"connection indices: {missing[:20]}"
            + (" ..." if len(missing) > 20 else "")
        )

    # T2WELL's connection sequence is:
    #   index 1: cell 2 -> cell 1
    #   index 2: cell 3 -> cell 2
    #   ...
    # Therefore each index k corresponds to the midpoint of geometry[k-1]
    # and geometry[k], after the geometry has been sorted by physical depth.
    rows = []
    for index in range(1, n_expected_connections + 1):
        cell_up = geometry.iloc[index - 1]
        cell_down = geometry.iloc[index]

        depth_mid = 0.5 * (
            float(cell_up["depth_m"]) + float(cell_down["depth_m"])
        )

        rows.append(
            {
                "index": index,
                "elem1": str(cell_down["block"]),
                "elem2": str(cell_up["block"]),
                "depth_m": depth_mid,
                "friction_gradient_Pam": rows_by_index[index],
            }
        )

    df = pd.DataFrame(rows).sort_values("depth_m").reset_index(drop=True)

    return df



# =============================================================================
# 1B. T2WELL / ECO2N: READ THE COMPLETE TEMPERATURE/PRESSURE HISTORY
# =============================================================================

def read_t2well_history(listing_path: Path, geometry_df: pd.DataFrame):
    """
    Read the complete T2WELL pressure and temperature history for the
    wellbore cells using PyTOUGH/t2listing.

    T2WELL/ECO2N listing values are used in their native units:
        pressure = Pa
        temperature = degC

    The well-cell depths come from the already parsed listing geometry, so
    the 1 m + 19 m top split and all subsequent cell centres are preserved.
    """
    geometry = geometry_df.sort_values("depth_m").reset_index(drop=True)
    depths_m = geometry["depth_m"].to_numpy(dtype=float)

    # Use the same constructor form as the previously working validation
    # script. The history reader needs the element records from the listing.
    lst = t2listing(str(listing_path))

    times_s = None
    temperature_columns = []
    pressure_columns = []

    for _, geom in geometry.iterrows():
        block = str(geom["block"])

        # T2WELL uses fixed-width 5-character element names. The final-profile
        # parser strips the leading padding when storing the block name, but
        # t2listing.history() expects the padded element name. Restore it here
        # without changing the dataframe or element ordering. For example:
        #   *wa 1 -> *wa 1
        #   a85  ->   a85
        #   a 2  ->   a 2
        history_block = block if block.startswith("*") else block.rjust(5)

        history = lst.history(
            [
                ("e", history_block, "P"),
                ("e", history_block, "T"),
            ]
        )

        if history is None:
            raise RuntimeError(
                f"t2listing.history() returned None for T2WELL element "
                f"'{block}' requested as '{history_block}'."
            )

        (time_p, pressure_pa), (time_t, temperature_c) = history

        time_p = np.asarray(time_p, dtype=float)
        time_t = np.asarray(time_t, dtype=float)
        pressure_pa = np.asarray(pressure_pa, dtype=float)
        temperature_c = np.asarray(temperature_c, dtype=float)

        if times_s is None:
            times_s = time_p.copy()
        else:
            if (
                len(time_p) != len(times_s)
                or not np.allclose(time_p, times_s, rtol=0.0, atol=1.0e-6)
            ):
                raise RuntimeError(
                    f"T2WELL pressure-history times for cell '{block}' "
                    "do not match the reference time grid."
                )

        if (
            len(time_t) != len(times_s)
            or not np.allclose(time_t, times_s, rtol=0.0, atol=1.0e-6)
        ):
            raise RuntimeError(
                f"T2WELL temperature-history times for cell '{block}' "
                "do not match the reference time grid."
            )

        pressure_columns.append(pressure_pa)
        temperature_columns.append(temperature_c)

    if times_s is None or not pressure_columns:
        raise RuntimeError("No T2WELL wellbore history was extracted.")

    times_s = np.asarray(times_s, dtype=float)
    pressure_pa = np.column_stack(pressure_columns)
    temperature_c = np.column_stack(temperature_columns)

    order = np.argsort(times_s)

    times_days = times_s[order] / 86400.0
    temperature_c = temperature_c[order]
    pressure_pa = pressure_pa[order]

    # Sanity checks on the native units used by the plotting code.
    if not np.all(np.isfinite(pressure_pa)):
        raise RuntimeError("T2WELL pressure history contains non-finite values.")
    if not np.all(np.isfinite(temperature_c)):
        raise RuntimeError("T2WELL temperature history contains non-finite values.")

    return times_days, depths_m, temperature_c, pressure_pa


# =============================================================================
# 2. CO2LINK / LEDAFLOW: READ ONLY THE FINAL TSV PROFILE BLOCK
# =============================================================================

def find_last_marker_position(binary_file, marker: bytes, chunk_size=8 * 1024 * 1024):
    """
    Search a large file backwards for the final occurrence of marker.
    This avoids reading a multi-GB TSV from the beginning.
    """
    binary_file.seek(0, 2)
    file_size = binary_file.tell()

    pos = file_size
    overlap = b""
    overlap_len = max(len(marker) - 1, 0)

    while pos > 0:
        read_size = min(chunk_size, pos)
        pos -= read_size
        binary_file.seek(pos)
        chunk = binary_file.read(read_size)
        data = chunk + overlap

        idx = data.rfind(marker)
        if idx != -1:
            return pos + idx

        overlap = data[:overlap_len] if overlap_len else b""

    return None


def _find_column_index(headers, exact_name):
    clean = [h.strip() for h in headers]
    try:
        return clean.index(exact_name)
    except ValueError as exc:
        raise KeyError(
            f"Column '{exact_name}' was not found in the TSV profile.\n"
            f"Available columns are:\n  " + "\n  ".join(clean)
        ) from exc


def read_co2link_final_profile(tsv_path: Path):
    """
    Read the LAST 'Time [s]:' block from the CO2LINK/LedaFlow TSV.

    The file is opened in binary mode and searched backwards, which makes
    this practical even when the TSV is several GB.

    IMPORTANT
    ---------
    In the supplied CO2LINK export, ``Mesh centers`` increases from the
    bottom of the well toward the wellhead. It is therefore a profile/elevation
    coordinate, not depth below wellhead. The physical depth used in every
    comparison is converted ONCE here, while pressure and temperature remain
    attached to the same row/state. This automatically puts pressure,
    temperature, enthalpy, and JT coefficient on the same corrected depth axis.

    Returns
    -------
    df : pandas.DataFrame
        mesh_center_raw_m, depth_m, pressure_bar, pressure_Pa, temperature_C
    final_time_s : float
    """

    marker = b"Time [s]:"

    with open(tsv_path, "rb") as f:
        marker_pos = find_last_marker_position(f, marker)
        if marker_pos is None:
            raise RuntimeError("Could not find any 'Time [s]:' block in the TSV file.")

        f.seek(marker_pos)

        time_line = f.readline().decode("utf-8", errors="replace").strip()
        mt = re.search(
            r"Time\s*\[s\]\s*:\s*\t?\s*([+-]?[0-9.]+(?:[EeDd][+-]?\d+)?)",
            time_line,
            flags=re.IGNORECASE,
        )
        if not mt:
            raise RuntimeError(f"Could not parse final TSV time from line:\n{time_line}")

        final_time_s = float(mt.group(1).replace("D", "E").replace("d", "e"))

        header_line = f.readline().decode("utf-8", errors="replace").rstrip("\r\n")
        unit_line = f.readline().decode("utf-8", errors="replace").rstrip("\r\n")

        headers = header_line.split("\t")
        units = unit_line.split("\t")

        i_depth = _find_column_index(headers, TSV_DEPTH_COLUMN)
        i_pressure = _find_column_index(headers, TSV_PRESSURE_COLUMN)
        i_temperature = _find_column_index(headers, TSV_TEMPERATURE_COLUMN)
        i_velocity = _find_column_index(headers, "Velocity - gas")
        i_gravity_gradient = _find_column_index(
            headers, TSV_GRAVITY_GRADIENT_COLUMN
        )
        i_friction_gradient = _find_column_index(
            headers, TSV_FRICTION_GRADIENT_COLUMN
        )

        required_max_index = max(
            i_depth,
            i_pressure,
            i_temperature,
            i_velocity,
            i_gravity_gradient,
            i_friction_gradient,
        )

        mesh_center_raw = []
        pressure_bar = []
        temperature_C = []
        velocity_gas = []
        gravity_gradient_Pam = []
        friction_gradient_Pam = []

        for raw_line in f:
            text = raw_line.decode("utf-8", errors="replace").rstrip("\r\n")

            if not text.strip():
                continue

            # Defensive stop in case something follows the final profile block.
            if text.startswith("Time [s]:"):
                break

            fields = text.split("\t")
            if len(fields) <= required_max_index:
                continue

            try:
                z_raw = float(fields[i_depth])
                p = float(fields[i_pressure])
                t = float(fields[i_temperature])
                vgas = float(fields[i_velocity])
                gravity_gradient = float(fields[i_gravity_gradient])
                friction_gradient = float(fields[i_friction_gradient])
            except ValueError:
                continue

            mesh_center_raw.append(z_raw)
            pressure_bar.append(p)
            temperature_C.append(t)
            velocity_gas.append(vgas)
            gravity_gradient_Pam.append(gravity_gradient)
            friction_gradient_Pam.append(friction_gradient)

    if not mesh_center_raw:
        raise RuntimeError("No numeric rows were found in the final TSV profile block.")

    if len(velocity_gas) != len(mesh_center_raw):
        raise RuntimeError(
            f"CO2LINK velocity length mismatch: velocity={len(velocity_gas)}, "
            f"profile rows={len(mesh_center_raw)}."
        )

    if len(gravity_gradient_Pam) != len(mesh_center_raw):
        raise RuntimeError(
            "CO2LINK gravity-gradient length mismatch: "
            f"gravity={len(gravity_gradient_Pam)}, "
            f"profile rows={len(mesh_center_raw)}."
        )

    if len(friction_gradient_Pam) != len(mesh_center_raw):
        raise RuntimeError(
            "CO2LINK friction-gradient length mismatch: "
            f"friction={len(friction_gradient_Pam)}, "
            f"profile rows={len(mesh_center_raw)}."
        )

    raw = np.asarray(mesh_center_raw, dtype=float)

    df = pd.DataFrame(
        {
            "mesh_center_raw_m": raw,
            "pressure_bar": np.asarray(pressure_bar, dtype=float),
            "temperature_C": np.asarray(temperature_C, dtype=float),
            "velocity_gas_ms": np.asarray(velocity_gas, dtype=float),
            "gravity_gradient_Pam": np.asarray(
                gravity_gradient_Pam, dtype=float
            ),
            "friction_gradient_Pam": np.asarray(
                friction_gradient_Pam, dtype=float
            ),
        }
    )

    direction = CO2LINK_MESH_DIRECTION.strip().lower()

    if direction == "bottom_to_top":
        # The TSV lists mesh centers from bottom to top. Use the mesh-center
        # coordinate for every plotted point and reverse it to physical depth.
        top_coordinate_m = (
            CO2LINK_WELLHEAD_DEPTH_M + CO2LINK_DEPTH_OFFSET_M
        )
        df["depth_m"] = top_coordinate_m - df["mesh_center_raw_m"]

    elif direction == "top_to_bottom":
        # Future files may already increase downward from the wellhead.
        # Use the mesh centers directly; no boundary values are required.
        top_coordinate_m = float(np.nanmin(df["mesh_center_raw_m"]))
        df["depth_m"] = df["mesh_center_raw_m"] - top_coordinate_m

    else:
        raise ValueError(
            "CO2LINK_MESH_DIRECTION must be either 'bottom_to_top' "
            "or 'top_to_bottom'."
        )

    df["pressure_Pa"] = df["pressure_bar"] * 1.0e5

    # Sort only AFTER converting the coordinate. Pressure and temperature are
    # never independently reversed; each P-T state stays attached to its cell.
    # Any later derived property (h, mu_JT, etc.) therefore uses the corrected
    # physical depth automatically.
    df = (
        df.sort_values("depth_m")
          .drop_duplicates("depth_m")
          .reset_index(drop=True)
    )

    return df, final_time_s



# =============================================================================
# 2B. CO2LINK / LEDAFLOW: READ THE COMPLETE TEMPERATURE/PRESSURE HISTORY
# =============================================================================

def read_co2link_history(tsv_path: Path, target_depths_m):
    """
    Stream through the CO2LINK/LedaFlow TSV and retain:
      - temperature at every mesh center for every time step;
      - pressure and temperature at selected physical depths using the
        corresponding mesh centers directly.

    No interpolation is performed. The CO2LINK mesh coordinate is corrected
    exactly as in read_co2link_final_profile(): raw Mesh centers increase from
    bottom to top, so depth below the wellhead is the offset-corrected
    top-coordinate minus raw_mesh_center.
    """
    target_depths_m = np.asarray(target_depths_m, dtype=float)

    times_s = []
    temperature_rows = []
    pressure_target_rows = []
    reference_depths = None
    top_coordinate_m = (
        CO2LINK_WELLHEAD_DEPTH_M + CO2LINK_DEPTH_OFFSET_M
        if CO2LINK_MESH_DIRECTION.strip().lower() == "bottom_to_top"
        else None
    )

    current_time_s = None
    headers = None
    i_depth = i_pressure = i_temperature = None
    raw_depth = []
    pressure = []
    temperature = []

    def store_current_block():
        nonlocal reference_depths, top_coordinate_m

        if current_time_s is None or not raw_depth:
            return

        raw = np.asarray(raw_depth, dtype=float)
        p = np.asarray(pressure, dtype=float)
        t = np.asarray(temperature, dtype=float)

        if top_coordinate_m is None:
            # For a future top-to-bottom export, the shallowest mesh center is
            # the zero-depth reference. Mesh boundaries are not used.
            top_coordinate_m = float(np.nanmin(raw))

        if CO2LINK_MESH_DIRECTION.strip().lower() == "bottom_to_top":
            depth = top_coordinate_m - raw
        else:
            depth = raw - top_coordinate_m

        valid = np.isfinite(depth) & np.isfinite(p) & np.isfinite(t)
        depth = depth[valid]
        p = p[valid]
        t = t[valid]

        order = np.argsort(depth)
        depth = depth[order]
        p = p[order]
        t = t[order]

        if reference_depths is None:
            reference_depths = depth.copy()
        elif (
            len(depth) != len(reference_depths)
            or not np.allclose(depth, reference_depths, rtol=0.0, atol=1.0e-8)
        ):
            raise RuntimeError(
                "CO2LINK mesh centers changed between time steps. "
                "A fixed depth-time contour grid cannot be constructed safely."
            )

        # The updated CO2LINK mesh uses the same physical mesh centers as
        # T2WELL. Select the requested cells directly; do not interpolate.
        p_target = np.full(len(target_depths_m), np.nan, dtype=float)
        t_target = np.full(len(target_depths_m), np.nan, dtype=float)

        for j, z in enumerate(target_depths_m):
            matches = np.where(
                np.isclose(depth, z, rtol=0.0, atol=1.0e-3)
            )[0]

            if len(matches) != 1:
                raise RuntimeError(
                    f"CO2LINK target depth {z:.6f} m does not correspond "
                    f"to exactly one mesh center; found {len(matches)} matches."
                )

            idx = matches[0]
            p_target[j] = p[idx]
            t_target[j] = t[idx]

        times_s.append(current_time_s)
        temperature_rows.append(t)
        pressure_target_rows.append(p_target)
    with open(tsv_path, "rb") as f:
        while True:
            raw_line = f.readline()
            if not raw_line:
                store_current_block()
                break

            line = raw_line.decode("utf-8", errors="replace").rstrip("\r\n")

            if line.strip().startswith("Time [s]:"):
                store_current_block()

                current_time_s = None
                headers = None
                raw_depth = []
                pressure = []
                temperature = []

                mt = re.search(
                    r"Time\s*\[s\]\s*:\s*\t?\s*"
                    r"([+-]?[0-9.]+(?:[EeDd][+-]?\d+)?)",
                    line,
                    flags=re.IGNORECASE,
                )
                if not mt:
                    raise RuntimeError(f"Could not parse CO2LINK time from:\n{line}")

                current_time_s = float(
                    mt.group(1).replace("D", "E").replace("d", "e")
                )

                header_line = f.readline().decode(
                    "utf-8", errors="replace"
                ).rstrip("\r\n")
                unit_line = f.readline().decode(
                    "utf-8", errors="replace"
                ).rstrip("\r\n")

                headers = header_line.split("\t")
                _ = unit_line  # units are retained by the file format but not needed here

                i_depth = _find_column_index(headers, TSV_DEPTH_COLUMN)
                i_pressure = _find_column_index(headers, TSV_PRESSURE_COLUMN)
                i_temperature = _find_column_index(headers, TSV_TEMPERATURE_COLUMN)

                continue

            if headers is None or not line.strip():
                continue

            fields = line.split("\t")
            max_index = max(i_depth, i_pressure, i_temperature)
            if len(fields) <= max_index:
                continue

            try:
                z_raw = float(fields[i_depth])
                p = float(fields[i_pressure])
                t = float(fields[i_temperature])
            except ValueError:
                continue

            raw_depth.append(z_raw)
            pressure.append(p)
            temperature.append(t)

    if not times_s:
        raise RuntimeError("No CO2LINK time-history blocks were found.")

    times_s = np.asarray(times_s, dtype=float)
    temperature_C = np.asarray(temperature_rows, dtype=float)
    pressure_target_bar = np.asarray(pressure_target_rows, dtype=float)

    order = np.argsort(times_s)

    return (
        times_s[order],
        reference_depths,
        temperature_C[order],
        pressure_target_bar[order] * 1.0e5,
    )


# =============================================================================
# 3. ECO2N CO2TAB READER AND THERMODYNAMIC INTERPOLATION
# =============================================================================

def read_co2tab(co2tab_path: Path):
    """
    Read the main ECO2N CO2TAB P-T property table.

    File structure used here:
        nP  nT
        nP pressure values [Pa]
        nT temperature values [degC]
        for each pressure:
            nT density values   [kg/m3]
            nT viscosity values [Pa s]
            nT enthalpy values  [J/kg]

    CO2TAB may contain an additional saturation-related table after the main
    P-T table. It is intentionally not used here.
    """

    with open(co2tab_path, "r", encoding="latin1", errors="replace") as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise RuntimeError("CO2TAB appears to be incomplete.")

    dims = lines[1].split()
    if len(dims) < 2:
        raise RuntimeError("Could not read nP and nT from CO2TAB line 2.")

    nP = int(dims[0])
    nT = int(dims[1])

    values = []
    for line in lines[2:]:
        for token in line.split():
            try:
                values.append(float(token.replace("D", "E").replace("d", "e")))
            except ValueError:
                pass

    minimum_needed = nP + nT + 3 * nP * nT
    if len(values) < minimum_needed:
        raise RuntimeError(
            f"CO2TAB has too few numeric values. Need at least {minimum_needed}, "
            f"found {len(values)}."
        )

    cursor = 0
    pressure_grid_Pa = np.asarray(values[cursor:cursor + nP], dtype=float)
    cursor += nP

    temperature_grid_C = np.asarray(values[cursor:cursor + nT], dtype=float)
    cursor += nT

    main_count = 3 * nP * nT
    main = np.asarray(values[cursor:cursor + main_count], dtype=float)
    main = main.reshape(nP, 3, nT)

    density = main[:, 0, :]
    viscosity = main[:, 1, :]
    enthalpy = main[:, 2, :]

    if not np.all(np.diff(pressure_grid_Pa) > 0):
        raise RuntimeError("CO2TAB pressure grid is not strictly increasing.")
    if not np.all(np.diff(temperature_grid_C) > 0):
        raise RuntimeError("CO2TAB temperature grid is not strictly increasing.")

    return {
        "pressure_Pa": pressure_grid_Pa,
        "temperature_C": temperature_grid_C,
        "density_kgm3": density,
        "viscosity_Pas": viscosity,
        "enthalpy_Jkg": enthalpy,
    }


def build_co2tab_interpolators(co2tab):
    """
    Build interpolators for h(P,T), viscosity(P,T), dh/dP|T and dh/dT|P.
    The dh/dT|P interpolator is also used as the ECO2N specific heat capacity.

    The derivatives are first evaluated on the native non-uniform CO2TAB grid
    using numpy.gradient, then linearly interpolated along the actual P-T path.
    """

    P = co2tab["pressure_Pa"]
    T = co2tab["temperature_C"]
    density = co2tab["density_kgm3"]
    h = co2tab["enthalpy_Jkg"]
    viscosity = co2tab["viscosity_Pas"]

    # Units:
    # dh_dP : (J/kg)/Pa
    # dh_dT : (J/kg)/K   [temperature increments in degC and K are identical]
    dh_dP = np.gradient(h, P, axis=0, edge_order=2)
    dh_dT = np.gradient(h, T, axis=1, edge_order=2)

    density_interp = RegularGridInterpolator(
        (P, T), density, method="linear", bounds_error=True
    )

    h_interp = RegularGridInterpolator(
        (P, T), h, method="linear", bounds_error=True
    )
    viscosity_interp = RegularGridInterpolator(
        (P, T), viscosity, method="linear", bounds_error=True
    )
    dhdP_interp = RegularGridInterpolator(
        (P, T), dh_dP, method="linear", bounds_error=True
    )
    dhdT_interp = RegularGridInterpolator(
        (P, T), dh_dT, method="linear", bounds_error=True
    )

    return density_interp, h_interp, viscosity_interp, dhdP_interp, dhdT_interp


def evaluate_eco2n_properties(P_Pa, T_C, co2tab):
    """
    Evaluate ECO2N CO2TAB enthalpy, viscosity, specific heat capacity,
    and Joule-Thomson coefficient.

    The specific heat capacity is obtained from the thermodynamic relation

        cp = (dh/dT)_P

    using the ECO2N CO2TAB enthalpy surface.

    mu_JT = -(dh/dP)_T / (dh/dT)_P

    Returns
    -------
    h_Jkg : ndarray
    viscosity_Pas : ndarray
    cp_JkgK : ndarray
    mu_JT_K_per_MPa : ndarray
    """

    P_Pa = validate_numeric_array("ECO2N pressure", P_Pa)
    T_C = validate_numeric_array("ECO2N temperature", T_C)

    Pmin, Pmax = co2tab["pressure_Pa"][[0, -1]]
    Tmin, Tmax = co2tab["temperature_C"][[0, -1]]

    bad = (P_Pa < Pmin) | (P_Pa > Pmax) | (T_C < Tmin) | (T_C > Tmax)
    if np.any(bad):
        idx = np.where(bad)[0]
        details = [
            f"i={i}, P={P_Pa[i]/1e5:.4f} bar, T={T_C[i]:.4f} C"
            for i in idx[:10]
        ]
        raise ValueError(
            "Some ECO2N final-profile states lie outside the CO2TAB range.\n"
            f"CO2TAB range: P={Pmin/1e5:.3f}-{Pmax/1e5:.3f} bar, "
            f"T={Tmin:.3f}-{Tmax:.3f} C.\n"
            "First offending states:\n  " + "\n  ".join(details)
        )

    density_interp, h_interp, viscosity_interp, dhdP_interp, dhdT_interp = (
        build_co2tab_interpolators(co2tab)
    )

    points = np.column_stack((P_Pa, T_C))
    density_kgm3 = np.asarray(density_interp(points), dtype=float)
    h_Jkg = np.asarray(h_interp(points), dtype=float)
    viscosity_Pas = np.asarray(viscosity_interp(points), dtype=float)
    dhdP = np.asarray(dhdP_interp(points), dtype=float)
    dhdT = np.asarray(dhdT_interp(points), dtype=float)

    if np.any(np.isclose(dhdT, 0.0, atol=1e-12)):
        warnings.warn("Very small dh/dT encountered while calculating ECO2N mu_JT.")

    mu_JT_K_per_Pa = -dhdP / dhdT
    mu_JT_K_per_MPa = mu_JT_K_per_Pa * 1.0e6

    cp_JkgK = dhdT

    return density_kgm3, h_Jkg, viscosity_Pas, cp_JkgK, mu_JT_K_per_MPa


# =============================================================================
# 4. CO2LINK: COOLPROP / SPAN-WAGNER PROPERTIES FROM P AND T
# =============================================================================

def evaluate_coolprop_properties(P_Pa, T_C, fluid=COOLPROP_FLUID):
    """
    Calculate h, dynamic viscosity, specific heat capacity, and mu_JT
    for CO2LINK states with CoolProp.

    The Joule-Thomson coefficient is evaluated using

        mu_JT = (T*alpha_p - 1) / (rho*cp)

    where alpha_p is the isobaric thermal-expansion coefficient.

    Returns
    -------
    h_Jkg : ndarray
    viscosity_Pas : ndarray
    cp_JkgK : ndarray
    mu_JT_K_per_MPa : ndarray
    """

    if CP is None:
        raise ImportError(
            "CoolProp is required for the CO2LINK thermodynamic properties. "
            "Install it with: pip install CoolProp"
        )

    P_Pa = validate_numeric_array("CO2LINK pressure", P_Pa)
    T_C = validate_numeric_array("CO2LINK temperature", T_C)

    density = np.full(P_Pa.shape, np.nan, dtype=float)
    h = np.full(P_Pa.shape, np.nan, dtype=float)
    viscosity = np.full(P_Pa.shape, np.nan, dtype=float)
    cp = np.full(P_Pa.shape, np.nan, dtype=float)
    mu_jt = np.full(P_Pa.shape, np.nan, dtype=float)

    failed = []

    for i, (p, tc) in enumerate(zip(P_Pa, T_C)):
        tk = tc + 273.15
        try:
            h_i = CP.PropsSI("Hmass", "P", p, "T", tk, fluid)
            viscosity_i = CP.PropsSI("VISCOSITY", "P", p, "T", tk, fluid)
            rho_i = CP.PropsSI("Dmass", "P", p, "T", tk, fluid)
            cp_i = CP.PropsSI("Cpmass", "P", p, "T", tk, fluid)
            alpha_i = CP.PropsSI(
                "ISOBARIC_EXPANSION_COEFFICIENT", "P", p, "T", tk, fluid
            )

            # K/Pa
            mu_i = (tk * alpha_i - 1.0) / (rho_i * cp_i)

            density[i] = rho_i
            h[i] = h_i
            viscosity[i] = viscosity_i
            cp[i] = cp_i
            mu_jt[i] = mu_i * 1.0e6  # K/MPa

        except Exception as exc:
            failed.append((i, p, tc, str(exc)))

    if failed:
        msg = [
            f"i={i}, P={p/1e5:.4f} bar, T={tc:.4f} C -> {err}"
            for i, p, tc, err in failed[:10]
        ]
        warnings.warn(
            f"CoolProp failed for {len(failed)} CO2LINK state(s). "
            "Those values were set to NaN.\n  " + "\n  ".join(msg)
        )

    return density, h, viscosity, cp, mu_jt



# =============================================================================
# 5. FFlow: READ FINAL T2WELL VELOCITY PROFILE
# =============================================================================

def read_fflow_final_profile(fflow_path: Path):
    """
    Read the final FFlow time block.

    FFlow columns:
        Time, Depth, Fliq, Fgas, VLiq, VGas, Umix

    The comparison uses VGas because the injected stream is pure CO2.
    """
    rows = []

    with open(fflow_path, "r", encoding="latin1", errors="replace") as f:
        for line in f:
            if "," not in line:
                continue

            parts = [x.strip() for x in line.split(",")]
            if len(parts) < 7:
                continue

            try:
                values = [float(x.replace("E", "e")) for x in parts[:7]]
            except ValueError:
                continue

            rows.append(values)

    if not rows:
        raise RuntimeError("No numerical data found in FFlow file.")

    arr = np.asarray(rows)

    # Keep the last printed time step only.
    final_time = arr[-1, 0]
    arr = arr[np.isclose(arr[:, 0], final_time)]

    df = pd.DataFrame(
        {
            "time_s": arr[:, 0],
            "depth_m": arr[:, 1],
            "velocity_gas_ms": arr[:, 5],
            "velocity_mixture_ms": arr[:, 6],
        }
    )

    df = df.sort_values("depth_m").reset_index(drop=True)
    return df, final_time


# =============================================================================
# 5. PLOTTING
# =============================================================================


def plot_combined_thermo_profiles(
    eco2n_df,
    co2link_df,
    co2tab,
    output_path,
):
    """
    Create a compact 2x2 comparison of the final temperature, enthalpy-change,
    specific-heat-capacity, and Joule-Thomson coefficient profiles.
    """

    # Use the same common P-T reference state as plot_delta_h_profile().
    p_ref = float(eco2n_df.iloc[0]["pressure_Pa"])
    t_ref = float(eco2n_df.iloc[0]["temperature_C"])

    _, h_eco_ref, _, _, _ = evaluate_eco2n_properties(
        np.asarray([p_ref], dtype=float),
        np.asarray([t_ref], dtype=float),
        co2tab,
    )
    _, h_link_ref, _, _, _ = evaluate_coolprop_properties(
        np.asarray([p_ref], dtype=float),
        np.asarray([t_ref], dtype=float),
    )

    delta_h_eco_kJkg = (
        eco2n_df["enthalpy_Jkg"].to_numpy(dtype=float) - h_eco_ref[0]
    ) / 1000.0

    delta_h_link_kJkg = (
        co2link_df["enthalpy_Jkg"].to_numpy(dtype=float) - h_link_ref[0]
    ) / 1000.0

    fig, axes = plt.subplots(
        2, 3,
        figsize=(9.0, 6.5),
        sharey=True,
    )

    profiles = [
        (
            axes[0, 0],
            eco2n_df["temperature_C"],
            co2link_df["temperature_C"],
            "Temperature [°C]",
            "Temperature",
        ),
        (
            axes[0, 1],
            delta_h_eco_kJkg,
            delta_h_link_kJkg,
            "Specific enthalpy change, Δh [kJ/kg]",
            "Enthalpy change",
        ),
        (
            axes[0, 2],
            eco2n_df["cp_JkgK"] / 1000.0,
            co2link_df["cp_JkgK"] / 1000.0,
            "Specific heat capacity, c$_p$ [kJ/(kg·K)]",
            "Specific heat capacity",
        ),
        (
            axes[1, 0],
            eco2n_df["JT_K_per_MPa"],
            co2link_df["JT_K_per_MPa"],
            "Joule–Thomson coefficient [K/MPa]",
            "Joule–Thomson coefficient",
        ),
        (
            axes[1, 1],
            eco2n_df["density_kgm3"],
            co2link_df["density_kgm3"],
            "Density [kg/m³]",
            "Density",
        ),
        (
            axes[1, 2],
            eco2n_df["viscosity_mPas"],
            co2link_df["viscosity_mPas"],
            "Dynamic viscosity [mPa·s]",
            "Viscosity",
        ),
    ]

    for ax, x_eco, x_link, xlabel, title in profiles:
        ax.plot(
            x_eco,
            eco2n_df["depth_m"],
            linewidth=1.5,
            linestyle="--",
            label="T2WELL / ECO2N",
        )
        ax.plot(
            x_link,
            co2link_df["depth_m"],
            linewidth=1.5,
            linestyle="-",
            label="CO2LINK",
        )

        # Shade the perforated interval (1450-1500 m).
        ax.axhspan(1450.0, 1500.0, alpha=0.3, zorder=0)

        ax.set_xlabel(xlabel, fontsize=10)
        # No subplot titles in the combined 2x2 figure.
        ax.tick_params(axis="both", labelsize=8)
        ax.grid(True, alpha=0.30)

    # The four subplots share the same y-axis; invert it only once.
    # Calling invert_yaxis() on every shared axis would toggle the direction
    # repeatedly and return it to the original orientation.
    axes[0, 0].invert_yaxis()

    axes[0, 0].set_ylabel("Depth [m]", fontsize=10)
    axes[1, 0].set_ylabel("Depth [m]", fontsize=10)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=2,
        fontsize=8,
        frameon=True,
        bbox_to_anchor=(0.5, 0.03),
    )

    fig.tight_layout(rect=[0, 0.08, 1, 1])
    fig.savefig(output_path, dpi=SAVE_DPI, bbox_inches="tight")
    fig.savefig(
        Path(output_path).with_suffix(".eps"),
        format="eps",
        bbox_inches="tight",
    )

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)


# =============================================================================
# Figure 9 is the only retained profile figure.
# =============================================================================


# =============================================================================
# 6. MAIN
# =============================================================================

def main():
    ensure_file(LISTING_FILE, "T2WELL/ECO2N listing file")
    ensure_file(TSV_FILE, "CO2LINK/LedaFlow TSV profile file")
    ensure_file(CO2TAB_FILE, "ECO2N CO2TAB file")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Reading final T2WELL/ECO2N profile...")
    eco2n_full_df, eco2n_time_days = read_t2well_final_profile(LISTING_FILE)

    # Match the profile comparison to the same bottomhole/midpoint definition
    # used in the trend comparison. For this model, the selected reservoir-layer center is T2WELL
    # cell a85. CO2LINK is intentionally not truncated here.
    eco2n_df = eco2n_full_df.copy()
    eco2n_df, t2well_endpoint_depth_m = truncate_t2well_profile_at_cell(
        eco2n_df, T2WELL_LAST_PROFILE_CELL
    )
    print(
        f"T2WELL profile truncated at {T2WELL_LAST_PROFILE_CELL} "
        f"(depth = {t2well_endpoint_depth_m:.3f} m)."
    )

    print("Reading final CO2LINK TSV profile...")
    co2link_df, co2link_time_s = read_co2link_final_profile(TSV_FILE)
    co2link_time_days = co2link_time_s / 86400.0

    print("Reading ECO2N CO2TAB...")
    co2tab = read_co2tab(CO2TAB_FILE)

    print("Evaluating ECO2N enthalpy, viscosity, specific heat capacity and Joule-Thomson coefficient from CO2TAB...")
    eco_density, eco_h, eco_viscosity, eco_cp, eco_jt = evaluate_eco2n_properties(
        eco2n_df["pressure_Pa"].to_numpy(),
        eco2n_df["temperature_C"].to_numpy(),
        co2tab,
    )
    eco2n_df["density_kgm3"] = eco_density
    eco2n_df["enthalpy_Jkg"] = eco_h
    eco2n_df["enthalpy_kJkg"] = eco_h / 1000.0
    eco2n_df["viscosity_Pas"] = eco_viscosity
    eco2n_df["viscosity_mPas"] = eco_viscosity * 1000.0
    eco2n_df["cp_JkgK"] = eco_cp
    eco2n_df["JT_K_per_MPa"] = eco_jt

    print("Evaluating CO2LINK enthalpy, viscosity, specific heat capacity and Joule-Thomson coefficient with CoolProp...")
    link_density, link_h, link_viscosity, link_cp, link_jt = evaluate_coolprop_properties(
        co2link_df["pressure_Pa"].to_numpy(),
        co2link_df["temperature_C"].to_numpy(),
    )
    co2link_df["density_kgm3"] = link_density
    co2link_df["enthalpy_Jkg"] = link_h
    co2link_df["enthalpy_kJkg"] = link_h / 1000.0
    co2link_df["viscosity_Pas"] = link_viscosity
    co2link_df["viscosity_mPas"] = link_viscosity * 1000.0
    co2link_df["cp_JkgK"] = link_cp
    co2link_df["JT_K_per_MPa"] = link_jt

    # -------------------------------------------------------------------------
    # Diagnostics
    # -------------------------------------------------------------------------
    print("\n---------------- FINAL PROFILE SUMMARY ----------------")
    print(f"T2WELL/ECO2N final time : {eco2n_time_days:.10g} days")
    print(f"CO2LINK final time       : {co2link_time_days:.10g} days")
    print(f"T2WELL/ECO2N well cells  : {len(eco2n_df)}")
    print(
        f"T2WELL profile endpoint  : {T2WELL_LAST_PROFILE_CELL} "
        f"at {t2well_endpoint_depth_m:.3f} m"
    )
    print(f"CO2LINK profile points   : {len(co2link_df)}")
    print(
        f"T2WELL depth range       : {eco2n_df['depth_m'].min():.3f} "
        f"to {eco2n_df['depth_m'].max():.3f} m"
    )
    print(
        f"CO2LINK depth range      : {co2link_df['depth_m'].min():.3f} "
        f"to {co2link_df['depth_m'].max():.3f} m"
    )
    print(
        "CO2LINK shallowest cell : "
        f"z={co2link_df.iloc[0]['depth_m']:.3f} m, "
        f"P={co2link_df.iloc[0]['pressure_bar']:.3f} bar, "
        f"T={co2link_df.iloc[0]['temperature_C']:.3f} C"
    )
    print(
        "CO2LINK deepest cell    : "
        f"z={co2link_df.iloc[-1]['depth_m']:.3f} m, "
        f"P={co2link_df.iloc[-1]['pressure_bar']:.3f} bar, "
        f"T={co2link_df.iloc[-1]['temperature_C']:.3f} C"
    )
    if eco2n_time_days is not None:
        time_difference_s = abs(eco2n_time_days * 86400.0 - co2link_time_s)
        relative = time_difference_s / max(abs(co2link_time_s), 1.0)
        if time_difference_s > 1.0 and relative > 1.0e-6:
            warnings.warn(
                "The final ECO2N and CO2LINK output times are not identical: "
                f"difference = {time_difference_s:.6g} s. "
                "The plots still compare each file's last available profile."
            )

    # Save numerical data
    # -------------------------------------------------------------------------
    eco2n_csv = OUTPUT_DIR / "ECO2N_final_profile_with_CO2TAB_properties.csv"
    co2link_csv = OUTPUT_DIR / "CO2LINK_final_profile_with_CoolProp_properties.csv"

    eco2n_df.to_csv(eco2n_csv, index=False)
    co2link_df.to_csv(co2link_csv, index=False)

    # -------------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------------
    combined_thermo_png = OUTPUT_DIR / "depth_vs_thermo_profiles_2x3.png"
    plot_combined_thermo_profiles(
        eco2n_df,
        co2link_df,
        co2tab,
        combined_thermo_png,
    )

    print("\nSaved:")
    print(f"  {eco2n_csv}")
    print(f"  {co2link_csv}")
    print(f"  {combined_thermo_png}")
    print(f"  {combined_thermo_png.with_suffix('.eps')}")

if __name__ == "__main__":
    main()
