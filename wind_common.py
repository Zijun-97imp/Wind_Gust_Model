# -*- coding: utf-8 -*-
"""
wind_common.py
==============
Shared utilities for the multi-scale gust wind model.

This module centralises:

* project-root / data-file resolution (so every sub-program works
  regardless of the current working directory);
* unit conversion constants (metre <-> foot, knot <-> m/s);
* altitude-band classification;
* helper functions for saving figures and CSV results.

Reference
---------
MIL-HDBK-1797 / MIL-F-8785C (Dryden and Von Karman turbulence models)
Fang Z., Zhang G., "Characteristics Research of Multi-scale Wind Model by
the Space and Time".
"""

import os
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Project / data paths
# --------------------------------------------------------------------------- #
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

DATA_XLSX = os.path.join(PROJECT_ROOT, "wind_info_over_altitude.xlsx")
WHITE_NOISE_CSV = os.path.join(PROJECT_ROOT, "white_noise.csv")
WIND_DATA_CSV = os.path.join(PROJECT_ROOT, "wind_data.csv")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "output")


def ensure_output_dir():
    """Create (if necessary) and return the output directory."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    return OUTPUT_DIR


# --------------------------------------------------------------------------- #
# Unit conversions
# --------------------------------------------------------------------------- #
M_TO_FT = 3.280839895
FT_TO_M = 1.0 / M_TO_FT
KNOT_TO_MS = 0.514444444
MS_TO_KNOT = 1.0 / KNOT_TO_MS


# --------------------------------------------------------------------------- #
# Altitude classification
# --------------------------------------------------------------------------- #
# Altitude boundaries in metres (MIL-HDBK-1797 uses feet).
LOW_ALT_LIMIT = 304.8        # 1000 ft  -> low-altitude model
HIGH_ALT_LIMIT = 609.6       # 2000 ft  -> medium/high-altitude model
OVER_ALT_LIMIT = 1000.0      # 1000 m   -> above this use the lookup table


def classify_altitude(height):
    """
    Classify the flight altitude into one of the four model bands.

    Returns
    -------
    str
        One of ``"low"`` (h < 1000 ft), ``"transition"`` (1000..2000 ft),
        ``"mid-high"`` (2000 ft .. 1000 m), ``"over"`` (>= 1000 m).
    """
    if height < LOW_ALT_LIMIT:
        return "low"
    if height <= HIGH_ALT_LIMIT:
        return "transition"
    if height < OVER_ALT_LIMIT:
        return "mid-high"
    return "over"


# --------------------------------------------------------------------------- #
# Turbulence strength
# --------------------------------------------------------------------------- #
STRENGTH_LABELS = {1: "light", 2: "moderate", 3: "severe"}


def validate_strength(strength):
    if int(strength) not in STRENGTH_LABELS:
        raise ValueError(
            "strength must be 1 (light), 2 (moderate) or 3 (severe); "
            f"got {strength!r}"
        )
    return int(strength)


# --------------------------------------------------------------------------- #
# High-altitude lookup table
# --------------------------------------------------------------------------- #
def _load_over_altitude_table():
    """Load ``wind_info_over_altitude.xlsx`` once and cache the result."""
    global _OVER_ALTITUDE_CACHE
    if "_OVER_ALTITUDE_CACHE" not in globals():
        df = pd.read_excel(DATA_XLSX)
        # Columns:
        # altitude, light_sgima_u, light_sigma_vw, moderate_sigma_u,
        # moderate_sigma_vw, severe_sigma_u, severe_sigma_vw, Lu, Lvw
        _OVER_ALTITUDE_CACHE = df
    return _OVER_ALTITUDE_CACHE


def over_altitude_values(height):
    """
    Return the tabulated scale lengths (Lu, Lv, Lw) and turbulence speeds
    (sigma_u, sigma_v, sigma_w) at ``height`` (>= 1000 m).

    The spreadsheet stores Lu/Lvw in metres and the sigma values in knots;
    everything is returned in SI (m / m/s).
    """
    df = _load_over_altitude_table()
    altitudes = df["altitude"].astype(float).tolist()

    # Snap to the nearest tabulated altitude if the exact value is not found.
    idx = _nearest_index(altitudes, float(height))

    Lu = float(df["Lu"].iloc[idx])             # metre
    Lvw = float(df["Lvw"].iloc[idx])           # metre
    return Lu, Lvw, Lvw


def over_altitude_sigma(height, strength):
    """Return (sigma_u, sigma_v, sigma_w) in m/s for h >= 1000 m."""
    strength = validate_strength(strength)
    df = _load_over_altitude_table()
    altitudes = df["altitude"].astype(float).tolist()
    idx = _nearest_index(altitudes, float(height))

    col_u = {
        1: "light_sgima_u",
        2: "moderate_sigma_u",
        3: "severe_sigma_u",
    }[strength]
    col_vw = {
        1: "light_sigma_vw",
        2: "moderate_sigma_vw",
        3: "severe_sigma_vw",
    }[strength]

    sigma_u = float(df[col_u].iloc[idx]) * KNOT_TO_MS
    sigma_vw = float(df[col_vw].iloc[idx]) * KNOT_TO_MS
    return sigma_u, sigma_vw, sigma_vw


def _nearest_index(values, target):
    """Index of the element of ``values`` closest to ``target``."""
    arr = np.asarray(values, dtype=float)
    return int(np.argmin(np.abs(arr - target)))


# --------------------------------------------------------------------------- #
# Plot / file helpers
# --------------------------------------------------------------------------- #
def save_figure(fig, name, out_dir=None, dpi=200):
    """Save a matplotlib figure into the output directory."""
    out_dir = out_dir or ensure_output_dir()
    path = os.path.join(out_dir, name)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


def save_csv(columns, out_dir=None, name="result.csv"):
    """
    Save a dict-of-arrays as a CSV file.

    Parameters
    ----------
    columns : dict
        Mapping of column name -> 1-D array.
    """
    out_dir = out_dir or ensure_output_dir()
    path = os.path.join(out_dir, name)
    df = pd.DataFrame(columns)
    df.to_csv(path, index=False)
    return path
