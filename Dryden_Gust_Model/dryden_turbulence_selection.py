# -*- coding: utf-8 -*-
"""
dryden_turbulence_selection.py
==============================
Dryden turbulence *scale length* and *intensity* selection versus altitude.

The Dryden model describes **small-scale** atmospheric turbulence.  This
module provides the turbulence scale lengths ``Lu, Lv, Lw`` and the gust
intensities ``sigma_u, sigma_v, sigma_w`` (all returned in SI units:
metres and m/s) for a given flight altitude and turbulence strength.

Altitude bands (MIL-HDBK-1797 / MIL-F-8785C)
--------------------------------------------
* ``low``        : h < 1000 ft  (304.8 m)   -> analytical formulas
* ``transition`` : 1000..2000 ft            -> linear interpolation
* ``mid-high``   : 2000 ft .. 1000 m        -> constant 2500 ft scales
* ``over``       : h >= 1000 m              -> lookup table

Turbulence strength:  1 = light,  2 = moderate,  3 = severe.

References
----------
MIL-HDBK-1797, "Flying Qualities of Piloted Aircraft".
"""

import os
import sys

# Allow this module to import the shared helpers located in the project root.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import wind_common as wc


# --------------------------------------------------------------------------- #
# Reference wind speeds at 20 ft and 2000 ft for each strength (ft/s)
# --------------------------------------------------------------------------- #
_W20 = {1: 15.0, 2: 30.0, 3: 45.0}      # low-altitude reference (ft/s)
_W2000 = {1: 25.0, 2: 40.0, 3: 55.0}    # mid/high-altitude reference (ft/s)

# Mid/high-altitude constant scale lengths (ft)
_LU_MIDHIGH = 2500.0
_LV_MIDHIGH = _LU_MIDHIGH / 2.0
_LW_MIDHIGH = _LV_MIDHIGH


def _low_altitude(h_ft, strength):
    """Analytical low-altitude scales/intensities (input/output in ft, ft/s)."""
    denom_scale = (0.177 + 0.000823 * h_ft) ** 1.2
    denom_inten = (0.177 + 0.000823 * h_ft) ** 0.4

    Lw = h_ft / 2.0
    Lu = h_ft / denom_scale
    Lv = Lu / 2.0

    sigma_w = 0.1 * _W20[strength]
    sigma_u = sigma_w / denom_inten
    sigma_v = sigma_u
    return Lu, Lv, Lw, sigma_u, sigma_v, sigma_w


def _midhigh_altitude(strength):
    """Constant mid/high-altitude scales/intensities (ft, ft/s)."""
    sigma_w = 0.1 * _W2000[strength]
    return (
        _LU_MIDHIGH, _LV_MIDHIGH, _LW_MIDHIGH,
        sigma_w, sigma_w, sigma_w,
    )


def _linear_interp(h_ft, v1000, v2000):
    """Linear interpolation between the value at 1000 ft and 2000 ft."""
    frac = (h_ft - 1000.0) / (2000.0 - 1000.0)
    return v1000 + (v2000 - v1000) * frac


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def dryden_length_scale(height):
    """
    Dryden turbulence scale lengths ``(Lu, Lv, Lw)`` in metres.

    Parameters
    ----------
    height : float
        Flight altitude in metres.
    """
    band = wc.classify_altitude(height)

    if band == "low":
        h_ft = height * wc.M_TO_FT
        Lu, Lv, Lw, _, _, _ = _low_altitude(h_ft, 1)
    elif band == "mid-high":
        Lu, Lv, Lw = _LU_MIDHIGH, _LV_MIDHIGH, _LW_MIDHIGH
    elif band == "transition":
        h_ft = height * wc.M_TO_FT
        Lu1000, Lv1000, Lw1000, _, _, _ = _low_altitude(1000.0, 1)
        Lu = _linear_interp(h_ft, Lu1000, _LU_MIDHIGH)
        Lv = _linear_interp(h_ft, Lv1000, _LV_MIDHIGH)
        Lw = _linear_interp(h_ft, Lw1000, _LW_MIDHIGH)
    else:  # "over" -> table is already in metres, convert to ft for uniformity
        Lu_m, Lv_m, Lw_m = wc.over_altitude_values(height)
        Lu, Lv, Lw = Lu_m * wc.M_TO_FT, Lv_m * wc.M_TO_FT, Lw_m * wc.M_TO_FT

    return Lu * wc.FT_TO_M, Lv * wc.FT_TO_M, Lw * wc.FT_TO_M


def dryden_turbulence_speed(height, strength):
    """
    Dryden turbulence intensities ``(sigma_u, sigma_v, sigma_w)`` in m/s.

    Parameters
    ----------
    height : float
        Flight altitude in metres.
    strength : int
        1 (light), 2 (moderate) or 3 (severe).
    """
    strength = wc.validate_strength(strength)
    band = wc.classify_altitude(height)

    if band == "low":
        h_ft = height * wc.M_TO_FT
        _, _, _, sigma_u, sigma_v, sigma_w = _low_altitude(h_ft, strength)
    elif band == "mid-high":
        _, _, _, sigma_u, sigma_v, sigma_w = _midhigh_altitude(strength)
    elif band == "transition":
        h_ft = height * wc.M_TO_FT
        *_, su1000, sv1000, sw1000 = _low_altitude(1000.0, strength)
        *_, su2000, sv2000, sw2000 = _midhigh_altitude(strength)
        sigma_u = _linear_interp(h_ft, su1000, su2000)
        sigma_v = _linear_interp(h_ft, sv1000, sv2000)
        sigma_w = _linear_interp(h_ft, sw1000, sw2000)
    else:  # "over" -> table sigma is in m/s, convert to ft/s for uniformity
        su_m, sv_m, sw_m = wc.over_altitude_sigma(height, strength)
        sigma_u, sigma_v, sigma_w = su_m * wc.M_TO_FT, sv_m * wc.M_TO_FT, sw_m * wc.M_TO_FT

    return sigma_u * wc.FT_TO_M, sigma_v * wc.FT_TO_M, sigma_w * wc.FT_TO_M
