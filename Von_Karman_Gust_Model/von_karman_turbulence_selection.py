# -*- coding: utf-8 -*-
"""
von_karman_turbulence_selection.py
==================================
Von Karman turbulence *scale length* and *intensity* selection versus altitude.

The Von Karman and Dryden models share the same turbulence scale lengths and
intensities (both derive them from MIL-HDBK-1797 / MIL-F-8785C); they differ
only in the spectral shaping filter.  This module therefore re-exports the
shared selection logic from :mod:`dryden_turbulence_selection` under
Von-Karman-specific names, keeping the existing call sites unchanged while
avoiding duplicated (and previously diverging) code.
"""

import os
import sys

# Project root (for wind_common) and the Dryden module (for the shared logic).
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, "Dryden_Gust_Model"))

import dryden_turbulence_selection as _dts  # noqa: E402


def vonkarman_length_scale(height):
    """Von Karman turbulence scale lengths ``(Lu, Lv, Lw)`` in metres."""
    return _dts.dryden_length_scale(height)


def vonkarman_turbulence_speed(height, strength):
    """Von Karman turbulence intensities ``(sigma_u, sigma_v, sigma_w)`` in m/s."""
    return _dts.dryden_turbulence_speed(height, strength)
