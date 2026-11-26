# -*- coding: utf-8 -*-
"""
normal_distribution_prediction.py
=================================
Gust amplitude *correction / prediction* via the normal-distribution method.

This implements the paper's amplitude-correction calculation.  Given the
auto-correlation function of the wind field between two locations separated
by ``d``, the gust amplitude at the target location is obtained for any
desired confidence level.

Theory
------
Auto-correlation functions (Dryden)::

    R_u(d) = exp(-d/Lu)
    R_v(d) = exp(-d/Lv) (1 - d/(2 Lv))
    R_w(d) = exp(-d/Lw) (1 - d/(2 Lw))

The *equivalent* (conditional) strength is ``sigma_e = sigma sqrt(1 - R^2)``.
The amplitude at a location separated by ``d``, for a two-sided confidence
level ``C``, follows a normal distribution whose percentile is
``P = (1 + C) / 2``::

    U(d) = z_P * sigma_e(d)                      (initial velocity V1 = 0)
    U(d) = R(d) * V1 + z_P * sigma_e(d)          (initial velocity V1 != 0)

where ``z_P = norm.ppf(P)``.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt

import wind_common as wc
import dryden_turbulence_selection as dts


# --------------------------------------------------------------------------- #
# Correlation / strength helpers
# --------------------------------------------------------------------------- #
def autocorrelation(height, strength, d):
    """Dryden auto-correlation ``(Ru, Rv, Rw)`` at separation ``d`` (m)."""
    Lu, Lv, Lw = dts.dryden_length_scale(height)
    d = np.asarray(d, dtype=float)
    Ru = np.exp(-d / Lu)
    Rv = np.exp(-d / Lv) * (1.0 - d / (2.0 * Lv))
    Rw = np.exp(-d / Lw) * (1.0 - d / (2.0 * Lw))
    return Ru, Rv, Rw


def equivalent_strength(height, strength, d):
    """Conditional strength ``sigma sqrt(1 - R^2)`` for each component."""
    sigma_u, sigma_v, sigma_w = dts.dryden_turbulence_speed(height, strength)
    Ru, Rv, Rw = autocorrelation(height, strength, d)
    e_u = sigma_u * np.sqrt(np.clip(1.0 - Ru ** 2, 0.0, None))
    e_v = sigma_v * np.sqrt(np.clip(1.0 - Rv ** 2, 0.0, None))
    e_w = sigma_w * np.sqrt(np.clip(1.0 - Rw ** 2, 0.0, None))
    return e_u, e_v, e_w


def separation_grid(height, num=500, span=6.0):
    """Log-spaced separation distances ``d`` from 0.01 m up to ``span * Lmax``."""
    Lu, Lv, Lw = dts.dryden_length_scale(height)
    Lmax = max(Lu, Lv, Lw)
    return np.logspace(np.log10(0.01), np.log10(span * Lmax), num)


def _percentile(confidence):
    """Two-sided confidence ``C`` -> cumulative percentile ``(1 + C) / 2``."""
    return (1.0 + confidence) / 2.0


# --------------------------------------------------------------------------- #
# Prediction functions
# --------------------------------------------------------------------------- #
def gust_amplitude(confidence, height, strength, d):
    """
    Gust amplitude (m/s) versus separation ``d``.

    Returns ``(U0, V0, W0)`` for zero initial velocity and ``(U2, V2, W2)``
    for a non-zero initial velocity (conditioned on the zero-velocity value).
    """
    P = _percentile(confidence)
    e_u, e_v, e_w = equivalent_strength(height, strength, d)
    Ru, Rv, Rw = autocorrelation(height, strength, d)

    U0 = norm.ppf(P, loc=0.0, scale=e_u)
    V0 = norm.ppf(P, loc=0.0, scale=e_v)
    W0 = norm.ppf(P, loc=0.0, scale=e_w)

    U2 = norm.ppf(P, loc=Ru * U0, scale=e_u)
    V2 = norm.ppf(P, loc=Rv * V0, scale=e_v)
    W2 = norm.ppf(P, loc=Rw * W0, scale=e_w)

    return U0, V0, W0, U2, V2, W2


def gust_length_ratio(height, strength, d,
                      confidences=(0.99, 0.98, 0.95, 0.90, 0.80)):
    """
    Non-dimensional amplitude ``U/sigma`` versus ``d/L`` for several
    confidence levels.  Returns dicts keyed by confidence level.
    """
    Lu, Lv, Lw = dts.dryden_length_scale(height)
    sigma_u, sigma_v, sigma_w = dts.dryden_turbulence_speed(height, strength)
    e_u, e_v, e_w = equivalent_strength(height, strength, d)

    dLu, dLv, dLw = d / Lu, d / Lv, d / Lw

    result = {}
    for C in confidences:
        z = norm.ppf(_percentile(C))
        result[C] = (
            z * e_u / sigma_u,
            z * e_v / sigma_v,
            z * e_w / sigma_w,
        )
    return (dLu, dLv, dLw, result)


# --------------------------------------------------------------------------- #
# Plotting + export
# --------------------------------------------------------------------------- #
def wind_model_prection_output(percentage, height, strength, airspeed, num,
                               out_dir=None):
    """
    Generate and save the three prediction figures.

    Parameters
    ----------
    percentage : float
        Confidence coefficient in [0, 1] (e.g. 0.99 for 99 %).
    num : int
        Number of separation points.
    """
    out_dir = out_dir or wc.ensure_output_dir()
    confidence = float(percentage)

    d = separation_grid(height, num=int(num))
    Lu, Lv, Lw = dts.dryden_length_scale(height)
    dLu, dLv, dLw = d / Lu, d / Lv, d / Lw
    Ru, Rv, Rw = autocorrelation(height, strength, d)

    # --- Fig 1: relativity (correlation) function ------------------------ #
    fig1, ax = plt.subplots(figsize=(8, 6))
    ax.plot(dLu, Ru, "-r", label="Longitudinal wind")
    ax.plot(dLv, Rv, "--b", label="Lateral wind")
    ax.plot(dLw, Rw, "--g", label="Vertical wind")
    ax.set_xlabel("d / L")
    ax.set_ylabel("R (correlation)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig1.suptitle("Relativity Function Diagram")
    fig1.tight_layout()
    wc.save_figure(fig1, "Relativity_Function_Diagram.png", out_dir, dpi=200)

    # --- Fig 2: amplitude prediction curve ------------------------------- #
    U0, V0, W0, U2, V2, W2 = gust_amplitude(confidence, height, strength, d)
    fig2, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, (A0, A2, ylabel) in zip(
            axes,
            ((U0, U2, "Um (m/s)"), (V0, V2, "Vm (m/s)"), (W0, W2, "Wm (m/s)"))):
        ax.semilogx(d, A0, "--r", label="Initial velocity = 0")
        ax.semilogx(d, A2, "-b", label="Initial velocity != 0")
        ax.set_xlabel("dm (m)")
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(alpha=0.3, which="both")
    fig2.suptitle("Gust Wind Amplitude Prediction Curve")
    fig2.tight_layout()
    wc.save_figure(fig2, "Gust_Wind_Amplitude_Prediction_Curve.png",
                   out_dir, dpi=200)

    # --- Fig 3: non-dimensional amplitude vs confidence ------------------ #
    dLu, dLv, dLw, ratios = gust_length_ratio(height, strength, d)
    labels = {0.99: "99%", 0.98: "98%", 0.95: "95%",
              0.90: "90%", 0.80: "80%"}
    styles = ["--r", "--b", "--g", "-k", "-m"]
    fig3, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, (dL, comp, ylabel) in zip(
            axes,
            ((dLu, 0, "Um / sigma"), (dLv, 1, "Vm / sigma"),
             (dLw, 2, "Wm / sigma"))):
        for (C, data), st in zip(sorted(ratios.items(), reverse=True), styles):
            ax.semilogx(dL, data[comp], st, label=labels.get(C, f"{C:.0%}"))
        ax.set_xlabel("dm / L")
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(alpha=0.3, which="both")
    fig3.suptitle("Non-dimensional Wind Amplitude in Different Confidence Levels")
    fig3.tight_layout()
    wc.save_figure(fig3, "Gust_Wind_Amplitude_Ratio.png", out_dir, dpi=200)

    # --- export ----------------------------------------------------------- #
    wc.save_csv({"d_m": d, "Ru": Ru, "Rv": Rv, "Rw": Rw},
                out_dir, "correlation_functions.csv")

    return [fig1, fig2, fig3]
