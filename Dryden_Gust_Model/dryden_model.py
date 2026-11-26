# -*- coding: utf-8 -*-
"""
dryden_model.py
===============
Discrete/continuous **Dryden** small-scale turbulence model (MIL-HDBK-1797).

This module builds the Dryden *shaping filters* (transfer functions) for the
three wind components, generates turbulence time series by filtering white
noise, evaluates the theoretical power spectral densities (PSD), and produces
plots + CSV export.

Wind components
---------------
``u`` : longitudinal   ``v`` : lateral   ``w`` : vertical

Shaping filters (SI, ``V`` = airspeed, ``L`` = scale length, ``sigma`` =
turbulence intensity)::

    H_u(s) = sigma_u sqrt(2 L_u/(pi V)) / (1 + (L_u/V) s)
    H_v(s) = sigma_v sqrt(  L_v/(pi V)) (1 + sqrt(3)(L_v/V) s) / (1 + (L_v/V) s)^2
    H_w(s) = sigma_w sqrt(  L_w/(pi V)) (1 + sqrt(3)(L_w/V) s) / (1 + (L_w/V) s)^2

Power spectral densities (one-sided, ``omega`` in rad/s)::

    Phi_u(omega) = sigma_u^2 (2 L_u/(pi V)) / (1 + (L_u omega/V)^2)
    Phi_v(omega) = sigma_v^2 (2 L_v/(pi V)) (1 + 3 (L_v omega/V)^2) / (1 + (L_v omega/V)^2)^2
    Phi_w(omega) = sigma_w^2 (2 L_w/(pi V)) (1 + 3 (L_w omega/V)^2) / (1 + (L_w omega/V)^2)^2
"""

import os
import sys

# Make `wind_common` importable from the project root.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy import signal
from matplotlib import pyplot as plt

import wind_common as wc
import dryden_turbulence_selection as dts


# --------------------------------------------------------------------------- #
# Transfer functions (shaping filters)
# --------------------------------------------------------------------------- #
def dryden_transfer_function(height, strength, airspeed):
    """
    Return the Dryden shaping filters ``(H_u, H_v, H_w)`` as
    :class:`scipy.signal.TransferFunction` objects (all in SI units).
    """
    Lu, Lv, Lw = dts.dryden_length_scale(height)
    sigma_u, sigma_v, sigma_w = dts.dryden_turbulence_speed(height, strength)
    V = float(airspeed)

    # H_u = sigma_u sqrt(2 Lu/(pi V)) / (1 + (Lu/V) s)
    a_u = sigma_u * np.sqrt(2.0 * Lu / (np.pi * V))
    H_u = signal.TransferFunction([a_u], [Lu / V, 1.0])

    # H_v = sigma_v sqrt(Lv/(pi V)) (1 + sqrt(3)(Lv/V) s) / (1 + (Lv/V) s)^2
    b = sigma_v * np.sqrt(Lv / (np.pi * V))
    tau_v = Lv / V
    H_v = signal.TransferFunction(
        [np.sqrt(3.0) * tau_v * b, b],
        [tau_v ** 2, 2.0 * tau_v, 1.0],
    )

    # H_w = sigma_w sqrt(Lw/(pi V)) (1 + sqrt(3)(Lw/V) s) / (1 + (Lw/V) s)^2
    c = sigma_w * np.sqrt(Lw / (np.pi * V))
    tau_w = Lw / V
    H_w = signal.TransferFunction(
        [np.sqrt(3.0) * tau_w * c, c],
        [tau_w ** 2, 2.0 * tau_w, 1.0],
    )

    return H_u, H_v, H_w


# --------------------------------------------------------------------------- #
# Theoretical power spectral density
# --------------------------------------------------------------------------- #
def dryden_power_spectral_density(height, strength, airspeed, omega=None,
                                  num_freq=1000):
    """
    Theoretical one-sided Dryden PSD.

    Parameters
    ----------
    height, strength, airspeed :
        Turbulence scenario.
    omega : array_like, optional
        Angular frequencies (rad/s).  If ``None`` a logarithmic grid from
        1e-3 to 1e3 rad/s is used.
    num_freq : int
        Number of frequency points when ``omega`` is not supplied.

    Returns
    -------
    omega : ndarray
    phi_u, phi_v, phi_w : ndarray
        PSD values in (m/s)^2 per (rad/s).
    """
    Lu, Lv, Lw = dts.dryden_length_scale(height)
    sigma_u, sigma_v, sigma_w = dts.dryden_turbulence_speed(height, strength)
    V = float(airspeed)

    if omega is None:
        omega = np.logspace(-3, 3, num_freq)

    omega = np.asarray(omega, dtype=float)
    x_u = Lu * omega / V
    x_v = Lv * omega / V
    x_w = Lw * omega / V

    phi_u = (2.0 * sigma_u ** 2 * Lu / (np.pi * V)) / (1.0 + x_u ** 2)
    phi_v = ((2.0 * sigma_v ** 2 * Lv / (np.pi * V))
             * (1.0 + 3.0 * x_v ** 2) / (1.0 + x_v ** 2) ** 2)
    phi_w = ((2.0 * sigma_w ** 2 * Lw / (np.pi * V))
             * (1.0 + 3.0 * x_w ** 2) / (1.0 + x_w ** 2) ** 2)

    return omega, phi_u, phi_v, phi_w


# --------------------------------------------------------------------------- #
# Time-domain turbulence generation
# --------------------------------------------------------------------------- #
def _normalize_to_intensity(x, sigma):
    """Remove the mean and rescale so that ``std(x) == sigma``."""
    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)
    std = np.std(x)
    if std > 0:
        x = x * (sigma / std)
    return x


def dryden_wind_velocities(time, height, strength, airspeed, num_sample=2000,
                           seed=0):
    """
    Generate Dryden turbulence velocity time series.

    White noise is shaped by the Dryden filters and then normalised so that
    each component's standard deviation equals the target turbulence intensity
    ``sigma`` (this guarantees a correct overall intensity while preserving the
    Dryden spectral shape).

    Returns
    -------
    t : ndarray  (num_sample,)
    u, v, w : ndarray  (num_sample,)
        Wind components in m/s.
    """
    num_sample = int(num_sample)
    t = np.linspace(0.0, float(time), num_sample)
    H_u, H_v, H_w = dryden_transfer_function(height, strength, airspeed)
    sigma_u, sigma_v, sigma_w = dts.dryden_turbulence_speed(height, strength)

    rng = np.random.default_rng(seed)
    u_in = rng.standard_normal(num_sample)
    v_in = rng.standard_normal(num_sample)
    w_in = rng.standard_normal(num_sample)

    _, u, _ = signal.lsim(H_u, U=u_in, T=t)
    _, v, _ = signal.lsim(H_v, U=v_in, T=t)
    _, w, _ = signal.lsim(H_w, U=w_in, T=t)

    u = _normalize_to_intensity(u, sigma_u)
    v = _normalize_to_intensity(v, sigma_v)
    w = _normalize_to_intensity(w, sigma_w)

    return t, u, v, w


# --------------------------------------------------------------------------- #
# Plotting + export (keeps the historical `dryden_plottings_output` signature)
# --------------------------------------------------------------------------- #
def dryden_plottings_output(frequency, time, height, strength, airspeed,
                            num_sample=2000, seed=0, out_dir=None):
    """
    Generate and save the Dryden model figures.

    Returns
    -------
    list of matplotlib.figure.Figure
    """
    out_dir = out_dir or wc.ensure_output_dir()

    t, u, v, w = dryden_wind_velocities(time, height, strength, airspeed,
                                        num_sample=num_sample, seed=seed)
    omega, phi_u, phi_v, phi_w = dryden_power_spectral_density(
        height, strength, airspeed)

    # --- time domain ----------------------------------------------------- #
    fig1, axes = plt.subplots(3, 1, figsize=(12, 15))
    for ax, y, name in zip(axes, (u, v, w), ("Longitudinal", "Lateral", "Vertical")):
        ax.plot(t, y, lw=0.8)
        ax.set_ylabel(f"{name} Wind Speed (m/s)")
        ax.set_xlabel("Time (s)")
        ax.set_title(f"{name} wind")
        ax.grid(alpha=0.3)
    fig1.suptitle("Dryden Gust Wind Velocity", fontsize=16)
    fig1.tight_layout(rect=[0, 0, 1, 0.97])
    wc.save_figure(fig1, "Dryden_Gust_Wind_Model.png", out_dir, dpi=200)

    # --- frequency domain (PSD) ------------------------------------------ #
    fig2, axes = plt.subplots(3, 1, figsize=(12, 15))
    for ax, phi, name in zip(axes, (phi_u, phi_v, phi_w),
                             ("Longitudinal", "Lateral", "Vertical")):
        ax.loglog(omega, phi, lw=1.2)
        ax.set_ylabel(r"$\Phi$ (m$^2$/s per rad/s)")
        ax.set_xlabel(r"$\omega$ (rad/s)")
        ax.set_title(f"{name} wind")
        ax.grid(alpha=0.3, which="both")
    fig2.suptitle("Dryden Gust Wind Power Spectral Density", fontsize=16)
    fig2.tight_layout(rect=[0, 0, 1, 0.97])
    wc.save_figure(fig2, "Dryden_Gust_Wind_Model_PSD.png", out_dir, dpi=200)

    # --- export CSV ------------------------------------------------------- #
    wc.save_csv({"time_s": t, "u_mps": u, "v_mps": v, "w_mps": w},
                out_dir, "dryden_turbulence.csv")

    return [fig1, fig2]
