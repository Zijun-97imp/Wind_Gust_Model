# -*- coding: utf-8 -*-
"""
von_karman_model.py
===================
Continuous **Von Karman** small-scale turbulence model (MIL-HDBK-1797).

Identical structure to the Dryden model but using the Von Karman rational
shaping filters and the Von Karman power spectral densities.

Shaping filters (SI; ``tau = L/V``)::

    H_u(s) = sigma_u sqrt(2 tau_u/pi) (1 + 0.25 tau_u s)
             / (1 + 1.357 tau_u s + 0.1987 tau_u^2 s^2)

    H_v(s) = sigma_v sqrt(tau_v/pi)
             (1 + 2.7478 tau_v s + 0.3398 tau_v^2 s^2)
             / (1 + 2.9958 tau_v s + 1.9754 tau_v^2 s^2 + 0.1539 tau_v^3 s^3)

    H_w(s) = sigma_w sqrt(tau_w/pi)
             (1 + 2.7478 tau_w s + 0.3398 tau_w^2 s^2)
             / (1 + 2.9958 tau_w s + 1.9754 tau_w^2 s^2 + 0.1539 tau_w^3 s^3)

Power spectral densities (one-sided, ``omega`` in rad/s)::

    Phi_u = sigma_u^2 (2 L_u/(pi V)) / (1 + (1.339 L_u omega/V)^2)^(5/6)
    Phi_v = sigma_v^2 (2 L_v/(pi V)) (1 + (8/3)(1.339 L_v omega/V)^2)
            / (1 + (1.339 L_v omega/V)^2)^(11/6)
    Phi_w = sigma_w^2 (2 L_w/(pi V)) (1 + (8/3)(1.339 L_w omega/V)^2)
            / (1 + (1.339 L_w omega/V)^2)^(11/6)
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy import signal
from matplotlib import pyplot as plt

import wind_common as wc
import von_karman_turbulence_selection as vts


# --------------------------------------------------------------------------- #
# Transfer functions (shaping filters)
# --------------------------------------------------------------------------- #
def vk_transfer_function(height, strength, airspeed):
    """Return the Von Karman shaping filters ``(H_u, H_v, H_w)`` (SI units)."""
    Lu, Lv, Lw = vts.vonkarman_length_scale(height)
    sigma_u, sigma_v, sigma_w = vts.vonkarman_turbulence_speed(height, strength)
    V = float(airspeed)

    tau_u = Lu / V
    tau_v = Lv / V
    tau_w = Lw / V

    # H_u
    a = sigma_u * np.sqrt(2.0 * tau_u / np.pi)
    H_u = signal.TransferFunction(
        [0.25 * a * tau_u, a],
        [0.1987 * tau_u ** 2, 1.357 * tau_u, 1.0],
    )

    # H_v
    b = sigma_v * np.sqrt(tau_v / np.pi)
    H_v = signal.TransferFunction(
        [0.3398 * b * tau_v ** 2, 2.7478 * b * tau_v, b],
        [0.1539 * tau_v ** 3, 1.9754 * tau_v ** 2, 2.9958 * tau_v, 1.0],
    )

    # H_w
    c = sigma_w * np.sqrt(tau_w / np.pi)
    H_w = signal.TransferFunction(
        [0.3398 * c * tau_w ** 2, 2.7478 * c * tau_w, c],
        [0.1539 * tau_w ** 3, 1.9754 * tau_w ** 2, 2.9958 * tau_w, 1.0],
    )

    return H_u, H_v, H_w


# --------------------------------------------------------------------------- #
# Theoretical power spectral density
# --------------------------------------------------------------------------- #
def vk_power_spectral_density(height, strength, airspeed, omega=None,
                              num_freq=1000):
    """Theoretical one-sided Von Karman PSD (see module docstring)."""
    Lu, Lv, Lw = vts.vonkarman_length_scale(height)
    sigma_u, sigma_v, sigma_w = vts.vonkarman_turbulence_speed(height, strength)
    V = float(airspeed)

    if omega is None:
        omega = np.logspace(-3, 3, num_freq)
    omega = np.asarray(omega, dtype=float)

    x_u = 1.339 * Lu * omega / V
    x_v = 1.339 * Lv * omega / V
    x_w = 1.339 * Lw * omega / V

    phi_u = (2.0 * sigma_u ** 2 * Lu / (np.pi * V)) / (1.0 + x_u ** 2) ** (5.0 / 6.0)
    phi_v = ((2.0 * sigma_v ** 2 * Lv / (np.pi * V))
             * (1.0 + (8.0 / 3.0) * x_v ** 2)
             / (1.0 + x_v ** 2) ** (11.0 / 6.0))
    phi_w = ((2.0 * sigma_w ** 2 * Lw / (np.pi * V))
             * (1.0 + (8.0 / 3.0) * x_w ** 2)
             / (1.0 + x_w ** 2) ** (11.0 / 6.0))

    return omega, phi_u, phi_v, phi_w


# --------------------------------------------------------------------------- #
# Time-domain turbulence generation
# --------------------------------------------------------------------------- #
def _normalize_to_intensity(x, sigma):
    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)
    std = np.std(x)
    if std > 0:
        x = x * (sigma / std)
    return x


def von_karman_wind_velocities(time, height, strength, airspeed,
                               num_sample=2000, seed=0):
    """Generate Von Karman turbulence velocity time series (u, v, w in m/s)."""
    num_sample = int(num_sample)
    t = np.linspace(0.0, float(time), num_sample)
    H_u, H_v, H_w = vk_transfer_function(height, strength, airspeed)
    sigma_u, sigma_v, sigma_w = vts.vonkarman_turbulence_speed(height, strength)

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
# Plotting + export
# --------------------------------------------------------------------------- #
def von_karman_plottings_output(frequency, time, height, strength, airspeed,
                                num_sample=2000, seed=0, out_dir=None):
    """Generate and save the Von Karman model figures.  Returns list of figs."""
    out_dir = out_dir or wc.ensure_output_dir()

    t, u, v, w = von_karman_wind_velocities(time, height, strength, airspeed,
                                            num_sample=num_sample, seed=seed)
    omega, phi_u, phi_v, phi_w = vk_power_spectral_density(height, strength,
                                                            airspeed)

    fig1, axes = plt.subplots(3, 1, figsize=(12, 15))
    for ax, y, name in zip(axes, (u, v, w),
                           ("Longitudinal", "Lateral", "Vertical")):
        ax.plot(t, y, lw=0.8)
        ax.set_ylabel(f"{name} Wind Speed (m/s)")
        ax.set_xlabel("Time (s)")
        ax.set_title(f"{name} wind")
        ax.grid(alpha=0.3)
    fig1.suptitle("Von Karman Gust Wind Velocity", fontsize=16)
    fig1.tight_layout(rect=[0, 0, 1, 0.97])
    wc.save_figure(fig1, "Von_Karman_Gust_Wind_Model.png", out_dir, dpi=200)

    fig2, axes = plt.subplots(3, 1, figsize=(12, 15))
    for ax, phi, name in zip(axes, (phi_u, phi_v, phi_w),
                             ("Longitudinal", "Lateral", "Vertical")):
        ax.loglog(omega, phi, lw=1.2)
        ax.set_ylabel(r"$\Phi$ (m$^2$/s per rad/s)")
        ax.set_xlabel(r"$\omega$ (rad/s)")
        ax.set_title(f"{name} wind")
        ax.grid(alpha=0.3, which="both")
    fig2.suptitle("Von Karman Gust Wind Power Spectral Density", fontsize=16)
    fig2.tight_layout(rect=[0, 0, 1, 0.97])
    wc.save_figure(fig2, "Von_Karman_Gust_Wind_Model_PSD.png", out_dir, dpi=200)

    wc.save_csv({"time_s": t, "u_mps": u, "v_mps": v, "w_mps": w},
                out_dir, "von_karman_turbulence.csv")

    return [fig1, fig2]
