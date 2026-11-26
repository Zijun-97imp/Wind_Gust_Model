# -*- coding: utf-8 -*-
"""
large_model.py
==============
Large-scale (mean / background) wind model built from real meteorological
data **plus a mean-wind vertical profile (wind shear)**.

The large-scale wind field is the slowly-varying background wind ``Vm`` of the
multi-scale model (the ``Vm`` and ``HV`` terms in Eq. (7)-(9) of the paper).
Its horizontal spatial characteristic scale is of the order of hundreds of
metres, and its vertical characteristic scale of the order of tens of metres,
so over a short flight segment it can be treated as a slowly-varying mean plus
a low-frequency trend; the paper explicitly notes that large-scale winds
"exhibit directional and vertical wind shear".

This module provides:

* the raw large-scale time series (``wind_data.csv``);
* the large-scale **mean** wind ``(Vm_u, Vm_v, Vm_w)``;
* an interpolation of the full large-scale time series onto an arbitrary grid;
* a **mean-wind vertical profile** (power law and logarithmic law) that makes
  the mean wind altitude-dependent, i.e. ``Vm = Vm(z)``;
* the associated **vertical wind shear** (the wind-speed difference over a
  vertical extent ``dz`` around the flight altitude).

Mean-wind profile (vertical wind shear)
---------------------------------------
The planetary boundary layer mean wind increases with altitude.  Two standard
profiles are implemented:

* Power law:  ``V(z) = V_ref (z / z_ref)^alpha``  (``alpha`` ~ 0.1-0.4).
* Log law  :  ``V(z) = (u* / kappa) ln(z / z0)``   (``kappa`` = 0.4).

The vertical wind shear across a vertical extent ``dz`` (e.g. an aircraft's
reference height) is the difference of the profile between ``z - dz/2`` and
``z + dz/2``.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from matplotlib import pyplot as plt

import wind_common as wc


# --------------------------------------------------------------------------- #
# Raw large-scale data
# --------------------------------------------------------------------------- #
def load_large_scale_data(csv_path=None):
    """
    Load the large-scale wind time series.

    Returns
    -------
    t, u, v, w : ndarray
        Time (s) and wind components (m/s).
    """
    csv_path = csv_path or wc.WIND_DATA_CSV
    data = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
    t = data[:, 0]
    u = data[:, 1]
    v = data[:, 2]
    w = data[:, 3]
    return t, u, v, w


def large_scale_mean_wind(csv_path=None):
    """Large-scale mean wind ``(Vm_u, Vm_v, Vm_w)`` in m/s (from data)."""
    t, u, v, w = load_large_scale_data(csv_path)
    return float(np.mean(u)), float(np.mean(v)), float(np.mean(w))


def large_scale_time_series(time, num_sample=2000, csv_path=None):
    """
    Large-scale wind interpolated onto a target time grid.

    Uses linear interpolation of the loaded data and holds the last value for
    times beyond the data record.
    """
    t_data, u, v, w = load_large_scale_data(csv_path)
    t_target = np.linspace(0.0, float(time), int(num_sample))

    u_out = np.interp(t_target, t_data, u)
    v_out = np.interp(t_target, t_data, v)
    w_out = np.interp(t_target, t_data, w)
    return t_target, u_out, v_out, w_out


# --------------------------------------------------------------------------- #
# Mean-wind vertical profile (wind shear)
# --------------------------------------------------------------------------- #
def wind_profile_power_law(z, ref_height=10.0, ref_speed=8.0, alpha=0.15):
    """
    Power-law mean-wind vertical profile.

    .. math:: V(z) = V_ref (z / z_ref)^alpha

    Parameters
    ----------
    z : array_like
        Altitude(s) in metres.
    ref_height : float
        Reference altitude ``z_ref`` (m).
    ref_speed : float
        Mean wind speed at ``z_ref`` (m/s).
    alpha : float
        Wind-shear exponent (0.1 = smooth terrain, 0.2 = open farmland,
        0.3+ = urban / forest).
    """
    z = np.asarray(z, dtype=float)
    return ref_speed * (np.maximum(z, 1e-3) / ref_height) ** alpha


def wind_profile_log_law(z, z0=0.03, ustar=0.30):
    """
    Logarithmic mean-wind vertical profile.

    .. math:: V(z) = (u_* / kappa) ln(z / z0)

    Parameters
    ----------
    z : array_like
        Altitude(s) in metres.
    z0 : float
        Aerodynamic roughness length (m): ~0.0002 open water, ~0.03 grass,
        ~1.0 urban.
    ustar : float
        Friction velocity (m/s).
    """
    kappa = 0.4
    z = np.asarray(z, dtype=float)
    return (ustar / kappa) * np.log(np.maximum(z, 1e-3) / z0)


def wind_shear_rate(z, ref_height=10.0, ref_speed=8.0, alpha=0.15):
    """
    Vertical wind-shear rate ``dV/dz`` (1/s) of the power-law profile at ``z``.

    .. math:: dV/dz = V_ref * alpha * z^(alpha-1) / z_ref^alpha
    """
    z = np.asarray(z, dtype=float)
    z = np.maximum(z, 1e-3)
    return ref_speed * alpha * (z / ref_height) ** (alpha - 1.0) / ref_height


def vertical_wind_shear(height, extent=10.0, ref_height=10.0, ref_speed=8.0,
                        alpha=0.15, law="power", direction="u"):
    """
    Vertical wind-shear velocity difference across a vertical extent.

    The wind-speed change experienced by an aircraft of vertical extent ``dz``
    (reference height / wingspan) centred at altitude ``height`` is the
    difference of the mean-wind profile between ``height - dz/2`` and
    ``height + dz/2``.  It is returned as a slowly-varying (quasi-steady)
    component along one Cartesian direction.

    Parameters
    ----------
    height : float
        Flight altitude (m).
    extent : float
        Vertical extent ``dz`` (m) over which the shear is evaluated.
    ref_height, ref_speed, alpha :
        Power-law profile parameters (see :func:`wind_profile_power_law`).
    law : {'power', 'log'}
        Profile law.  For ``'log'`` the arguments ``ref_height/ref_speed`` are
        ignored and ``z0``/``ustar`` (passed via ``alpha`` slot) are unused;
        use :func:`wind_profile_log_law` directly for full control.
    direction : {'u', 'v', 'w'}
        Cartesian axis along which the shear-induced wind difference acts.

    Returns
    -------
    (dV_u, dV_v, dV_w) : tuple of float
        Shear velocity components in m/s.
    """
    z_top = height + 0.5 * extent
    z_bot = height - 0.5 * extent

    if law == "log":
        V_top = wind_profile_log_law(z_top, z0=0.03, ustar=ref_speed)
        V_bot = wind_profile_log_law(z_bot, z0=0.03, ustar=ref_speed)
    else:
        V_top = wind_profile_power_law(z_top, ref_height, ref_speed, alpha)
        V_bot = wind_profile_power_law(z_bot, ref_height, ref_speed, alpha)

    dV = float(V_top - V_bot)
    return {"u": (dV, 0.0, 0.0), "v": (0.0, dV, 0.0),
            "w": (0.0, 0.0, dV)}[direction]


# --------------------------------------------------------------------------- #
# Full large-scale field (mean + fluctuation)
# --------------------------------------------------------------------------- #
def large_scale_field(time, height=None, num_sample=2000, csv_path=None,
                      ref_height=10.0, ref_speed=None, alpha=0.15,
                      use_shear=True):
    """
    Large-scale wind field as a function of time at altitude ``height``.

    Returns the full large-scale contribution used in the multi-scale
    superposition: the altitude-dependent mean wind (from the measured data,
    optionally scaled onto a power-law profile) plus the measured low-frequency
    fluctuation, i.e. ``Vm(z) + [V_data(t) - mean(V_data)]``.

    Parameters
    ----------
    time : float
        Target duration (s).
    height : float, optional
        Altitude (m).  If ``None`` the measured mean is used unchanged.
    ref_height, ref_speed, alpha :
        Power-law profile reference.  ``ref_speed`` defaults to the magnitude of
        the measured mean wind at 10 m so that the profile passes through the
        measured mean at the reference altitude.
    use_shear : bool
        If True, scale the mean wind onto the altitude profile; otherwise keep
        the measured mean constant.

    Returns
    -------
    t : ndarray
    (u, v, w) : ndarray
    """
    t_data, u, v, w = load_large_scale_data(csv_path)
    Vm_u, Vm_v, Vm_w = (float(np.mean(u)), float(np.mean(v)), float(np.mean(w)))
    mag = float(np.hypot(Vm_u, np.hypot(Vm_v, Vm_w)))

    t_target = np.linspace(0.0, float(time), int(num_sample))

    if height is not None and use_shear:
        ref_speed = ref_speed if ref_speed is not None else mag
        # Altitude-dependent magnitude of the mean wind.
        factor = float(wind_profile_power_law(height, ref_height, ref_speed,
                                              alpha) / max(ref_speed, 1e-9))
        mean_u, mean_v, mean_w = factor * Vm_u, factor * Vm_v, factor * Vm_w
    else:
        mean_u, mean_v, mean_w = Vm_u, Vm_v, Vm_w

    # Low-frequency fluctuation around the (possibly rescaled) mean.
    fluct_u = np.interp(t_target, t_data, u) - Vm_u
    fluct_v = np.interp(t_target, t_data, v) - Vm_v
    fluct_w = np.interp(t_target, t_data, w) - Vm_w

    return (t_target,
            mean_u + fluct_u, mean_v + fluct_v, mean_w + fluct_w)


# --------------------------------------------------------------------------- #
# Plotting + export
# --------------------------------------------------------------------------- #
def large_scale_plottings_output(time=None, num_sample=2000, out_dir=None,
                                 height=None, ref_height=10.0, ref_speed=8.0,
                                 alpha=0.15):
    """Plot and export the large-scale wind.  Returns a list of figs."""
    out_dir = out_dir or wc.ensure_output_dir()

    t, u, v, w = load_large_scale_data()
    Vm_u, Vm_v, Vm_w = large_scale_mean_wind()

    fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    for ax, y, name, mean in zip(axes, (u, v, w),
                                 ("Longitudinal", "Lateral", "Vertical"),
                                 (Vm_u, Vm_v, Vm_w)):
        ax.plot(t, y, lw=0.8)
        ax.axhline(mean, color="r", ls="--", label=f"mean = {mean:.3f} m/s")
        ax.set_ylabel(f"{name} wind (m/s)")
        ax.legend()
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle("Large-scale Wind (measured data)")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    wc.save_figure(fig, "Large_scale_Wind.png", out_dir, dpi=200)

    # --- wind profile / shear (if a height scan is requested) ------------- #
    if height is not None:
        z = np.linspace(1.0, max(float(height) * 1.5, 30.0), 200)
        V = wind_profile_power_law(z, ref_height, ref_speed, alpha)
        fig2, ax = plt.subplots(figsize=(7, 6))
        ax.plot(V, z, lw=1.6)
        ax.axhline(height, color="r", ls="--", label=f"flight h = {height:g} m")
        ax.set_xlabel("Mean wind speed V(z) (m/s)")
        ax.set_ylabel("Altitude z (m)")
        ax.set_title(r"Mean-wind vertical profile (power law, $\alpha$="
                     f"{alpha:g})")
        ax.legend()
        ax.grid(alpha=0.3)
        fig2.tight_layout()
        wc.save_figure(fig2, "Wind_Profile_Shear.png", out_dir, dpi=200)

        wc.save_csv({"altitude_m": z, "mean_wind_mps": V},
                    out_dir, "wind_profile.csv")
        return [fig, fig2]

    wc.save_csv({"time_s": t, "u_mps": u, "v_mps": v, "w_mps": w},
                out_dir, "large_scale_wind.csv")

    return [fig]
