# -*- coding: utf-8 -*-
"""
multi_scale_model.py
====================
Multi-scale wind-field superposition (Eq. (9) of the paper).

The total wind field is the linear superposition, in Cartesian coordinates, of

* large-scale mean wind ``Vm`` (altitude-dependent via a mean-wind profile),
* mid-scale vertical wind shear ``V_shear`` (wind-speed change across a
  vertical extent around the flight altitude),
* the measured large-scale low-frequency fluctuation ``V_fluc(t)``,
* mesoscale discrete (``1-cos``) gust ``Vg``,
* small-scale random turbulence ``Vturb`` (Dryden or Von Karman).

::

    V_w(t) = Vm(z) + V_shear(z) + V_fluc(t) + Vg(t) + Vturb(t)

Each component is returned separately so the contributions can be inspected.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from matplotlib import pyplot as plt

import wind_common as wc

from Dryden_Gust_Model import dryden_model as dm
from Von_Karman_Gust_Model import von_karman_model as vkm
from Mesoscale_Gust import meso_model as meso
from Largescale_Gust import large_model as large


def multi_scale_wind(time, height, strength, airspeed, model="dryden",
                     num_sample=2000, seed=0,
                     meso_um=0.0, meso_dm=0.0, meso_dir="u", t_start=0.0,
                     vm=None, use_large_scale_data=True,
                     use_shear=False, shear_extent=10.0,
                     ref_height=10.0, ref_speed=None, alpha=0.15,
                     shear_dir="u", use_large_scale_series=False):
    """
    Build the multi-scale wind field.

    Parameters
    ----------
    time : float
        Simulation duration (s).
    height : float
        Flight altitude (m).
    strength : int
        1 = light, 2 = moderate, 3 = severe.
    airspeed : float
        Flight airspeed (m/s).
    model : {'dryden', 'von_karman'}
        Small-scale turbulence model.
    num_sample : int
        Number of time points.
    seed : int
        RNG seed for the turbulence component.
    meso_um : float
        Mesoscale gust amplitude (m/s); 0 disables the gust.
    meso_dm : float
        Mesoscale gust length (m).
    meso_dir : {'u', 'v', 'w'}
        Direction along which the mesoscale gust acts.
    t_start : float
        Entry time of the mesoscale gust (s).
    vm : tuple (Vm_u, Vm_v, Vm_w), optional
        Explicit large-scale mean wind.  Overrides the data/profile logic.
    use_large_scale_data : bool
        Read the large-scale mean wind from ``wind_data.csv``.
    use_shear : bool
        Scale the large-scale mean onto a power-law vertical profile and add
        the vertical wind-shear component.
    shear_extent : float
        Vertical extent (m) over which the wind shear is evaluated.
    ref_height, ref_speed, alpha :
        Power-law profile ``V(z) = V_ref (z/z_ref)^alpha``.  ``ref_speed``
        defaults to the magnitude of the measured mean wind at 10 m.
    shear_dir : {'u', 'v', 'w'}
        Axis along which the wind shear acts.
    use_large_scale_series : bool
        Add the measured low-frequency large-scale fluctuation time series.

    Returns
    -------
    dict with keys ``t, u, v, w`` (totals) and ``components`` (a dict of the
    individual contributions keyed by name).
    """
    num_sample = int(num_sample)
    t = np.linspace(0.0, float(time), num_sample)

    # 1) small-scale turbulence ------------------------------------------- #
    if model.lower() == "von_karman":
        _, tu, tv, tw = vkm.von_karman_wind_velocities(
            time, height, strength, airspeed, num_sample=num_sample, seed=seed)
    else:
        _, tu, tv, tw = dm.dryden_wind_velocities(
            time, height, strength, airspeed, num_sample=num_sample, seed=seed)

    # 2) mesoscale 1-cos gust --------------------------------------------- #
    gu = gv = gw = np.zeros(num_sample)
    if meso_um > 0 and meso_dm > 0:
        gust = meso.meso_gust_time_series(t, meso_um, meso_dm, airspeed,
                                          t_start=t_start)
        if meso_dir == "u":
            gu = gust
        elif meso_dir == "v":
            gv = gust
        else:
            gw = gust

    # 3) large-scale mean + wind shear + fluctuation ---------------------- #
    # Start from the measured mean wind (constant).
    if vm is not None:
        Vm_u, Vm_v, Vm_w = (float(vm[0]), float(vm[1]), float(vm[2]))
    elif use_large_scale_data:
        Vm_u, Vm_v, Vm_w = large.large_scale_mean_wind()
    else:
        Vm_u = Vm_v = Vm_w = 0.0

    # Wind shear / altitude profile scaling of the mean wind.
    dV_u = dV_v = dV_w = 0.0
    if use_shear:
        mag = float(np.hypot(Vm_u, np.hypot(Vm_v, Vm_w)))
        ref = ref_speed if ref_speed is not None else mag
        if ref > 1e-9 and mag > 1e-9:
            factor = float(large.wind_profile_power_law(
                height, ref_height, ref, alpha) / ref)
            Vm_u, Vm_v, Vm_w = factor * Vm_u, factor * Vm_v, factor * Vm_w
        # Shear-induced wind-speed difference across the vertical extent.
        dV_u, dV_v, dV_w = large.vertical_wind_shear(
            height, extent=shear_extent, ref_height=ref_height,
            ref_speed=(ref if ref is not None else 1.0),
            alpha=alpha, direction=shear_dir)

    # Measured low-frequency fluctuation (zero-mean by construction).
    fluc_u = fluc_v = fluc_w = np.zeros(num_sample)
    if use_large_scale_series and use_large_scale_data:
        t_data, u_d, v_d, w_d = large.load_large_scale_data()
        fluc_u = np.interp(t, t_data, u_d) - float(np.mean(u_d))
        fluc_v = np.interp(t, t_data, v_d) - float(np.mean(v_d))
        fluc_w = np.interp(t, t_data, w_d) - float(np.mean(w_d))

    # 4) superposition ------------------------------------------------------ #
    u = Vm_u + dV_u + fluc_u + gu + tu
    v = Vm_v + dV_v + fluc_v + gv + tv
    w = Vm_w + dV_w + fluc_w + gw + tw

    return {
        "t": t,
        "u": u, "v": v, "w": w,
        "components": {
            "large_scale_mean": (Vm_u, Vm_v, Vm_w),
            "wind_shear": (dV_u, dV_v, dV_w),
            "large_scale_fluctuation": (fluc_u, fluc_v, fluc_w),
            "turbulence": (tu, tv, tw),
            "mesoscale_gust": (gu, gv, gw),
        },
    }


def multi_scale_plottings_output(time, height, strength, airspeed,
                                 model="dryden", num_sample=2000, seed=0,
                                 meso_um=0.0, meso_dm=0.0, meso_dir="u",
                                 t_start=0.0, vm=None, out_dir=None,
                                 use_shear=False, shear_extent=10.0,
                                 ref_height=10.0, ref_speed=None, alpha=0.15,
                                 shear_dir="u", use_large_scale_series=False):
    """Plot and export the multi-scale wind field.  Returns a list of figs."""
    out_dir = out_dir or wc.ensure_output_dir()
    res = multi_scale_wind(time, height, strength, airspeed, model=model,
                           num_sample=num_sample, seed=seed,
                           meso_um=meso_um, meso_dm=meso_dm, meso_dir=meso_dir,
                           t_start=t_start, vm=vm,
                           use_large_scale_data=(vm is None),
                           use_shear=use_shear, shear_extent=shear_extent,
                           ref_height=ref_height, ref_speed=ref_speed,
                           alpha=alpha, shear_dir=shear_dir,
                           use_large_scale_series=use_large_scale_series)
    t = res["t"]
    comp = res["components"]
    idx = {"u": 0, "v": 1, "w": 2}

    fig, axes = plt.subplots(3, 1, figsize=(12, 15), sharex=True)
    for ax, key, name in zip(axes, ("u", "v", "w"),
                             ("Longitudinal", "Lateral", "Vertical")):
        i = idx[key]
        ax.plot(t, comp["turbulence"][i], alpha=0.6, lw=0.7,
                label="Small-scale turbulence")
        if np.any(np.abs(comp["large_scale_fluctuation"][i]) > 1e-12):
            ax.plot(t, comp["large_scale_mean"][i]
                    + comp["large_scale_fluctuation"][i],
                    "b", alpha=0.7, lw=0.8, label="Large-scale field")
        ax.plot(t, res[key], "k", lw=1.0, label="Total wind")
        ax.axhline(comp["large_scale_mean"][i], color="r", ls="--",
                   label="Large-scale mean")
        if abs(comp["wind_shear"][i]) > 1e-12:
            ax.axhline(comp["large_scale_mean"][i] + comp["wind_shear"][i],
                       color="m", ls=":", lw=1.2, label="Mean + wind shear")
        ax.set_ylabel(f"{name} wind (m/s)")
        ax.legend(loc="upper right", fontsize="small")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle(f"Multi-scale Wind Field (h={height} m, {model}, "
                 f"strength={strength})")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    wc.save_figure(fig, "Multi_scale_Wind_Field.png", out_dir, dpi=200)

    wc.save_csv({"time_s": t, "u_mps": res["u"], "v_mps": res["v"],
                 "w_mps": res["w"]},
                out_dir, "multi_scale_wind.csv")

    return [fig]
