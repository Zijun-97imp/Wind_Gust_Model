# -*- coding: utf-8 -*-
"""
meso_model.py
=============
Mesoscale **discrete gust** model (``1-cos`` profile).

The mesoscale gust is the medium-scale discrete gust used to represent wind
shear / high-intensity gusts for UAV disturbance studies.  Its spatial shape
is the ``1-cos`` profile (MIL-HDBK-1797, "Discrete Gust")::

    U(x) = 0,                                   x < 0
    U(x) = (Um/2) (1 - cos(pi x / dm)),         0 <= x <= 2 dm
    U(x) = 0,                                   x > 2 dm

where ``Um`` is the gust amplitude and ``dm`` the gust length (the gust rises
from 0 to ``Um`` and returns to 0 over a distance ``2 dm``).  In the time
domain ``x = V t``::

    Vg(t) = (Um/2) (1 - cos(pi V t / dm)),      0 <= t <= 2 dm / V
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from matplotlib import pyplot as plt

import wind_common as wc


def one_cosine_gust(x, Um, dm):
    """Spatial ``1-cos`` gust profile at distances ``x`` (m)."""
    x = np.asarray(x, dtype=float)
    gust = np.zeros_like(x)
    inside = (x >= 0) & (x <= 2.0 * dm)
    gust[inside] = (Um / 2.0) * (1.0 - np.cos(np.pi * x[inside] / dm))
    return gust


def meso_gust_time_series(t, Um, dm, V, t_start=0.0):
    """Mesoscale ``1-cos`` gust velocity (m/s) along the flight track.

    Parameters
    ----------
    t : array_like
        Time vector (s).
    Um : float
        Gust amplitude (m/s).
    dm : float
        Gust length (m).
    V : float
        Flight airspeed (m/s).
    t_start : float
        Time at which the aircraft enters the gust (s).
    """
    t = np.asarray(t, dtype=float)
    x = V * (t - t_start)
    return one_cosine_gust(x, Um, dm)


def mesoscale_plottings_output(Um, dm, V, time=None, t_start=0.0,
                               num_sample=2000, out_dir=None):
    """Plot and export the mesoscale ``1-cos`` gust.  Returns a list of figs."""
    out_dir = out_dir or wc.ensure_output_dir()
    if time is None:
        time = 4.0 * dm / V + t_start

    t = np.linspace(0.0, float(time), int(num_sample))
    vg = meso_gust_time_series(t, Um, dm, V, t_start=t_start)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(t, vg, "r", lw=1.2)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Gust velocity (m/s)")
    ax1.set_title("Mesoscale 1-cos gust (time domain)")
    ax1.grid(alpha=0.3)

    x = V * (t - t_start)
    ax2.plot(x, vg, "b", lw=1.2)
    ax2.set_xlabel("Distance (m)")
    ax2.set_ylabel("Gust velocity (m/s)")
    ax2.set_title("Mesoscale 1-cos gust (spatial domain)")
    ax2.grid(alpha=0.3)

    fig.suptitle(f"Mesoscale Gust (Um={Um} m/s, dm={dm} m)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    wc.save_figure(fig, "Mesoscale_1cos_Gust.png", out_dir, dpi=200)

    wc.save_csv({"time_s": t, "distance_m": x, "gust_mps": vg},
                out_dir, "mesoscale_gust.csv")

    return [fig]
