# -*- coding: utf-8 -*-
"""
Main.py
=======
Command-line driver for the multi-scale gust wind model.

Run without arguments to be prompted for the model parameters interactively
(no need to edit the source code), or use ``--help`` to see all options.  The
program generates:

* Dryden and Von Karman small-scale turbulence (time series + PSD),
* the normal-distribution amplitude-correction (prediction) figures,
* the mesoscale ``1-cos`` gust,
* the large-scale (measured) wind,
* the multi-scale superposition (Eq. 9 of the paper),

and writes every figure and CSV table into the ``output`` directory.

Example
-------
    python Main.py --height 10 --airspeed 10 --time 100 --strength 3 \
        --confidence 0.99 --model both
"""

import argparse
import os
import sys

# Make the project root and every sub-package importable regardless of CWD.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _sub in ("", "Dryden_Gust_Model", "Von_Karman_Gust_Model",
             "Normal_Distribution_Model_Prediction", "Mesoscale_Gust",
             "Largescale_Gust"):
    _p = os.path.join(_ROOT, _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import wind_common as wc
import dryden_model as dm
import dryden_turbulence_selection as dts
import von_karman_model as vkm
import von_karman_turbulence_selection as vts
import normal_distribution_prediction as ndp
import meso_model as meso
import large_model as large
import multi_scale_model as msm


def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Multi-scale gust wind model (Dryden / Von Karman / "
                    "1-cos gust / large-scale + amplitude correction).",
    )
    p.add_argument("--height", type=float, default=500.0,
                   help="flight altitude in metres (default 500)")
    p.add_argument("--airspeed", type=float, default=10.0,
                   help="airspeed in m/s (default 10)")
    p.add_argument("--time", type=float, default=100.0,
                   help="simulation duration in s (default 100)")
    p.add_argument("--frequency", type=float, default=1.0,
                   help="reference gust frequency in Hz (default 1)")
    p.add_argument("--strength", type=int, default=2,
                   choices=[1, 2, 3],
                   help="turbulence strength 1=light 2=moderate 3=severe")
    p.add_argument("--confidence", type=float, default=0.99,
                   help="confidence coefficient in [0,1] (default 0.99)")
    p.add_argument("--num", type=int, default=500,
                   help="number of separation points for prediction (default 500)")
    p.add_argument("--num-sample", type=int, default=2000,
                   help="time-domain sample count (default 2000)")
    p.add_argument("--model", choices=["dryden", "von_karman", "both"],
                   default="both")
    p.add_argument("--meso-um", type=float, default=5.0,
                   help="mesoscale gust amplitude m/s; 0 disables (default 5)")
    p.add_argument("--meso-dm", type=float, default=100.0,
                   help="mesoscale gust length m (default 100)")
    p.add_argument("--meso-dir", choices=["u", "v", "w"], default="u")
    p.add_argument("--t-start", type=float, default=10.0,
                   help="mesoscale gust entry time s (default 10)")
    p.add_argument("--no-large-scale", action="store_true",
                   help="do not read large-scale mean wind from data")
    p.add_argument("--shear", action="store_true",
                   help="enable vertical wind shear (power-law mean-wind profile)")
    p.add_argument("--shear-extent", type=float, default=10.0,
                   help="vertical extent (m) for wind shear (default 10)")
    p.add_argument("--ref-height", type=float, default=10.0,
                   help="reference altitude z_ref (m) of the wind profile")
    p.add_argument("--ref-speed", type=float, default=None,
                   help="mean wind speed at z_ref (m/s); default = measured mean")
    p.add_argument("--alpha", type=float, default=0.15,
                   help="wind-shear exponent of the power law (default 0.15)")
    p.add_argument("--shear-dir", choices=["u", "v", "w"], default="u",
                   help="axis along which the wind shear acts")
    p.add_argument("--large-scale-series", action="store_true",
                   help="add the measured large-scale fluctuation time series")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--table", action="store_true",
                   help="also print the paper Table 1/2 parameter sweep")
    p.add_argument("--self-check", action="store_true",
                   help="run the automated physics validation and exit")
    p.add_argument("-i", "--interactive", action="store_true",
                   help="prompt interactively for the model parameters at runtime")
    p.add_argument("--out-dir", type=str, default=None)
    return p.parse_args(argv)


# Ordered spec of the parameters that may be entered at runtime.
# (attribute name, prompt label, type cast, allowed choices or None)
_PROMPT_SPECS = [
    ("height", "Flight altitude (m)", float, None),
    ("airspeed", "Airspeed (m/s)", float, None),
    ("time", "Simulation duration (s)", float, None),
    ("frequency", "Reference gust frequency (Hz)", float, None),
    ("strength", "Turbulence strength (1=light 2=moderate 3=severe)", int,
     [1, 2, 3]),
    ("confidence", "Confidence coefficient (0-1)", float, None),
    ("num", "Separation points for prediction", int, None),
    ("num_sample", "Time-domain sample count", int, None),
    ("model", "Small-scale model", str, ["dryden", "von_karman", "both"]),
    ("meso_um", "Mesoscale gust amplitude (m/s, 0=off)", float, None),
    ("meso_dm", "Mesoscale gust length (m)", float, None),
    ("meso_dir", "Mesoscale gust direction", str, ["u", "v", "w"]),
]


def _prompt(label, default, cast=str, choices=None):
    """Prompt for one value, showing the default; validate and cast the reply."""
    hint = f" [{default}]"
    if choices:
        hint += f" ({'/'.join(str(c) for c in choices)})"
    while True:
        raw = input(f"  {label}{hint}: ").strip()
        if raw == "":
            return default
        try:
            value = cast(raw)
        except ValueError:
            print(f"    invalid value '{raw}', please re-enter")
            continue
        if choices and value not in choices:
            print(f"    must be one of {list(choices)}, please re-enter")
            continue
        return value


def _interactive_input(args):
    """
    Prompt for the model parameters at runtime and update ``args`` in place.

    Pressing Enter accepts the current default (which is either the value given
    on the command line or the program default), so no parameter value needs to
    be edited inside the source code.
    """
    print("\n=== Multi-scale wind model — parameter input ===")
    print("    (press Enter to accept the default shown in brackets)\n")
    for attr, label, cast, choices in _PROMPT_SPECS:
        setattr(args, attr, _prompt(label, getattr(args, attr), cast, choices))
    print("")
    return args


def _report(height, strength, confidence, model, meso_um, meso_dm, airspeed):
    """Print a short summary of scales, intensities and predicted amplitudes."""
    Lu, Lv, Lw = dts.dryden_length_scale(height)
    su, sv, sw = dts.dryden_turbulence_speed(height, strength)
    from scipy.stats import norm
    z = norm.ppf((1.0 + confidence) / 2.0)

    print("\n================ SUMMARY ================")
    print(f"Altitude            : {height:g} m  ({wc.classify_altitude(height)})")
    print(f"Airspeed            : {airspeed:g} m/s")
    print(f"Strength            : {wc.STRENGTH_LABELS[strength]}")
    print(f"Confidence          : {confidence:.1%}")
    print(f"Small-scale model   : {model}")
    print("-" * 42)
    print(f"Length scales (m)   : Lu={Lu:.3f}, Lv={Lv:.3f}, Lw={Lw:.3f}")
    print(f"Intensity sigma(m/s): u={su:.4f}, v={sv:.4f}, w={sw:.4f}")
    print(f"Amplitude at {confidence:.0%} conf. (z={z:.3f}): "
          f"u={z*su:.4f}, v={z*sv:.4f}, w={z*sw:.4f} m/s")
    if meso_um > 0:
        print(f"Mesoscale gust      : Um={meso_um:g} m/s, dm={meso_dm:g} m")
    print("=========================================\n")


def _sweep_table(confidence=0.99, strength=3, airspeed=10.0):
    """Reproduce the paper Table 1/2 altitude sweep with the corrected model."""
    from scipy.stats import norm
    z = norm.ppf((1.0 + confidence) / 2.0)
    heights = [10.0, 30.0, 100.0, 300.0, 500.0, 1000.0]

    print("\nTable sweep (V=%.0f m/s, %s turbulence, %.0f%% confidence)"
          % (airspeed, wc.STRENGTH_LABELS[strength], confidence * 100))
    print(" H(m)   Lu(m)     Lv(m)     sigma_u   sigma_w   W_99%%(m/s)")
    print("-" * 62)
    for h in heights:
        Lu, Lv, Lw = dts.dryden_length_scale(h)
        su, sv, sw = dts.dryden_turbulence_speed(h, strength)
        w_99 = z * sw
        print(f"{h:6g} {Lu:9.2f} {Lv:9.2f} {su:9.4f} {sw:9.4f} {w_99:9.3f}")
    print("-" * 62)


def run(args, out_dir=None):
    """
    Run the full multi-scale wind pipeline for the given parsed arguments.

    Returns the output directory path.  This function is shared by the
    command-line driver and the graphical interface (``gui.py``).
    """
    out_dir = out_dir or args.out_dir or wc.ensure_output_dir()
    print(f"Output directory: {out_dir}")

    _report(args.height, args.strength, args.confidence, args.model,
            args.meso_um, args.meso_dm, args.airspeed)

    # --- small-scale turbulence ------------------------------------------- #
    if args.model in ("dryden", "both"):
        dm.dryden_plottings_output(args.frequency, args.time, args.height,
                                   args.strength, args.airspeed,
                                   num_sample=args.num_sample, seed=args.seed,
                                   out_dir=out_dir)
    if args.model in ("von_karman", "both"):
        vkm.von_karman_plottings_output(args.frequency, args.time, args.height,
                                        args.strength, args.airspeed,
                                        num_sample=args.num_sample,
                                        seed=args.seed, out_dir=out_dir)

    # --- amplitude correction / prediction -------------------------------- #
    ndp.wind_model_prection_output(args.confidence, args.height, args.strength,
                                   args.airspeed, args.num, out_dir=out_dir)

    # --- mesoscale 1-cos gust --------------------------------------------- #
    if args.meso_um > 0:
        meso.mesoscale_plottings_output(args.meso_um, args.meso_dm,
                                        args.airspeed, time=args.time,
                                        t_start=args.t_start,
                                        num_sample=args.num_sample,
                                        out_dir=out_dir)

    # --- large-scale wind -------------------------------------------------- #
    large.large_scale_plottings_output(
        time=args.time, num_sample=args.num_sample, out_dir=out_dir,
        height=args.height if args.shear else None,
        ref_height=args.ref_height,
        ref_speed=(args.ref_speed if args.ref_speed is not None else 8.0),
        alpha=args.alpha)

    # --- multi-scale superposition ----------------------------------------- #
    vm = None if (not args.no_large_scale) else (0.0, 0.0, 0.0)
    msm.multi_scale_plottings_output(
        args.time, args.height, args.strength, args.airspeed,
        model=args.model if args.model != "both" else "dryden",
        num_sample=args.num_sample, seed=args.seed,
        meso_um=args.meso_um, meso_dm=args.meso_dm, meso_dir=args.meso_dir,
        t_start=args.t_start, vm=vm, out_dir=out_dir,
        use_shear=args.shear, shear_extent=args.shear_extent,
        ref_height=args.ref_height, ref_speed=args.ref_speed, alpha=args.alpha,
        shear_dir=args.shear_dir, use_large_scale_series=args.large_scale_series)

    if args.table:
        _sweep_table(args.confidence, args.strength, args.airspeed)

    print(f"\nDone. Results written to: {out_dir}")
    return out_dir


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = _parse_args(argv)

    if args.self_check:
        import self_check
        sys.exit(self_check.main())

    # Interactive runtime input: requested explicitly with --interactive, or
    # automatically when the program is launched with no command-line arguments.
    # This lets the user enter the model parameters at runtime instead of
    # editing the values in the source code.
    if args.interactive or not argv:
        _interactive_input(args)

    run(args)


if __name__ == "__main__":
    main()
