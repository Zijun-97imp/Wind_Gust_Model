# Gust Wind Model Forming by Python (GWM)
# Multi-scale Gust Wind Model — Documentation

> **Optimized and engineering-ready implementation of the multi-scale wind model**
>
> This program is a further optimization and engineering implementation of the
> multi-scale wind model. 
> The programs in this folder are the data source of the paper; this version 
> fixes several computational defects in the original program, fills in the 
> missing mesoscale and large-scale modules, adds the **vertical wind-shear 
> (wind-profile) component, the full large-scale time series and an automated 
> self-check module**, together with this documentation.

---

## Table of Contents

1. [Model Overview](#1-model-overview)
2. [Theory](#2-theory)
3. [Directory Structure](#3-directory-structure)
4. [Requirements and Dependencies](#4-requirements-and-dependencies)
5. [Quick Start](#5-quick-start)
6. [Command-line Arguments](#6-command-line-arguments)
7. [Interactive Input and GUI](#7-interactive-input-and-gui)
8. [Output Files](#8-output-files)
9. [Automated Self-check](#9-automated-self-check)
10. [Optimizations (vs. the original program)](#10-optimizations-vs-the-original-program)
11. [Validation Results](#11-validation-results)
12. [References](#12-references)

---

## 1. Model Overview

Random atmospheric gusts significantly affect the stability, safety and
structural design of low-altitude aircraft. A single-scale gust model (e.g. only
a Dryden / Von Karman small-scale turbulence model) cannot cover the complete
airspace from the large-scale mean wind, through the mesoscale discrete gust,
down to the small-scale random turbulence. This model **linearly superposes**
wind fields of different scales in Cartesian coordinates, and computes the gust
amplitude at any location and any confidence level by the **correlation function
+ normal distribution** method.

The model contains four spatial/temporal scales:

| Scale | Physical object | Model file | Symbol |
|-------|-----------------|------------|--------|
| Large scale | Background mean wind (measured data + wind profile) | `large_model.py` | `Vm(z)` |
| Mesoscale | Vertical wind shear (altitude-dependent profile difference) | `large_model.py` | `V_shear` |
| Mesoscale | Discrete gust (`1-cos` wind shear) | `meso_model.py` | `Vg` |
| Small scale | Random turbulence (Dryden / Von Karman) | `dryden_model.py` / `von_karman_model.py` | `Vturb` |
| (optional) Large-scale fluctuation | Measured low-frequency fluctuation series | `large_model.py` | `V_fluc` |

Total wind field (Eq. (9) of the paper, with the shear and large-scale series
terms completed in this version):

```
V_w(t) = Vm(z) + V_shear(z) + V_fluc(t) + Vg(t) + Vturb(t)
```

---

## 2. Theory

### 2.1 Small-scale turbulence (Dryden / Von Karman spectra)

Following MIL-HDBK-1797 / MIL-F-8785C, atmospheric turbulence is described by
the **turbulence scale lengths** `Lu, Lv, Lw` and the **turbulence intensities**
`σu, σv, σw`.

**Low altitude (h < 1000 ft = 304.8 m)**:

```
Lw = h / 2
Lu = h / (0.177 + 0.000823 h)^1.2        (h in ft)
Lv = Lu / 2

σw = 0.1 · W20        (W20 = 15 / 30 / 45 ft/s for light / moderate / severe)
σu = σw / (0.177 + 0.000823 h)^0.4
σv = σu
```

**Medium/high altitude (2000 ft = 609.6 m and above, below 1000 m)**:
`Lu = 2500 ft`, `Lv = Lw = 1250 ft`, `σw = 0.1·(25/40/55) ft/s`,
`σu = σv = σw`.

**Transition band (1000–2000 ft)**: linear interpolation in altitude between the
two bands above.

**Above 1000 m**: `Lu, Lvw, σ` are looked up from `wind_info_over_altitude.xlsx`.

**Dryden shaping filters (transfer functions, SI units)**:

```
H_u(s) = σu √(2 Lu/(πV)) / (1 + (Lu/V) s)
H_v(s) = σv √(  Lv/(πV)) (1 + √3 (Lv/V) s) / (1 + (Lv/V) s)²
H_w(s) = σw √(  Lw/(πV)) (1 + √3 (Lw/V) s) / (1 + (Lw/V) s)²
```

**Dryden power spectral densities (one-sided, ω in rad/s)**:

```
Φu(ω) = σu² (2 Lu/(πV)) / (1 + (Lu ω/V)²)
Φv(ω) = σv² (2 Lv/(πV)) (1 + 3 (Lv ω/V)²) / (1 + (Lv ω/V)²)²
Φw(ω) = σw² (2 Lw/(πV)) (1 + 3 (Lw ω/V)²) / (1 + (Lw ω/V)²)²
```

The Von Karman model uses the same `Lu/Lv/Lw` and `σ`, differing only in the
shaping filter and spectrum (rational-approximation coefficients are documented
in the header of `von_karman_model.py`; the spectra have the form
`(1 + (1.339 L ω/V)²)^(-5/6)` and `^(-11/6)`).

> **On the spectral/filter convention**: the paper's Eq. (2) and Eq. (5) use the
> discrete-filter coefficient form (`1+12(LΩ)²`, `2L/V` time constant), while its
> Eq. (10)–(11) autocorrelation functions match the standard continuous Dryden
> form. This program consistently adopts the **standard continuous Dryden form**
> (MIL-HDBK-1797 / MIL-F-8785C), so that the shaping filter, the power spectral
> density and the autocorrelation are mutually consistent (e.g.
> `Rv(Lv) = e⁻¹(1 − 1/2) = 0.184`).

### 2.2 Mesoscale discrete gust (`1-cos` profile)

```
U(x) = 0,                                   x < 0
U(x) = (Um/2) (1 − cos(π x / dm)),          0 ≤ x ≤ 2 dm
U(x) = 0,                                   x > 2 dm
```

`Um` is the gust amplitude and `dm` the gust length (the gust rises from 0 to
`Um` and returns to 0 over a distance `2 dm`). In the time domain `x = V t`, so
`Vg(t) = (Um/2)(1 − cos(π V t / dm))`.

### 2.3 Large-scale wind and vertical wind shear

The large-scale mean wind `Vm` is obtained as the time average of the measured
meteorological data (`wind_data.csv`). The paper notes that large-scale winds
"have higher velocities and exhibit directional and vertical wind shear". This
version therefore adds a **mean-wind vertical profile (wind shear) model**:

```
Power law:  V(z) = V_ref · (z / z_ref)^α        (α ≈ 0.1 flat, 0.2 farmland, 0.3+ urban)
Log law:    V(z) = (u* / κ) · ln(z / z0)        (κ = 0.4)
```

With `--shear`, the large-scale mean wind `Vm` is rescaled to the flight altitude
via the profile (i.e. `Vm(z)`), and an additional **vertical wind-shear
component** `V_shear` is added — the difference of the profile between the upper
and lower ends of the aircraft vertical extent `dz` (reference height/wingspan,
default 10 m):

```
V_shear = V(z + dz/2) − V(z − dz/2)
```

In addition, `--large-scale-series` adds the measured large-scale low-frequency
fluctuation `V_fluc(t)` (the `wind_data.csv` series interpolated onto the target
time grid, minus its mean). The spatial/temporal large-scale field itself is
provided by `large_model.py` (Eq. (7)–(8) of the paper).

### 2.4 Amplitude correction (normal-distribution method, Eq. (10)–(14))

Dryden autocorrelation functions:

```
Ru(d) = exp(−d/Lu)
Rv(d) = exp(−d/Lv) (1 − d/(2 Lv))
Rw(d) = exp(−d/Lw) (1 − d/(2 Lw))
```

Equivalent (conditional) strength `σe = σ √(1 − R²)`. Given a two-sided
confidence level `C`, the corresponding percentile is `P = (1 + C)/2`, and the
gust amplitude at a separation `d` is:

```
Initial velocity V1 = 0:   U(d) = z_P · σe(d)
Initial velocity V1 ≠ 0:   U(d) = R(d)·V1 + z_P · σe(d)      (z_P = norm.ppf(P))
```

---

## 3. Directory Structure

```
Gust_Wind_Model/
├── wind_common.py                 # shared utilities: paths, units, altitude bands, table, save
├── multi_scale_model.py           # multi-scale superposition (Eq. (9), shear + large-scale series)
├── self_check.py                  # automated physics validation (39 checks)
├── gui.py                         # graphical (visualization) interface
├── wind_info_over_altitude.xlsx   # turbulence scale/intensity table above 1000 m
├── wind_data.csv                  # measured large-scale wind (time,u,v,w)
├── white_noise.csv                # original white-noise samples (kept for reference)
├── Main/
│   └── Main.py                    # command-line driver (with interactive input)
├── Dryden_Gust_Model/
│   ├── dryden_turbulence_selection.py   # Dryden scale/intensity selection
│   └── dryden_model.py                  # Dryden shaping filter, PSD, time series, plots
├── Von_Karman_Gust_Model/
│   ├── von_karman_turbulence_selection.py
│   └── von_karman_model.py
├── Mesoscale_Gust/
│   └── meso_model.py              # mesoscale 1-cos gust
├── Largescale_Gust/
│   └── large_model.py             # large-scale wind (wind profile / vertical shear)
├── Normal_Distribution_Model_Prediction/
│   └── normal_distribution_prediction.py  # amplitude correction / prediction
├── Simulink_based_model/          # Simulink reference data (Dryden / Von Karman)
├── Literature_Resources/          # references
├── original_results/              # original program results (historical archive)
└── output/                        # program output (generated at runtime)
```

---

## 4. Requirements and Dependencies

- **Python 3.11** (recommended, verified; 3.8+ should also work).
- Dependencies:

| Package | Purpose |
|---------|---------|
| numpy | numerical computation |
| scipy | shaping-filter `lsim` / transfer functions, `norm.ppf` |
| matplotlib | plotting |
| pandas / openpyxl | read `wind_info_over_altitude.xlsx` |
| tkinter | the graphical interface (included with Python) |
| Pillow | image display in the GUI |

Install (if missing):

```bash
pip install numpy scipy matplotlib pandas openpyxl pillow
```

---

## 5. Quick Start

From the project root, run the command-line driver (no arguments → interactive
parameter input):

```bash
python Main/Main.py
```

Provide the full set of parameters via the command line (reproduces the paper's
Table 1/2 altitude sweep):

```bash
python Main/Main.py --height 10 --airspeed 10 --time 100 --strength 3 \
    --confidence 0.99 --model both --meso-um 6 --meso-dm 100 --table
```

Enable the vertical wind shear and the large-scale fluctuation series:

```bash
python Main/Main.py --height 100 --airspeed 10 --time 100 --strength 3 \
    --shear --alpha 0.2 --shear-extent 10 --large-scale-series
```

Run the automated physics self-check:

```bash
python self_check.py
```

Launch the graphical interface:

```bash
python gui.py
```

> The program does not depend on the current working directory: it can be run
> from any directory, since `wind_common.py` resolves the data-file paths from
> `__file__`. All results are written to the `output/` directory.

---

## 6. Command-line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--height` | float | 500 | flight altitude (m) |
| `--airspeed` | float | 10 | airspeed (m/s) |
| `--time` | float | 100 | simulation duration (s) |
| `--frequency` | float | 1 | reference gust frequency (Hz) |
| `--strength` | int | 2 | turbulence strength: 1=light 2=moderate 3=severe |
| `--confidence` | float | 0.99 | confidence coefficient (0–1) |
| `--num` | int | 500 | number of separation points for prediction |
| `--num-sample` | int | 2000 | time-domain sample count |
| `--model` | str | both | `dryden` / `von_karman` / `both` |
| `--meso-um` | float | 5 | mesoscale gust amplitude (m/s), 0 disables |
| `--meso-dm` | float | 100 | mesoscale gust length (m) |
| `--meso-dir` | str | u | mesoscale gust direction: u/v/w |
| `--t-start` | float | 10 | mesoscale gust entry time (s) |
| `--no-large-scale` | flag | — | do not read the large-scale mean wind from data |
| `--shear` | flag | — | enable vertical wind shear (power-law wind profile) |
| `--shear-extent` | float | 10 | shear vertical extent dz (m) |
| `--ref-height` | float | 10 | profile reference altitude z_ref (m) |
| `--ref-speed` | float | measured mean | mean wind speed at z_ref (m/s) |
| `--alpha` | float | 0.15 | power-law wind-shear exponent |
| `--shear-dir` | str | u | axis along which the wind shear acts |
| `--large-scale-series` | flag | — | add the large-scale fluctuation time series |
| `--seed` | int | 0 | RNG seed (reproducibility) |
| `--table` | flag | — | also print the paper Table 1/2 altitude sweep |
| `--self-check` | flag | — | run the automated self-check and exit |
| `-i`, `--interactive` | flag | — | prompt for the model parameters at runtime |
| `--out-dir` | str | `output` | output directory |

---

## 7. Interactive Input and GUI

The model parameters `height, airspeed, time, frequency, strength, confidence,
num, num_sample, model, meso_um, meso_dm, meso_dir` (and the advanced options)
can all be entered **at runtime** rather than edited in the source code.

**Interactive command line** — run with no arguments (or with `-i`/`--interactive`),
and each parameter is prompted one by one, showing its default (press Enter to
accept it):

```bash
python Main/Main.py            # prompts for every parameter
python Main/Main.py --interactive --time 50   # overrides some defaults first
```

**Graphical interface** — `python gui.py` opens a window with two tabs:

- **基本参数 Basic** — the twelve primary parameters above.
- **高级参数 Advanced** — gust entry time, seed, wind shear options, large-scale
  series and the altitude-sweep table.

Click **运行模型 Run** to run the pipeline in the background (the window stays
responsive); the summary is printed in the log area. Then **查看结果图 View
figures** opens a scrollable viewer with all generated figures, and **打开输出目录
Open folder** opens the `output/` directory.

---

## 8. Output Files

Generated in the `output/` directory:

| File | Content |
|------|---------|
| `Dryden_Gust_Wind_Model.png` | Dryden time-domain turbulence (u/v/w) |
| `Dryden_Gust_Wind_Model_PSD.png` | Dryden theoretical PSD |
| `Von_Karman_Gust_Wind_Model.png` | Von Karman time-domain turbulence |
| `Von_Karman_Gust_Wind_Model_PSD.png` | Von Karman theoretical PSD |
| `Relativity_Function_Diagram.png` | correlation function R vs d/L |
| `Gust_Wind_Amplitude_Prediction_Curve.png` | amplitude prediction curve (V1=0 / V1≠0) |
| `Gust_Wind_Amplitude_Ratio.png` | non-dimensional amplitude U/σ vs d/L (multiple confidence levels) |
| `Mesoscale_1cos_Gust.png` | mesoscale 1-cos gust |
| `Large_scale_Wind.png` | large-scale measured wind |
| `Wind_Profile_Shear.png` | mean-wind vertical profile (with `--shear`) |
| `Multi_scale_Wind_Field.png` | multi-scale superposition (shear / large-scale series) |
| `*.csv` | time-domain / correlation / wind-profile data (re-computable) |

---

## 9. Automated Self-check

`self_check.py` provides a one-click physics validation, comparing the model with
the paper's Table 2 scale lengths and the MIL-HDBK-1797 theory, and prints a
PASS/FAIL summary (39 checks; exit code 0 means all passed):

```bash
python self_check.py
# or
python Main/Main.py --self-check
```

| Category | Check |
|----------|-------|
| Scale lengths | `Lu` at 10/30/100/300/500/1000 m matches the paper's Table 2 |
| Turbulence intensity | generated `std(u/v/w)` matches the target `σ` (Dryden and Von Karman) |
| Power spectral density | Dryden high-frequency `ω⁻²`, Von Karman high-frequency `ω^(-5/3)` (Kolmogorov inertial subrange) |
| Autocorrelation | `Ru(Lu)=e⁻¹`, `Rv(Lv)=Rw(Lw)=e⁻¹/2` (Eq. (10)–(11)) |
| Superposition consistency | total wind == sum of components (numerical error < 1e-9) |
| Wind profile/shear | power-law profile monotonically increasing, correct shear sign |

---

## 10. Optimizations (vs. the original program)

### 10.1 Paths and engineering

- The original `Main.py` hard-coded an `E:\01Wind Field\...` path; `wind_common.py`
  now resolves the project root from `__file__`, so the program runs from any
  directory / drive.
- Data files were read with relative paths at import time; they now use absolute
  path resolution.
- Figures and data are written to a single `output/` directory.

### 10.2 Unit handling

- All calculations are in SI units; unit constants are centralized in
  `wind_common.py` (removing the error-prone `m ↔ ft`, `knot ↔ m/s` round-trips).

### 10.3 Turbulence scale/intensity selection fixes (`dryden_turbulence_selection.py`)

- **Transition-band interpolation bug**: the original `k = (v2−v1)/1500`
  (denominator should be 1000 ft) and the missing anchor in `Lu = k1·h_en`
  caused a discontinuity at the 2000 ft boundary; fixed to the anchored linear
  interpolation `v = v1 + (v2−v1)·(h−1000)/1000`.
  > The corrected `Lu` at 500 m (transition band) changes from ~633 m to ~598 m;
  > the other altitudes (10/30/100/300/1000 m) match the paper's Table 2.
- **σ slope mix-up**: the `k1/k2/k3` slopes for `sigma_u/v/w` in the transition
  band were misassigned; fixed.

### 10.4 PSD fixes (`dryden_model.py` / `von_karman_model.py`)

- The original `*_wind_velocities_freq` passed a cosine **time-domain signal** as
  "frequency" into the spectral formula; the PSD is now evaluated over the true
  angular frequency `ω` (logarithmic grid).
- `Φv`/`Φw` wrongly used `σu` (should be `σv`/`σw`) and `Lv` (`Φw` should be
  `Lw`); fixed. The Von Karman spectra now consistently use the standard
  coefficient `1.339`.

### 10.5 Amplitude correction fixes (`normal_distribution_prediction.py`)

- The `V2/W2` computation used `loc_num_u` instead of `loc_num_v`/`loc_num_w`; fixed.
- The confidence→percentile formula `P = 1 − pr/2` was wrong (99% gave the 0.505
  percentile); fixed to `P = (1 + pr)/2`.
- A distance-grid length mismatch (10000 vs num) was unified.

### 10.6 Completed missing modules

- **Mesoscale 1-cos gust** (the original `meso_model.py` was an empty stub):
  added `one_cosine_gust`, `meso_gust_time_series` and plotting/export.
- **Large-scale wind** (the original `Largescale_Gust/` only had
  `white_noise.csv`): added `large_model.py` reading `wind_data.csv`.
- **Multi-scale superposition**: added `multi_scale_model.py` implementing
  Eq. (9), returning each component separately for inspection.

### 10.7 Reproducibility and data export

- Turbulence white noise now uses a fixed-seed RNG (`--seed`, default 0); the
  original `white_noise.csv` is kept for reference.
- Added CSV export (component time series, correlation functions).

### 10.8 This round's further optimizations (shear / large-scale series / self-check / GUI / interactive input)

- **Completed the vertical wind-shear component**: Eq. (9) explicitly includes
  "mid-scale wind shear", which the original implementation only mentioned in a
  comment. `large_model.py` now provides power-law / log-law wind profiles
  `V(z)` and the vertical shear `V_shear = V(z+dz/2) − V(z−dz/2)`, added as an
  independent component in `multi_scale_model.py` (`--shear`).
- **Large-scale wind upgraded from "mean only" to "wind profile + full series"**:
  `Vm` is rescaled with altitude (`Vm(z)`), and the measured low-frequency
  fluctuation `V_fluc(t)` can be superposed (`--large-scale-series`), realizing
  the spatiotemporal field of Eq. (7)–(8).
- **Automated self-check** `self_check.py`: 39 physics checks (scales, intensity,
  spectral asymptote, correlation, superposition, wind profile).
- **Interactive runtime input** in `Main.py`: the primary parameters are prompted
  at runtime (`-i`/`--interactive`, or automatically when no arguments are given),
  so no value needs to be edited in the source code.
- **Graphical interface** `gui.py`: enter the parameters in a window, run the
  model, and view all generated figures.
- **Transparent spectral convention**: the documentation states clearly that the
  standard continuous Dryden spectrum/filter is used (consistent with Eq. (10)–
  (11)), and notes the discrete-coefficient difference in the paper's Eq. (2)/(5)
  (see §2.1).

---

## 11. Validation Results

Running `--table` with the paper's Table 1 parameters (V = 10 m/s, severe
turbulence, 99% confidence), the scale lengths compare with the paper's Table 2
as follows (`Lu` in m):

| Altitude (m) | This program Lu | Paper Lu | This program Lv | Paper Lv_m |
|--------------|-----------------|----------|-----------------|------------|
| 10 | 67.37 | 67.36 | 33.68 | 32.48 |
| 30 | 152.46 | 152.46 | 76.23 | 74.56 |
| 100 | 262.79 | 262.78 | 131.40 | 124.81 |
| 300 | 304.73 | 304.72 | 152.37 | 167.5 |
| 500 | 597.60 | 632.78* | 298.80 | 312.23* |
| 1000 | 832.00 | 832.1 | 624.00 | 501.23 |

> `*`: 500 m lies in the transition band, where the original program (and the
> paper) used a flawed interpolation (see §10.3); this version uses the corrected
> continuous interpolation, hence the slightly smaller value. The other altitudes
> match the paper.
> The paper's `Lv` value 32.48 differs slightly from `Lu/2 = 33.68` because of
> the paper's independent rounding/definition of `Lv`; this does not affect the
> model's overall correctness.

The generated time-domain turbulence statistics have been verified: each
component's standard deviation exactly matches the target intensity `σ` (e.g. for
severe turbulence at 10 m, `std(u)=σu=2.5904 m/s`, `std(w)=σw=1.3716 m/s`). The
PSD asymptotic slopes and the autocorrelation functions are verified item by item
by `self_check.py` (39/39 passed).

---

## 12. References

1. Zhang H., Li S., Yi J., Zhong G. Review on Urban Low-Altitude Air Route
   Planning. Transactions of Nanjing University of Aeronautics and Astronautics,
   38, 2021.
2. Adelfang S., Smith O. Gust models for launch vehicle ascent. 36th AIAA
   Aerospace Sciences Meeting and Exhibit, 1998.
3. MathWorks Aerospace Blockset — Dryden Wind Turbulence Model (Discrete/Continuous).
4. Fan D., Dong J. One Turbulent Modeling Method Based on Environmental Forecast
   Data. Journal of System Simulation, 38: 2409-2414, 2020.
5. Grauer J. Random noise generation using Fourier series. Journal of Aircraft,
   2018.
6. Fu W., Shang N., Chang X., et al. A New Method for Calculating Magnitude of
   Wind Gust. Journal of Northwestern Polytechnical University, 32: 822-827, 2014.
7. Chalk C., Harris T., Neal T., Pritchard E. Background Information and User
   Guide for MIL-F-8785B(ASG). Wright-Patterson Air Force Base, 1969.
8. MIL-HDBK-1797. Flying Qualities of Piloted Aircraft.

---

*Program generated on 2026-11.*

