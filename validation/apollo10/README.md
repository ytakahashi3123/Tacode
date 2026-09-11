# Apollo 10 (AS-505) entry

Tacode against the Apollo 10 entry of 1969-05-26.

**This is the case where nothing about the lift is assumed.** NASA TN D-6725 reports the
state vector at the entry interface (table II), the roll angle the guidance actually held
through the entry (figure 12) and the altitude and velocity the vehicle actually flew
(figure 13). The roll angle is the direction of the lift, so it goes in as an **input**,
and the altitude and velocity are the comparison. Nothing is left to tune except the
entry mass, which is not in the reference.

The companion case `validation/apollo4` has better flight data for the atmosphere and the
heating, but no bank history, so there a constant bank had to be assumed. Read the two
together.

## Running it

```bash
cd validation/apollo10
./run_tacode.sh                  # 2.9 s
python3 compare_apollo10.py      # the numbers below, and the figures
python3 entry_state.py           # the entry state of table II converted for config.yml
python3 digitize_figure.py       # rebuild reference/ from the scanned report
```

`digitize_figure.py` needs the report itself (in `../../../references`, outside the
repository) plus `pdftoppm` and Pillow; the digitised result is committed in `reference/`,
so the case runs without any of that.

## Where every number comes from

| Quantity | Value | Source |
|---|---|---|
| Entry state | 11067.06 m/s inertial, -6.616 deg. flight-path angle, 71.928 deg. azimuth, 174.244 E, 23.652 S, 123.883 km | [1] table II (best-estimated trajectory) |
| Initial velocity for Tacode | [10016.49, 3413.34, -1266.84] m/s | the same, converted by `entry_state.py` (see below) |
| **Bank angle** | `bank_angle_table` in `config.yml`, 89 nodes from 0 to 440 s | **[1] figure 12, the "Roll angle" trace**, digitised |
| Flight altitude and velocity | `reference/apollo10_flight.dat` | [1] figure 13, the "Actual" curves, digitised |
| Drag and lift coefficients | `CD = 1.2891`, `CL = 0.38773` | [1] table I, trim values at M > 29.5 (`L/D = 0.30076`). Figure 14 of the same report shows the L/D measured in flight within the predicted 0.03 of the preflight value |
| Reference area, length | `12.0171 m^2`, `3.9116 m` | 154 in. maximum diameter [3] |
| Entry mass | `5500 kg` | **not in the references**; the figure commonly quoted for a lunar-return command module. The sensitivity is measured below |
| Atmosphere | `database/atmosphere/atmospheremodel_apollo10.txt` | NRLMSISE-00 for 1969-05-26 16:37:52 UTC at the entry point, written by `database/atmosphere/generate_atmosphere_table.py` (F10.7 = F10.7A = 150, Ap = 4) |
| Epoch | 1969-05-26T16:37:52Z | lift-off 1969-05-18 16:49:00 UTC plus the 191:48:52.16 of table II |

### The entry state has to be converted

Table II gives the **inertial** speed, flight-path angle and azimuth at a geodetic
position. Tacode is given the **ECEF** velocity in the geocentric local horizon, so
`entry_state.py` does three things: builds the inertial velocity in the local horizon,
subtracts the rotation of the Earth (434.54 m/s eastward at that point), and rotates the
result into the geocentric frame. The outcome:

| | |
|---|---|
| Relative (ECEF) speed | 10657.67 m/s |
| Relative flight-path angle | -6.827 deg. (geocentric horizon) |
| Relative azimuth | 71.182 deg. |

The report does not say whether its flight-path angle is measured from the geodetic or the
geocentric horizon; the two differ by up to 0.19 deg. in the direction of "up", which is
30 m/s of vertical speed here. `entry_state.py --horizon geocentric` gives the other
reading, and the effect on the result is in the sensitivity table.

### The bank angle is the flight trace

Figure 12 plots the roll command (solid) and the roll angle the vehicle held (dashed);
they lie on top of each other except during the first reversal, where the vehicle takes
about 10 s to swing round while the command steps. The digitiser follows the continuous
branch, which is the command during those 10 s. **That is the same attitude either way**:
the plot wraps at +-180 deg., the reversal passes through +-180, and both the vertical and
the lateral component of the lift are the same for +180 and -180.

The guidance flew: lift up (about +2 deg.) for the first 85 s, a roll to lift-down
(+-180 deg.) between 95 and 135 s to cap the first pull-up, back to 0 and then to about
+50 deg. until 210 s, a reversal to -75 deg. until 320 s, and two more to +100 deg. at
360 s and -80 deg. at 400 s while the range was trimmed out.

## Results

Measured on 2026-09-11 with `~/venvs/myenv/bin/python` (Python 3.12.3, numpy 2.2.6,
scipy 1.15.3). The run takes 2.9 s. The digitised flight data spans 12 to 414 s after the
entry interface, 311 points.

| | Altitude | Inertial velocity |
|---|---|---|
| **Difference from the flight, all 311 points** | max **+5.5 km**, rms **3.2 km** | max +360 m/s, rms **189 m/s** |

The computed trajectory follows the flight through the first dip to 54 km, the long
shallow phase with its two shoulders, and the final descent; see
`output_comparison/profile_time-altitude.png`. The velocity is the **inertial** speed, as
plotted in the report: `compare_apollo10.py` adds the rotation of the Earth back to the
ECEF speed Tacode writes.

**The 3.2 km is the whole error budget of a 3-DOF entry simulation** — the atmosphere
table, the constant drag and lift coefficients, the assumed mass, the digitisation of both
the input and the reference, and the missing degrees of freedom, all together.

### Sensitivity

Altitude rms over the same 311 points.

| Variant | Altitude rms | Altitude max | Velocity rms |
|---|---|---|---|
| **Baseline** | **3.19 km** | +5.47 km | 189 m/s |
| Mass -10 % (4950 kg) | 2.33 km | +5.47 km | 149 m/s |
| Mass +10 % (6050 kg) | 4.33 km | -7.46 km | 292 m/s |
| `CL` -0.03 (0.35773) | 8.71 km | -15.20 km | 910 m/s |
| `CL` +0.03 (0.41773) | 5.56 km | +10.99 km | 968 m/s |
| Entry angles read from the geocentric horizon | 3.49 km | -6.37 km | 223 m/s |
| dt = 0.2 s | 3.19 km | +5.47 km | 189 m/s |
| **No lift at all** (`CL = 0`) | **32.94 km** | -48.50 km | 3632 m/s |

Three things to read out of this table:

- **The lift is what the case is about.** Turning it off costs 33 km rms; the trajectory
  then looks like the drag-only Apollo 4 run.
- **The lift coefficient is the limiting input**, not the solver: a change of 0.03 in `CL`
  — the accuracy the report itself claims for the preflight L/D — moves the altitude by
  two to three times the baseline error.
- The assumed mass is worth about 1 km per 10 %, and the time step nothing at all. A
  smaller mass fits slightly better, which is not a reason to change it: fitting the mass
  would turn the comparison into a calibration.

## What this case does not show

- **Nothing about the atmosphere or the heating.** Figure 13 carries no density, and the
  report no heating. That is what `validation/apollo4` is for.
- **Nothing about the attitude.** The bank angle is prescribed; the vehicle is a point
  mass. A 6-DOF run would have to reproduce the roll manoeuvres themselves, which needs
  the reaction-control system and the guidance law, neither of which Tacode has.
- **The reference is a digitised figure**, not a table of numbers. The curves are read
  from a 600 dpi scan by following the dark pixels, and the axis calibration comes from
  the tick marks, so the reading error is of the order of the line width — about 1 km in
  altitude and 100 m/s in velocity. `digitize_figure.py --check figure.png` draws the
  digitised curves back over the scan.

## References

1. C. A. Graves, J. C. Harpold, "Apollo Experience Report - Mission Planning for Apollo
   Entry", NASA TN D-6725, March 1972. (Table I, table II, figures 12, 13 and 14.)
2. R. C. Ried, Jr., W. C. Rochelle, J. D. Milhoan, "Radiative Heating to the Apollo Command
   Module: Engineering Prediction and Flight Measurement", NASA TM X-58091, April 1972.
3. W. C. Moseley, Jr., R. H. Moore, Jr., J. E. Hughes, "Stability Characteristics of the
   Apollo Command Module", NASA TN D-3890, March 1967.

Reference [1] is `Apollo_FlightPathAngle/19720013191.pdf` in
`../../../references/20260423_Apollo_radiation`.
