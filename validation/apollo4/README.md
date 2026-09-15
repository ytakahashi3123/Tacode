# Apollo 4 (AS-501) entry

Tacode against the Apollo 4 entry of 1967-11-09, the flight that flew a lunar-return
entry with an unmanned command module.

The case is split into three comparisons on purpose, because they fail in different ways:

| | What is compared | What it exercises |
|---|---|---|
| **(1) Trajectory** | altitude and relative velocity against the flight values | the whole solver: gravity, the rotating frame, the drag and the lift |
| **(2) Atmosphere** | the density of the table at the **flight** altitudes | the atmosphere table and its interpolator only — no trajectory integration, so the trajectory cannot hide in it |
| **(3) Heating** | stagnation-point correlations evaluated with the flight density and with the table density | whether the density and velocity Tacode produces can carry an aeroheating estimate |

Apollo 4 flew a **lifting** entry: it dipped to 55.6 km, pulled back up to 73.5 km and only
then came down. Drag alone cannot produce that, and the direction of the lift — the bank
angle — is what shapes it. **That bank angle is flight data here**, not an assumption:
NASA TN D-5399 plots the bank angle Apollo 4 actually held through the entry, and this case
flies it.

| Configuration | |
|---|---|
| `config.yml` | drag only (`lift_coefficient` absent, which is the default) |
| `config_lift.yml` | three degrees of freedom with the lift, at the **measured bank angle** |
| `config_6dof.yml` | six degrees of freedom: the capsule is flown at an angle of attack, trimmed by its **measured centre-of-gravity offset**, and the lift comes out of a coefficient table |

## Running it

```bash
cd validation/apollo4
./run_tacode.sh                            # drag only,   1.0 s
./run_tacode.sh -file config_lift.yml      # measured bank, 3.1 s
./run_tacode.sh -file config_6dof.yml      # 6-DOF,       9.8 s
python3 compare_apollo4.py                 # the numbers below, and the figures
python3 compare_apollo4.py --no-figure
python3 compare_apollo4.py --nose-radius 3.048
```

The runs write `output_result/`, `output_result_lift/` and `output_result_6dof/`, the
comparison writes `output_comparison/`; all are ignored by git (the numbers are recorded
here instead).

Two scripts rebuild the reference data from the reports, which live outside the repository
in `../../../references/20260423_Apollo_radiation`:

```bash
python3 make_reference.py      # the trajectory and the heating, from TM X-58091 table I
python3 digitize_bank.py       # the bank angle, from TN D-5399 figure 6(b)
python3 digitize_bank.py --check bank.png    # the reading drawn over the scan
```

`attitude_for_bank.py` turns a trim angle of attack and a bank angle into the three Euler
angles `config_6dof.yml` wants. **Rolling about the body x axis is not a bank**: the body
axis is 25 deg. off the velocity, so a roll about it tilts the crossflow plane out of the
plane the centre of gravity is offset in, a rolling moment appears, and the vehicle turns
out of the bank it was given. The script rotates the whole vehicle about the *velocity*
instead.

## Where every number comes from

| Quantity | Value | Source |
|---|---|---|
| Flight trajectory and heating | `reference/apollo4_flight.dat` | [1] table I, "Apollo 4 trajectory (45-day best estimated trajectory, 12/20/67)" |
| Radiometer (radiative) heating | last column of the same file | [1] figure 22, radiometer CA3363K, digitised |
| **Vertical lift-to-drag ratio** | `reference/apollo4_liftvertical.dat`, turned into `bank_angle_table` in `config_lift.yml` | **[2] figure 16(a), "Vertical lift-to-drag ratio"**, digitised |
| Bank angle (for the record) | `reference/apollo4_bank.dat` | [2] figure 6(b), "Bank-angle time history", digitised |
| Entry mass | `5424.6 kg` | [2] table II: 11 959 lb at the entry interface (11 659 lb at drogue deployment, a 2.5 % loss this case does not model) |
| Centre-of-gravity offset | `0.16688 m` | [2] table II: 6.57 in. off the axis at entry |
| Drag coefficient | `CD = 1.2891` | [3] table I, trim value at M > 29.5. [2] finds the wind-tunnel `CD` within 5 % of the value derived from the Apollo 4 flight |
| Lift coefficient | `CL = 0.47439` | `0.368 x CD`, where **0.368 is the lift-to-drag ratio [2] derives from the Apollo 4 flight itself** (the average between trim and the first peak g). The trim table of [3] gives 0.30076, but that is the Apollo 11 command module, which trimmed at a smaller angle of attack |
| Reference area, length | `12.0171 m^2`, `3.9116 m` | pi/4 x 3.9116^2; maximum diameter 154 in. [4] |
| Entry state | 123.4794 km, 10733.84 m/s, -7.129 deg. | [1] table I, first row (29 968 s from lift-off); the flight-path angle is a quadratic through the first three rows, differentiated at the first point |
| Entry point, azimuth | 153.932 E, 18.467 N, 62.65 deg. | **an assumption**; see below |
| Atmosphere | `database/atmosphere/atmospheremodel_apollo4.txt` | NRLMSISE-00 for 1967-11-09 20:19:29 UTC at the entry point, written by `database/atmosphere/generate_atmosphere_table.py` (F10.7 = F10.7A = 150, Ap = 4) |
| Epoch | 1967-11-09T20:19:29Z | lift-off 12:00:01 UTC plus the 29 968 s of table I; [2] table I puts the entry interface at 29 968.54 s |

### The direction of the lift, and how the reading was checked

[2] carries the same information twice, and the case uses the second form:

- **figure 6(b)**, the bank angle of the entry-control programs, and
- **figure 16(a)**, the vertical lift-to-drag ratio derived from the flight.

`digitize_bank.py` reads both. The data are open circles over a dense grid, so the symbols
are separated from the rules by the **local density** of dark pixels, and the axis
calibration is taken from the positions of the tick labels measured in the image rather
than by eye (reading the labels by eye put the time axis 20 s out).

The bank angle of figure 6(b) reproduces the entry as [2] describes it in words, event by
event:

| What the report says | What the reading gives |
|---|---|
| entered "in the planned full-positive lift attitude", bank 0 = lift up | -14 to 0 deg. until GET 30 050 |
| "a period of full-negative lift had to be flown shortly after the first peak g (t = 30 045)" | rolls to -169 deg. at GET 30 055-30 070 |
| "the UPCONTROL phase was entered at the same time that the lift vector was rolled back to lift vector up" (UPCONTROL at GET 30 085, table I) | back through 0 at GET 30 085-30 096 |
| "in the UPCONTROL phase, the bank angle was controlled between 40 and 90 deg." | +27 to +87 deg. over GET 30 100-30 200 |

**The case is nevertheless driven from figure 16(a)**, through
`bank = arccos((L/D)v / 0.368)`. The two figures agree to 0.03 in `(L/D)v` over the entry,
but they differ where it costs most: the roll reversal at the first peak g, where the
dynamic pressure is at its maximum and a few seconds either way changes the vertical lift
by a factor of two. Reading the vertical lift directly avoids having to resolve a fast
roll from a bank-angle plot.

The reading of figure 16(a) has one trap worth recording: the legend box of the figure sits
at `(L/D)v = 0.43`, and the symbol detector takes it for data unless values beyond the
**physical** limit are rejected. The report's own maximum flight-derived `L/D` is 0.41, so
that is the cut. With the legend read as data the reversal appears 10 s late, and the
computed skip apex comes out 60 km high — the whole case hinged on it.

The sign of the bank (left or right) is lost in the `arccos`; it does not enter the
vertical motion, only the ground track.

**The ground track is still an assumption.** Table I of [1] carries no latitude or
longitude, and figure 7(a) of [2] has them but only as curves; they were not digitised. The
entry point here is the reported splashdown point (30 deg. 06' N, 172 deg. 32' W) walked
back 3620 km along the great circle of 32.6 deg. inclination that passes through it. Figure
7(a) of [2] puts the entry at roughly 155 deg. E and 22 deg. N, about 3 deg. north of that;
turning the azimuth by +-20 deg. moves the altitude at peak heating by less than 0.3 km
(see the sensitivity table), so this does not enter the result.

The second row of table I is printed as `28 970` in the scan. The altitude drops 8679 ft
from the first row, which at 35 220 ft/s is 2.0 s at -6.9 deg., so the row belongs at
**29 970 s** and is read that way. It is the row that fixes the entry flight-path angle.

## Results

Measured on 2026-09-11 with `~/venvs/myenv/bin/python` (Python 3.12.3, numpy 2.2.6,
scipy 1.15.3). The flight table spans 552 s, from 123.5 km down to 37.3 km.

### (1) Trajectory

| Run | Points compared | Altitude difference | Relative-velocity difference |
|---|---|---|---|
| Drag only, whole run | 34 of 45 | max -70.5 km, rms 27.1 km | max -6504 m/s, rms 3445 m/s |
| Drag only, to the peak heating (72 s) | | max -3.2 km, rms 1.2 km | max -194 m/s, rms 65 m/s |
| **Measured lift, whole run** | **45 of 45** | max -14.6 km, rms **5.4 km** | max -2085 m/s, rms 788 m/s |
| Measured lift, to the peak heating | | max **+0.5 km**, rms **0.4 km** | max +46 m/s, rms 21 m/s |
| 6-DOF (constant bank), whole run | 45 of 45 | max -11.8 km, rms 4.7 km | max -1105 m/s, rms 474 m/s |
| 6-DOF, to the peak heating | | max -1.1 km, rms 0.4 km | max +57 m/s, rms 26 m/s |

The drag-only run descends monotonically and reaches the ground at 30 303 s, 725 s before
the actual splashdown; **its 27 km rms is the size of the missing lift, not an accuracy**.

With the measured vertical lift and nothing tuned, the run holds the flight to half a
kilometre through the dive, the dip (56.0 km at 80 s against 55.6 km at 80 s) and the
pull-up, and to 1.8 km out to 200 s. It then falls behind steadily: -3.8 km at 250 s,
-6.0 km at 300 s, -14.6 km at the end of the table. That drift has the sign of too little
lift or too much drag late in the entry, which is what holding `CD` and `L/D` constant
does: [2] reports the flight `L/D` rising from 0.36 to 0.41 by M = 6 (reached at 537 s),
and the wind-tunnel `CD` running 10 to 15 % above the flight value below M = 3.
See `output_comparison/profile_time-altitude.png`.

#### Why the vertical lift is the right thing to prescribe

Prescribing `CL` and a bank angle separately makes this trajectory very stiff: Apollo 4
skipped close to the skip-out boundary of the Apollo corridor, and with the bank angle of
figure 6(b) the apex of the skip runs from 68 km at `L/D` = 0.295 to 138 km at 0.368 —
7 km for every 5 % of `L/D`. Prescribing the **vertical** lift removes that: the quantity
which shapes the trajectory is then the measured one, and the remaining freedom (`CD`,
and the reference `L/D` used to split the lift into a magnitude and a bank angle) is
second order. Sweeping `CD` over the range the flight data itself implies:

| `CD` | Altitude rms | Altitude max | Skip apex |
|---|---|---|---|
| 1.19 | 8.0 km | -20.4 km | 66.7 km |
| 1.22 | 7.3 km | -18.8 km | 67.5 km |
| 1.25 | 6.5 km | -17.1 km | 68.3 km |
| **1.2891** (shipped) | **5.4 km** | -14.6 km | 69.4 km |

The flight apex is 73.5 km at 282 s. An inverse calculation on the flight table itself —
differentiating the tabulated altitude and velocity and solving for the forces — gives an
effective `CD` of 1.19 to 1.25 against the 1.2891 used here, and a vertical `L/D` of 0.39
during the pull-up, which is what figure 16(a) shows. **The flight data are consistent
with themselves and with the solver.**

#### Six degrees of freedom

Here nothing is told about the lift or the trim. The capsule carries a coefficient table of
its own shape (modified Newtonian,
`database/aerodynamic/generate_aerodynamic_table.py --shape capsule`), its centre of gravity
is offset from the axis by the **measured** 6.57 in. = 0.16688 m, and the aerodynamic moment
about that point is what holds it at an angle of attack.

| | Computed | Flight |
|---|---|---|
| Trim angle of attack (continuum) | 25.65 deg. | 24.4 deg. ([1]) |
| Angle of attack held, 50 to 120 s | mean 25.4 deg. (19.4 to 32.6) | |
| `CD` at trim | 1.234 | 1.19 to 1.25 derived from the flight table; 1.2891 wind tunnel [3] |
| `L/D` at trim | 0.394 | 0.368 flight-derived [2] |

**The trim angle is not tuned and comes out 1.3 deg. from the reported one**, and the
Newtonian coefficients land within a few per cent of the flight-derived ones.

The bank angle is a constant (55 deg.) in this run: Tacode has no roll control, so a
six-degree-of-freedom vehicle cannot fly the roll history `config_lift.yml` uses. Its 4.7 km
rms is therefore a fitted number, and the meaningful outputs of this configuration are the
trim angle and the coefficients above.

The time step is set by the trim oscillation, not by the trajectory: its period falls to
0.97 s at the peak dynamic pressure. The trajectory is identical to three digits at
dt = 0.05, 0.02 and 0.01 s; 0.02 s is shipped because it is the largest of those which
passes the `T/20` check of `solver.check_timestep_attitude` without a warning.

### (2) Atmosphere: the table at the flight altitudes

| | table / flight |
|---|---|
| Density, 42 points from 85.7 km down to 37.3 km | mean **1.000**, min 0.911, max 1.114, rms of (ratio - 1) **0.043** |

The flight density in table I is derived from the measured stagnation pressure; the three
rows above 100 km carry no pressure reading and are left out. A 4.3 % rms agreement over
two and a half decades of density is as close as this comparison can get: it says the
NRLMSISE-00 table, the reader and the cubic interpolator are all sound, and that using
F10.7 = 150 for a day whose indices were not available costs nothing at these altitudes.
**This part does not pass through the trajectory at all**, so it holds whatever the lift
does.

### (3) Stagnation-point convective heating

The flight column is the **cold-wall heating at S/R = 0.732** — the radiometer location in the
stagnation region at 24.4 deg. angle of attack — not the stagnation-point value, so a ratio
of exactly 1 is not the target.

| Correlation | Nose radius 4.694 m (heat-shield radius) | Nose radius 3.048 m (the 10-ft equivalent sphere of [1]) |
|---|---|---|
| Detra-Kemp-Riddell, flight density | mean 0.840, rms 0.162 | mean **1.042**, rms **0.055** |
| Detra-Kemp-Riddell, table density | mean 0.840, rms 0.164 | mean 1.042, rms 0.060 |
| Sutton-Graves, flight density | mean 0.743, rms 0.262 | mean 0.922, rms 0.100 |

At the peak (30 040 s, 56.7 km, 9798 m/s) the flight value is 2.430e6 W/m^2, against
2.010e6 W/m^2 (DKR, 4.694 m) and 2.495e6 W/m^2 (DKR, 3.048 m).

Reference [1] states that the normal-shock conditions on the command module at 24.4 deg.
angle of attack correspond to a **10-foot-radius sphere**, and the heating goes as
1/sqrt(R_n); with that radius the DKR correlation lands within 6 % of the flight value over
the whole entry. **The line that matters for Tacode is the second one:** replacing the flight
density with the table density moves the heating by less than a percent, so the density the
solver would hand to an aeroheating model is good enough for it.

The radiative heating (radiometer CA3363K, 30 000 to 30 090 s, peaking at 2.15e6 W/m^2 at
30 032 s, of the same order as the convective peak) is
plotted in `output_comparison/profile_time-radiation.png` but is not compared against a
correlation: the Tauber-Sutton f(V) table is not among the references collected here.

### Sensitivity

Altitude and velocity at 72 s after entry (the peak-heating point), and the time at which the
computed trajectory reaches the ground. Drag only, so that the bank angle does not mix in.

| Variant | Altitude at 72 s | Velocity at 72 s | Ground |
|---|---|---|---|
| **Baseline** (dt = 0.05 s, 5424.6 kg, azimuth 62.65 deg.) | 53.485 km | 9604.3 m/s | 335 s |
| dt = 0.2 s | 53.485 km | 9604.3 m/s | 336 s |
| dt = 0.01 s | 53.485 km | 9604.3 m/s | 335 s |
| Mass -10 % | 53.541 km | 9492.2 m/s | 345 s |
| Mass +10 % | 53.439 km | 9698.6 m/s | 328 s |
| Azimuth 42.65 deg. | 53.201 km | 9579.9 m/s | 333 s |
| Azimuth 82.65 deg. | 53.328 km | 9579.6 m/s | 335 s |
| `atmospheremodel.txt` (the shipped 2015, 55 N table) | 53.332 km | 9903.9 m/s | 334 s |

The time step does not matter at all at this resolution, which is why 0.05 s is shipped.
The mass is now the measured one, and the azimuth moves the velocity by a quarter of a
percent. Swapping in the atmosphere table of the tutorials costs 300 m/s (3 %) at peak
heating — that is the value of generating the table for the date and place of the flight.
None of these approach the lift coefficient, which is what limits the case.

## What this case does not show

- **The lift coefficient is the limit, and it is not settled by the references.** The
  wind-tunnel trim value and the value derived from this very flight differ by 22 %, and
  the skip apex moves 68 km between them. A case built on a less marginal entry is a
  sharper test of the solver: `validation/apollo10` flies a measured bank history too and
  lands at 3.2 km rms over the whole entry, because that trajectory is not sitting on the
  skip-out boundary.
- **No roll control.** The six-degree-of-freedom run cannot follow the measured roll
  history: an unguided capsule keeps the bank it starts with. Reproducing the manoeuvres
  themselves would need the reaction-control system and the guidance law, neither of which
  Tacode has.
- **The moments of inertia and the damping derivatives are assumed.** They are not in the
  references; they set the period and the decay of the trim oscillation, not the trim angle,
  and the trajectory is insensitive to the oscillation at this amplitude.
- **The ground track is assumed**, as described above.
- **Nothing about the heating models of Tacode**, which has none. The correlations live in
  `compare_apollo4.py`.

## References

1. R. C. Ried, Jr., W. C. Rochelle, J. D. Milhoan, "Radiative Heating to the Apollo Command
   Module: Engineering Prediction and Flight Measurement", NASA TM X-58091, April 1972.
   (Table I: the trajectory and the heating. Figure 22: the radiometer.)
2. E. R. Hillje, "Entry Aerodynamics at Lunar Return Conditions Obtained From the Flight of
   Apollo 4 (AS-501)", NASA TN D-5399, 1969. (Figure 6(b): the bank angle. Table I: the
   event times. Table II: the mass and the centre of gravity. The text: the flight-derived
   lift-to-drag ratio.)
3. C. A. Graves, J. C. Harpold, "Apollo Experience Report - Mission Planning for Apollo
   Entry", NASA TN D-6725, March 1972. (Table I: the trim coefficients.)
4. W. C. Moseley, Jr., R. H. Moore, Jr., J. E. Hughes, "Stability Characteristics of the
   Apollo Command Module", NASA TN D-3890, March 1967.

In `../../../references/20260423_Apollo_radiation`: [1] is `19720017314.pdf`, [2] is
`Apollo_EntryAerodynamics/19690029435.pdf` (public, from
https://ntrs.nasa.gov/citations/19690029435), [3] is
`Apollo_FlightPathAngle/19720013191.pdf` and [4] is
`Paper_STABILITY_CHARACTERISTICS/19670010564.pdf`.
