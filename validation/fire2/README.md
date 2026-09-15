# Project Fire flight II entry

Tacode against the Project Fire flight II reentry of 1965-05-22, in three and in
six degrees of freedom.

**This is the case where nothing at all is assumed about the attitude.** The Fire
reentry package was completely uncontrolled: no reaction jets, no guidance, no
roll programme. It was spin-stabilised, axisymmetric, and entered at 11 351 m/s
on a −14.7° path. There is no bank angle to guess (as `validation/apollo4` had to)
and no roll history to supply (as `validation/apollo10` does). The inputs are the
entry state, the mass, the inertias, the reference diameter, the spin rate and the
measured aerodynamics — every one of them a number printed in the two NASA reports.

It was chosen as the **six-degree-of-freedom** validation case: it is the one
uncontrolled, axisymmetric, hypersonic entry with a published attitude
measurement (rate gyros on all three axes).

## What this case can and cannot check

Decided before any number was measured, so that the write-up is not a selection
of the intervals that happen to agree.

**It can check:**

- the **trajectory** — altitude, velocity and flight-path angle over the whole
  312 s from the 121 920 m entry interface to 79 m above the sea
- the **atmosphere table**, against the density measured by two Nike-Apache
  sounding rockets fired from Ascension Island 3 and 15 hours after impact
- the **period of the pitch-yaw oscillation**, which is set by the static moment,
  the inertia, the spin and the dynamic pressure together

**It cannot check** — and the reasons are in the physics, not in the data:

- the **growth of the angle-of-attack envelope** from 3° at entry to 7.7°, 13° and
  finally 19.5°. Every step of that growth was made by a discrete disturbance: the
  separation spring at t = 1610.43 s, the asymmetric melting of the first beryllium
  calorimeter at 1640.38 s, the ejection of the second phenolic heat shield at
  1648.18 s, and a series of events from 1652.7 s that the report could not
  identify. Tacode has no disturbance events and assumes an axisymmetric body, so
  it can reproduce none of them
- the **decay of the spin**, which fell from 17.9 rad/s at separation to 13.6 rad/s
  by t = 1658 s and then to 6.6 rad/s by 1662 s through a disturbance the report
  could not explain either. With no roll damping derivative in the references,
  Tacode holds the spin constant
- **anything below about Mach 3** (t > 1700 s). Neither the Tacode aerodynamic
  table nor NASA's own table V carries a Mach axis, so the hypersonic drag is flown
  all the way to the sea. The trajectory numbers below are split at t = 1700 s for
  that reason

## Running it

```bash
cd validation/fire2
./run_tacode.sh                             # 3-DOF, drag only, 2.0 s
./run_tacode.sh -file config_matched.yml    # the same at table V's own drag
./run_tacode.sh -file config_6dof.yml       # 6-DOF, 12 s
python3 run_segments.py                     # 6-DOF as a chain over the shield ejections, 14 s
python3 compare_fire2.py                    # the trajectory numbers, and figures
python3 compare_attitude.py                 # the attitude numbers, and a figure
```

Rebuilding the inputs from the scanned reports (they live in `../../../references/pdf`,
outside the repository, and need `pdftotext`, `pdftoppm` and Pillow):

```bash
python3 entry_state.py            # table V converted for config.yml
python3 initial_attitude.py       # the initial Euler angles and spin for config_6dof.yml
python3 extract_trajectory.py     # reference/fire2_trajectory.dat from table V
python3 digitize_aerodynamics.py  # the aerodynamic table from figure 4
python3 digitize_figure.py        # reference/ from figures 17 and 19
```

The results of all of those are committed under `reference/` and `database/`, so
the case runs without the reports. The two tables under `database/` are copies of
the masters at the top of the repository, as every case is; `digitize_aerodynamics.py`
writes both, and the atmosphere table came from

```bash
python3 ../../database/atmosphere/generate_atmosphere_table.py \
        --datetime 1965-05-22T22:21:57Z --longitude -16.9 --latitude -7.4 \
        --f107 75 --f107a 75 --ap 4 --altitude 0 200 --altitude-step 1 \
        -o ../../database/atmosphere/atmospheremodel_fire2.txt
```

## The trajectory reference is a simulation, not a measurement

**NASA TN D-3569 table V is titled "REENTRY TRAJECTORY PARAMETERS OBTAINED FROM
COMPUTER SIMULATION".** It is NASA's own three-degree-of-freedom particle
trajectory, started from the Ascension Island TPQ-18 radar at t = 1608 s
(150 578 m, 11 326.6 m/s, −15.30°, azimuth 122.987°) and flown through the
sounding-rocket atmosphere with the US Standard Atmosphere 1962 above it. The
report shows it merging with the radar beacon track on both sides of the telemetry
blackout and reaching the measured impact point within 500 m.

So agreement with table V is a **check of the solver**, not of the physics of the
flight: same entry state, same measured atmosphere, does Tacode fly the same
trajectory. That is worth having — it is the sharpest such check in `validation/`,
because table V is a numerical table at 0.5 s intervals rather than a figure to be
digitised — but it is a different claim from the Apollo cases, and it is labelled
as such in `reference/fire2_trajectory.dat` itself.

The **attitude** comparison is against a measurement: the rate gyros.

### Table V was flown at a constant ballistic coefficient, and that was fair

Table V prints both the dynamic pressure and the acceleration, so the ballistic
coefficient it was flown with can be read straight back out of it:

    CD S / m = -a g / q

and it comes to **6.1490e-3 m²/kg at every one of the 500 usable records, to four
figures**. NASA therefore flew a single constant value: no mass loss, no Mach
dependence, no incidence dependence. With m = 86.568 kg and S = 0.354657 m² that
is **CD = 1.5009**, 1.6 % above the 1.4766 that figure 4 of TN D-4183 gives.

Holding the ballistic coefficient constant was not an approximation they had to
apologise for. The vehicle **sheds mass and diameter together**: over the five
configurations of table I the mass falls 24 % and the reference area 24 %, so

| Configuration | Mass, kg | d, m | m/(π d²/4), kg/m² | Iy/d³, kg/m |
|---|---|---|---|---|
| Complete reentry package | 86.568 | 0.672 | 244.1 | 9.25 |
| Less first calorimeter | 83.189 | 0.651 | 249.9 | 9.58 |
| Less first phenolic layer | 76.022 | 0.630 | 243.9 | 9.22 |
| Less second calorimeter | 72.166 | 0.607 | 249.4 | 9.52 |
| Less second phenolic layer | 66.179 | 0.587 | 244.5 | 9.18 |

— the ballistic coefficient varies by **±2.4 %** and Iy/d³, which sets the attitude
frequency, by **±3.6 %**. `run_segments.py` flies the steps anyway and the result
is below; it changes almost nothing, which is the point.

## Where every number comes from

| Quantity | Value | Source |
|---|---|---|
| Entry state | −16.9054° E, −7.3966° N, 120.478 km, 11 351.53 m/s, −14.708° (geodetic horizon), azimuth 122.861° | [1] table V at t = 1618.25 s |
| Initial velocity for Tacode | [9222.71, −5955.14, −2887.09] m/s | the same, converted by `entry_state.py` |
| Mass | `86.568 kg` | [1] table III. A measured mass: [1] table I builds it up from the 2004.88 kg at lift-off by subtracting each jettisoned item, ending at 86.57 kg |
| Reference diameter, area | `0.672 m`, `0.354657 m²` | [2] table I |
| Moments of inertia | `Ix 3.511`, `Iy 2.806`, `Iz 2.874 kg m²` | [2] table I — **measured, not estimated** |
| Centre of gravity | 0.029 m forward of the aerodynamic reference point | [2] table I, x_mc − x_cg = 0.306 − 0.277 |
| Spin at separation | `171 rpm = 17.907 rad/s` | [1], from the telemetry signal strength |
| Initial angle of attack | `1.52°` | [2], the angle between the body axis and the velocity at separation |
| **Drag coefficient** | `1.4766` | **[2] figure 4**, \|Cx\| at α = 0, digitised by `digitize_aerodynamics.py` |
| Aerodynamic table | `database/aerodynamic/aerodynamic_fire2.txt` | [2] figure 4, all three panels |
| Trajectory reference | `reference/fire2_trajectory.dat` | [1] table V, extracted by `extract_trajectory.py` |
| Frequency reference | `reference/fire2_dynamic_pressure.dat` | [2] figure 17, the solid segments |
| Angle-of-attack envelope | `reference/fire2_angle_attack.dat` | [2] figure 19 |
| Atmosphere | `database/atmosphere/atmospheremodel_fire2.txt` | NRLMSISE-00 for 1965-05-22 22:21:57 UTC at the entry point. F10.7 = F10.7A = 75, Ap = 4 |
| Epoch | 1965-05-22T22:21:57.953Z | lift-off 21:54:59.703 GMT [1] plus 1618.25 s |

Two of those need a word.

**The mass is 86.568 kg, not the 86.586 kg that [2] table I prints.** [1] table I
builds the entry mass up from lift-off and ends at 94.73 − 8.16 = 86.57 kg, which
86.568 rounds to and 86.586 does not, and [1] table III prints 86.568. The
difference is 0.02 %, but the two reports disagree and this is why one was chosen.

**The run starts at t = 1618.25 s**, half a second after the 121 920 m entry
interface, because the scan lost the latitude and longitude of the first record.
That is 1.4 km of altitude and 1 m/s of speed; the comparison is made in elapsed
flight time either way.

### The aerodynamics are measured, and they are the data the report itself used

Figure 4 of TN D-4183 gives Cx, CR and Cm against angle of attack from 0 to 80°,
measured in the **Ames hypersonic free-flight facility at Mach 35** on a basic
Apollo body shape. TN D-4183 used exactly this figure to turn the flight rate-gyro
records into angles of attack, so the coefficients Tacode flies and the
coefficients behind the "measured" angles of attack are the same numbers.

`digitize_aerodynamics.py` reads all three panels from a 600 dpi scan. The
calibration is checked against the figure's own grid lines, which land on their
nominal values to **0.005 in Cx, 0.002 in CR and 0.0003 in Cm**; that is the
accuracy of the digitisation. At α = 0 the value of CR and Cm is not read but taken
as zero, which is exact for an axisymmetric body and is what the figure draws — the
left edge of the panel is where the axis and the Cm = 0 rule lie on top of the curve.

**The origin of the angle-of-attack axis is not the left edge of the frame.** The
sixteen vertical grid lines are 5° apart and the leftmost is α = 5°, not 0; α = 0 is
144 pixels further left. Taking the frame for the origin shifts every coefficient by
5°, which was caught by the physics rather than by eye: it put CD at 1.467 instead
of 1.477 and, far worse, made the small-angle slope of Cm twice as steep, because
the first 5° of that curve is where most of its curvature is. The calibration used
now puts α = 80° within one pixel of where the PDF's own text layer puts the "80"
label.

The resulting numbers: **CD(0) = 1.4766**, Cm(5°) = −0.01126 so Cm_alpha = −0.129
per radian about the figure's own reference point, Cm minimum −0.064 at 40°, Cm back
through zero near 75°. The table has **no Knudsen dependence** — figure 4 is one test
condition — so it is written at the same value for every Knudsen number. The
consequence is confined to the first 20 s, where the flow is rarefied but the
dynamic pressure is below 2 N/m².

### Reading the scanned table

`extract_trajectory.py` recovers table V from the 1966 scan, whose embedded OCR is
badly damaged. It is documented in the file itself; the three things worth knowing
here are that **nothing is filled in by interpolation** (a value that cannot be read
is written as `nan` and counted), that three of the columns are printed twice, in SI
and in US customary units, which gives a second independent reading of the same
number, and that the row a value belongs to is decided by **clustering the vertical
positions within each record** rather than by fixed thresholds — some pages are
tilted enough that the right-hand columns of a row sit higher than its left-hand
ones by most of a row spacing, and a threshold there puts the density into the
flight-path-angle column.

Measured on 2026-09-14: **597 of the 625 records** are recovered (28 have an
unreadable time column and are dropped rather than guessed at) and **918 of the 7761
values are `nan`**. Per column the fraction recovered runs from 82 % (Reynolds
number) to 97 % (pressure); altitude 95 %, velocity 91 %, density 93 %.

The identities the table has to satisfy are reported as a diagnostic, not enforced:
q = ρV²/2 holds to a median of 4e-7 and M = V/a to 5e-6. p = ρRT holds to 2e-5
except above 95 km, where the composition of the atmosphere changes and R is no
longer 287.05 — that is the table, not the scan.

## Results: the trajectory

Measured on 2026-09-14 with `~/venvs/myenv/bin/python` (Python 3.12.3, numpy 2.2.6,
scipy 1.15.3). The three-degree-of-freedom run takes 3.0 s at dt = 0.02 s.

Difference from table V, over the records where the scan preserved the value.

| | | Altitude | Relative velocity | Flight-path angle |
|---|---|---|---|---|
| **CD = 1.4766** (figure 4) | whole entry, 312 s | max +231 m, **rms 126 m** | max +65 m/s, rms 12.2 m/s | max −0.50°, rms 0.16° |
| | hypersonic, to t = 1700 s | max −179 m, rms 117 m | max +65 m/s, rms 23.1 m/s | max +0.24°, rms 0.11° |
| **CD = 1.5009** (table V's own) | whole entry | max +333 m, rms 182 m | max −59 m/s, rms 13.0 m/s | max −0.50°, rms 0.16° |
| | hypersonic, to t = 1700 s | max −79 m, **rms 48 m** | max −59 m/s, rms 25.0 m/s | max +0.21°, **rms 0.09°** |
| Mass stepped by segment | whole entry | max +217 m, rms 124 m | max +65 m/s, rms 12.0 m/s | max −0.51°, rms 0.17° |

Four things to read out of this table.

- **With the input made identical, Tacode reproduces NASA's integration of the same
  entry to 48 m rms in altitude and 0.09° in flight-path angle** over the 82 s from
  120 km to 20 km, through a peak deceleration of 73 g. That is 4e-4 of the altitude
  and is the sharpest solver check in `validation/`.
- The 117 m of the figure-4 run is **not solver error, it is the 1.6 % difference
  between the drag coefficient measured at Mach 35 and the one NASA flew.** Which of
  the two is right is a question about the 1965 aerodynamic estimate, not about
  Tacode.
- **Stepping the mass and the diameter through the five configurations changes
  nothing** (126 m against 124 m), for the reason in the table above: this vehicle
  sheds mass and frontal area in the same proportion.
- Below Mach 3 the runs drift by a few hundred metres, which is the missing Mach
  dependence of the drag showing up. The velocity difference peaks at 65 m/s on
  1800 m/s at t = 1660 s, 3.6 %.

The time step does nothing: from dt = 0.005 s to dt = 0.1 s the altitude rms moves
from 122.3 m to 120.6 m. The difference from table V is in the physics of the two
runs, not in their arithmetic.

### The atmosphere table against the density measured in flight

Table V carries the ambient density, pressure and temperature of the atmosphere the
trajectory was flown through, which came from the two Nike-Apache soundings.
Comparing it with the NRLMSISE-00 table Tacode flies does not involve the trajectory
at all, so it is an independent check of the table.

| Altitude | Points | Mean ratio, flight / table |
|---|---|---|
| all, 0 to 122 km | 535 | **0.9938** |
| 0 to 20 km | 392 | 0.999 |
| 20 to 40 km | 91 | 0.999 |
| 40 to 60 km | 15 | 0.999 |
| 60 to 80 km | 11 | 1.093 |
| 80 to 100 km | 12 | 0.965 |
| 100 to 122 km | 14 | 0.750 |

The mean ratio is 0.994 with an rms scatter of 7.9 %. Below 60 km the two agree to
0.1 %, which is less a triumph than a statement that both are close to a standard
atmosphere there. The 25 % at 100 to 122 km is the thermosphere, where NRLMSISE-00
depends on the solar indices: F10.7 = 75 was used for the May 1965 solar minimum
because the indices of the day were not available offline. The trajectory does not
care — the dynamic pressure is below 2 N/m² above 100 km.

## Results: the attitude

The six-degree-of-freedom run takes 24 s at dt = 0.005 s, which is what the
oscillation needs: its period is 0.14 s at peak dynamic pressure and the spin period
is 0.35 s.

### What is compared, and why it is the frequency

Figure 17 of TN D-4183 plots two things: the dynamic pressure of the trajectory
(dashed) and **the dynamic pressure the report's own six-degree-of-freedom
simulations needed in order to match the frequency of the measured pitch and yaw
rates** (solid, in short segments). For a spinning axisymmetric body the transverse
rate seen in body axes oscillates at

    lambda = p - a + sqrt(a² + k q),    a = Ix p / (2 Iy),    k = -Cm_alpha S d / Iy

so at a fixed spin, inertia and Cm_alpha those solid segments **are** the measured
frequency, written as a dynamic pressure. Comparing them with Tacode's own dynamic
pressure is therefore comparing the predicted oscillation period with the measured
one.

`compare_attitude.py` does not assume that relation, it **fits it to the Tacode run**
and reports how well it holds — which is a check in its own right, of the
six-degree-of-freedom dynamics against the analytic tricyclic solution:

- the relation fits 75 windows of the run to **0.05 rad/s rms** on frequencies of 27
  to 46 rad/s, that is 0.15 %
- the fitted k corresponds to **Cm_alpha = −0.1418 per radian**, against the −0.1423
  that the digitised table gives once its moment is transferred from the figure's
  reference point to the centre of gravity (the table itself gives −0.1290 about its
  own reference point). **0.4 %** — so the aerodynamic table, the centre-of-gravity
  transfer and the rigid-body solver are mutually consistent

### The oscillation frequency

Four of the six simulation windows of figure 17 come out of the scan cleanly; the
other two lie on the steep flank of the dashed curve and cannot be separated from it.

| Window, s | q̄ Tacode, N/m² | q̄ required, N/m² | ratio | ω computed, rad/s | ω measured, rad/s | ratio |
|---|---|---|---|---|---|---|
| 1644.9 | 66 652 | 67 569 | 1.014 | 37.17 | 37.35 | 1.005 |
| 1649.3 | 117 344 | 133 438 | 1.137 | 45.93 | 48.33 | **1.052** |
| 1652.6 | 108 241 | 112 185 | 1.036 | 44.51 | 45.13 | 1.014 |
| 1662.4 | 26 337 | 26 156 | 0.993 | 27.74 | 27.69 | 0.998 |

**The oscillation Tacode predicts is within 5 % of the one the vehicle flew** — 1.7 %
on average over the four windows, and the largest error is at the peak of the
dynamic pressure. The run through the segmented configurations gives the same answer
(0 to 5.0 %, 1.6 % on average).

That 5 % is close to the resolution of the test. The report attributes the gap
between the two curves of its own figure 17 to "uncertainties in the atmospheric
measurements or the wind-tunnel data, or both", and the digitisation of figure 17 is
itself worth about 600 N/m² rms — which is how well the dashed curve of that figure
reproduces the dynamic pressure of table V, and therefore how well the calibration
is known.

### The angle-of-attack envelope, which does not compare

Over the interval figure 19 covers, Tacode's total angle of attack stays between
0.64° and 5.29°, while the measured envelope climbs from 0.36° to 19.59°. **These are
not comparable and no error is quoted for them.** The measured envelope is a
staircase and every step is a disturbance: the first beryllium calorimeter melting
asymmetrically at 1640.38 s (an angular impulse of 6.440 N·m·s), the second phenolic
heat shield coming off at 1648.18 s (9.49 N·m·s), and a series of events from 1652.7 s
that the report itself could not identify. Tacode has no way to represent any of
them. `output_comparison/profile_attitude.png` draws the two together so that the
size of what is missing is visible.

The angular impulses **could** have been put in — the report measures their magnitude
— but not their direction, so a run that reproduced the envelope would have had its
direction chosen to make it do so. That is fitting, not validation, and it is why
`run_segments.py` steps only the mass, the inertia and the diameter.

## References

1. Lewis, John H., Jr.; and Scallion, William I.: *Flight Parameters and Vehicle
   Performance for Project Fire Flight II, Launched May 22, 1965.* NASA TN D-3569,
   1966. (NTRS 19660025567)
2. Scallion, William I.; and Lewis, John H., Jr.: *Body Motions and Angles of Attack
   During Project Fire Flight II Reentry.* NASA TN D-4183, 1967. (NTRS 19670029017)
3. Murphy, Charles H.: *Free Flight Motion of Symmetric Missiles.* BRL Report 1216,
   1963 — the tricyclic solution the frequency relation above comes from. A copy is
   in `../../../references/pdf`.
