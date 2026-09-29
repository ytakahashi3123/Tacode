# Mars Pathfinder entry

Tacode against the Mars Pathfinder entry of 1997-07-04, the first validation case at Mars.

Pathfinder flew a **ballistic** entry: a 70° sphere-cone of 2.65 m and 585 kg, spinning at
2 rpm, with no lift, no guidance and a small angle of attack. Nothing about the lift has to
be assumed, which is what makes it the natural first case at Mars. It is also the first
case that needs **the drag coefficient to depend on the Mach number**: the vehicle is
hypersonic at the top and at Mach 1.8 when the parachute mortar fires, and the drag
coefficient falls from 1.71 to 1.55 on the way. Tacode's aerodynamic table got a Mach axis
for it (`database/aerodynamic/README.md`, "Mach number").

The comparison runs from the entry interface (r = 3522.2 km) to the mortar fire, 169.6 s
later at 1.7 km. The parachute phase is not modelled.

## Running it

```bash
cd validation/pathfinder
./run_tacode.sh                  # the three configs, 0.6 s each
python3 compare_pathfinder.py    # the numbers below, and the figures
python3 convert_pds.py           # rebuild reference/ and the two tables from the PDS files
python3 entry_state.py           # the entry state of Spencer et al. table 1 for config.yml
```

`config_mcd.yml` reads the Mars Climate Database table, which is not in the repository (it is
MCD data, see `database/atmosphere/README.md`); `run_tacode.sh` stops and prints the command
that fetches it. It is the same table as `tutorial/work_reentry_mars`, so copying it from there
does as well. `convert_pds.py` needs the PDS files themselves (in
`../../../references/20260929_Mars_EDL/pathfinder_pds`, outside the repository); its outputs
are committed, so the case runs without them.

## The three configs

Each changes one thing, so that each result can be put down to one cause:

| Config | Atmosphere | Drag coefficient | What it tests |
|---|---|---|---|
| `config.yml` | measured in flight | the flight's, **CD(Mach, Kn)** | the solver |
| `config_mcd.yml` | Mars Climate Database | the flight's, CD(Mach, Kn) | the atmosphere model |
| `config_panel.yml` | measured in flight | **Tacode's own panel method, CD(Kn)** | the aerodynamic model |

## What the reference data are

As with `fire2`, the case keeps apart three things that all get called flight data:

- **A — measured.** The axial acceleration from the science accelerometer, 32 Hz, and the
  entry state from orbit determination. These are the only measurements.
- **B — reconstructed.** The trajectory and the atmosphere of the PDS archive (`edl_ddr.tab`).
  NASA integrated the acceleration from the entry state for the trajectory, and took the
  density from the drag equation, `ρ = 2 m a/(CD A V²)`.
- **C — the aerodynamic database.** The `CD` in that equation came from the preflight
  database of Braun et al. (1995), as a function of the Knudsen number, the Mach number and
  the angle of attack. **It is not published as numbers.**

Since B was made from A and C, **C can be traced back out of A and B**:
`CD = 2 m a/(ρ A V²)`, at every point of the flight. That is what `aerodynamic_pathfinder.txt`
is. `edl_ddr.tab` does not carry the speed, so `convert_pds.py` integrates the measured
acceleration from the documented entry state again — the same procedure NASA describes, written
independently. It reproduces the archived positions to 0.2 km horizontally, and the speed it
gives does not depend on which of the two published entry states it starts from (0.02 %).

**`config.yml` is therefore a consistency check, not an independent test of the drag.** The
atmosphere and the drag coefficient were both taken from the flight. What it does test is the
rest of the solver: `CD(Mach, Kn)` and `ρ(h)` are functions of the state, not of the time, so
the measured deceleration comes back only if Tacode, starting from the entry state, flies the
same altitudes at the same times with the same speeds, and the Mach numbers it forms from them
match. The two other configs are real tests, one of the atmosphere model and one of the
aerodynamic model.

## Where every number comes from

| Quantity | Value | Source |
|---|---|---|
| Entry state (Tacode) | 3522.2 km, 22.6303 N, 338.1691 E, inertial 7264.2 m/s at -14.0610 deg., azimuth 253.1479 deg., 1997-07-04 16:51:50.482 UTC | [2] table 1 (orbit determination) |
| Initial velocity for Tacode | [-6974.38, -2042.81, -1764.87] m/s | the same, less the rotation of Mars (`entry_state.py`): 7478.63 m/s, -13.6499 deg., 253.6746 deg. |
| Entry state (reconstruction) | 3597.2 km, 23.0 N, 343.67 E, 7444.7 m/s, -16.85 deg., 255.41 deg., Mars-fixed, 16:51:12.28 UTC | [1] `edlddrds.htm`; used by `convert_pds.py` only |
| Mass, area | 585.3 kg, 5.526 m² | [1] `edlddrds.htm` |
| Acceleration | `reference/pathfinder_acceleration.dat` | [1] `r_sacc_s.tab`, the z axis, 1 g = 9.795433 m/s² |
| Trajectory and atmosphere | `reference/pathfinder_ddr.dat`, `database/atmosphere/atmospheremodel_pathfinder_ddr.txt` | [1] `edl_ddr.tab` |
| Drag coefficient | `database/aerodynamic/aerodynamic_pathfinder.txt` | traced back from [1], see above |
| Speed of sound | γ = 1.335, M = 43.49 g/mol | the ratio of specific heats of the MSL reconstruction [5, figure 8]; the mean molecular weight of `edl_ddr.tab` below 100 km (Viking [6]) |
| Mortar fire | 169.6 s after the entry interface | [2] table 2 |
| Gravity, planet | JGMRO_120F, a sphere of 3396.0 km, IAU 2015 rotation | as `tutorial/work_reentry_mars` |
| MCD atmosphere | the landing site, Ls 142.7°, 03:00 local true solar time, climatology | MCD v6.2 [4], fetched 2026-09-29 |
| Panel-method table | continuum 1.627, free-molecular 2.000, no Mach number | `database/aerodynamic/aerodynamic_mars_spherecone70.txt` |

## What the data turned out to hold

Three things had to be dealt with before any number could be trusted. Each would have gone
unnoticed in a comparison that only looked at the end result.

**1. The accelerometer changes range, and the samples after a change are wrong.** The science
accelerometer switches from 0.016 g to 0.8 g to 40 g full scale as the drag grows. At each
switch the gain column changes at once, but the raw counts take 0.3 to 0.5 s to settle from the
old range, so the calibrated value jumps: at 82.28 s (0.8 g → 40 g) it reads **36 g for a tenth
of a second** where the true value is 0.7 g. Integrated as it is, that is a false deceleration
of 35 m/s, and it made the traced-back drag coefficient 2–4 % too large everywhere below.
`convert_pds.py` replaces the 16 samples (0.5 s) after each of the six switches by linear
interpolation. The sign that this was the right fix: the along-track drift between the
integrated trajectory and the archived one, 4 km by the end, disappeared with it.

After the 40 g switch the reading still sits about 0.06 g above the continuation of the
0.8 g range for a while. That is 7 % of the signal at 60 km and 0.4 % at the peak; the drag
coefficient read between 59 and 62 km is the least certain part of the table.

**2. The entry state printed in the archive does not lie on the archived trajectory.** With its
speed, flight-path angle and azimuth, the position that makes the integration follow
`edl_ddr.tab` is 54 km north and 49 km west of the one printed (23.0 N, 343.67 E). The speed
history does not care (the position moves the trajectory sideways, not along), which is why
`convert_pds.py` can still use it for the speed. Tacode instead starts from the orbit-determination
state of Spencer et al. [2], which is printed to four decimals and runs **parallel to the
archived trajectory at a constant 6.9 km** (6.1 km south, 3.1 km west) through the whole entry.

**3. The archived altitude is 1.2 km above the radius the entry state gives.** Integrating the
measured acceleration from either published entry state puts the vehicle 1.202 km lower than
`edl_ddr.tab` does, with a spread of only 18 m over the entry. A constant with that little
spread is not dynamics — a difference in the entry state or in the drag would grow or shrink
through the flight — it is a difference of datum (the archive gives the altitude above the
landing site, radius 3389.72 km). It is within the 1.7 km uncertainty of the entry radius.
Tacode flies in the frame of the entry state, so **the flight atmosphere table and the
reference altitude are both put in that frame**, 1.202 km below the archive's altitude. That
choice matters for one thing only, the comparison with the MCD, and the numbers are given both
ways below.

One thing was checked and turned out fine: the time tags of the acceleration and of the archive
agree. Between 35 and 55 km, where the density grows by a fifth every second, the traced-back
drag coefficient scatters by 0.05 % about a smooth curve as the two files stand, and by
0.16–0.21 % when one is shifted against the other by 0.1 s (0.3 % at 0.25 s).

## The drag coefficient

![CD against the Mach number](output_comparison/profile_mach-cd.png)

Traced back along the flight (the dots; `reference/pathfinder_reconstruction.dat`), the drag
coefficient falls from 1.706 at Mach 18–22 to 1.667 at Mach 37 and 1.658 at Mach 43, and on
the other side to 1.612 at Mach 6–7 and 1.548 at Mach 1.8. Two effects move it, and in this
flight they act on **different stretches**: the Knudsen number only above 62 km, where the
vehicle is at Mach 41–48, and the Mach number only below that, in the continuum (Kn < 5e-3).
So the table is built in the separable form `CD(Mach, Kn) = CD_c(Mach) + ΔCD(Kn)`, each part read
from its own stretch. It reproduces the traced-back values to **0.12 % rms** (1.7 % at worst,
at the 40 g range switch). Above 77 km the acceleration is too small to read (below 0.5 m/s²);
there the increment is bridged to the free-molecular 2.000 of the panel method. That part of the
entry takes 1.85 m/s off the speed in all, so its drag coefficient does not matter.

The panel method of the tutorial (1.627 throughout the continuum) is 2–5 % low in the hypersonic
part and up to 5 % high below Mach 8.

## Results

Measured on 2026-09-29 with `~/venvs/myenv/bin/python` (Python 3.12.3, numpy 2.2.6,
scipy 1.15.3), `compare_pathfinder.py`. Time from the entry interface; altitude above the
3396 km sphere in the frame of the entry state; the speed against the one integrated from the
measured acceleration; the deceleration against the measured one wherever it exceeds 1 g.

| | `config.yml` | `config_mcd.yml` | `config_panel.yml` |
|---|---|---|---|
| | flight atmosphere, flight CD | MCD, flight CD | flight atmosphere, panel CD |
| Altitude, rms (max) | **0.047 km** (+0.065) | 0.95 km (+1.50) | 0.095 km (-0.17) |
| Speed, rms (max) | **4.5 m/s** (+10.4) | 145 m/s (-333) | 29 m/s (+76) |
| Deceleration / measured, mean | **1.0002** | 0.959 | 1.006 |
| Deceleration / measured, rms of the difference | **0.39 %** | 11 % | 2.5 % |
| Deceleration / measured, range | 0.985 – 1.011 | 0.85 – 1.20 | 0.965 – 1.051 |
| Peak deceleration | **155.36 m/s²** at 77.70 s | 161.50 m/s² at 76.96 s | 152.48 m/s² at 78.16 s |
| At the mortar fire | 1.722 km, 388.8 m/s | 3.148 km, 374.1 m/s | 1.680 km, 377.3 m/s |

The measured peak is **155.38 m/s² (15.84 g₀) at 77.52 s**, and the reconstruction has the
vehicle at 1.677 km and 387.7 m/s at the mortar fire, at Mach 1.78.

![Deceleration](output_comparison/profile_time-deceleration.png)

**The solver (`config.yml`).** Given the atmosphere and the drag coefficient of the flight, Tacode
reproduces the measured deceleration to 0.39 % rms, the peak to 0.02 m/s² and 0.18 s, the
altitude to 47 m and the speed to 4.5 m/s. The largest speed difference, 10 m/s (0.2 %), sits at
the peak deceleration and is what puts the peak 0.18 s late. The 47 m is mostly the 60 m by which
the two entry states differ in radius (Spencer's state against the archive's, both carried into
the frame of the latter). As said above, this is a check of consistency, and it is passed.

**The aerodynamic model (`config_panel.yml`).** Tacode's own drag coefficient, with no Mach
number in it, gives the deceleration 3 % low in the hypersonic part and up to 5 % high at the low
Mach numbers (the lower panel of the figure), which is exactly its difference from the flight's drag
coefficient. The speed is 29 m/s off on average and 76 m/s at worst. **The altitude at the mortar
fire agrees (1.680 against 1.677 km) only because the two errors cancel**; the speed there is
10 m/s low. This is the size of the missing Mach dependence, not an accuracy.

**The atmosphere model (`config_mcd.yml`).** The Mars Climate Database climatology for the site,
season and local time is a different atmosphere from the one Pathfinder flew through:

| Altitude | MCD / flight density (entry-state frame) | (the archive's altitude as it is) | MCD - flight temperature |
|---|---|---|---|
| 100–130 km | 2.39 (1.98 – 3.99) | 2.06 | -2 K |
| 60–100 km | 1.88 (1.04 – 2.60) | 1.59 | **+15 K** |
| 40–60 km | 1.06 (0.91 – 1.21) | 0.91 | 0 K |
| 20–40 km | **1.25** (1.21 – 1.29) | 1.10 | -6 K |
| 10–20 km | **1.26** (1.24 – 1.28) | 1.12 | 0 K |
| 1.5–10 km | 1.16 (0.93 – 1.27) | 1.06 | +11 K |

![Density ratio](output_comparison/profile_density-ratio.png)

The flight was 15 K colder than the climatology between 60 and 100 km, and its density above
60 km was half of the MCD's. Below 40 km the MCD is 25 % denser in the frame Tacode flies in,
and 10–12 % denser on the archive's own altitude: **the 1.2 km datum question above is as large
as the model error itself there**, and this case does not settle it. What does not depend on it
is the effect on the flight: with the MCD, the drag builds up early (20 % above the measured
deceleration at 55–70 s), the peak comes 0.56 s early and 4 % high, the vehicle then decelerates
too little (10–15 % below the measured deceleration after the peak) and arrives at the mortar
fire 1.5 km too high and 14 m/s too slow. One profile of a climatology has no day-to-day
variability in it; the MCD's own `rmsrho` is not used here.

## Time step

`config.yml` at 0.2, 0.1, 0.05 and 0.025 s: the state at the mortar fire moves by 1 cm from 0.2 s
to 0.025 s. At 0.5 s the run stops at 170.0 s instead of 169.6 s, because 169.6 is not a whole
number of steps (the loop runs while the time is below the end), which is why the case uses
0.1 s.

## What is not compared, and why

- **The horizontal position.** The two published entry states differ from the archived trajectory
  by 6.9 and 73 km sideways, far more than anything the solver could contribute, so a
  cross-range number would measure the entry states.
- **The angle of attack.** The entry is flown in three degrees of freedom. The vehicle did
  oscillate, growing in the low supersonic range [2, figure 13], and whatever that did to the drag
  is folded into the traced-back coefficient.
- **The wind.** Neither the reconstruction nor Tacode includes one here.
- **The parachute phase**, after the mortar fire.
- **The MCD's height datum.** The Mars tutorial takes the MCD height to be above the 3396.0 km
  sphere (`zkey = 5`); that value was not confirmed against the MCD documentation for this case.
  An error there would be 0.2 km at most, small against the 1.2 km above.

## References

1. PDS, Mars Pathfinder ASI/MET archive, `MPFL-M-ASIMET-4-DDR-EDL-V1.0` (volume MPAM_0001):
   `edl_ddr.tab`, `edlddrds.htm`, `edl_erdr/r_sacc_s.tab`.
   <https://pds.nasa.gov/data/mpfl-m-asimet-4-ddr-edl-v1.0/mpam_0001/>
2. D. A. Spencer, R. C. Blanchard, S. W. Thurman, R. D. Braun, C.-Y. Peng and P. H. Kallemeyn,
   *Mars Pathfinder atmospheric entry reconstruction*, AAS 98-146 (1998).
3. J. A. Magalhães, J. T. Schofield and A. Seiff, *Results of the Mars Pathfinder atmospheric
   structure investigation*, J. Geophys. Res. 104, 8943 (1999).
4. E. Millour et al., *The Mars Climate Database (version 6.1)*, EPSC 2022; F. Forget et al.,
   J. Geophys. Res. 104, 24155 (1999).
5. C. D. Karlgaard, P. Kutty, M. Schoenenberger and J. Shidner, *Mars Entry Atmospheric Data
   System trajectory reconstruction algorithms and flight results*, AIAA Aerospace Sciences
   Meeting, 2013 (NTRS 20130003188).
6. T. Owen et al., *The composition of the atmosphere at the surface of Mars*, J. Geophys. Res. 82,
   4635 (1977).
7. R. D. Braun, R. W. Powell, W. C. Engelund, P. A. Gnoffo, K. J. Weilmuenster and
   R. A. Mitcheltree, *Mars Pathfinder six-degree-of-freedom entry analysis*, J. Spacecraft
   Rockets 32, 993 (1995) — the aerodynamic database, not consulted directly.

The PDS files and the papers are kept in `../../../references/20260929_Mars_EDL/`, outside the
repository.
