# Atmosphere tables

All three files are NRLMSISE-00 output. `atmosphere.filename_atmosphere` in `config.yml`
selects which one is used; the file format is detected from the file itself, so no
other setting has to change.

| File | Altitude | Format | Profile | Solar activity |
|---|---|---|---|---|
| `atmospheremodel.txt` | 0–400 km, 1 km | CCMC / VITMO | single point, 55°N 45°E | 2015-01-01 01:30 UT, F10.7 left at the site default |
| `atmospheremodel_2015_700km.txt` | 0–700 km, 1 km | pymsis | single point, 55°N 45°E | 2015-01-01 01:30 UT, F10.7 = 140, F10.7A = 146, Ap = 8 |
| `atmospheremodel_700km.txt` | 0–700 km, 1 km | NRLMSISE-00 Fortran | mean of 88 points, 50°S–50°N × 0–315°E | IDAY 0, 00:00 UT, F10.7 = F10.7A = 76.6, Ap = 8.4 7.0 5.0 4.0 3.0 3.9 6.8 |

**The first two are the same atmosphere; the second just goes higher.**
`atmospheremodel_2015_700km.txt` was generated to reproduce `atmospheremodel.txt` and
carry on above its top, and it does: over 100–400 km the densities agree to **0.37 % rms**
and at 400 km itself to 0.07 %. The cases that fly above 400 km — `tutorial/work`,
`tutorial/work_montecarlo` and `tutorial/template` — use it; the re-entry cases, which
start at 150 km and descend, keep `atmospheremodel.txt`.

`atmospheremodel.txt` is kept unchanged as the reference table.

This directory is the master copy. Each case directory keeps the table it actually uses
in `<case>/database/atmosphere` and reads it from there, so to run a case against
`atmospheremodel_700km.txt` copy the file into that directory first and then name it in
`atmosphere.filename_atmosphere`.

**The two are not interchangeable.** They were generated for different conditions, and
the difference is large enough to matter:

| Altitude | `atmospheremodel.txt` | `atmospheremodel_700km.txt` | Ratio |
|---|---|---|---|
| 200 km | 2.484e-10 kg/m³ | 1.866e-10 kg/m³ | 1.33 |
| 300 km | 1.693e-11 kg/m³ | 9.197e-12 kg/m³ | 1.84 |
| 400 km | 1.938e-12 kg/m³ | 8.950e-13 kg/m³ | 2.17 |

Temperature at 400 km is 823.3 K against 783.0 K. Running the tutorial case with the
700 km table instead of the default moves the position by about 3.5 km after 5000 s —
that is the change in solar conditions, not the change in altitude range.

## Choosing between them

Use `atmospheremodel_2015_700km.txt` for a trajectory that goes above 400 km and
`atmospheremodel.txt` for one that does not; they are the same atmosphere, so the choice
costs nothing either way.

Use `atmospheremodel_700km.txt` when the trajectory spends time well above 400 km and
low solar activity is the condition you want. Above the top of whichever table is
loaded, `Tacode` extrapolates the density exponentially with a scale height fitted to
the table (see the Requirements and Notes sections of the top-level `README.md`), so a
shorter table is not a hard limit — but real data is better than extrapolated data.

To use it:

```yaml
atmosphere:
  filename_atmosphere: atmospheremodel_700km.txt
```

## How `atmospheremodel_2015_700km.txt` was produced

CCMC does not record which solar indices it used for `atmospheremodel.txt` — its header
says only `not specified` for all three — so they were recovered by fitting. Running
NRLMSISE-00 through `pymsis` at the same date and point over a grid of F10.7, F10.7A and
Ap, and minimising the root-mean-square difference of `log(rho)` above 100 km, gives

    F10.7 = 140, F10.7A = 146, Ap = 8

with a residual of 0.37 %, which is small enough that the two tables can be treated as the
same atmosphere. The file is then:

```console
python3 generate_atmosphere_table.py --datetime 2015-01-01T01:30:00Z \
        --longitude 45.0 --latitude 55.0 --f107 140 --f107a 146 --ap 8 \
        --altitude 0 700 --altitude-step 1 -o atmospheremodel_2015_700km.txt
```

## How `atmospheremodel_700km.txt` was produced

`NRLMSISE-00_readctl_2fequal.FOR` (kept outside this repository, alongside the
NRLMSISE-00 source). For each altitude from 0 to 700 km in 1 km steps it calls `GTD7`
at 11 latitudes (50°S to 50°N, every 10°) and 8 longitudes (0° to 315°, every 45°),
computes the local solar time for each longitude, and writes the arithmetic mean of the
88 results. The `_2fequal` variant sets `F107 = F107A`, which is why both appear as
76.6 in the file header even though the control file lists 70.4 for `F107`.

The columns are the same seven in both files: height, O, N2, O2, total mass density,
neutral temperature, N.

## Knudsen number

`Tacode` builds the Knudsen number from N2, O2, N and O only, since those are the
species in these tables. They account for 96% of the mass at 400 km, but only 66% at
600 km and 35% at 700 km, where helium dominates. The drag coefficient in
`database/aerodynamic/aerodynamic.txt` is flat over that Knudsen range, so the effect
on a trajectory is negligible, but the reported Knudsen number is increasingly
approximate above roughly 500 km.
