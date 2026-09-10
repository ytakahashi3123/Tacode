# Wind tables

Tables of the wind read by `wind.kind_wind_model: fileread`. The solver itself is not
built to reach out to a data service, and its dependencies are kept to
numpy / scipy / PyYAML / simplekml, so the tables are prepared offline —
`generate_wind_table.py` here fetches and converts, and `src/wind/wind.py` only reads
text. This is the same division as the atmosphere, where NRLMSISE-00 is sampled into
`database/atmosphere/*.txt` rather than being called during a run.

## File format

One format, whatever the dimensionality. Comment lines start with `#`; every other line
holds six numbers, or seven when the table carries a time axis:

```
# Longitude[deg.] Latitude[deg.] Altitude[km] East[m/s] North[m/s] Up[m/s]
     235.00000   -10.00000     0.13440     -5.68424      1.18740      0.00000
```

```
# Time[s] Longitude[deg.] Latitude[deg.] Altitude[km] East[m/s] North[m/s] Up[m/s]
      0.000     235.00000   -10.00000     0.13967     -5.68424      1.18740      0.00000
```

The time is measured in seconds from the `# Epoch (UTC):` line of the table. Every data
row must have the same width; a table that mixes six- and seven-number rows is rejected.

The grid is reconstructed from the rows, so their order does not matter, but it must be
filled completely: a table whose rows do not add up to
(longitudes × latitudes × altitudes) is rejected rather than padded, because a missing
point and a calm point would otherwise be indistinguishable.

**A one-dimensional vertical profile is a field with one longitude node and one latitude
node.** Axes carrying a single node are dropped when the interpolator is built, so the
profile becomes a function of altitude alone and travels the same code path as a full
field. The same applies to a single-altitude table (constant in the vertical), to a
single-time table (one snapshot) and to a one-point table (a uniform wind).

Longitudes are folded into `[-180, 180)`, so a table written in the meteorological
`0-360` convention is read correctly, and `-180` and `180` are the same meridian and must
not both appear. A table that covers the whole globe is detected from its spacing and a
wrapped node is added at the seam, so interpolation stays continuous across it.

A line of the form

```
# Epoch (UTC): 2024-01-01T00:00:00.000Z
```

is read as the time the data belongs to. For a single snapshot, the run warns and
continues if the `epoch` section of `config.yml` names a different time. For a table with
a time axis the epoch is not optional: it is what places the run on that axis, so a table
with a time axis and no epoch — in the table or in the configuration — stops the run.
Outside the time range of the table the nearest snapshot is used, and the run warns once.

### Sign and frame

`East`, `North` and `Up` are the components in the **local geocentric** horizon, the same
frame `initial_settings.velocity` is given in. It is geocentric rather than geodetic, so
it is tilted from the true horizontal by up to 0.19 deg; see the note in the top-level
`README.md`. Meteorological `u` and `v` are defined in the geodetic horizon and inherit
that tilt when converted here, which mixes about 0.2 m/s of vertical component into a
horizontal wind of 60 m/s.

The vertical wind is written as zero. Synoptic-scale vertical motion is of the order of
cm/s and does not matter for a trajectory, and NCEP carries it as `omega` in Pa/s, whose
conversion to m/s needs the density.

## `generate_wind_table.py`

```console
# A region of NCEP/NCAR Reanalysis 1 at a given time
python3 generate_wind_table.py --source ncep --datetime 2024-01-01T00:00:00Z \
        --longitude -130 -110 --latitude -25 -5 -o wind_ncep_20240101_pacific.txt

# The same region and time from HWM14, 20 to 500 km, on the axes of the table above
python3 generate_wind_table.py --source hwm14 --datetime 2024-01-01T00:00:00Z \
        --like wind_ncep_20240101_pacific.txt --altitude 20 500 --altitude-step 20 \
        --hwm14-data ~/hwm14 -o wind_hwm14_20240101_pacific.txt

# Join them with a transition layer between 20 and 30 km
python3 generate_wind_table.py --source merge \
        -i wind_ncep_20240101_pacific.txt wind_hwm14_20240101_pacific.txt \
        --transition 20 30 -o wind_merged_20240101_pacific.txt

# Twelve hours of NCEP, which gives the table a time axis
python3 generate_wind_table.py --source ncep --datetime 2024-01-01T00:00:00Z \
        --duration 43200 --longitude -130 -110 --latitude -25 -5 -o wind_series.txt

# A vertical profile from a text file of altitude[km] east[m/s] north[m/s]
python3 generate_wind_table.py --source profile -i profile.txt --centre 140 35 \
        -o wind_profile.txt

# Zero wind everywhere, to check the file path and the reading
python3 generate_wind_table.py --source calm -o wind_calm.txt
```

`--source ncep` reads NOAA PSL's OPeNDAP service and **needs nothing beyond the standard
library and numpy**: it requests the `.ascii` response and parses that, so neither
`netCDF4` nor `cfgrib` has to be installed. Only the range asked for is transferred, so a
region of a few degrees is a handful of small requests. The data are the reanalysis, that
is, past observations assimilated into a model — not a forecast.

Its vertical extent is the limit worth knowing: the pressure levels run 1000 to 10 hPa,
which is **the ground to about 31 km**. Above that the table has nothing, and
`wind.kind_extrapolation: zero` (the default) leaves the air co-rotating with the Earth
there. Clamping the top value upwards instead would stretch the stratospheric wind into
the thermosphere, which is why `zero` rather than `clamp` is the default.

The wind sits on pressure levels, so the altitudes come from the geopotential height
(`hgt`) of the same time and grid. That height varies from place to place, so the
region-mean height of each pressure level is taken as the common altitude node and every
vertical column is interpolated onto those nodes. Geopotential height is converted to
geometric altitude by `z = R H / (R - H)`, a correction of 0.14 km at 30 km.

`--duration` asks for a time axis: NCEP is six-hourly, so every sample covering
`[--datetime, --datetime + duration]` is fetched and written with its offset in seconds in
front. The altitude nodes are then the region- **and time**-mean geopotential height of
each pressure level, so that the grid stays rectangular.

### HWM14

`--source hwm14` covers the thermosphere, where NCEP has nothing: HWM14 (Drob et al.,
*Earth and Space Science* **2**, 301-319, 2015) reaches 500 km. It is Fortran plus the NRL
coefficient files `hwm123114.bin`, `dwm07b104i.dat` and `gd2qd.dat`, none of which is a
pip install, so **it is not a dependency of this script**: the module is imported only
when `--source hwm14` is asked for, and if it is absent the script says how to obtain and
build it and stops.

To build it, put `hwm14.f90` and the three data files in one directory and run

```console
python3 -m numpy.f2py -c hwm14.f90 -m hwm14
```

which needs a Fortran compiler and, from numpy 1.26 on, `meson` and `ninja`. Then run
this script with that directory on `PYTHONPATH` and named by `--hwm14-data`, since the
model opens its data files from the current directory.

HWM14 returns a meridional and a zonal wind and no vertical component. `--ap` sets the
three-hour ap index the disturbance model uses; `-1`, the default, leaves the quiet
component alone. `--altitude`, `--altitude-step`, `--longitude-step` and `--latitude-step`
set the grid, and `--like` copies the time and horizontal axes from an existing table,
which is how a table is made that can be merged with an NCEP one.

### Merging a lower and an upper table

`--source merge` joins two tables along the altitude with a transition layer. Switching
at a single altitude would make the wind jump — the two models disagree by tens of m/s
where they meet, since NCEP's top level and HWM14's lower boundary are each at the edge of
their validity — so the weight moves linearly across `--transition LOW HIGH`. Below `LOW`
the lower table is used alone, above `HIGH` the upper one, and where only one of the two
has data that one is used whatever the weight says.

The two tables must share their time, longitude and latitude axes, which is what `--like`
is for; a mismatch is refused rather than silently reinterpreted. The altitudes of the
result are the union of both, and the transition layer must lie inside both ranges.

## Files here

| File | Contents |
|---|---|
| `wind_ncep_20240101_pacific.txt` | NCEP/NCAR Reanalysis 1 at 2024-01-01 00 UTC, 130-110°W and 25-5°S, 17 levels to 31 km (9 × 9 × 17 points). The region around where `tutorial/work_reentry` comes down; public-domain NOAA data |
| `wind_merged_20240101_pacific.txt` | The table above joined with HWM14 over the same region and time, 0 to 500 km with a 20-30 km transition layer (9 × 9 × 42 points). The one to use for a re-entry that starts above the stratosphere |
| `wind_profile_sample.txt` | A one-dimensional profile with a mid-latitude jet peaking at 10 km, to show the degenerate case. Illustrative, not measured |

None of them is used by any tutorial case: every shipped `config.yml` has
`wind.flag_wind: False`. Copy the one you want into the case's own `database/wind/`, as
the atmosphere and aerodynamic tables are copied, and switch the flag on. All three cover
the region `tutorial/work_reentry` comes down in and nothing else, so the horizontal
clamp holds the edge value over the rest of that trajectory; regenerate them for another
region when it matters.

The merged table was produced with HWM14 built as described above. It is model output
rather than the NRL coefficient files, so it is a plain table here, but reproducing it
needs HWM14.

## Sources

- NCEP/NCAR Reanalysis 1 — https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.pressure.html
- NOAA PSL OPeNDAP — https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis/pressure/catalog.html
- HWM14 — https://map.nrl.navy.mil/map/pub/nrl/HWM/HWM14/HWM14_ess224.pdf
