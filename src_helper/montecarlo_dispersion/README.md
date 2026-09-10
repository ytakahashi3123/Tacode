# src_helper/montecarlo_dispersion

Summarises a Monte-Carlo run made by `src/tacode-montecarlo.py`: it collects the last
point of every case (the impact point, for an entry) and reports how far each one lies
from the others. **Nothing here runs the solver** — the script only reads the files the
cases have already written.

| File | What it is |
|---|---|
| `montecarlo_dispersion.py` | reads every case of a Monte-Carlo run and reports the dispersion of the terminal points |

The Tecplot output is read with `../animate_trajectory/tecplot_reader.py`, so the format
is not written down twice.

## Requirements

Nothing beyond what Tacode itself needs (numpy, PyYAML). `matplotlib` is required only
for `--plot`, and the script says so and carries on without it.

## Usage

```console
cd tutorial/work_montecarlo_wind
./run_tacode-mc.sh
python3 ../../src_helper/montecarlo_dispersion/montecarlo_dispersion.py work_montecarlo_wind
```

The argument is the working directory the Monte-Carlo driver created
(`montecarlo.work_dir` in its `config.yml`), the one holding `case0001`, `case0002`, …
and `case_template`.

```
    Case  Time[s]  Long[deg.]  Lati[deg.]  Alti[km]  East[km]  North[km]  Range[km]  wind.velocity[0]  wind.velocity[1]
-----------------------------------------------------------------------------------------------------------------------
case0001   2775.0  -119.66347   -13.54402   -0.0086    +9.502     -6.099     11.291           27.4989           7.78745
case0002   2775.0  -119.65517   -13.52426   -0.0104   +10.401     -3.913     11.113           28.3371           9.58371
...

Number of cases:  16
Origin: the mean of the cases
Standard deviation:  East 9.047 km,  North 4.605 km,  Up 0.0085 km
Horizontal distance from the origin: median 9.621 km,  maximum 15.929 km
```

The cases are drawn at random, so the numbers differ from run to run; these are from one
run of the tutorial.

`East` and `North` are the offsets of each terminal point from the origin, resolved in
the local horizon at the origin. The frame is **geocentric**, the same one the initial
velocity, the wind and the `Upl/Vpl/Wpl` columns already use.

The columns on the right are the inputs that were actually dispersed. They are not
configured here: the script compares each case's control file with the one in
`case_template`, which the driver leaves next to the cases, and lists whatever differs.
A run that scatters the density instead of the wind therefore reports
`initial_settings.density_factor[0]` without being told to.

### Measuring from a reference case

With `--reference` the offsets are measured from the terminal point of another run
instead of from the mean of the cases. Pointing it at the same entry without wind
separates the two things the wind does — where it moves the mean, and how far it
scatters the cases:

```console
python3 ../../src_helper/montecarlo_dispersion/montecarlo_dispersion.py work_montecarlo_wind \
    --reference ../work_reentry/output_result/tecplot.dat
```

```
Mean offset from the origin:  East +28.808 km,  North +9.664 km
Standard deviation:  East 9.043 km,  North 4.612 km,  Up 0.0399 km
```

The mean offset is the wind itself (20 m/s east, 10 m/s north here), and the standard
deviation is what the +-50 % scatter on it costs. Against the mean of the cases the first
of those is invisible, because it is common to every case.

### The dispersion plot

```console
python3 ../../src_helper/montecarlo_dispersion/montecarlo_dispersion.py work_montecarlo_wind \
    --reference ../work_reentry/output_result/tecplot.dat \
    --plot work_montecarlo_wind/dispersion.png
```

The plot puts the cases in the same east-north plane as the table, and adds

| Drawn | What it is |
|---|---|
| `+` at the origin | the reference case, or the mean of the cases when `--reference` is not given |
| `x` | the mean of the cases; the vector from `+` to `x` is the bias, here the nominal wind |
| shaded ellipses | the 1, 2 and 3 sigma covariance ellipses of the cases. Ellipses rather than circles, because the scatter is not isotropic — a wind uncertainty spreads the impact point mostly along the direction it blows |
| red circle | the CEP 50 %: the radius, measured from the mean, that half the cases fall inside |

`--mark LABEL=PATH` puts a star on another run that is **not** part of the Monte-Carlo
set and leaves it out of every statistic. Pointing it at the same entry flown through a
real wind table shows how far a uniform-wind estimate sits from the measured field:

```console
python3 ../../src_helper/montecarlo_dispersion/montecarlo_dispersion.py work_montecarlo_wind \
    --reference ../work_reentry/output_result/tecplot.dat \
    --mark "NCEP+HWM14 table=../work_reentry_wind_table/output_result/tecplot.dat" \
    --plot work_montecarlo_wind/dispersion.png
```

### Settings file

Everything below can be written in a settings file instead of being typed on the command
line. The tools read **`config_helper.yml` in the current directory** by default, in the
same way the solver reads `config.yml`, and take their own section of it:

```yaml
montecarlo_dispersion:
  directory: work_montecarlo_wind
  reference: ../work_reentry/output_result/tecplot.dat
  mark:
    - NCEP+HWM14 table=../work_reentry_wind_table/output_result/tecplot.dat
  plot: dispersion.png
```

```console
cd tutorial/work_montecarlo_wind
python3 ../../src_helper/montecarlo_dispersion/montecarlo_dispersion.py
```

The keys are the long options without their dashes, with `-` written as `_`
(`--altitude-max` becomes `altitude_max`), and an unknown key stops the run rather than
being quietly ignored. The order of
precedence is **command line > settings file > default**, so a one-off change is still
made on the command line. `-file` points at another file, and `--save-config` writes the
settings in effect back out (comments are not kept, so save a handwritten file under
another name to compare).

### Options

| Option | Meaning |
|---|---|
| `--reference` | Tecplot output of a reference case; the offsets are measured from its terminal point |
| `--mark` | `LABEL=PATH` of another run to mark on the plot; drawn as a star and left out of the statistics. May be repeated |
| `-o, --output` | write the table to this CSV file as well |
| `--plot` | write a scatter plot (east against north) to this file; needs matplotlib |
| `--dpi` | resolution of the plot (default 140) |
| `--tecplot` | path of the Tecplot output inside a case (default `output_result/tecplot.dat`) |
| `--control` | name of the control file inside a case (default `config.yml`) |
| `-file` | settings file to read (default `config_helper.yml` in the current directory) |
| `--save-config` | write the settings in effect to that file (or to the path given) and exit |

## Notes

- The **last row** of each case is used, whatever it is. For an entry that reaches the
  ground it is the impact point; for a run that stops on `time_elapsed_maximum` it is
  simply where the vehicle was at that time. The reported range of terminal altitudes
  tells the two apart.
- Cases whose output is missing are reported and skipped, so a partially finished run can
  still be summarised.
