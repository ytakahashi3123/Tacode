# src_helper/montecarlo_animation

Animates a Monte-Carlo run made by `src/tacode-montecarlo.py`: the trajectories of the
cases and, next to them, the **dispersion ellipses** growing as the cases separate.
**Nothing here runs the solver** — the script only reads the files the cases have already
written.

| File | What it is |
|---|---|
| `montecarlo_animation.py` | draws the cases of a Monte-Carlo run as an animation (or a single frame) |

The Tecplot output is read with `../animate_trajectory/tecplot_reader.py`, and the local
horizon and the covariance ellipse come from `../montecarlo_dispersion/`, so neither the
file format nor the convention is written down twice.

## Requirements

`matplotlib` (an optional dependency of Tacode; `setup_env.sh --with-optional` installs
it). An `.mp4` output additionally needs `ffmpeg`; `.gif` and `.html` do not.

## Usage

```console
cd tutorial/work_montecarlo_wind
./run_tacode-mc.sh
python3 ../../src_helper/montecarlo_animation/montecarlo_animation.py work_montecarlo_wind \
    --reference ../work_reentry/output_result/tecplot.dat \
    --mark "NCEP+HWM14 table=../work_reentry_wind_table/output_result/tecplot.dat" \
    -o montecarlo_dispersion.mp4
```

The argument is the working directory the Monte-Carlo driver created
(`montecarlo.work_dir` in its `config.yml`), the one holding `case0001`, `case0002`, …
The **output format follows the extension** of `-o`, exactly as in `animate_trajectory`:
`.mp4` (needs ffmpeg), `.gif`, or `.html` (one self-contained page, no ffmpeg needed).
With `--snapshot TIME` a single frame is written instead, in whatever image format
matplotlib infers from the extension; a negative time takes the last frame, which is the
figure to put in a report.

## What is drawn

There are five views: `--view flat` (the default), `3d`, `3d-relative`, `globe` and
`follow`. All of them keep the dispersion plane and the history; what changes is the panel
next to them. `3d` and `globe` draw the trajectories **as they are** — the absolute path
the run computed — while `3d-relative` and `follow` measure everything from the reference
case at the same time, which is what shows the dispersion growing.

| Panel | What it shows |
|---|---|
| Ground track (`flat`, upper left) | where the vehicle is. At the scale of a whole entry the cases lie on top of each other, so this panel is the map, not the result; an arrow marks the area the dispersion panel looks at |
| Trajectories (`3d`, left) | the paths themselves, in a box of **longitude, latitude and altitude** — the absolute trajectory, not a difference. The dispersion ellipses lie on the floor at the impact point, and a grey shadow follows each case. At this scale the 100 cases lie on one another, so they are drawn as a single bundle |
| Trajectories (`3d-relative`, left) | the same box, but measured **from the reference case at the same time**, in east and north from it. This is the one that shows the bundle coming down tight and unravelling into the ellipse; the cases are coloured by their wind |
| Trajectories in ECEF (`globe`, left) | the **absolute** trajectories, in ECEF, drawn with the Earth and a graticule every 30 deg. The 100 cases lie within tens of kilometres of each other against a radius of 6378 km, so they are drawn as one bundle in a single colour; the dispersion ellipses are placed on the tangent plane at the impact point, where they are the size of a dot. `--exaggerate` stretches the altitude about the surface (a 150 km descent is 2 % of the radius), and the camera looks at the middle of the track unless `--azimuth` and `--elevation` say otherwise |
| Trajectories (`follow`, left) | the same ECEF axes as `globe`, but with the **camera following the reference case**: the window is a cube around it, so the tens of kilometres of dispersion are not lost against a radius of 6378 km. What is drawn inside it is each case's **departure from the reference at the same time**, carried to the reference's current position — in one frame the vehicle covers some 90 km, so an absolute trail would leave the window at once, while a departure trail keeps the whole history of how the case came away, and its head is the case's true position. The window widens with the dispersion (`--window` fixes it), the cases are coloured by their wind, the ellipses lie on the tangent plane at the reference, the ground appears once it is inside the window, and the bar is the scale, since the axes are switched off. Unlike `3d-relative` the three directions are at the same scale and the ground is there |
| Altitude and dispersion (lower) | the altitude against time in grey, and on the right-hand axis the 1 sigma spread of the cases, east and north. This is where the wind is read off: the spread grows only while the vehicle is in the air the wind is blowing |
| Dispersion | each case as a point, with the path it has taken, measured **from the reference case at the same time** in the local horizon (east, north, km). The 1, 2 and 3 sigma covariance ellipses and the CEP 50 % circle are redrawn every frame, so the cloud is seen growing from a point into the final ellipse |

The horizontal spread of an entry is tens of kilometres against a descent of 150 km, so
the 3D box is **not to scale**: its proportions are set for legibility (`set_box_aspect`)
and its axes are labelled in kilometres so the numbers can still be read off.
`--altitude-max` cuts the box down to the part where anything happens — 40 km for the
wind tutorial, where the dispersion is built during the long terminal descent. Points
above the limit are dropped rather than drawn outside the box, which is what matplotlib
would otherwise do with a 3D line.

The points are coloured by the east component of the wind of each case (the `WindE`
column), so the fan is ordered by the input that made it. A run that scatters something
else — the density, say — leaves every case at the same colour.

`--reference` names the run the offsets are measured from, typically the same entry
without wind: the centre of the cloud is then the mean effect of the wind and the size of
the ellipse is what the uncertainty on it costs. Without it the mean of the cases is used
and the bias, being common to every case, is invisible.

`--mark LABEL=PATH` follows another run that is **not** part of the Monte-Carlo set — for
instance the same entry flown through a real wind table — as a star with its own trail,
and leaves it out of every statistic.

## Settings file

Everything below can be written in a settings file instead of being typed on the command
line. The tools read **`config_helper.yml` in the current directory** by default, in the
same way the solver reads `config.yml`, and take their own section of it:

```yaml
montecarlo_animation:
  directory: work_montecarlo_wind
  reference: ../work_reentry/output_result/tecplot.dat
  view: 3d
  spin: 0.0
  output: montecarlo_dispersion_3d.mp4
```

```console
cd tutorial/work_montecarlo_wind
python3 ../../src_helper/montecarlo_animation/montecarlo_animation.py
```

The keys are the long options without their dashes, with `-` written as `_`
(`--altitude-max` becomes `altitude_max`), and an unknown key stops the run rather than
being quietly ignored. The order of
precedence is **command line > settings file > default**, so a one-off change is still
made on the command line. `-file` points at another file, and `--save-config` writes the
settings in effect back out (comments are not kept, so save a handwritten file under
another name to compare).

## Options

| Option | Meaning |
|---|---|
| `-o, --output` | output file; the extension chooses the format (`.mp4`, `.gif`, `.html`, or an image with `--snapshot`) |
| `--reference` | Tecplot output of a reference case; the dispersion is measured from it |
| `--mark` | `LABEL=PATH` of another run to follow as a star. May be repeated |
| `--snapshot` | write a single frame at this time (s) instead of an animation; a negative value takes the last frame |
| `--frames` | number of frames the run is resampled to (default 180) |
| `--fps` | frames per second (default 20) |
| `--tail` | length of the trail behind each case in seconds; 0 (the default) keeps the whole path |
| `--view` | `flat` (default), `3d`, `3d-relative`, `globe` or `follow` |
| `--window` | half width of the `follow` window, km; 0 (the default) widens it with the dispersion itself |
| `--exaggerate` | stretch the altitude by this factor in the globe view (1.0, the default, is the true scale) |
| `--altitude-max` | upper limit of the altitude axis of the 3D view, km; 0 (default) covers the whole descent |
| `--elevation`, `--azimuth` | camera of the 3D view at the first frame, deg. Left out, the box uses 20 and -62, and the globe looks at the middle of the track. In `follow` the elevation is 20 and `--azimuth` is the offset from the longitude being tracked (-55) |
| `--spin` | how far the 3D camera turns over the animation, deg (default 35) |
| `--dpi` | resolution (default 110) |
| `--embed-limit` | size limit of the frames embedded in an `.html` output, MB (default 512) |
| `--tecplot` | path of the Tecplot output inside a case (default `output_result/tecplot.dat`) |
| `-file` | settings file to read (default `config_helper.yml` in the current directory) |
| `--save-config` | write the settings in effect to that file (or to the path given) and exit |

## Notes

- The cases are resampled onto one time axis with `numpy.interp`, which holds the ends,
  so a case that has already landed stays at its impact point while the others come down.
- The frame is **geocentric**, the same one the initial velocity, the wind and the
  `Upl/Vpl/Wpl` columns already use.
- Longitudes are unwrapped when they are read, so an entry crossing the 180th meridian
  averages correctly; the axis is labelled back in `[-180, 180)`.
- 100 cases at 180 frames take about 30 s to write as an `.mp4` of some 1.3 MB, about
  40 s for `--view 3d`, about a minute for `--view globe` (the whole Earth is redrawn
  every frame) and about 20 s for `--view follow` (3.8 MB; only the patch of ground inside
  the window is drawn, and only once it is inside it).
- The globe and follow views borrow the Earth's surface, the local ground patch and the
  equal-scale box from `../animate_trajectory/`, which is imported only when one of those
  views is asked for, so the tool still runs without matplotlib being importable at
  start-up. The follow camera is the same idea as `animate_trajectory --window`.
- The follow window is a **cube**, so a point is inside it only while its distance from
  the centre is below the half width; the scale bar is placed with that in mind (it was
  disappearing when it was put at 0.8 of the half width in two directions at once).
- A single case can also be animated in three dimensions over the globe, with its
  attitude, by `../animate_trajectory/`; this tool is for what the cases do *relative to
  each other*.
