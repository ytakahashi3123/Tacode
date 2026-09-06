# src_helper/animate_trajectory

Animates a Tacode trajectory, with the attitude if the output has it. **Nothing here runs
the solver** — the scripts only read the files that `src/tacode.py` has already written.
Keeping them out of `src/` means the computation has no dependency on a plotting library.

Each tool under `src_helper/` lives in its own directory, the same way each module under
`src/` does.

| File | What it is |
|---|---|
| `animate_trajectory.py` | animates a trajectory, with the attitude if the file has it |
| `tecplot_reader.py` | reads the Tecplot output into named columns |
| `vehicle_shape.py` | the vehicle shapes that are drawn (display only, no aerodynamics) |

## Requirements

`matplotlib` in addition to what Tacode itself needs. It is optional for the solver, so
install it separately:

```console
python3 -m pip install matplotlib
```

The output format follows the extension of `-o`:

| Extension | Writer | Needs |
|---|---|---|
| `.gif` | `PillowWriter` | nothing else |
| `.html`, `.htm` | `HTMLWriter` | nothing else |
| `.mp4` | `FFMpegWriter` | `ffmpeg` on the PATH |

`.html` is a single self-contained page: the frames are embedded as base64 PNG, and the
page has play, pause and a frame slider. It needs no `ffmpeg`, and unlike a GIF it does
not quantise the colours, so it is the format to reach for when `ffmpeg` is missing. Only
the play and pause icons come from a CDN, so they are blank when the page is opened
offline; the buttons themselves still work.

`--snapshot` writes a still instead of an animation, in whatever format matplotlib infers
from the extension (`.png`, `.pdf`, `.eps`).

## Animating a trajectory

```console
cd tutorial/work_reentry_6dof
./run_tacode.sh
python3 ../../src_helper/animate_trajectory/animate_trajectory.py output_result/tecplot.dat -o attitude.gif
```

The window is split into four:

| Panel | What it shows |
|---|---|
| left | the ground under the vehicle, the trajectory, and the vehicle drawn far larger than life so its attitude is visible; a globe inset gives the position on the planet |
| top right | the vehicle in the local horizon frame, with the velocity vector (dashed) and the body axes (x red, y green, z blue) — the angle between the body x axis and the dashed line is the angle of attack |
| middle right | altitude against time, with a cursor at the current step |
| bottom right | angle of attack and sideslip against time (velocity, for a 3-DOF file) |

A file from a three-degree-of-freedom run has no attitude columns; the script says so and
draws the trajectory alone, without orienting the vehicle.

### Options

| Option | Meaning |
|---|---|
| `-o, --output` | output file. `.gif`, `.html`, or `.mp4` when ffmpeg is installed |
| `-s, --shape` | `capsule` (sphere-cone, the default), `satellite` (bus with solar panels), `aircraft` (the most legible attitude) |
| `--frames` | number of frames; the trajectory is subsampled to this (default 150) |
| `--fps` | frames per second of the animation (default 20) |
| `--scale` | how large the vehicle is drawn on the ground view, in km (default 500) |
| `--window` | half width of the ground view, in km (default 1800). `0` fits the whole trajectory instead of following the vehicle |
| `--fixed-view` | keep the camera still instead of turning it with the vehicle |
| `--elevation` | camera elevation of the ground view, deg |
| `--dpi` | resolution. Lower it for a smaller file |
| `--embed-limit` | how many MB of frames an `.html` output may embed (default 512). Beyond it the animation is silently cut short, so the script says when the limit was reached |
| `--snapshot` | write a single frame at this time (s) as an image, instead of an animation |

### Examples

```console
# a light animation
python3 src_helper/animate_trajectory/animate_trajectory.py tecplot.dat -o attitude.gif --frames 100 --dpi 80

# one frame for a report, at t = 300 s
python3 src_helper/animate_trajectory/animate_trajectory.py tecplot.dat -o attitude.png --snapshot 300 --shape aircraft

# the whole trajectory in view, no following
python3 src_helper/animate_trajectory/animate_trajectory.py tecplot.dat -o orbit.gif --window 0 --fixed-view

# a page to open in a browser, no ffmpeg needed
python3 src_helper/animate_trajectory/animate_trajectory.py tecplot.dat -o attitude.html --dpi 80

# the smallest file, if ffmpeg is installed
python3 src_helper/animate_trajectory/animate_trajectory.py tecplot.dat -o attitude.mp4 --dpi 80
```

### File size

Every frame is a full picture in a `.gif` and in an `.html`; only the `.mp4` compresses
across frames, and since the animation is mostly a still background it compresses very
well. The `.html` is the largest of the three because it embeds the frames as PNG, which
keeps every colour, and base64 adds a third on top; the GIF is smaller only because it
quantises to 256 colours. The 6-DOF tutorial at 80 frames and `--dpi 70` comes to:

| Format | Size |
|---|---|
| `.mp4` | 1.1 MB |
| `.gif` | 4.6 MB |
| `.html` | 16.1 MB |

Write an `.mp4` if ffmpeg is available, or lower `--dpi` and `--frames` — a GIF runs about
70 kB a frame at `--dpi 80`, so 150 frames come to roughly 10 MB. The tutorial cases do not
keep an animation in the repository for that reason; generate one when you need it.

The script prints the size of what it wrote, so a file that came out larger than expected
is visible without looking it up.

## Vehicle shapes

The shapes exist so that the attitude can be read at a glance. **They have nothing to do
with the aerodynamics** — the coefficients come from `database/aerodynamic`, whatever
shape is drawn. `capsule` matches the geometry the supplied sphere-cone table was
generated for; `aircraft` is the easiest to read because its nose, wings and fin
distinguish all three axes at once.

To add a shape, write a function in `vehicle_shape.py` that returns a list of
`(vertices, colour)` in body axes — [forward, right, down], sized to about 1 — and add it
to `LIST_KIND`.

## Conventions

The scripts import `attitude.attitude` from `src/` for the quaternion and local-horizon
transformations, rather than repeating them, so the drawing always follows the same
convention as the solver. The quaternion columns are the transformation from ECEF to body
axes; see the Coordinate Transformation section of the top-level `README.md`.
