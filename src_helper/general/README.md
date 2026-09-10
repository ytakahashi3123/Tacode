# src_helper/general

Shared code for the post-processing tools, in the same way `src/general/` holds what the
solver's modules share. **This is not a tool** — there is nothing here to run.

| File | What it is |
|---|---|
| `helper_config.py` | reads the settings of a tool from `config_helper.yml`, merges them with the command line, and writes them back out for `--save-config` |
| `coastline.py` | the coastline of the 3-D views: reads the shipped text, cuts it to a window, turns it into segments in ECEF, and drops the far side of the globe |
| `coastline_110m.txt` | the coastline itself, 134 polylines and 5128 points (86 kB) |
| `generate_coastline.py` | writes that text from the Natural Earth shapefile. **Run once**; it needs the network, the tools do not |

## The settings file

The tools read **`config_helper.yml` in the current directory** by default, exactly as
the solver reads `config.yml`, and `-file` points at another one. A single file carries
one section per tool:

```yaml
animate_trajectory:
  filename: output_result/tecplot.dat
  output: attitude.mp4

montecarlo_dispersion:
  directory: work_montecarlo_wind
  reference: ../work_reentry/output_result/tecplot.dat

montecarlo_animation:
  directory: work_montecarlo_wind
  view: 3d
  spin: 0.0
```

- The keys are the long options without their dashes, with `-` written as `_`
  (`--altitude-max` becomes `altitude_max`). A repeated option such as `--mark` is a list.
- The order of precedence is **command line > settings file > default**. Whether an
  option was given on the command line is decided by it differing from the default, so
  writing the default value out explicitly simply leaves the file in charge.
- An **unknown key stops the run**, listing what is available. A key that is ignored
  because of a typo is the worst outcome: the setting looks made and is not.
- A missing file is not an error, so a tool still works with the command line alone.
- Paths are relative to the current directory, as everywhere else in Tacode.
- `--save-config` writes the settings in effect into the file (other sections are kept,
  **comments are not**), and exits. It is the quickest way to start a file, and a way to
  record exactly what produced a figure.

## Adding it to a tool

```python
import helper_config as helper_config

NAME_SECTION = 'my_tool'

def argument():
  parser = argparse.ArgumentParser(...)
  parser.add_argument('directory', nargs='?', default=None, help=...)
  ...
  helper_config.add_argument(parser)
  return helper_config.get_setting(parser, NAME_SECTION)

def main():
  args = argument()
  if helper_config.save_file(args, NAME_SECTION) :
    return
  helper_config.require(args, 'directory', NAME_SECTION)
```

A positional argument has to be declared with `nargs='?'` and a default of `None`, since
it may come from the file instead; `require` then reports it in terms of both places.

## The coastline

`coastline_110m.txt` is Natural Earth's 1:110 m physical coastline (public domain, no
attribution required), reduced to plain text: one polyline per block, `longitude latitude`
in degrees, blocks separated by a blank line. The reader needs **numpy and nothing else**.

That is the whole point of shipping it. Drawing a real map with `cartopy` would pull in
GEOS and PROJ, and the machine that makes the figures is often the machine that runs the
solver, whose dependencies are deliberately four packages. A coastline as text costs
86 kB and works offline.

```python
import coastline as coastline

segment = coastline.read_coastline()                        # 折れ線（経度・緯度, deg）
segment = coastline.select_range(segment, (135.0, 145.0), (30.0, 40.0))
line    = coastline.get_segment_ecef(segment)               # (M, 2, 3), Line3DCollection 用
line    = coastline.select_visible(line, coastline.get_view_direction(20.0, -60.0))
```

Three things are worth knowing:

- `select_range` **cuts a polyline where it leaves the window**, and extends each piece by
  one point on either side. Keeping "every point inside plus its neighbours" instead joins
  a coast that only passes through, and draws a straight line across the window.
- A window may cross the 180th meridian, since it is built as a centre plus a margin, so
  the longitudes are folded by ±360° onto the window rather than onto `[-180, 180)`.
- `Line3DCollection` only takes segments that all hold the same number of points, so
  `get_segment_ecef` splits the polylines into pairs. And matplotlib's 3-D does not hide
  lines behind a surface: without `select_visible` the coastline of the far hemisphere
  shows through an opaque Earth and the continents appear mirrored on top of each other.

To rebuild the file:

```console
python3 src_helper/general/generate_coastline.py
```

It downloads the zipped shapefile, reads the polylines itself (the `.shp` format is simple
enough that no library is needed for it) and writes the text. `--zip` reads an archive
already on disk, and `--url` takes another Natural Earth file — the 1:50 m coastline is
the same format, some ten times as large.
