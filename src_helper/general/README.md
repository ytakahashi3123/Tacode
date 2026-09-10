# src_helper/general

Shared code for the post-processing tools, in the same way `src/general/` holds what the
solver's modules share. **This is not a tool** — there is nothing here to run.

| File | What it is |
|---|---|
| `helper_config.py` | reads the settings of a tool from `config_helper.yml`, merges them with the command line, and writes them back out for `--save-config` |

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
