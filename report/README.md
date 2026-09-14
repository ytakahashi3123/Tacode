# report

The Tacode version 2 report: theory, verification and validation.

It is the successor to `../../manuscript/manuscript_ver1.2_en`, which described the
three-degree-of-freedom formulation of version 1 and the operation of the code. This one
describes the version 2 formulation (six degrees of freedom, lift and bank angle, wind,
absolute time, divergence detection) and, at greater length, what has been done to
establish that the code is correct.

| File | |
|---|---|
| `tacode_v2_en.tex` | **The authoritative version.** English, builds with `pdflatex` |
| `tacode_v2_ja.tex` | Japanese, builds with `platex` + `dvipdfmx`, as the version 1 manuscript did |
| `figure/` | Every figure, as PDF |
| `generate_figure.py` | Rebuilds the figures from the data in this repository |
| `check_tex.py` | Checks the sources without typesetting them |
| `Makefile` | `make`, `make en`, `make ja`, `make figure`, `make check`, `make clean` |

```console
make            # both PDFs
make check      # structure of the sources, no TeX needed
make figure     # rebuild the figures (needs matplotlib; some of them run the solver)
```

## The figures come from the repository

Every figure except the four coordinate-system sketches is drawn by
`generate_figure.py` from data that lives in this repository: the reference outputs
committed under `tutorial/`, the reference data and run outputs of `validation/`, and,
for the convergence study, runs of the solver at a sequence of time steps. Nothing is
drawn by hand and nothing is imported from outside.

The four sketches (`cartesian_coordinate`, `geodetic_coordinate`, `polar_coordinate`,
`potential_schematic`) are copies of the version 1 manuscript's figures, converted from
EPS to PDF.

A figure whose inputs are missing is **not drawn**; the script says which files are
missing and which command produces them, and carries on. So after a fresh checkout:

```console
cd ../validation/fire2      && ./run_tacode.sh && ./run_tacode.sh -file config_matched.yml \
                                              && ./run_tacode.sh -file config_6dof.yml
cd ../validation/apollo4    && ./run_tacode.sh -file config_lift.yml
cd ../validation/apollo10   && ./run_tacode.sh
cd ../../report && make figure
```

`generate_figure.py --skip-run` leaves out the convergence study, which is the only
figure that runs the solver itself (about 20 s).

## `check_tex.py` is not a compiler

The machine this report was written on has no TeX system, so the sources could not be
typeset while they were being written. `check_tex.py` was written to catch what can be
caught without one: unbalanced environments, braces and `$`, a `\ref` with no `\label`, a
`\cite` with no `\bibitem`, a duplicated `\label`, an `\includegraphics` pointing at a
file that is not there, and Markdown notation that has leaked into the LaTeX (it caught a
`**bold**` in a caption).

**It is not a substitute for typesetting.** Passing it does not mean the document
compiles. The first build on a machine with TeX should be treated as the real check.

## Keeping the numbers honest

Every measured number in the report is also written down in the README of the case it
came from, with the date and the interpreter it was measured with:

- `../validation/apollo4/README.md`
- `../validation/apollo10/README.md`
- `../validation/fire2/README.md`

If a case is re-run and a number moves, both places have to change. The report is not the
primary record; those READMEs are.
