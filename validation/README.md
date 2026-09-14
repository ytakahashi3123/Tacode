# validation

Cases that put Tacode against **flight data**, as opposed to `tutorial/`, which demonstrates
how to run it, and `test/`, which checks the code against exact solutions, identities and its
own reference outputs.

| Case | Flight | What it compares |
|---|---|---|
| `apollo4/` | Apollo 4 (AS-501), 1967-11-09 | the entry trajectory at the measured bank angle, the atmosphere table against the flight-derived density, and stagnation-point heating correlations against the flight heating |
| `apollo10/` | Apollo 10 (AS-505), 1969-05-26 | the entry trajectory at the measured bank angle: altitude to 3.2 km rms and inertial velocity to 189 m/s rms over the whole entry |
| `fire2/` | Project Fire flight II, 1965-05-22 | an **uncontrolled, axisymmetric** entry at 11.35 km/s, where nothing about the attitude has to be assumed: the trajectory against NASA's own 0.5 s numerical table (48 m rms in altitude at the same drag), the atmosphere table against the sounding-rocket density, and **the six-degree-of-freedom oscillation period against the measured rate gyros, to 5 %** |

The two Apollo cases fly **the bank angle measured in flight**, so neither assumes the
direction of the lift. They are complementary: Apollo 4 also has the density and the
heating, but its trajectory skipped close to the skip-out boundary and is therefore
hypersensitive to the lift coefficient; Apollo 10 has only altitude and velocity to compare
against, but its trajectory is well conditioned, which makes it the sharper test of the
solver.

Fire II is of a different kind again. Its vehicle was uncontrolled and axisymmetric, so
there is no lift direction to supply at all, and its trajectory reference is a **numerical
table at 0.5 s intervals** rather than a digitised figure. That table is NASA's own
three-degree-of-freedom simulation of the same entry rather than a measurement, which is
stated wherever its numbers are used: matching it checks the solver. The case's measured
comparisons are the atmospheric density and, above all, **the attitude** — it is the only
case here that exercises the six-degree-of-freedom side against flight data.

Each case directory holds everything it needs: `config.yml`, its own `database/`, the flight
data converted to SI in `reference/`, the script that converted it, and the script that draws
the comparison. The run outputs (`output_result/`, `output_restart/`, `output_comparison/`)
are ignored by git; the numbers measured from them are written down in the README of the
case, with the date and the interpreter they were measured with.

**These cases are not part of `run_tests.sh`.** They need matplotlib for the figures, they
are compared against measurements rather than against an exact answer, and their value is in
the write-up rather than in a pass or fail. Re-run one when the physics it covers is touched
(the atmosphere table and its interpolation, the drag, the time integration) and check the
numbers in its README still hold.

Raw sources — scanned reports, digitised figures, notes — live outside the repository in
`../../references/`, which is deliberately not tracked. A case must run without it.
