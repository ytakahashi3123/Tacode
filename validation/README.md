# validation

Cases that put Tacode against **flight data**, as opposed to `tutorial/`, which demonstrates
how to run it, and `test/`, which checks the code against exact solutions, identities and its
own reference outputs.

| Case | Flight | What it compares |
|---|---|---|
| `apollo4/` | Apollo 4 (AS-501), 1967-11-09 | the entry trajectory at the measured bank angle, the atmosphere table against the flight-derived density, and stagnation-point heating correlations against the flight heating |
| `apollo10/` | Apollo 10 (AS-505), 1969-05-26 | the entry trajectory at the measured bank angle: altitude to 3.2 km rms and inertial velocity to 189 m/s rms over the whole entry |

Both cases fly **the bank angle measured in flight**, so neither assumes the direction of
the lift. They are complementary: Apollo 4 also has the density and the heating, but its
trajectory skipped close to the skip-out boundary and is therefore hypersensitive to the
lift coefficient; Apollo 10 has only altitude and velocity to compare against, but its
trajectory is well conditioned, which makes it the sharper test of the solver.

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
