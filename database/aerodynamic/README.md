# Aerodynamic tables

`satellite.filename_aerodynamic` in `config.yml` selects the file. The 3-DOF
computation only reads the drag coefficient (`CFx` at the smallest angle of attack in
the table); the 6-DOF computation reads the whole set of force and moment coefficients
and interpolates them in the angle of attack as well.

| File | Angles of attack | Knudsen range | Source |
|---|---|---|---|
| `aerodynamic.txt` | 0° only | 0.203 – 4.35e4 | DSMC + CFD of the EGG re-entry capsule |
| `aerodynamic_spherecone_aoa.txt` | 0–180°, every 10° | 1e-4 – 1e5 | analytic sphere-cone model, see below |
| `aerodynamic_apollo_aoa.txt` | 0–180°, every 5° | 1e-4 – 1e5 | the same model over the Apollo command module, see below |

`aerodynamic.txt` is the default used by `tutorial/work`, `tutorial/work_reentry` and
`tutorial/template`, and it is what the committed reference outputs were produced with.
It stays unchanged.

This directory is the master copy. Each case directory keeps the table it actually uses
in `<case>/database/aerodynamic` and reads it from there — `aerodynamic.txt` for the
3-DOF cases, `aerodynamic_spherecone_aoa.txt` for `tutorial/work_reentry_6dof`. Copy a
table into that directory before naming it in `satellite.filename_aerodynamic`.

**A table with a single angle of attack gives no restoring moment.** The coefficients
are then the same whatever the attitude, so a 6-DOF run with `aerodynamic.txt` will
tumble rather than oscillate about trim. The code prints a warning when this happens.

## File format

```
<any header line that is not a row of numbers>
variables = Kn, CFx, CFy, CFz, CMx, CMy, CMz, SDV_CFx ... SDV_CMz, Altitude
AOA 0
<14 numbers per row, one row per Knudsen number>
AOA 10
<the same Knudsen numbers, in the same order>
...
```

* A line whose first word is `AOA` starts a block and gives the angle of attack in
  degrees. A file with no `AOA` line at all is read as a single block at 0°, which is
  how `aerodynamic.txt` (written before the 6-DOF extension) is still read unchanged.
* A data row is a line of exactly 14 numbers. Anything else is skipped, so headers,
  blank lines and trailing junk are harmless.
* Every block must list the same Knudsen numbers, in ascending order — the blocks form
  the (angle of attack, Knudsen number) grid that is interpolated. The code stops with
  an error otherwise.
* Outside the table the values are clamped to the edge, in both the angle of attack and
  the Knudsen number.

## Sign convention

The body axes are [forward, right, down], and the table is given for a velocity in the
body x–z plane, that is, the crossflow on the +z side. Tacode applies

    F_body = -q S [CFx, CFy, CFz]        M_body = q S L [CMx, CMy, CMz]

so `CFx` at 0° is the ordinary drag coefficient, and a statically stable body has
`CMy < 0` for a positive angle of attack. For a general attitude the coefficient
vectors are rotated about the body x axis onto the actual crossflow plane, which is
exact for a body of revolution.

**The 6-DOF computation therefore assumes an axisymmetric body.** The table is looked
up with the total angle of attack alone, so it carries no sideslip dependence, and an
axisymmetric body has `CFy = CMx = CMz = 0` over the whole table: no side force, no
rolling and no yawing moment inside the crossflow plane. If those columns are not zero,
the rotation moves them as if they lay in that plane, and the side force and the rolling
and yawing moments end up pointing in the wrong direction — the table does not hold
enough information to recover the right one.

Tacode checks this when the table is read for a 6-DOF run and prints a warning naming
the offending columns; it does not stop, since the in-plane coefficients are still
usable. A column counts as non-zero when it exceeds 1e-3 of the in-plane coefficients
it is compared against (`CFy` against `max(|CFx|, |CFz|)`, `CMx` and `CMz` against
`max(|CMy|)`), which passes the round-off of an analytically generated table. The
warning is only raised for `kind_aerodynamic_model: fileread` with
`attitude.flag_attitude: True`; a 3-DOF run reads `CFx` alone and is unaffected.

`aerodynamic.txt` triggers it: `CMx` and `CMz` reach 2% and 7% of `CMy`, which is
measurement scatter of a nominally axisymmetric capsule rather than a real asymmetry.
It is harmless in the 3-DOF runs the table is meant for.

The moments are about the reference point the table was generated for. If the centre of
gravity is elsewhere, give the offset in `attitude.center_of_gravity` and the code moves
the moments for you.

## `aerodynamic_spherecone_aoa.txt`

Produced by `generate_aerodynamic_table.py`, which is committed next to it:

```
python3 generate_aerodynamic_table.py -o aerodynamic_spherecone_aoa.txt
```

**This is analytic sample data, not a measurement and not DSMC.** It exists so that the
6-DOF tutorial has an attitude-dependent table to run against. The geometry is a
sphere-cone of 1 m base diameter, 45° half angle, 0.25 m nose radius, with the moment
reference point on the axis 0.20 m ahead of the base; the reference area and length are
π/4 m² and 0.5 m, matching `tutorial/work_reentry_6dof/config.yml`.

The coefficients come from a panel integration of two limits, blended by Knudsen number:

| Regime | Pressure | Shear |
|---|---|---|
| continuum (Kn ≤ 1e-3) | modified Newtonian, `Cp = 1.839 sin²δ` | none |
| free molecular (Kn ≥ 10) | `Cp = 2 sin²δ` | `Cτ = 2 sinδ cosδ` |

with `f = sin²(π/2 · (log₁₀Kn + 3)/4)` bridging between them, δ being the local
inclination of the surface to the flow. Shadowing is decided by `sinδ > 0`, which is
exact for this convex shape.

Known limitations of the model, which matter if you compare against real data:

* The free-molecular limit is the hyperthermal one with complete momentum transfer and
  no re-emission, so the resulting force there is pure drag: the normal force is zero
  and only the moment about the centre of gravity survives. A real diffuse-reflecting
  surface produces some normal force.
* Modified Newtonian carries no viscous drag and no base pressure, so the continuum
  drag is the pressure contribution only.
* The dynamic damping derivatives are not part of this model at all. They are given
  separately in `config.yml` (`attitude.damping_coefficient`).

## `aerodynamic_apollo_aoa.txt`

The same panel model over the profile of the Apollo command module, written for the
validation case in `validation/apollo4`:

```
python3 generate_aerodynamic_table.py --shape capsule --angle-step 5 -o aerodynamic_apollo_aoa.txt
```

The geometry is the published one: 154 in. (3.9116 m) maximum diameter, 4.6939 m
spherical heat shield, 0.1956 m corner radius and a 32.5° afterbody half-angle. **The
vehicle flies heat-shield forward**, so the body x axis points out of the heat shield and
the angle of attack is measured from it. The moment reference point is the apex of the
heat shield, **on the axis**, and the reference area and length are π/4 · 3.9116² m² and
3.9116 m.

The table is therefore axisymmetric and has no trim of its own. The command module trims
because its centre of gravity is off the axis; that offset lives in the configuration
(`attitude.center_of_gravity`), and the solver moves the moment to it. In
`validation/apollo4/config_6dof.yml` the offset is the **measured** 0.16688 m (6.57 in.,
NASA TN D-5399 table II), which trims this table at 25.65° against the 24.4° of the
flight report. Solving the offset from the table instead would give 0.1583 m; the
measured value is used so that nothing in the case is fitted to the answer.

How good the model is, measured against NASA TN D-6725 table I (trim values at M > 29.5):

| | this table at 24.4° | flight report |
|---|---|---|
| `CD` | 1.266 | 1.2891 |
| `L/D` | 0.373 | 0.301 |

The drag is within 2 %; the lift-to-drag ratio is 24 % high, which is what modified
Newtonian usually does on a blunt body — it puts no base pressure and no viscous force on
the afterbody. The consequence for the trajectory is worked out in
`validation/apollo4/README.md`.
