# Aerodynamic tables

`satellite.filename_aerodynamic` in `config.yml` selects the file. The 3-DOF
computation only reads the drag coefficient (`CFx` at the smallest angle of attack in
the table); the 6-DOF computation reads the whole set of force and moment coefficients
and interpolates them in the angle of attack as well.

| File | Angles of attack | Knudsen range | Source |
|---|---|---|---|
| `aerodynamic.txt` | 0° only | 0.203 – 4.35e4 | DSMC + CFD of the EGG re-entry capsule |
| `aerodynamic_spherecone_aoa.txt` | 0–180°, every 10° | 1e-4 – 1e5 | analytic sphere-cone model, see below |

`aerodynamic.txt` is the default used by `tutorial/work`, `tutorial/work_reentry` and
`tutorial_template`, and it is what the committed reference outputs were produced with.
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
