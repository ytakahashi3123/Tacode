# Tacode
Trajectory analysis code

[![tests](https://github.com/ytakahashi3123/Tacode/actions/workflows/tests.yml/badge.svg)](https://github.com/ytakahashi3123/Tacode/actions/workflows/tests.yml)


# Code description

`Tacode` solves the equation of motion of an object on an Earth-centered, Earth-fixed (ECEF) noninertial frame by a python script.
The planetary gravity, Coriolis, centrifugal, and aerodynamic forces act on the object.
The gravity force is obtained by differentiating the gravitational potential considering J20, J22, J30, and J40.
The atmospheric data is given by NRLMSISE-00 Atmosphere Model.
The equation of motion is numerically solved using fourth-order Runge-Kutta method in four stages.

By default the object is a point mass with **three degrees of freedom**, and the aerodynamic force is the drag along the direction of motion, given by atmospheric density, drag coefficient, characteristic (projection) area, and velocity.

Since v2.3.0 the attitude can be solved as well, giving **six degrees of freedom**. The rigid body then rotates under the aerodynamic, damping and gravity-gradient moments, and the aerodynamic force follows the attitude through the angle of attack, so that the trajectory and the attitude are coupled in both directions. The attitude is carried by a quaternion, and the aerodynamic coefficients are interpolated from a table in the angle of attack as well as in the Knudsen number.
The attitude is switched off unless `attitude.flag_attitude` is set, and a configuration without an `attitude` section reproduces the earlier three-degree-of-freedom results bit for bit.

![Atmospheric-entry trajectories for initial velocities of 7250, 7450, and 7650 m/s](figure/trajectory.jpg)

The trajectories above were written to `geodetic.kml` by the KML output and rendered
in Google Earth. They are the `tutorial/work_reentry` case (7450 m/s) together with
two runs of the same configuration at 7250 and 7650 m/s.

## Governing equation

`Tacode` advances a state vector in the Earth-Centered Earth-Fixed (ECEF) rotating frame.
In three degrees of freedom the state is the position and the velocity of a point mass.
Setting `attitude.flag_attitude: True` adds the attitude quaternion and the angular
velocity of a rigid body, giving six degrees of freedom:

```math
	{\boldsymbol y} = \left[ {\boldsymbol x}, \; \dot{\boldsymbol x} \right] \quad \textrm{(3-DOF)} ,
	\qquad
	{\boldsymbol y} = \left[ {\boldsymbol x}, \; \dot{\boldsymbol x}, \; {\boldsymbol q}, \; {\boldsymbol \Omega} \right] \quad \textrm{(6-DOF)} .
```

Both sets are advanced by the same time integrator (four-stage fourth-order Runge-Kutta,
or explicit Euler) with the same timestep. The translational and the rotational equations
are coupled in both directions: the moments turn the body, and the aerodynamic force
follows the attitude. Of the four forces below, only the aerodynamic one depends on the
attitude, so with the attitude switched off the translational set reduces exactly to the
three-degree-of-freedom equations that earlier versions solved.

### Translational motion

The equations of motion for a mass point in the non-inertial ECEF frame are

```math
	m \frac { { \partial  }^{ 2 } {\boldsymbol x} }{ \partial { t }^{ 2 } } = 
	{\boldsymbol F}_{\rm grav} 
	- 2m{ \boldsymbol \omega }\times \frac { { \partial  } {\boldsymbol x} }{ \partial { t } } 
	- m {\boldsymbol  \omega  }\times \left( { \boldsymbol \omega  }\times { \boldsymbol x } \right) 
	+ {\boldsymbol F}_{\rm aero} ,
```

where, $t$ is time, ${\boldsymbol x}$ and $m$ are the position vector of the mass point and its mass, respectively, and ${\boldsymbol \omega}$ denotes the angular velocity vector.
For the Earth where the coordinate system rotates around the $z$-axis, ${\boldsymbol \omega}=(0, 0, 7.292115\times10^{-5})$ rad/s.
Here, $`{\boldsymbol F}_{\rm grav}`$ is the force of gravity and $`{\boldsymbol F}_{\rm aero}`$ is the aerodynamic force.
The second and third terms on the right-hand side are the Coriolis and the centrifugal forces that appear because the frame rotates.

### Rotational motion

The rotation of the rigid body about its centre of gravity is governed by Euler's equation,
written in the body frame where the inertia tensor is constant:

```math
	{\boldsymbol I} \frac{ \partial {\boldsymbol \Omega} }{ \partial t }
	+ {\boldsymbol \Omega} \times \left( {\boldsymbol I} {\boldsymbol \Omega} \right)
	= {\boldsymbol M}_{\rm aero} + {\boldsymbol M}_{\rm damp} + {\boldsymbol M}_{\rm grav} .
```

Here ${\boldsymbol I}$ is the inertia tensor about the centre of gravity and
${\boldsymbol \Omega}$ is the angular velocity **with respect to the inertial frame**,
expressed in body components; $`{\boldsymbol M}_{\rm aero}`$, $`{\boldsymbol M}_{\rm damp}`$ and
$`{\boldsymbol M}_{\rm grav}`$ are the aerodynamic, the dynamic damping and the
gravity-gradient moments, each of which can be switched off individually. The three
equations are written out in [Aerodynamic Moment](#aerodynamic-moment) and
[Gravity-gradient Moment](#gravity-gradient-moment) below.

### Gravity Force

Because the Earth is a spheroid shape where the radii at the equator and the poles are different due to its rotation, it is convienient to describe the potential with spherical harmonics in polar coordinates [@CWagner1966], which can be represented as a summation of $J$ terms as follows:

```math
	U = - \frac{G M}{r} 
	\left[ 1
	- \left( \frac{a_{e}}{r} \right)^2 J_{20} \frac{ 3 \sin^2 \beta - 1}{2} 
	- \left( \frac{a_{e}}{r} \right)^2 J_{22} 3 \cos^2 \beta \cos 2 \left( \alpha + \alpha_{22} \right) 
	- \left( \frac{a_{e}}{r} \right)^3 J_{30} \frac{ 5 \sin^3 \beta - 3 \sin \beta}{2}  
	- \left( \frac{a_{e}}{r} \right)^4 J_{40} \frac{1}{8} \left( 35 \sin^4 \beta - 30 \sin^2 \beta + 3 \right) \right] .
```

Here, $G$ is the universal gravitational constant, $M$ is the planetary mass, and $r$ represents the radial distance. 
For Earth, $`G=6.67408\times10^{-11}`$ $`{\rm m}^3`$/{kg$`\cdot`$s$`^2`$} and $M=5.9722\times10^{24}$ kg are given.
In addition, $a_e$ represents the planet's equatorial radius, while $J_{20}$, $J_{22}$, $J_{30}$, $J_{40}$, and $\alpha_{22}$ are the dynamic form factors associated with the planet's shape given by the reference [@liu2019guidance]. 
By taking partial derivatives of this potential with respect to each coordinate components, the gravity in the polar coordinate system, denoted as $`\boldsymbol{F}_{grav}^{'}=(F_{g_r}^{'}, F_{g_{\alpha}}^{'}, F_{g_{\beta}}^{'})`$, are expressed as follows.

$$
\begin{aligned}
	F_{g_r}^{'}
	=& \frac{GMm}{r^2} 
	\left[ 
	- 1
	+ \frac{3}{2} \left( \frac{a_{e}}{r} \right)^2 J_{20} \left( 3 \sin^2 \beta - 1 \right)
	+ 9 \left( \frac{a_{e}}{r} \right)^2 J_{22} \left( \cos^2 \beta \right) \cos 2 \left( \alpha + \alpha_{22} \right)
	+ 2 \left( \frac{a_{e}}{r} \right)^3 J_{30} \left( 5 \sin^3 \beta - 3 \sin \beta \right)
	+ \frac{5}{8} \left( \frac{a_{e}}{r} \right)^4 J_{40} \left( 35 \sin^4 \beta -30 \sin^2 \beta + 3  \right)
	\right] , \\
	F_{g_{\alpha}}^{'} 
	=& \frac{GMm}{r^2} 
	\left[ 
	6 \left( \frac{a_{e}}{r} \right)^2 J_{22}  \cos \beta \sin 2 \left( \alpha + \alpha_{22} \right)
	\right] , \\
	F_{g_{\beta}}^{'} 
	=& \frac{GMm}{r^2} 
	\left[ 
	- 3 \left( \frac{a_{e}}{r} \right)^2 J_{20} \left( \sin \beta \cos \beta \right)
	+ 6 \left( \frac{a_{e}}{r} \right)^2 J_{22} \left( \sin \beta \cos \beta \right) \cos 2 \left( \alpha + \alpha_{22} \right)
	-  \frac{3}{2} \left( \frac{a_{e}}{r} \right)^3 J_{30} \left( 5 \sin^2 \beta - 1 \right) \cos \beta
	- \frac{5}{2} \left( \frac{a_{e}}{r} \right)^4 J_{40} \left( 7 \sin^2 \beta - 3 \right) \sin \beta \cos \beta
	\right] .
\end{aligned}
$$

By transformation from polar coordinates to the Cartesian coordinate system, the gravitational force $\boldsymbol{F}_{grav}$ in the Cartesian coordinate system can be obtained.

### Aerodynamic Force

The atmosphere is assumed to co-rotate with the Earth, so the ECEF velocity
${\boldsymbol v} = \partial {\boldsymbol x} / \partial t$ is also the velocity relative to
the air, and the dynamic pressure is

```math
	q_{\infty} = \frac{1}{2} \rho \left| {\boldsymbol v} \right|^2 .
```

The atmospheric density $\rho$ is interpolated from the NRLMSISE-00 table, and it is
multiplied by `initial_settings.density_factor`, which is what the Monte-Carlo driver
disperses. Winds are not modelled.

**In three degrees of freedom** the attitude is unknown, so the aerodynamic force can only
be a drag along the direction of motion:

```math
	{\boldsymbol F}_{\rm aero} = - q_{\infty} C_D S \frac{ \boldsymbol{v} }{\left| \boldsymbol{v} \right|} ,
```

where $C_D$ is the drag coefficient and $S$ is the characteristic (projected) area. The
drag coefficient is either a constant or interpolated from the aerodynamic table against
the Knudsen number.

**In six degrees of freedom** the force is read from the aerodynamic table in body axes as
a function of the total angle of attack as well as of the Knudsen number. Writing
$`[u, v, w] = {\boldsymbol C}_{b/e} {\boldsymbol v}`$ for the velocity in body components,
the total angle of attack and the aerodynamic roll angle are

```math
	\alpha_t = \arccos \frac{u}{\left| {\boldsymbol v} \right|} , \qquad
	\varphi_a = \arctan \frac{v}{w} ,
```

$\alpha_t$ being the angle between the body $x$ axis and the velocity, and $\varphi_a$ the
angle between the plane the coefficients were tabulated in and the plane the crossflow
actually lies in. The table is entered at $\alpha_t$ and rotated about the body $x$ axis
onto that plane:

```math
	{\boldsymbol F}_{\rm aero} = - q_{\infty} S \, {\boldsymbol C}_{F} ,
	\qquad
	{\boldsymbol C}_{F} = {\boldsymbol R}_{x} \left( \varphi_a \right)
	\begin{bmatrix} C_{Fx} \\ C_{Fy} \\ C_{Fz} \end{bmatrix} \left( \alpha_t, Kn \right) ,
	\qquad
	{\boldsymbol R}_{x} \left( \varphi_a \right) = \begin{bmatrix}
	1 & 0 & 0 \\
	0 & \cos \varphi_a & \sin \varphi_a \\
	0 & -\sin \varphi_a & \cos \varphi_a \end{bmatrix} .
```

This rotation is exact for a body of revolution, which is what a table indexed by
$\alpha_t$ alone describes. The force is transformed back to the ECEF frame with
$`{\boldsymbol C}_{b/e}^{T}`$ before it enters the translational equation.

The sign is chosen so that a positive $C_{Fx}$ pushes the body along $-x$: at zero
incidence the table gives $`{\boldsymbol C}_F = [C_D, 0, 0]`$ and the force becomes
$-q_{\infty} C_D S \hat{\boldsymbol x}_b$, which is the three-degree-of-freedom drag above
whenever the nose is on the velocity vector. The two formulations therefore agree exactly
at $\alpha_t = 0$, and this is checked by the test suite.

### Attitude kinematics

The attitude is carried by a quaternion ${\boldsymbol q} = [q_0, q_1, q_2, q_3]$ (scalar
first) that represents the transformation from ECEF to body axes,
$`{\boldsymbol z}_b = {\boldsymbol C}_{b/e} {\boldsymbol z}_e`$, with

```math
	{\boldsymbol C}_{b/e} \left( {\boldsymbol q} \right) = \begin{bmatrix}
	q_0^2 + q_1^2 - q_2^2 - q_3^2 & 2 \left( q_1 q_2 + q_0 q_3 \right) & 2 \left( q_1 q_3 - q_0 q_2 \right) \\
	2 \left( q_1 q_2 - q_0 q_3 \right) & q_0^2 - q_1^2 + q_2^2 - q_3^2 & 2 \left( q_2 q_3 + q_0 q_1 \right) \\
	2 \left( q_1 q_3 + q_0 q_2 \right) & 2 \left( q_2 q_3 - q_0 q_1 \right) & q_0^2 - q_1^2 - q_2^2 + q_3^2
	\end{bmatrix} .
```

A quaternion is used rather than the Euler angles themselves because it has no
singularity: a re-entry body passing through a pitch of 90 deg would hit the gimbal lock
of a 3-2-1 sequence. Its evolution is

```math
	\frac{ \partial {\boldsymbol q} }{ \partial t } = \frac{1}{2} {\boldsymbol Q} \left( {\boldsymbol \Omega}_{\rm rel} \right) {\boldsymbol q} ,
	\qquad
	{\boldsymbol Q} = \begin{bmatrix}
	0 & -p & -q & -r \\
	p & 0 & r & -q \\
	q & -r & 0 & p \\
	r & q & -p & 0 \end{bmatrix} .
```

**Two angular velocities appear, and they must not be confused.** Euler's equation is
written for the rate with respect to the inertial frame, ${\boldsymbol \Omega}$, whereas
the quaternion is referred to the rotating ECEF frame and therefore evolves with the
relative rate

```math
	{\boldsymbol \Omega}_{\rm rel} = \left[ p, q, r \right] = {\boldsymbol \Omega} - {\boldsymbol C}_{b/e} {\boldsymbol \omega} .
```

Leaving the Earth rate ${\boldsymbol \omega}$ in would drift the attitude by 0.0042 deg/s.
The same relative rate is the one the aerodynamic damping sees, because the atmosphere
co-rotates with the Earth. The state that is integrated is ${\boldsymbol \Omega}$; the
input in `config.yml` and the output in the Tecplot file are both
${\boldsymbol \Omega}_{\rm rel}$, which is the more intuitive of the two, so the code
converts at the boundaries.

The quaternion is renormalised at every stage of the integration to keep
$`\left| {\boldsymbol q} \right| = 1`$ and the transformation matrix orthonormal.

### Aerodynamic Moment

The moment coefficients are read from the same table, at the same $\alpha_t$, and rotated
by the same ${\boldsymbol R}_x(\varphi_a)$ — a proper rotation, so the pseudo-vector of the
moment transforms exactly as the force vector does. The coefficients are referred to the
moment reference point of the table, and are moved to the centre of gravity by the usual
transfer term:

```math
	{\boldsymbol M}_{\rm aero} = q_{\infty} S L \, {\boldsymbol C}_{M}
	+ \left( {\boldsymbol r}_{\rm ref} - {\boldsymbol r}_{\rm cg} \right) \times {\boldsymbol F}_{\rm aero} ,
	\qquad
	{\boldsymbol C}_{M} = {\boldsymbol R}_{x} \left( \varphi_a \right)
	\begin{bmatrix} C_{Mx} \\ C_{My} \\ C_{Mz} \end{bmatrix} \left( \alpha_t, Kn \right) ,
```

where $L$ is the characteristic length and
$`{\boldsymbol r}_{\rm cg} - {\boldsymbol r}_{\rm ref}`$ is `attitude.center_of_gravity`.
With the body axes taken as [forward, right, down], a positive $C_{My}$ pitches the nose
up, so a statically stable body has $C_{My} < 0$ at a positive angle of attack and the
moment restores it towards $\alpha_t = 0$.

The table gives the static coefficients only, that is, the moment of a body held at a
fixed attitude. A rotating body also feels a moment opposing its rotation, which is
modelled with the dynamic damping derivatives:

```math
	{\boldsymbol M}_{\rm damp} = q_{\infty} S L \frac{L}{2 \left| {\boldsymbol v} \right|}
	\begin{bmatrix} C_{lp} \, p \\ C_{mq} \, q \\ C_{nr} \, r \end{bmatrix} .
```

The derivatives are given in `config.yml` (`attitude.damping_coefficient`) and have to be
negative to damp. Without them the oscillation of a statically stable body never decays,
because the static coefficients alone are conservative.

If `satellite.kind_aerodynamic_model` is `constant`, there is no table to read: the axial
force is taken from `drag_coefficient`, the normal force is zero, and the pitching moment
is the linear model $C_m = C_{m\alpha} \alpha_t$ with $C_{m\alpha}$ given by
`attitude.static_stability_derivative`. This is the configuration whose oscillation period
can be written in closed form, and the test suite uses it for exactly that reason.

### Gravity-gradient Moment

Gravity pulls harder on the near side of the body than on the far side, and the resultant
of that difference about the centre of gravity is the gravity-gradient moment. Expanding
the attraction to first order in the body size over the orbital radius leaves

```math
	{\boldsymbol M}_{\rm grav} = \frac{3 G M}{r^3} \, \hat{\boldsymbol u} \times \left( {\boldsymbol I} \hat{\boldsymbol u} \right) ,
	\qquad
	\hat{\boldsymbol u} = {\boldsymbol C}_{b/e} \frac{ {\boldsymbol x} }{ \left| {\boldsymbol x} \right| } ,
```

with $\hat{\boldsymbol u}$ the unit vector along the line joining the vehicle and the planet
centre, in body components. It vanishes for a body whose inertia tensor is isotropic, and
whenever a principal axis is aligned with $\hat{\boldsymbol u}$; the sign of
$\hat{\boldsymbol u}$ does not matter, since it appears twice. For a small vehicle in low
Earth orbit the term is of the order of $10^{-6}$ N m, negligible against the aerodynamic
moment during an entry but the only moment left in vacuum.

## Coordinate Transformation

The equations of motion are solved in a Cartesian coordinate system. 
However, for initial conditions such as initial position and velocity, atmospheric density model, and output formats, it is convenient to use the geodetic coordinate system, which represents latitude, longitude, and altitude coordinates. 
In addition, the gravity forces are calculated through the use of Polar coordinates. 
`Tacode` includes the coordinate transformation modules among the Cartesian coordinate,  WGS84 geodetic [@WGS84], and polar coordinate systems.

The attitude computation adds two more frames. The following are used throughout:

| Frame | Definition | Used for |
|---|---|---|
| ECEF Cartesian | rotating with the planet, origin at its centre | the equations of motion; the quaternion is referred to it |
| WGS84 geodetic | longitude, latitude, altitude above the reference ellipsoid | initial position, atmosphere lookup, KML output |
| Polar (geocentric) | radius, geocentric latitude $\beta$, geocentric longitude $\alpha$ | the gravity potential; the frame of the initial velocity |
| Local horizon (NED) | geocentric [north, east, down] at the current position | the Euler angles |
| Body | [forward, right, down], fixed to the vehicle | inertia tensor, aerodynamic coefficients, angular velocity |

The local horizon is deliberately built on the *geocentric* latitude and longitude, so
that it is the same frame the initial velocity is already given in (see
[Reference frame of the initial velocity](#reference-frame-of-the-initial-velocity)). It is
obtained from the ECEF frame by

```math
	{\boldsymbol C}_{n/e} \left( \alpha, \beta \right) = \begin{bmatrix}
	-\sin \beta \cos \alpha & -\sin \beta \sin \alpha & \cos \beta \\
	-\sin \alpha & \cos \alpha & 0 \\
	-\cos \beta \cos \alpha & -\cos \beta \sin \alpha & -\sin \beta \end{bmatrix} ,
```

whose rows are the north, east and down directions. The attitude with respect to that
frame is given by the 3-2-1 Euler sequence yaw $\psi$ → pitch $\theta$ → roll $\phi$,

```math
	{\boldsymbol C}_{b/n} = {\boldsymbol R}_{x} \left( \phi \right) {\boldsymbol R}_{y} \left( \theta \right) {\boldsymbol R}_{z} \left( \psi \right) ,
	\qquad
	{\boldsymbol C}_{b/e} = {\boldsymbol C}_{b/n} {\boldsymbol C}_{n/e} ,
```

where $`{\boldsymbol R}_x`$, $`{\boldsymbol R}_y`$ and $`{\boldsymbol R}_z`$ are the
elementary rotations of the frame about each axis. Yaw is measured from the north towards
the east, so an eastward flight with the nose on the velocity vector has $\psi = 90$ deg;
pitch is positive nose-up and roll positive right-wing-down.

The Euler angles are used for the input and the output only. Internally the attitude is
the quaternion of ${\boldsymbol C}_{b/e}$, which is referred to a frame that does not move
with the vehicle. The Euler angles written to the output file are therefore rebuilt at
each output point, since the local horizon they refer to turns as the vehicle travels.

Two further angles are derived from the velocity in body components,
$`[u, v, w] = {\boldsymbol C}_{b/e} {\boldsymbol v}`$, and written to the output:

```math
	\alpha = \arctan \frac{w}{u} , \qquad
	\beta_{\rm side} = \arcsin \frac{v}{\left| {\boldsymbol v} \right| } ,
```

the angle of attack and the sideslip angle. They are the components of the total angle of
attack $\alpha_t$ used to enter the aerodynamic table, and reduce to it when the motion is
confined to one plane.


# How to start calculation

`Tacode` reads `config.yml` from the **current directory**, and writes its output there.
Run it from a case directory, not from the repository root.

## Trajectory simulation

```console
cd tutorial/work
./run_tacode.sh
```

or directly:

```console
cd tutorial/work
python3 ../../src/tacode.py
```

Tutorial cases:

| Case | Description |
|---|---|
| `tutorial/work` | orbital flight, 5000 s (about 1 s of run time) |
| `tutorial/work_reentry` | atmospheric entry at 7450 m/s from 150 km |
| `tutorial/work_reentry_6dof` | the same entry with the attitude solved, 1000 s (about 13 s of run time) |
| `tutorial/work_montecarlo` | several cases run in parallel |

`tutorial_template` is the template used to create a new case (copy it to a new directory).

Each case directory carries the tables it uses in its own `database` subdirectory, and
`config.yml` points at them with a path relative to the current directory:

```
tutorial/work
|-- config.yml
|-- run_tacode.sh
`-- database
    |-- atmosphere
    |   `-- atmospheremodel.txt
    `-- aerodynamic
        `-- aerodynamic.txt
```

```yaml
satellite:
  directory_path_specify: manual
  directory_aerodynamic: database/aerodynamic
  filename_aerodynamic: aerodynamic.txt

atmosphere:
  directory_path_specify: manual
  directory_atmosphere: database/atmosphere
  filename_atmosphere: atmospheremodel.txt
```

A case is therefore self-contained: it can be copied elsewhere, and its tables can be
edited without touching any other case. The `database` directory at the top of the
repository is the master copy that the cases are seeded from; it holds every table,
including the ones a given case does not use. Setting `directory_path_specify` to
`default` (or `auto`) reads that master copy instead, resolved from the location of the
source file and independent of the current directory.

Output is written to the directories named in `config.yml`:

| Directory | File | Description |
|---|---|---|
| `output_result` | `tecplot.dat` | Trajectory in Tecplot point format |
| `output_result` | `geodetic.kml` | Ground track for Google Earth |
| `output_restart` | `restart.dat` | State in ECEF Cartesian coordinates |

## Six-degree-of-freedom (attitude) simulation

```console
cd tutorial/work_reentry_6dof
./run_tacode.sh
```

The attitude is solved when `attitude.flag_attitude` is `True`. A configuration without an
`attitude` section, or with the flag set to `False`, behaves exactly as the earlier
three-degree-of-freedom versions did, down to the last bit of the output files; that is
checked against the committed reference outputs by `test/test_regression_3dof.py`.

In addition to the three-degree-of-freedom input, the following settings are read:

| Setting | Meaning |
|---|---|
| `attitude.flag_attitude` | switches the whole 6-DOF computation on |
| `attitude.inertia_tensor` | `Ixx`, `Iyy`, `Izz` and the products of inertia `Ixy`, `Iyz`, `Izx`, about the centre of gravity, in body axes, kg m² |
| `attitude.center_of_gravity` | vector from the moment reference point of the aerodynamic table to the centre of gravity, body axes, m |
| `attitude.damping_coefficient` | `Clp`, `Cmq`, `Cnr`, 1/rad; negative values damp |
| `attitude.static_stability_derivative` | $C_{m\alpha}$, 1/rad, used only with `kind_aerodynamic_model: constant` |
| `attitude.flag_moment_aerodynamic` | switch for the static aerodynamic moment |
| `attitude.flag_moment_damping` | switch for the dynamic damping moment |
| `attitude.flag_moment_gravity_gradient` | switch for the gravity-gradient moment |
| `initial_settings.attitude` | initial [yaw, pitch, roll] in deg, referred to the local horizon |
| `initial_settings.angular_velocity` | initial body rates [p, q, r] in deg/s, relative to the ECEF frame |

The aerodynamic table has to depend on the angle of attack for the motion to be physical;
`database/aerodynamic/README.md` describes the file format, the sign convention and the
sphere-cone sample table supplied with the code.

The Tecplot file then carries 13 more columns, after the existing ones, and the restart
file carries the quaternion and the angular velocity after the position and the velocity:

| Column | Meaning |
|---|---|
| `q0`–`q3` | attitude quaternion, ECEF to body |
| `Yaw`, `Pitch`, `Roll` | Euler angles referred to the local horizon, deg |
| `P`, `Q`, `R` | body rates relative to the ECEF frame, deg/s |
| `AoA`, `Sideslip`, `AoAtotal` | angle of attack, sideslip angle and total angle of attack, deg |

Note that the timestep now has to resolve the attitude oscillation rather than the
trajectory; see [Timestep of the attitude computation](#timestep-of-the-attitude-computation).

## Visualising the result

`src_helper/` holds the post-processing tools, one directory per tool, the same way `src/`
holds one directory per module. **They never run the solver** — they read the files it has
already written, which is why matplotlib stays out of the requirements for the computation
itself.

```console
cd tutorial/work_reentry_6dof
./run_tacode.sh
python3 ../../src_helper/animate_trajectory/animate_trajectory.py output_result/tecplot.dat -o attitude.gif
```

The animation puts the vehicle on its trajectory and turns it with the computed attitude,
so that the motion can be read at a glance:

| Panel | What it shows |
|---|---|
| left | the ground under the vehicle, the trajectory, and the vehicle drawn far larger than life so its attitude is visible; a globe inset gives the position on the planet |
| top right | the vehicle in the local horizon frame, seen from ahead of the velocity vector, with the velocity (dashed) and the body axes (x red, y green, z blue). The angle between the body x axis and the dashed line is the angle of attack |
| middle right | altitude against time, with a cursor at the current step |
| bottom right | angle of attack and sideslip against time |

Three shapes are available with `--shape`: `capsule` (a sphere-cone, matching the geometry
of the supplied aerodynamic table), `satellite` and `aircraft`. They are drawn for the eye
only and have no effect on any coefficient. `--snapshot <time>` writes a single frame as an
image instead of an animation, for a report or a paper.

A file from a three-degree-of-freedom run has no attitude columns; the tool says so and
draws the trajectory without orienting the vehicle.

**The output format follows the extension of `-o`.** The three animated formats are the
same picture; they differ in what they need and what they cost:

| Extension | Writer | Needs | Size |
|---|---|---|---|
| `.gif` | `PillowWriter` | nothing beyond matplotlib | 4.6 MB |
| `.html` | `HTMLWriter` | nothing beyond matplotlib | 16.1 MB |
| `.mp4` | `FFMpegWriter` | `ffmpeg` on the PATH | 1.1 MB |

The sizes are for the 6-DOF tutorial at 80 frames and `--dpi 70`. `.mp4` is by far the
smallest; `.html` is a single self-contained page with play, pause and a frame slider,
which needs no `ffmpeg` and, unlike a GIF, does not quantise the colours. `--snapshot`
writes a still instead, in whatever format matplotlib infers from the extension (`.png`,
`.pdf`, `.eps`).

See `src_helper/animate_trajectory/README.md` for the full list of options.

## Monte-Carlo simulation

```console
cd tutorial/work_montecarlo
./run_tacode-mc.sh
```

Tutorial case: `tutorial/work_montecarlo`.
It copies `tutorial_template` for each case and runs them in parallel
(`number_iteration` cases, up to `maximum_number_execution` at a time).


## Configuration file

Trajectory simulation by `Tacode` is controled by the configuration file: `config.yml`.
`tutorial_template/config.yml` carries every setting with a comment on what it does,
including the `attitude` section, which is present but switched off there.

## Tests

```console
./run_tests.sh              # run everything
./run_tests.sh -v           # verbose
./run_tests.sh test_kepler  # one module
```

The tests live in `test/` and use only the standard library's `unittest`, so they
need nothing beyond what `Tacode` itself requires. They cover:

| Module | What it checks |
|---|---|
| `test_coordinate_system.py` | Round trips between the Cartesian, geodetic and polar systems; behaviour at longitude 180 deg, at the poles and on the equator |
| `test_force_term.py` | The gravity vector equals `-grad U` for a potential written independently from the README; the Coriolis, centrifugal and drag terms |
| `test_kepler.py` | Energy and angular momentum are conserved in the two-body limit (J terms, rotation and drag switched off); RK4 converges at fourth order |
| `test_atmosphere.py` | Table interpolation reproduces the nodes, clamps outside the table range, and scales the Knudsen number with the characteristic length; the aerodynamic table is read in both formats, and its axisymmetry is checked so that a table a 6-DOF run cannot represent is reported rather than used silently |
| `test_solver_invariants.py` | The time loop stops at the requested time, all history arrays share one length, the atmospheric values line up with the position of the same index, and the Tecplot header matches the number of rows |
| `test_attitude.py` | Quaternion, Euler-angle and local-horizon conversions round trip and agree with analytic rotations and with the frame the initial velocity already uses; the quaternion kinematics reproduce a constant-rate rotation and drop the Earth rate |
| `test_moment_term.py` | Torque-free motion conserves the angular momentum and the energy, an axisymmetric body precesses at the analytic rate, the gravity-gradient torque matches its closed form and vanishes for an isotropic body, and the damping term removes rotational energy |
| `test_solver_attitude.py` | A 6-DOF run keeps the state arrays aligned with the trajectory, the pitch oscillation matches its analytic period, planar motion stays planar, damping shrinks the amplitude, and the aerodynamic force at zero incidence equals the 3-DOF drag |
| `test_regression_3dof.py` | With the attitude switched off, both tutorial cases reproduce the reference outputs committed in `tutorial/`, and a configuration carrying an `attitude` section set to `False` gives exactly the same trajectory as one without the section |
| `test_helper_visualization.py` | The post-processing tools in `src_helper/`: the Tecplot reader on both a 3-DOF and a 6-DOF output, the vehicle shapes (front distinguishable from back, roll visible), and the animation script writing an actual still, an `.html` animation as one self-contained file, and an `.mp4`. The drawing tests are skipped when matplotlib is not installed, and the `.mp4` one when `ffmpeg` is not |
| `test_attitude_verification.py` | Problems whose answer is known in closed form, solved by the production solver: the order of convergence, the Jacobi-elliptic solution of the torque-free asymmetric body, conservation of the angular momentum vector in inertial space, the precession of an axisymmetric body, the logarithmic decrement of a damped oscillation, the gravity-gradient libration frequency in a circular orbit, and the axisymmetry of the tabulated aerodynamics |
| `test_regression_6dof.py` | The 6-DOF tutorial case, run for its full 1000 s, reproduces the `tecplot.dat`, `restart.dat` and `geodetic.kml` committed in `tutorial/work_reentry_6dof` |
| `test_timestep_attitude.py` | The timestep itself: the `T/20` criterion the solver warns at is the first one whose numerical damping disappears, the order of convergence survives the atmosphere table, the angle-of-attack table is only C0 and costs the fourth order, and the shipped `dt = 0.05 s` is converged |

`test/smoke_tutorial.py` is separate from the suite above: it runs the tutorial
case end to end in a temporary directory and inspects the three output files, then
does the same for a shortened 6-DOF case and checks the attitude columns.
It needs `simplekml`, which the unit tests do not.

```console
python3 test/smoke_tutorial.py
```

### Continuous integration

`.github/workflows/tests.yml` runs on every push to `main` and on every pull
request:

| Job | What it does |
|---|---|
| `unit` | The test suite on Python 3.9, 3.10, 3.11, 3.12 and 3.13 |
| `tutorial` | The tutorial case end to end, through `run_tacode.sh`, and a check that a bad configuration exits with a non-zero code |
| `setup-script` | `setup_env.sh` on a machine without the packages, then the suite using the `.venv` it built |

## Requirements

`Tacode` needs Python 3.9 or later and the following packages.

| Package | Version | Used for |
|---|---|---|
| numpy | >= 1.15.4 | Arrays and linear algebra throughout |
| scipy | >= 1.4.0 | Interpolation of the atmosphere and aerodynamic tables |
| PyYAML | >= 3.11 | Reading `config.yml` (imported as `yaml`) |
| simplekml | >= 1.3.2 | KML output (`post_process.kml.flag_output`) |

Optional:

| Package | Version | Used for |
|---|---|---|
| gpxpy | >= 1.3.5 | GPX output only, which is disabled in `output_routine()` |
| matplotlib | >= 3.3 | The post-processing tools in `src_helper/`. The solver never imports it |

`ffmpeg` is optional as well, and is a program rather than a Python package: the animation
tool needs it on the PATH to write an `.mp4`, but not for a `.gif` or an `.html`.

Beyond these, `Tacode` uses only the Python standard library:
`argparse`, `os`, `random`, `shutil`, `subprocess`, `sys`.

```console
python3 -m pip install -r requirements.txt
```

Versions verified by the developer (2026-09-05): Python 3.12.3 with
numpy 2.2.6 / scipy 1.15.3 / PyYAML 6.0.2 / simplekml 1.3.6, and also with
numpy 2.5.2 / scipy 1.18.1 / PyYAML 6.0.3 / simplekml 1.3.2.

### Optional: build a virtual environment

If those packages are not available, `setup_env.sh` creates a Python virtual
environment in `.venv` and installs them. This is optional — if you already have
a working interpreter, keep using it and ignore this script.

```console
./setup_env.sh              # build .venv and install
./setup_env.sh --check      # only report what is missing, install nothing
./setup_env.sh --help       # all options
```

Once `.venv` exists, `run_tacode.sh` uses it automatically. The interpreter is
resolved in this order:

1. `$TACODE_PYTHON`, if set — `TACODE_PYTHON=/path/to/python ./run_tacode.sh`
2. `.venv/bin/python` at the repository root, if it exists
3. `python3`

With no `.venv` and no `TACODE_PYTHON`, the scripts behave exactly as before.


# Notes and current limitations

## Reference frame of the initial velocity

The initial **position** (`initial_settings.coordinate`) is given in the **WGS84 geodetic**
system as longitude [deg], latitude [deg], and altitude [km].

The initial **velocity** (`initial_settings.velocity`) is, however, interpreted in the local
frame of the **geocentric** latitude and longitude, as [east, north, up] in m/s.
The "north" and "up" directions of the two systems differ by up to the flattening angle
(about 0.19 deg). The Tecplot output uses the same geocentric frame, so a run is
self-consistent, but take this distinction into account when comparing with external data.

## Atmosphere table format

Two formats are accepted, and the format is detected from the file itself:

| Format | Recognised by | Produced by |
|---|---|---|
| CCMC / VITMO | a line reading `Selected parameters are:` | the NRLMSISE-00 form at [CCMC](https://ccmc.gsfc.nasa.gov/) |
| NRLMSISE-00 Fortran | a line starting `1, ALTITUDE` | `NRLMSISE-00_readctl.FOR` writing the table directly |

Both give the same seven columns: height, O, N2, O2, total mass density, neutral
temperature, N. Rows that are not seven numbers are skipped, so trailing blank or
comment lines are harmless.

Two tables are supplied in the master `database/atmosphere`. Switch between them with
`atmosphere.filename_atmosphere` alone — the format is detected, so nothing else
changes. A case directory carries only the table it uses, so copy the other one into
`<case>/database/atmosphere` before naming it.

| File | Altitude | Profile | Solar activity |
|---|---|---|---|
| `atmospheremodel.txt` (default) | 0–400 km | single point, 55°N 45°E | 2015-01-01, site default F10.7 |
| `atmospheremodel_700km.txt` | 0–700 km | mean over 50°S–50°N × all longitudes | F10.7 = 76.6 (low activity) |

`database/atmosphere/README.md` records how each was generated and why they are not
interchangeable — the density at 400 km differs by a factor of 2.17 between them.

Two cautions when swapping tables:

- **The solar activity must match your case.** Two tables generated for different
  F10.7 differ by a factor of two in density at 400 km, which changes the drag by
  the same factor.
- **The Knudsen number is built from N2, O2, N and O only.** Those are 96% of the
  mass at 400 km but only 35% at 700 km, where helium dominates. The drag
  coefficient is flat over that range, so the effect on the trajectory is
  negligible, but the Knudsen number itself is increasingly approximate above
  roughly 500 km.

## Atmosphere table range

The atmosphere table (`database/atmosphere/atmospheremodel.txt`, NRLMSISE-00) covers
0–400 km, while the tutorial orbit reaches 522 km. Above the top of the table the
properties are extrapolated:

```math
\rho(z) = \rho_{\rm top} \exp \left( - \frac{z - z_{\rm top}}{H} \right)
```

The upper thermosphere is close to isothermal and in diffusive equilibrium, so this is the
physically appropriate form rather than an ad-hoc fit. The scale height $H$ is not a
hard-coded constant: it is fitted to the top 50 km of whichever table is loaded, so
replacing the table changes it automatically. For the table shipped here the fit gives
47.8 km, which agrees with $kT/(m_{\rm O} g) = 49.3$ km — atomic oxygen makes up 96% of
the mass at 400 km.

The temperature is held at its value at the top of the table (the region is isothermal),
and the Knudsen number is scaled inversely with the density.

Set `atmosphere.kind_extrapolation` to `clamp` to restore the earlier behaviour, in which
the values were held at those of the top of the table. That overestimates the density by a
factor of 11 at 522 km, which moves the position by 93 m after one orbit and by 14 km
after nine.

Extrapolation is a safety net, not a substitute for data: extending the table itself is
better whenever the trajectory spends a long time above it.

## Verification of the attitude computation

The attitude solver is checked against problems whose solution is known in closed form,
in `test/test_attitude_verification.py`. The production solver is what is run in every
case; nothing is reimplemented for the comparison except the reference solutions
themselves.

| Problem | Reference solution | Agreement |
|---|---|---|
| Rotation at a constant rate, isotropic inertia, no torque | rotation about the angular velocity by \|ω\|t | 4th-order convergence (Runge-Kutta) |
| The same, coupled with the trajectory and the aerodynamics | none; three timesteps compared with each other | 4th order (Runge-Kutta), 1st order (explicit Euler) |
| Torque-free asymmetric rigid body | Jacobi elliptic functions, `cn`, `sn`, `dn` | below 1e-9 in ω after 20 s |
| Torque-free motion, any body | the angular momentum vector is fixed in inertial space | below 1e-8 |
| Torque-free axisymmetric body | constant nutation angle, precession at \|H\|/I<sub>t</sub> | 6 decimal places |
| Damped pitch oscillation, constant model | ζ = c/(2√(kI)), envelope exp(−ζω<sub>n</sub>t) | within 5 % on the decrement, 1 % on the frequency |
| Gravity-gradient libration in a circular orbit | ω = n√(3(I<sub>xx</sub>−I<sub>zz</sub>)/I<sub>yy</sub>) | within 1 % |
| Tabulated aerodynamics at a general attitude | force in the plane of the axis and the velocity, moment normal to it | below 1e-12 |

Two caveats are worth stating. First, the closed-form solution of the torque-free
asymmetric body is itself checked, by substituting it back into Euler's equation and
differentiating numerically, so that an error in the reference cannot hide an error in the
code. Second, the explicit Euler scheme converges at **second** order for a pure rotation
rather than first, because renormalising the quaternion turns the update into a rotation
about the correct axis whose angle is wrong only at third order in the timestep; in the
coupled problem it is first order as expected.

All of this is verification, not validation: it shows that the equations in
[Governing equation](#governing-equation) are solved correctly, not that they describe any
particular vehicle. Validation would need flight or wind-tunnel data, and the sphere-cone
table supplied with the code is an analytic model rather than measured data.

## Timestep of the attitude computation

The attitude oscillation is far faster than the orbital motion, and the timestep has to
resolve it rather than the trajectory. For a statically stable body the period is

```math
	T = 2 \pi \sqrt{ \frac{ I_{yy} }{ q_\infty S L \left| C_{m\alpha} \right| } } ,
```

so it shrinks as the dynamic pressure rises: in `tutorial/work_reentry_6dof` it is about
13 s at 150 km, about 5 s at 120 km and of the order of 0.1 s near peak dynamic pressure.
The code estimates this period during the run and prints a warning the first time the
timestep exceeds one twentieth of it. There is no substepping — the attitude and the
trajectory share the timestep — so a full entry down to the ground needs a timestep of a
few milliseconds, which is why the 6-DOF tutorial stops at 1000 s.

## Aerodynamic table for the attitude computation

A table with a single angle of attack (such as the default `aerodynamic.txt`) leaves the
coefficients independent of the attitude, so no restoring moment appears and the body
tumbles. The code warns about this at the start of the run. The sphere-cone table
supplied for the tutorial is an analytic model, not measured data; its assumptions and
their consequences are listed in `database/aerodynamic/README.md`.

**The body is assumed to be axisymmetric.** The table is looked up with the total angle
of attack alone and the coefficient vectors are then rotated about the body x axis onto
the actual crossflow plane, which is exact only for a body of revolution. Such a body
has `CFy = CMx = CMz = 0` everywhere in the table; if those columns are not zero, they
are rotated as if they lay in the crossflow plane and the side force and the rolling and
yawing moments point in the wrong direction, with no way to recover the right one from a
table that has no sideslip dependence. Tacode therefore checks the columns when it reads
a table for a 6-DOF run and warns, naming the offending ones, rather than answering
silently; see `database/aerodynamic/README.md` for the threshold. A genuinely
non-axisymmetric vehicle needs both the table and the rotation extended to the sideslip
angle.

The angle of attack is measured against the ECEF velocity, that is, the atmosphere is
assumed to co-rotate with the Earth. Winds are not modelled.

## Restart

Restart (`computational_setup.flag_initial: False`) is **not implemented** in this version;
the code stops with a message, and that includes the attitude. `output_restart/restart.dat`
is written on every run, but it currently stores the whole history rather than the final
state only. `restart_process.frequency_output` thins that history out; it defaults to 1,
which writes every step as before, and the last step is always written whatever the
setting.


# Contact:

Yusuke Takahashi, Hokkaido University

ytakahashi@eng.hokudai.ac.jp


# References

- Yusuke Takahashi, Masahiro Saito, Nobuyuki Oshima, and Kazuhiko Yamada, “Trajectory Reconstruction for Nanosatellite in Very Low Earth Orbit Using Machine Learning.” Acta Astronautica 194: 301–8. 2022. https://doi.org/https://doi.org/10.1016/j.actaastro.2022.02.010.
- Wagner, Carl. A. 1971. “The Gravity Potential And Force Field of the Earth Through Fourth Order.” NASA TN D-3317, 1–60.
- Fucheng Liu, Shan Lu, and Yue Sun, Guidance and Control Technology of Spacecraft on Elliptical Orbit. Springer, 2019.
- NRLMSISE-00 Atmosphere Model
- Defense Mapping Agency, “Department of Defense World Geodetic System 1984: its definition and relationships with local geodetic system“, 8350, 1987
- James R. Wertz (ed.), Spacecraft Attitude Determination and Control. D. Reidel, 1978. (quaternion kinematics, Euler's equation, gravity-gradient torque)
- Frank J. Regan and Satya M. Anandakrishnan, Dynamics of Atmospheric Re-Entry. AIAA Education Series, 1993. (attitude motion of a re-entry body, static and dynamic stability derivatives)