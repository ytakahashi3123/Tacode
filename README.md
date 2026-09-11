# Tacode
Trajectory analysis code

[![tests](https://github.com/ytakahashi3123/Tacode/actions/workflows/tests.yml/badge.svg)](https://github.com/ytakahashi3123/Tacode/actions/workflows/tests.yml)


# Code description

`Tacode` solves the equation of motion of an object on an Earth-centered, Earth-fixed (ECEF) noninertial frame by a python script.
The planetary gravity, Coriolis, centrifugal, and aerodynamic forces act on the object.
The gravity force is obtained by differentiating the gravitational potential considering J20, J22, J30, and J40.
The atmospheric data is given by NRLMSISE-00 Atmosphere Model.
The equation of motion is numerically solved using fourth-order Runge-Kutta method in four stages.

By default the object is a point mass with **three degrees of freedom**, and the aerodynamic force is the drag along the direction of motion, given by atmospheric density, drag coefficient, characteristic (projection) area, and velocity. A lift perpendicular to it can be added with `satellite.lift_coefficient`, at a constant bank angle or at one read from a table in time; it is off by default.

Since v2.3.0 the attitude can be solved as well, giving **six degrees of freedom**. The rigid body then rotates under the aerodynamic, damping and gravity-gradient moments, and the aerodynamic force follows the attitude through the angle of attack, so that the trajectory and the attitude are coupled in both directions. The attitude is carried by a quaternion, and the aerodynamic coefficients are interpolated from a table in the angle of attack as well as in the Knudsen number.
The attitude is switched off unless `attitude.flag_attitude` is set, and a configuration without an `attitude` section reproduces the earlier three-degree-of-freedom results bit for bit.

Since v2.4.0 the **wind** can be taken into account, so that the atmosphere no longer has to co-rotate with the Earth as a rigid body: the aerodynamics is then evaluated at the air-relative velocity, which reaches the drag, the moments, the dynamic pressure and the angle of attack alike. The wind is either one uniform vector or a table of up to four dimensions — time, longitude, latitude and altitude — prepared offline from NCEP reanalysis for the lower atmosphere and HWM14 for the thermosphere. It is switched off unless `wind.flag_wind` is set, and a configuration without a `wind` section reproduces the earlier results bit for bit.

Since v2.5.0 the three-degree-of-freedom motion can **fly a lift force**: `satellite.lift_coefficient` with either a constant `satellite.bank_angle` or a `satellite.bank_angle_table` of bank angle against time, each Runge-Kutta stage reading the table at its own time. That is what makes an entry with a measured roll history reproducible, and `validation/` compares two such entries — Apollo 4 and Apollo 10 — with their flight data. The run is also guarded now: a solution which leaves the physical range stops with a non-zero exit code instead of writing a diverged trajectory, every error path returns a non-zero code, a Monte-Carlo case which does not complete is reported rather than dropped, and `montecarlo.random_seed` makes a Monte-Carlo set reproducible. A configuration written for an earlier version reproduces its results bit for bit, as before.

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

A **lift** may be added to it, which is what a blunt lifting body such as a re-entry capsule
flies with. It is perpendicular to the velocity, and its direction in that plane is set by
the bank angle $\sigma$:

```math
	{\boldsymbol F}_{\rm lift} = q_{\infty} C_L S \left( \cos\sigma \, \hat{\boldsymbol n}_u + \sin\sigma \, \hat{\boldsymbol n}_r \right) ,
	\qquad
	\hat{\boldsymbol n}_u = \frac{ \hat{\boldsymbol u} - ( \hat{\boldsymbol u} \cdot \hat{\boldsymbol v} ) \hat{\boldsymbol v} }{ \left| \cdots \right| } ,
	\qquad
	\hat{\boldsymbol n}_r = \hat{\boldsymbol v} \times \hat{\boldsymbol n}_u ,
```

where $`\hat{\boldsymbol u} = {\boldsymbol x}/|{\boldsymbol x}|`$ is the local vertical
(geocentric, the frame the initial velocity and the wind already use). So $\sigma = 0$ is
lift up, $\sigma = 90^\circ$ is lift to the right of the flight direction and
$\sigma = 180^\circ$ is lift down; the bank angle is not modulated during a run. Both
`satellite.lift_coefficient` and `satellite.bank_angle` default to zero, and with
$C_L = 0$ the term is not evaluated at all, so a configuration without them reproduces the
earlier results bit for bit. In six degrees of freedom the lift comes out of the attitude
and the table instead, and these two settings are not used.

The bank angle may also be given as a **table against time**,
`satellite.bank_angle_table`, a list of `[time, angle]` pairs which is interpolated
linearly and clamped outside its range; each Runge-Kutta stage reads it at its own time,
as the wind is read. That is how a measured roll history is flown — Tacode has no
guidance law of its own — and `validation/apollo10` does exactly that.

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
| `tutorial/work_montecarlo_wind` | the same entry with the wind switched on, scattered case by case |
| `tutorial/work_reentry_wind_table` | the same entry flown through a real wind table (NCEP + HWM14) |

`tutorial_template` is the template used to create a new case (copy it to a new directory),
and `tutorial_template_wind` is the one the wind Monte-Carlo tutorial copies.

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

Each tool takes its settings either on the command line or from a settings file,
**`config_helper.yml` in the current directory**, in the same way the solver reads
`config.yml`. One section per tool, keys named after the long options, and the command
line wins when both are given:

```yaml
animate_trajectory:
  filename: output_result/tecplot.dat
  output: attitude.mp4

montecarlo_animation:
  directory: work_montecarlo_wind
  view: 3d
  spin: 0.0
```

An unknown key stops the run rather than being ignored, a missing file is not an error,
and `-file` points at another one. `tutorial/work_reentry_6dof/config_helper.yml` and
`tutorial/work_montecarlo_wind/config_helper.yml` are shipped with those tutorials, so
the tools can be run there with no arguments at all; `--save-config` writes the settings
in effect back out as a starting point, or as a record of what produced a figure. See
`src_helper/general/README.md`.

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

### Dispersion of a Monte-Carlo run

```console
cd tutorial/work_montecarlo_wind
python3 ../../src_helper/montecarlo_dispersion/montecarlo_dispersion.py work_montecarlo_wind
```

`montecarlo_dispersion.py` reads the last point of every case of a Monte-Carlo run and
reports how far each one lies from the others, resolved east and north in the local
horizon. The inputs that were dispersed are listed alongside: the tool compares each
case's control file with the `case_template` the driver leaves next to the cases, so it
reports whatever actually differs without being configured. `--reference` measures the
offsets from another run instead of from the mean of the cases, `-o` writes the table as
CSV, and `--plot` draws the scatter (matplotlib, needed for that option only).

See `src_helper/montecarlo_dispersion/README.md` for the full list of options.

### Animation of a Monte-Carlo run

```console
cd tutorial/work_montecarlo_wind
python3 ../../src_helper/montecarlo_animation/montecarlo_animation.py work_montecarlo_wind \
    --reference ../work_reentry/output_result/tecplot.dat \
    --mark "NCEP+HWM14 table=../work_reentry_wind_table/output_result/tecplot.dat" \
    -o montecarlo_dispersion.mp4
```

`montecarlo_animation.py` animates the cases of a Monte-Carlo run together with their
dispersion ellipses. The right-hand panel carries the point of it: every case is drawn
where it is **relative to the reference case at the same time**, in the local horizon, so
the cloud is seen growing out of a single point into the final ellipse, with the 1, 2 and
3 sigma ellipses and the CEP 50 % circle redrawn at every frame. The panel on the lower
left puts the 1 sigma spread on the same time axis as the altitude, which is where one
reads off *when* — and so at what altitude — the wind did its work. The points are
coloured by the wind of each case, so the fan is ordered by the input that made it.

`--view 3d` swaps the ground track for the trajectories themselves, in a box of
**longitude, latitude and altitude**: the absolute path the run computed, with the
dispersion ellipses lying on the floor at the impact point. `--view globe` draws the same
absolute trajectories in ECEF with the Earth around them, and `--exaggerate` stretches
the altitude about the surface when the descent itself is what should be visible (150 km
is 2 % of the radius). In both, the 100 cases lie on one another at that scale, so they
are drawn as a single bundle. The Earth carries a coastline — Natural Earth 1:110 m,
shipped as text in `src_helper/general/`, read with numpy alone, so no `cartopy` and no
download are involved — with the far hemisphere dropped, since matplotlib's 3-D does not
hide lines behind a surface. `--no-coastline` leaves it out.

`--view 3d-relative` measures the same box **from the reference case at the same time**
instead. That is the view for the dispersion itself: the bundle comes down tight and
unravels into the ellipse, each case coloured by its wind. The box is not to scale — tens
of kilometres across against a descent of 150 km — and `--altitude-max` cuts it down to
the part where the spread is built (40 km for the wind tutorial). `--spin 0` keeps the
camera still in any of the three.

`--view follow` keeps the ECEF axes of `globe` but **puts the camera on the reference
case**: the window is a cube around it, which is what makes tens of kilometres of
dispersion readable against a radius of 6378 km. Inside that window each case is drawn as
its **departure from the reference at the same time**, carried to where the reference is
now. A trail in absolute coordinates would be useless there — the vehicle covers some
90 km between frames, so it would leave the window immediately — whereas a departure
trail keeps the whole history of how the case came away, and its head is still the case's
true position. The window is the larger of what the dispersion needs and `--ground` times
the altitude of the reference, so the surface is in the frame from the first breath of the
entry, the window closes in as the vehicle comes down, and it opens again as the scatter
grows (`--window` fixes it instead). The half width is written in the state box, since the
axes are switched off, and `--scale-bar` adds a bar with ticks and its length beside it.
Compared with `3d-relative`, the three directions are at the same scale and the ground is
there; compared with `globe`, the Earth is not what the frame is spent on.

The output format follows the extension of `-o` (`.mp4`, `.gif`, or a self-contained
`.html`), and `--snapshot -1` writes the last frame as a single image instead — the
figure to put in a report. `--mark` follows a run that is not part of the set, drawn as a
star with its own trail and left out of the statistics.

See `src_helper/montecarlo_animation/README.md` for the full list of options.

## Monte-Carlo simulation

```console
cd tutorial/work_montecarlo
./run_tacode-mc.sh
```

Tutorial case: `tutorial/work_montecarlo`.
It copies `tutorial_template` for each case and runs them in parallel
(`number_iteration` cases, up to `maximum_number_execution` at a time).

Each entry of `montecarlo.target_variable` is `[variable, section, dispersion]`. The
nominal value is read from the Monte-Carlo control file, and the matching lines of the
copied control file are overwritten case by case with

```
value * (1 + dispersion*(U - 0.5)),   U uniform in [0, 1)
```

that is, a uniform relative scatter of +-`dispersion`/2 on each element of the list. The
scatter is **relative**, so an element whose nominal value is `0.0` stays `0.0` in every
case. The **section** is what tells `wind.velocity` apart from `initial_settings.velocity`.

Once every case has run, the driver gathers their Tecplot outputs into a single file,
`montecarlo.result_dir/montecarlo.filename_tecplot`, with one zone per case, so that all
the trajectories can be loaded at once. The variables line is carried over from the cases
themselves, so the columns follow whatever the run wrote, and cases that disagree on their
columns stop the run rather than being mixed. `montecarlo.filename_trajectory` names the
file to gather inside each case, `output_result/tecplot.dat` by default; set
`flag_tecplot: False` to skip the step, which is worth doing for a long run with many
cases (the 100 cases of the wind tutorial below add up to some 90 MB). Statistics of the
impact points are not computed here — that is what `montecarlo_dispersion.py` above does.

### A case which does not complete

The driver waits for every case, looks at its exit code, and checks that it left a result
file behind. If any case did not complete it lists them and stops with a non-zero exit
code, so that a broken run is not mistaken for a finished one:

```
Cases that did not complete:  2
--Exit code 1 : work_montecarlo/case0003
--No result file: work_montecarlo/case0003/output_result/tecplot.dat
--Set montecarlo.flag_allow_failure: True to go on with the cases that did complete.
Program stopped.
```

Set `montecarlo.flag_allow_failure: True` to report the failed cases and carry on with the
ones that did complete. This matters most for a large run: the statistics are built from
the cases that wrote a result, so three cases failing quietly out of a hundred would
otherwise leave a plot that looks perfectly reasonable. `run_tacode-mc.sh` pipes the run
through `tee`, so it sets `pipefail` to return the exit code of the run itself rather than
that of `tee`.

### Scattering the wind

```console
cd tutorial/work_montecarlo_wind
./run_tacode-mc.sh
python3 ../../src_helper/montecarlo_dispersion/montecarlo_dispersion.py work_montecarlo_wind \
    --reference ../work_reentry/output_result/tecplot.dat
```

Tutorial case: `tutorial/work_montecarlo_wind`, which copies `tutorial_template_wind`.
The template is the entry of `tutorial/work_reentry` with the wind switched on: a uniform
20 m/s east and 10 m/s north, scattered by +-50 % case by case. The vehicle of that case
is light (`m/(CD A)` of about 9.8 kg/m2) and spends some 1000 s below 32 km, so the wind
is not a perturbation there but a leading term.

`montecarlo_dispersion.py` collects the impact point of every case. Measured from the
same entry **without** wind, one run of the 100 cases gives

```
Mean offset from the origin:  East +27.409 km,  North +7.835 km
Standard deviation:  East 8.407 km,  North 4.185 km
```

The mean offset is where the nominal wind puts the impact point, and the standard
deviation is what the uncertainty on that wind costs. The cases are drawn at random and
no seed is set, so the numbers still move by a few hundred metres from run to run at 100
cases; the whole run takes some 20 s.

`--plot` draws the cases with their 1, 2 and 3 sigma covariance ellipses and the CEP 50 %
circle. Ellipses rather than circles: a wind uncertainty spreads the impact point mostly
along the direction the wind blows, and in this case the scatter is twice as wide east to
west as it is north to south.

Scattering `initial_settings.density_factor` by +-10 % as well (it is in the control
file, commented out) spreads the impact point by some 80 km along the ground track, an
order of magnitude more than the wind, and the wind then hides inside it — which is the
reason the tutorial scatters one at a time.

### A real wind field

A uniform wind is only an estimate. `tutorial/work_reentry_wind_table` flies the same
entry through `database/wind/wind_merged_20240101_pacific.txt` — NCEP/NCAR Reanalysis 1
for 2024-01-01T00Z below 31 km, HWM14 above, merged through a 20-30 km transition layer:

```console
cd tutorial/work_reentry_wind_table
./run_tacode.sh
```

The impact point moves 12.8 km (10.1 km west, 7.9 km north) from the calm case, not the
29 km east a uniform 20 m/s suggests, and it moves the other way: most of the deflection
is picked up between 80 and 150 km, where the vehicle spends 1552 s and NCEP has no data
at all. Marking it on the Monte-Carlo plot puts the two side by side:

```console
cd tutorial/work_montecarlo_wind
python3 ../../src_helper/montecarlo_dispersion/montecarlo_dispersion.py work_montecarlo_wind \
    --reference ../work_reentry/output_result/tecplot.dat \
    --mark "NCEP+HWM14 table=../work_reentry_wind_table/output_result/tecplot.dat" \
    --plot work_montecarlo_wind/dispersion.png
```

The star lands outside the 3 sigma ellipse of the constant-wind cases, on the far side of
the calm impact point. Scattering a uniform wind measures the sensitivity to a wind; it
does not bracket the wind that was actually there.

A wind **table** cannot be scattered directly, since the driver varies numbers in the
control file and a table is selected by its file name. What can be scattered is
`wind.velocity_factor`, which scales the whole field:

```yaml
  target_variable:
    -
      - velocity_factor # Variable name
      - wind # Variable's root name
      - 0.4          # Dispersion in random
```

Copy `tutorial/work_reentry_wind_table` to serve as the template (its `run_tacode.sh`
takes the path of `src/` as `$1`, as `tutorial_template_wind/run_tacode.sh` does) and
point `montecarlo.template_path` at it. Scattering the merged NCEP+HWM14 field by +-20 %
in this way moves the impact point by 10.1 km west and 7.8 km north on average — the
nominal table, as it must be — with a standard deviation of 1.05 km east-west and 0.82 km
north-south. That measures the uncertainty on the *strength* of a known field, which is a
different question from the +-50 % on a uniform wind above; neither brackets the error of
the field itself.


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
| `test_force_term.py` | The gravity vector equals `-grad U` for a potential written independently from the README; the Coriolis, centrifugal and drag terms; the lift is perpendicular to the air-relative velocity, has the magnitude the coefficient asks for, turns about the velocity with the bank angle, and is not evaluated at all when `lift_coefficient` is absent or zero. The bank table is interpolated and clamped, wins over the constant, and stops the run when malformed; through the solver, a constant table reproduces the scalar bank bit for bit, lift up raises the trajectory and lift down lowers it, and a table which switches lands between the two |
| `test_kepler.py` | Energy and angular momentum are conserved in the two-body limit (J terms, rotation and drag switched off); RK4 converges at fourth order |
| `test_atmosphere.py` | Table interpolation reproduces the nodes, clamps outside the table range, and scales the Knudsen number with the characteristic length; the aerodynamic table is read in both formats, and its axisymmetry is checked so that a table a 6-DOF run cannot represent is reported rather than used silently; `kind_aerodynamic_model: constant` does not read the table at all, and a table written by `database/atmosphere/generate_atmosphere_table.py` is read back |
| `test_database_path.py` | The path of the atmosphere, aerodynamic and wind tables: `default` and `auto` resolve to the master under `database/`, `manual` is taken as given, a misspelt `directory_path_specify` stops the run instead of falling back, `manual` without its directory key stops as well, and the copy of every table in a case directory is identical to the master |
| `test_solver_invariants.py` | The time loop stops at the requested time, all history arrays share one length, the atmospheric values line up with the position of the same index, and the Tecplot header matches the number of rows |
| `test_attitude.py` | Quaternion, Euler-angle and local-horizon conversions round trip and agree with analytic rotations and with the frame the initial velocity already uses; the quaternion kinematics reproduce a constant-rate rotation and drop the Earth rate |
| `test_moment_term.py` | Torque-free motion conserves the angular momentum and the energy, an axisymmetric body precesses at the analytic rate, the gravity-gradient torque matches its closed form and vanishes for an isotropic body, and the damping term removes rotational energy |
| `test_solver_attitude.py` | A 6-DOF run keeps the state arrays aligned with the trajectory, the pitch oscillation matches its analytic period, planar motion stays planar, damping shrinks the amplitude, and the aerodynamic force at zero incidence equals the 3-DOF drag |
| `test_regression_3dof.py` | With the attitude switched off, both tutorial cases reproduce the reference outputs committed in `tutorial/`, and a configuration carrying an `attitude` section set to `False` gives exactly the same trajectory as one without the section |
| `test_montecarlo.py` | `montecarlo.random_seed` makes a run repeatable (the same seed gives the same dispersion, another seed does not, and no shipped configuration fixes a seed), the dispersion is multiplicative so a base value of zero stays put, and the case directories are found by their four-digit suffix. The driver rewrites the right line of the control file: the section tells `wind.velocity` apart from `initial_settings.velocity`, lines that happen to hold the same value are not rewritten together, the search is closed at the end of the section so that a key of another section is never rewritten, and a missing section, a missing key or a key which is not a list stops the run. A case which exits with a non-zero code or writes no result file is counted, and the run stops with a non-zero exit code unless `flag_allow_failure` is set. The wind tutorial and its template agree with each other, and two shortened cases actually run and come out different. The postprocess gathers the cases into one Tecplot file: one zone per case with its own point count, the template left out, a case without a result skipped, and cases whose columns disagree stopping the run. The animation helper is checked on its geometry — the offsets from the reference at the same time, the window that holds every point, the unwrapped longitude — and on actually writing a frame, a self-contained `.html` and the 3D view, which is skipped without matplotlib |
| `test_helper_config.py` | The settings file of the post-processing tools: the order of precedence (command line, then the file, then the default), a list-valued option, a missing file being no error, and an unknown key, a malformed section or a missing required value stopping the run. `--save-config` writes what can be read back and keeps the other sections, all three tools read the file the same way, and the `config_helper.yml` shipped with the tutorials is accepted by the tool it belongs to |
| `test_helper_visualization.py` | The post-processing tools in `src_helper/`: the Tecplot reader on both a 3-DOF and a 6-DOF output, a file of several zones being refused rather than joined into one trajectory (and taken apart by `read_tecplot_zone`), the vehicle shapes (front distinguishable from back, roll visible), and the animation script writing an actual still, an `.html` animation as one self-contained file, and an `.mp4`. The drawing tests are skipped when matplotlib is not installed, and the `.mp4` one when `ffmpeg` is not |
| `test_attitude_verification.py` | Problems whose answer is known in closed form, solved by the production solver: the order of convergence, the Jacobi-elliptic solution of the torque-free asymmetric body, conservation of the angular momentum vector in inertial space, the precession of an axisymmetric body, the logarithmic decrement of a damped oscillation, the gravity-gradient libration frequency in a circular orbit, and the axisymmetry of the tabulated aerodynamics |
| `test_regression_6dof.py` | The 6-DOF tutorial case, run for its full 1000 s, reproduces the `tecplot.dat`, `restart.dat` and `geodetic.kml` committed in `tutorial/work_reentry_6dof` |
| `test_timestep_attitude.py` | The timestep itself: the `T/20` criterion the solver warns at is the first one whose numerical damping disappears, the order of convergence survives the atmosphere table, the angle-of-attack table is only C0 and costs the fourth order, the shipped `dt = 0.05 s` is converged, and the check itself runs once per `INTERVAL_CHECK_TIMESTEP` rather than every step |
| `test_epoch.py` | The absolute time: an ISO 8601 epoch is read from the configuration in every form PyYAML can hand over, the derived day of year, universal time and Julian date are right across a leap year and a year boundary, the Julian date agrees with `astropy` where it is installed, and the epoch stays off in every configuration shipped with the code |
| `test_error_exit.py` | The exit code of an error path: no module in `src/` calls the built-in `exit()`, which returns 0 and hides a failure from the shell, the Monte-Carlo driver and CI, and a mistyped wind model, a wind velocity of the wrong length and a malformed epoch each stop with the code 1, both in process and when the solver is run as a child process |
| `test_divergence.py` | The sanity check on the state: the limits are built from the initial state (ten times the geocentric distance, ten times the escape velocity) and follow the configuration, a state which is not finite - position, velocity, quaternion or angular velocity - stops the run, so does a distance or a speed beyond the limit, a state just inside them does not, a coarse time step which used to write a diverged trajectory and return 0 now stops with the code 1, the same run goes through with the check off, and a sound run is untouched |
| `test_wind.py` | The wind, including the time axis, the merge of a lower and an upper table and the absence of HWM14: the frame of the given components agrees with the one the initial velocity already uses, the identity that the aerodynamics under a wind `w` at velocity `v` equals the aerodynamics in calm air at `v - w`, that gravity and the rotation terms and the angular velocity do not see the wind, the downwind drift of a full re-entry and the terminal speed relative to the air, the columns written to `tecplot.dat`, and that the wind stays off, and `velocity_factor` stays at 1.0, in every configuration shipped with the code. `velocity_factor` scales the constant wind and the table alike, 0.0 brings back the co-rotating atmosphere, and a factor which is not a single number stops the run. For a table: the nodes are reproduced and the interpolation between them is linear, the row order does not matter, a one-node axis degenerates so that a vertical profile needs no separate path, longitudes are folded and a global table joins across the seam, the horizontal is clamped and the vertical follows `kind_extrapolation`, an incomplete or duplicated grid is rejected, a mismatched epoch is reported, and the shipped tables and the offline paths of the generator are read back. With a time axis: the nodes and the linear interpolation between them, the shift by the epoch of the run, the clamp outside the range, the refusal to guess when there is no epoch, and that each Runge-Kutta stage reads its own time — shown by one step over the whole window, where the effective wind of a ramp is half its end value. For the merge: the blend across the transition layer, the alignment of the horizontal nodes, and the refusal of tables on different axes. The paths that fetch NCEP and that call HWM14 are not exercised, since they need the network and a Fortran build |

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

## Validation against flight data

`validation/` holds cases which put the code against measurements, as opposed to
`tutorial/`, which shows how to run it, and `test/`, which checks it against exact
solutions and its own reference outputs. They are **not** part of `run_tests.sh`: they
are compared against flight data rather than against an exact answer, and their result
is a write-up rather than a pass or a fail.

| Case | Flight | Result |
|---|---|---|
| `validation/apollo4` | Apollo 4 (AS-501) entry, 1967-11-09 | with the vertical lift measured in flight, altitude within 0.5 km to the peak heating and 5.4 km rms over the whole entry; the atmosphere table within 4.3 % of the flight-derived density |
| `validation/apollo10` | Apollo 10 (AS-505) entry, 1969-05-26 | **the bank angle measured in flight as the input**: altitude 3.2 km rms, inertial velocity 189 m/s rms over the whole entry |

```console
cd validation/apollo4
./run_tacode.sh                            # drag only,   1.4 s
./run_tacode.sh -file config_lift.yml      # with lift,   3.6 s
./run_tacode.sh -file config_6dof.yml      # 6-DOF,      19 s
python3 compare_apollo4.py                 # the numbers and the figures
```

Apollo 4 flew a lunar-return entry, and NASA TM X-58091 tabulates its altitude, relative
velocity, free-stream density and heating through the whole of it. The case compares
three things separately, because they fail in different ways:

- **The trajectory.** Apollo 4 flew with lift, roll-modulated: it dipped to 55.6 km,
  pulled back up to 73.5 km and only then came down. A drag-only run cannot produce that
  at all — it descends monotonically and is a comparison only up to the peak heating 72 s
  in, where it holds to 3.2 km and 194 m/s. **NASA TN D-5399 reports the vertical
  lift-to-drag ratio derived from the flight itself**, and driving the case with it (and
  nothing tuned) the run holds the flight to half a kilometre through the dive, the dip and
  the pull-up — **0.4 km and 46 m/s over the first 72 s** — and to 5.4 km rms over the
  whole 552 s, drifting low late in the entry where constant `CD` and `L/D` stop being
  a good description.
- **The same entry in six degrees of freedom**, where nothing is prescribed about the lift
  or the trim: the capsule carries a modified-Newtonian table of its own shape and trims
  at an angle of attack because its centre of gravity is offset from the axis by the
  measured 0.16688 m. It trims at 25.65 deg. against the 24.4 deg. reported for the
  flight, with `CD` and `L/D` within a few per cent of the flight-derived values, and
  nothing adjusted.
- **The atmosphere table**, evaluated at the flight altitudes so that the trajectory
  cannot hide in it: the NRLMSISE-00 table generated for the date and the place of the
  entry matches the flight-derived density to a mean ratio of 1.000 and an rms of 4.3 %
  over 42 points and two and a half decades of density.
- **The heating.** Stagnation-point correlations (Detra-Kemp-Riddell, Sutton-Graves) fed
  with the density of the table land within a percent of the same correlations fed with
  the flight density, which is what matters for the solver; the DKR value itself is
  within 6 % of the flight measurement when the equivalent nose radius of the reference
  (a 10-foot sphere at 24.4 deg. angle of attack) is used.

Every assumption, every source and the sensitivity to the assumed mass, azimuth, time
step and atmosphere table are written down in `validation/apollo4/README.md`.

**`validation/apollo10` closes the gap the bank angle leaves open.** NASA TN D-6725
reports the state vector at the entry interface, the roll angle the guidance held through
the entry, and the altitude and velocity flown. The roll angle is the direction of the
lift, so it goes in as an input (`satellite.bank_angle_table`, digitised from the figure
by the script in the case directory) and the trajectory is the comparison: **3.2 km rms in
altitude and 189 m/s rms in inertial velocity over the whole 414 s**, with nothing tuned.
Switching the lift off costs 33 km rms, and a change of 0.03 in `CL` — the accuracy the
report claims for its own L/D — costs two to three times the baseline error, which is
where this case runs out of resolution.

The atmosphere table of the case was written by
`database/atmosphere/generate_atmosphere_table.py`, which calls NRLMSISE-00 through
`pymsis` and writes the CCMC format the solver reads. `pymsis` is a dependency of that
generator only, in the same way that the wind-table generator owns its data sources.

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
| astropy | >= 5.0 | A reference for the Julian date in `test_epoch.py` only. That test is skipped when it is absent, and neither the solver nor the tools import it |

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
and the Knudsen number is scaled inversely with the density. The exponent is capped at 100
scale heights: beyond that the exponential underflows to a density of exactly zero and the
Knudsen number becomes `Inf`, which then reaches the output file. That far out the
atmosphere is a vacuum and the free-molecular coefficients no longer change, so holding
the values there costs nothing.

Set `atmosphere.kind_extrapolation` to `clamp` to restore the earlier behaviour, in which
the values were held at those of the top of the table. That overestimates the density by a
factor of 11 at 522 km, which moves the position by 93 m after one orbit and by 14 km
after nine.

Extrapolation is a safety net, not a substitute for data: extending the table itself is
better whenever the trajectory spends a long time above it.

## A run which diverges

Integration ends when the altitude goes below zero, or when `time_elapsed_maximum` is
reached. Neither of those catches a time step which is too coarse for the orbit: the
trajectory is then thrown outwards, and the numbers — a geocentric distance of 6.6e9 m, a
velocity of 1.4e6 m/s — were written to `tecplot.dat`, `restart.dat` and `geodetic.kml`
like any other result, with the exit code 0.

The state is therefore checked once a step:

- every component of the position, the velocity and, with the attitude solved, the
  quaternion and the angular velocity is a finite number,
- the geocentric distance is within `factor_radius_maximum` times its initial value
  (10 by default),
- the speed is within `factor_velocity_maximum` times the escape velocity at the initial
  position (10 by default).

Falling outside prints what happened and stops with a non-zero exit code, so nothing is
written:

```
The solution has diverged at 400.000 s: the geocentric distance 2.66051e+08 m is beyond the limit 6.56667e+07 m.
--A time step which is too coarse for the orbit is the usual cause;
--halve time_integration.timestep_constant and run it again.
```

The limits are taken from the initial state rather than given as absolute numbers, since
an absolute limit means something different for every orbit. A sound run never comes near
them — the tutorial cases reproduce their reference output byte for byte with the check in
place. `computational_setup.flag_check_divergence: False` switches it off, for a
deliberately hyperbolic trajectory or when the limits themselves are in the way.

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

A body of revolution can still fly at an angle of attack, and that is how a re-entry
capsule gets its lift: its centre of gravity sits off the axis, so the aerodynamic force
has a moment arm about it. The table is written about a point **on the axis** and stays
axisymmetric; `attitude.center_of_gravity` is the vector from that point to the centre of
gravity, and the solver adds `(-r_cg) x F` to the moment. The trim angle follows from the
table and the offset, and nothing else has to be told about the lift. The Apollo
validation case is built this way, and
`database/aerodynamic/generate_aerodynamic_table.py --shape capsule` writes a table for
that shape.

The angle of attack is measured against the ECEF velocity, that is, the atmosphere is
assumed to co-rotate with the Earth. Winds are not modelled.

## Wind

Without a wind the atmosphere is taken to co-rotate with the Earth as a rigid body, so
the ECEF velocity *is* the air-relative velocity. The `wind` section relaxes that: the
aerodynamics is then evaluated at

$$
\boldsymbol{v}_\mathrm{air} = \boldsymbol{v}_\mathrm{ECEF} - \boldsymbol{v}_\mathrm{wind}
$$

which changes the drag, the 6-DOF force and moment, the dynamic pressure, the angle of
attack and the sideslip angle, and the pitch period the timestep is checked against.
Gravity, the Coriolis and the centrifugal terms keep using the ECEF velocity, and the
angular velocity is untouched: a wind translates the air, it does not rotate it.

```yaml
wind:
  flag_wind: False
  kind_wind_model: constant
  velocity:
    - 0.0
    - 0.0
    - 0.0
```

`flag_wind: False` is the default and the shipped behaviour, and the outputs are then
exactly what they were before the section existed. `constant` applies one uniform wind
everywhere; `fileread` interpolates a table.

The components are `[East, North, Up]` in the **local geocentric** horizon in m/s, which
is the frame the initial velocity is already given in, converted through the same routine.
It is geocentric rather than geodetic, so it is tilted from the true horizontal by up to
0.19 deg, exactly as described under *Reference frame of the initial velocity* above; a
horizontal wind of 60 m/s therefore picks up about 0.2 m/s of spurious vertical
component. Meteorological `u` and `v` are defined in the geodetic horizon, so a table
converted from such a source inherits that tilt. Keeping one convention for the initial
velocity, the velocity output and the wind was preferred over mixing two.

With the wind switched on, `tecplot.dat` gains `WindE[m/s]`, `WindN[m/s]`, `WindU[m/s]`
and `VelairAbs[m/s]`. The ground-relative velocity `Upl/Vpl/Wpl` and its magnitude
`VelplAbs[m/s]` stay as they are, so both speeds can be read off the same row. The wind
is a function of the position alone, so it is not stored during the run and is looked up
again when the output is written, in the same way as the Euler angles and the aerodynamic
angles.

### Scaling the wind

```yaml
wind:
  velocity_factor:
    - 1.0
```

`velocity_factor` multiplies the whole wind field, whatever the model. It is written as a
one-element list exactly as `initial_settings.density_factor` is, and for the same
reason: that is what lets the Monte-Carlo driver scatter it, since the driver rewrites
numbers in the control file and a wind **table** is chosen by its file name rather than by
a number. The default of 1.0 is exact, so a configuration that leaves it alone gives the
same result bit for bit as one without the key.

### Wind tables

`kind_wind_model: fileread` reads a table of longitude × latitude × altitude, with the
path resolved exactly as the atmosphere and aerodynamic tables are:

```yaml
wind:
  flag_wind: True
  kind_wind_model: fileread
  directory_path_specify: manual
  directory_wind: database/wind
  filename_wind: wind_ncep_20240101_pacific.txt
  kind_extrapolation: zero
```

There is one file format regardless of dimensionality: comment lines, then one row per
grid point holding longitude, latitude, altitude and the three components — or seven
numbers, with the time in front, when the table carries a time axis. **A one-dimensional
vertical profile is written as a field with one longitude node and one latitude node** —
axes carrying a single node are dropped when the interpolator is built, so the profile
becomes a function of altitude alone without a second code path, and a single snapshot is
a table whose time axis has one node. The grid must be filled completely; a table whose
rows do not add up to the product of its axes is rejected rather than padded, since a
missing point and a calm point would otherwise look the same. Longitudes are folded into
`[-180, 180)` so a table in the meteorological `0-360` convention reads correctly, and a
table spanning the globe is detected and joined across the seam. The interpolator is built
once at start-up and only evaluated during the run, as the atmosphere and aerodynamic
interpolators are.

Outside the altitude range of the table, `kind_extrapolation` decides: `zero` (the
default) leaves the air co-rotating with the Earth, and `clamp` holds the value at the
end. `zero` is the default because NCEP reaches only about 31 km while a re-entry begins
far above it, and clamping would stretch the stratospheric wind into the thermosphere.
The horizontal directions are always clamped to the edge of the region. Leaving the range
of the table warns once, naming the range and the position.

### Interpolation in time

A table with a time axis carries its times in seconds from its own `# Epoch (UTC):` line,
and the `epoch` section places the run on that axis: the wind at elapsed time `t` is read
at `(epoch of the run - epoch of the table) + t`. The epoch is therefore not optional
here, and a table with a time axis but no epoch — in the table or in the configuration —
stops the run rather than guessing. Outside the range of the axis the nearest snapshot is
used and the run warns once. For a single snapshot the elapsed time is not consulted at
all, and a mismatch against the epoch is only reported.

Each Runge-Kutta stage reads the wind at its own time, `t`, `t + dt/2`, `t + dt/2` and
`t + dt`, since the force is now a function of time. It made no difference while every
term was time-independent.

### Where the tables come from

`database/wind/generate_wind_table.py` prepares the tables. It fetches a region of
NCEP/NCAR Reanalysis 1 — past observations, not a forecast — from NOAA PSL's OPeNDAP
service, needing nothing beyond the standard library and numpy, since it reads the
`.ascii` response rather than NetCDF; only the range asked for is transferred, and asking
for a duration rather than an instant gives the table a time axis. The wind sits on
pressure levels, so the altitudes come from the geopotential height of the same time and
grid.

Above about 31 km NCEP has nothing, and the model for the thermosphere is HWM14 (NRL),
which reaches 500 km. HWM14 is Fortran plus the NRL coefficient files, so it is not a
dependency: the script imports it only when asked for it, and otherwise says how to obtain
and build it and stops. Its output is written in the same format, and `--source merge`
joins a lower and an upper table across a transition layer, since the two models disagree
by tens of m/s where they meet and switching at a single altitude would make the wind
jump. Details, the sign convention and the shipped samples are in
`database/wind/README.md`.

### How much the wind matters

The effect is not a small perturbation in the re-entry case. The vehicle of
`tutorial/work_reentry` is light (`m/(C_D A)` of about 9.8 kg/m²) and reaches the ground
at 13 m/s after spending some 1000 s below 32 km, where a wind of tens of m/s is
comparable to its own speed. A uniform 20 m/s easterly moves the landing point by
about 29 km.

A real profile is not uniform, and its layers partly cancel. With
`wind_ncep_20240101_pacific.txt` — the actual analysis for 2024-01-01 00 UTC over the
region the vehicle comes down in, where the wind is easterly in the stratosphere and
westerly near 10 km — the landing point moves 2.1 km, almost all of it northward.

Adding the thermosphere changes that again. The vehicle spends 1552 s between 80 and
150 km, where HWM14 gives winds of tens of m/s, and with
`wind_merged_20240101_pacific.txt` covering 0 to 500 km the landing point moves 12.8 km
(10.1 km west, 7.9 km north) rather than 2.1 km. Most of the displacement is picked up
above the range NCEP alone can reach.

## Epoch (absolute time)

The solver marches an elapsed time that starts at zero; it needs no calendar. The models
that depend on the date do: the wind reads a NCEP time axis and hands HWM14 a day of year
and a universal time, and a lunar and solar ephemeris would need the same. The `epoch`
section carries that absolute time, once, for all of them.

```yaml
epoch:
  flag_epoch: False
  datetime: '2026-03-21T03:00:00Z'
```

`flag_epoch: False` is the default and the shipped behaviour: nothing reads the date and
the output files are exactly what they were before the section existed. With it switched
on, `datetime` is read as UTC in ISO 8601. Quoting it keeps it a string, but an unquoted
timestamp works as well, since PyYAML resolves it to a `datetime` of its own and the code
accepts both; a value with no time zone is taken as UTC, and one with an offset is
converted to it.

The epoch is written into the header of `tecplot.dat` and `restart.dat` as a comment
rather than as a column, so that the Tecplot point format keeps its purely numeric rows.
The UTC of a row is the epoch plus its `Time[s]`.

Leap seconds are not modelled: UTC is treated as a uniform time scale, which leaves the
epoch off a strict UTC count by the number of leap seconds since the date it was written
for. The models this feeds — a six-hourly wind grid, a lunar position good to a few
arcminutes — are insensitive to that at the level of a second.

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