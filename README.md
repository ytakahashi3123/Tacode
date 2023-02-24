# tacode
Trajectory analysis code


# Code description

This code solves the equation of motion of a point-mass object in three degrees of freedom on an Earth-centered, Earth-fixed (ECEF) noninertial frame by a python script.
The planetary gravity, Coriolis, centrifugal, and aerodynamic forces act on this point mass.
The gravity force is obtained by differentiating the gravitational potential considering J20, J22, J30, and J40.
The aerodynamic force was given by atmospheric density, drag coefficient, characteristic (projection) area, and velocity.
The atmospheric data is given by NRLMSISE-00 Atmosphere Model.
The equation of motion is numerically solved using fourth-order Runge\UTF{2013}Kutta method in four stages.


# How to start calculation

```
python3 tacode.py
```

# Contact:

Yusuke Takahashi, Hokkaido University

ytakahashi@eng.hokudai.ac.jp


# Date:

2023/02/21
