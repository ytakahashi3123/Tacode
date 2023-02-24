#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/05/23

import numpy as np
from orbital.orbital import orbital


def force_initialsettings(config):

  force = np.zeros(5*3).reshape(5,3)

  return force


def force_routine(config, coordinate, velocity, mass_satellite, area_satellite, cdmean_aerodynamic, density_factor, density, force):
  
  potential_factor     = config['planet']['potential_factor']
  radius_equat_planet  = config['planet']['radius']
  gravity_const_planet = config['planet']['gravitational_constant']
  mass_planet          = config['planet']['mass']
  rotation_rate_planet = config['planet']['rotation_rate']


  # Gravitational force
  radius_pmass    = np.sqrt( coordinate[0]**2 + coordinate[1]**2 + coordinate[2]**2 )
  radius_by_coord = radius_equat_planet/radius_pmass                   # a_e/r
  gme_by_radius2  = gravity_const_planet*mass_planet/(radius_pmass**2) # G*M_e/r2

  sin_beta = coordinate[2]/radius_pmass
  cos_beta = np.sqrt( coordinate[0]**2 + coordinate[1]**2 )/radius_pmass

  angle_long = np.arccos( coordinate[0]/np.sqrt( coordinate[0]**2 + coordinate[1]**2 ) ) * np.sign( coordinate[1] )
  sin_labd   = np.sin( angle_long )
  cos_labd   = np.cos( angle_long )

  J2  = potential_factor['J2']
  J22 = potential_factor['J22']
  J3  = potential_factor['J3']
  J4  = potential_factor['J4']
  labd22           = potential_factor['Lambda22']*orbital.deg2rad
  sin_2labd_labd22 = np.sin( 2.0*(angle_long + labd22 ) )
  cos_2labd_labd22 = np.cos( 2.0*(angle_long + labd22 ) )

  force_g_r = gme_by_radius2 * (- 1.0                                                                           \
                                + 1.5    *radius_by_coord**2*J2  *( 3.0*sin_beta**2 -  1.0 )                    \
                                + 9.0    *radius_by_coord**2*J22 *(     cos_beta**2         )*cos_2labd_labd22  \
                                + 2.0    *radius_by_coord**3*J3  *( 5.0*sin_beta**3 -  3.0*sin_beta )           \
                                + 5.0/8.0*radius_by_coord**4*J4  *(35.0*sin_beta**4 - 30.0*sin_beta**2 + 3.0 )  \
                                )
  force_g_a = gme_by_radius2 * (  6.0    *radius_by_coord**2*J22 *(     cos_beta            )*sin_2labd_labd22 )
  force_g_b = gme_by_radius2 * (                                                                                \
                                - 1.0    *radius_by_coord**2*J2  *( 3.0*sin_beta*cos_beta  )                    \
                                + 6.0    *radius_by_coord**2*J22 *(     sin_beta*cos_beta  )*cos_2labd_labd22   \
                                - 0.5    *radius_by_coord**3*J3  *(15.0*sin_beta**2 -  3.0 )*cos_beta           \
                                - 0.5    *radius_by_coord**4*J4  *(35.0*sin_beta**3 - 15.0*sin_beta )*cos_beta  \
                               ) 

  force[1,0] =  ( force_g_r*cos_beta - force_g_b*sin_beta )*cos_labd - force_g_a*sin_labd
  force[1,1] =  ( force_g_r*cos_beta - force_g_b*sin_beta )*sin_labd + force_g_a*cos_labd
  force[1,2] =    force_g_r*sin_beta + force_g_b*cos_beta


  # Coriolis
  force[2,0]  = 2.0*rotation_rate_planet*velocity[1]
  force[2,1]  =-2.0*rotation_rate_planet*velocity[0]
  force[2,2]  = 0.0


  # Centrifugal
  force[3,0] = rotation_rate_planet**2*coordinate[0]
  force[3,1] = rotation_rate_planet**2*coordinate[1]
  force[3,2] = 0.0


  # Aerodynamic (Fx = 1/2 rho U^2 * Ux/U)
  velocity_mag   = np.sqrt(velocity[0]**2+velocity[1]**2+velocity[2]**2)
  fact_aero      = 0.50*density_factor*density*velocity_mag*area_satellite*cdmean_aerodynamic/mass_satellite
  #"density factor" added by Tomoki Sakai 2023/2/3
  force[4,0] = -fact_aero * velocity[0]
  force[4,1] = -fact_aero * velocity[1]
  force[4,2] = -fact_aero * velocity[2]


  #print(force[1,:],force[2,:],force[3,:])

  # Total
  force[0,:] = force[1,:] + force[2,:] + force[3,:] + force[4,:]


  return force

