#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/05/31

import numpy as np
import atmosphere.atmosphere as atmosphere
import coordinate_system.coordinate_system as coordinate_system
import force_term.force_term as force_term
import satellite.satellite as satellite
from orbital.orbital import orbital

# Constants
one_sixth = 1.0/6.0
fact_rk   = [0.5, 0.5, 1.0, 0.0]
fact_up   = [1.0, 2.0, 2.0, 1.0]



def solve_equation_motion(config, iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict, atmosphere_dict, aerodynamic_dict):
  
  print( 'Start calculation of equation of motion...' )

  # Mass-point properties
  mass_satellite      = config['satellite']['mass']
  area_satellite      = config['satellite']['characteristic_area']
  length_satellite    = config['satellite']['characteristic_length']
  density_factor      = config['initial_settings']['density_factor'][0] #Added by Tomoki Sakai 2023/2/3

  # Atmosphere model
  kind_atmosphere_model = config['atmosphere']['kind_atmosphere_model']
  if kind_atmosphere_model == 'constant' :
    density_atmosphere     = config['atmosphere']['density']
    temperature_atmosphere = config['atmosphere']['temperature']
  elif kind_atmosphere_model == 'fileread' :
    altitude_atmosphere    = atmosphere_dict['Height']
    density_atmosphere     = atmosphere_dict['Mass_density']
    temperature_atmosphere = atmosphere_dict['Temperature_neutral']
    knudsen_atmosphere     = atmosphere_dict['Knudsen_number']
  else :
    print('kind_atmosphere_model in config is incorrect.')
    print('Program stopped.')
    exit()

  # Aerodynamic model
  kind_aerodynamic_model = config['satellite']['kind_aerodynamic_model']
  if kind_aerodynamic_model == 'constant' :
    cdmean_aerodynamic   = config['satellite']['drag_coefficient']
  elif kind_aerodynamic_model == 'fileread' :
    knudsen_aerodynamic  = aerodynamic_dict['Knudsen_number']
    cdmean_aerodynamic   = aerodynamic_dict['CD_mean']
    altitude_aerodynamic = aerodynamic_dict['Altitude']
  else :
    print('kind_aerodynamic_model in config is incorrect.')
    print('Program stopped.')
    exit()

  # Force setting
  force = force_term.force_initialsettings(config)

  # Calculation parameter settings
  kind_time_scheme = config['time_integration']['kind_time_scheme']
  delta_time       = config['time_integration']['timestep_constant']

  # Position and velocity
  coordinate_cart = coordinate_dict['cartesian']
  coordinate_geod = coordinate_dict['geodetic']
  velocity_cart   = velocity_dict['cartesian']
  velocity_pola   = velocity_dict['polar']

  # Trajectory properties
  density_traj     = trajectory_dict['density']
  temperature_traj = trajectory_dict['temperature']
  knudsen_traj     = trajectory_dict['knudsen']


  print('Time elapsed (s):, Longitude (deg.), Latitude (deg.), Altitude (km), Velocity Mag. (m/s)')

  # Main routine
  while time_elapsed <= config['computational_setup']['time_elapsed_maximum'] :

    coord_tmp = coordinate_cart[iteration]
    veloc_tmp = velocity_cart[iteration]

    v_res = 0.0
    r_res = 0.0

    if kind_time_scheme == 'explicit_euler' :
    # Euler explicit 

      # Cartesian --> Geodetic system
      coord_geod = coordinate_system.convert_cartesian_geodetic(config, coord_tmp)

      # Recalculate atmosphere status
      altitude_tmp = coord_geod[2] * orbital.m2km
      if kind_atmosphere_model == 'fileread' :
        density, temperature, knudsen = atmosphere.get_atmosphere_property(altitude_tmp, altitude_atmosphere, density_atmosphere, temperature_atmosphere, knudsen_atmosphere)

      # Aerodynamic coefficient
      if kind_aerodynamic_model == 'fileread' :
        cdmean = satellite.get_aerodynamic_coefficient(knudsen, knudsen_aerodynamic, cdmean_aerodynamic)

      # Calculate force
      force = force_term.force_routine(config, coord_tmp, veloc_tmp,   \
                                       mass_satellite, area_satellite, \
                                       cdmean, density_factor, density, force)
      #"density factor" added by Tomoki Sakai 2023/2/3
      force_total = force[0,:]
 
      # Update solution and Calculate residual
      coord_tmp, veloc_tmp, \
      r_res, v_res =  solve_eulerexplicit(delta_time, mass_satellite, force_total, \
                                          coord_tmp, veloc_tmp, \
                                          r_res, v_res)


    elif kind_time_scheme == 'runge_kutta' :
    # 4th stage Runge-Kutta method
    
      # Reinitialization
      r_virtual = coord_tmp
      v_virtual = veloc_tmp
      r_virprev = r_virtual
      v_virprev = v_virtual

      # R.K. 1st - 4th stages
      for m in range(0, 4):
        # Cartesian --> Geodetic system
        coord_geod = coordinate_system.convert_cartesian_geodetic(config, r_virtual)

        # Recalculate atmosphere status
        altitude_tmp = coord_geod[2] * orbital.m2km
        if kind_atmosphere_model == 'fileread' :
          density, temperature, knudsen = atmosphere.get_atmosphere_property(altitude_tmp, altitude_atmosphere, density_atmosphere, temperature_atmosphere, knudsen_atmosphere)

        # Aerodynamic coefficient
        if kind_aerodynamic_model == 'fileread' :
          cdmean = satellite.get_aerodynamic_coefficient(knudsen, knudsen_aerodynamic, cdmean_aerodynamic)

        # External force (but "Delta V" is given here)
        # --Not yet
        # --Convert cartesian coordinate
        #set_angle_longlat_cartesian(r_virtual, longitude_tmp, latitude_tmp)
        #exchanege_axis_xyz2zxy(externalforce_tmp)
        #convert_coordinate_ry( latitude_tmp ,externalforce_tmp)
        #convert_coordinate_rz(-longitude_tmp,externalforce_tmp)

        # Calculate force
        force = force_term.force_routine(config, r_virtual, v_virtual, \
                                         mass_satellite, area_satellite, \
                                         cdmean, density_factor, density, force)
        #"density factor" added by Tomoki Sakai 2023/2/3
        force_total = force[0,:]
 
        # Update solution and Calculate residual
        coord_tmp, veloc_tmp, \
        r_virtual, v_virtual, \
        r_res, v_res = solve_rungekutta(m, delta_time, mass_satellite, \
                                        force_total, \
                                        coord_tmp, veloc_tmp, \
                                        r_virtual, v_virtual, \
                                        r_virprev, v_virprev, \
                                        r_res, v_res)

    else :
      print('kind_time_scheme in config is incorrect.')
      print('Program stopped.')
      exit()

    # Residuals
    v_res = np.sqrt(v_res)
    r_res = np.sqrt(r_res)

    # Convert
    coord_geodetic = coordinate_system.convert_cartesian_geodetic(config, coord_tmp)
    coord_polar = coordinate_system.set_angle_polar(config, coord_tmp)
    angle_beta  = coord_polar[1]
    angle_alpha = coord_polar[2]
    veloc_polar = coordinate_system.convert_carteasian_polar(config, veloc_tmp, angle_alpha, angle_beta)

    # Update
    coordinate_cart.append( coord_tmp )
    coordinate_geod.append( coord_geodetic )
    #coordinate_pola.append( coord_polar )
    velocity_cart.append( veloc_tmp )
    velocity_pola.append( veloc_polar )

    density_traj.append( density )
    temperature_traj.append( temperature )
    knudsen_traj.append( knudsen )

    time_elapsed = time_elapsed + delta_time
    iteration    = iteration + 1

    coord_geodetic_display = np.multiply(coord_geodetic, orbital.unit_convert_geoditic)
    veloc_polar_mag = np.linalg.norm(veloc_polar)
    #print("Time elapsed:", time_elapsed, 'Longitude:', coord_geodetic_display[0], 'Latitude:', coord_geodetic_display[1], 'Altitude:', coord_geodetic_display[2], 'Velocity (Mag.)',veloc_polar_mag )
    print('{:.1f}'.format(time_elapsed)+', '+'{:.3f}'.format(coord_geodetic_display[0])+', '+'{:.3f}'.format(coord_geodetic_display[1])+', '+'{:.3f}'.format(coord_geodetic_display[2])+', '+'{:.3f}'.format(veloc_polar_mag) )

    if coord_geodetic_display[2] <= 0.0 :
      break

  return iteration, coordinate_dict, velocity_dict, trajectory_dict


def solve_eulerexplicit(dt, mass, force, coord_tmp, veloc_tmp, r_res, v_res):

  # Update solution and Calculate residual
  dv        = ( force )*dt
  dr        = ( veloc_tmp )*dt
  veloc_tmp = veloc_tmp + dv
  coord_tmp = coord_tmp + dr

  # Discrepancies from the previous step
  v_res = v_res + np.linalg.norm(dv)**2
  r_res = r_res + np.linalg.norm(dr)**2

  return coord_tmp, veloc_tmp, r_res, v_res


def solve_rungekutta(m, dt, mass, force, coord_tmp, veloc_tmp, r_virtual, v_virtual, r_virprev, v_virprev, r_res, v_res):

  # Update solution and Calculate residual
  kv_rk = ( force )*dt
  kr_rk = ( v_virtual )*dt

  # Update solution
  dv        = fact_up[m]*one_sixth*kv_rk
  dr        = fact_up[m]*one_sixth*kr_rk
  veloc_tmp = veloc_tmp + dv
  coord_tmp = coord_tmp + dr

  # Temporary variables
  v_virtual = v_virprev + fact_rk[m]*kv_rk
  r_virtual = r_virprev + fact_rk[m]*kr_rk

  # Discrepancies from the previous step
  v_res = v_res + np.linalg.norm(dv)**2
  r_res = r_res + np.linalg.norm(dr)**2
 
  return coord_tmp, veloc_tmp, r_virtual, v_virtual, r_res, v_res