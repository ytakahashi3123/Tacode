#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/05/31

import numpy as np
from general.general import general
import coordinate_system.coordinate_system as coordinate_system

class orbital(general):

  # Constants
  file_control_default = "config.yml"
  deg2rad = np.pi/180.0
  rad2deg = 180.0/np.pi
  m2km    = 1.e-3
  km2m    = 1.e+3

  unit_convert_geoditic     = [rad2deg,rad2deg,m2km]
  unit_convert_geoditic_inv = [deg2rad,deg2rad,km2m]

  KEY_COORD_GEODETIC  = 'geodetic'
  KEY_COORD_CARTESIAN = 'cartesian'
  KEY_COORD_POLAR     = 'polar'
  KEY_VELOC_CARTESIAN = 'cartesian'
  KEY_VELOC_POLAR     = 'polar'
  
  KEY_TRAJECTORY_DENSITY     = 'density'
  KEY_TRAJECTORY_TEMPERATURE = 'temperature'
  KEY_TRAJECTORY_KNUDSEN     = 'knudsen'

  newline_code='\n'
  blank_code=' '

  def __init__(self):
    print("Calling class: orbital")

    return


  def make_directory_output(self, config):
    # Make directory

    print('Making directories for output...')

    dir_restart  = config['restart_process']['directory_output']
    self.make_directory(dir_restart)

    dir_result  = config['post_process']['directory_output']
    self.make_directory(dir_result)

    return


  def initial_settings(self, config):

    print('Setting initial conditions')

    flag_initial = config['computational_setup']['flag_initial']
    timestep     = config['time_integration']['timestep_constant']

    # Initial start
    if flag_initial :
      print('--From initial condition set in control file')
      
      iteration    = 0
      time_elapsed = 0.0

      # Initial coordinate in the geodetic coordinate: 0:Long.(deg) 1: Lat.(deg), 2:Alt.(m)
      coord_init = np.array( config['initial_settings']['coordinate'] )
      # --Unit convert
      #coord_init[0] # Longitude, deg. --> rad.
      #coord_init[1] # Latitude, deg. --> rad.    
      #coord_init[2] # Altitude, km --> m
      coord_init = np.multiply(coord_init, self.unit_convert_geoditic_inv)
      # Initial velocity
      veloc_init = np.array( config['initial_settings']['velocity'] )

      coordinate_geodetic = [coord_init]
      velocity_polar      = [veloc_init]

    else :
      print('--from restart file')

      print('Restart routine is Not implemented in this version.')
      print('Please check: flag_initial.')
      print('Program stopped')
      exit()

      coordinate_geodetic = []
      velocity_polar      = []

      iteration, time_elapsed, coordinate_geodetic, velocity_geodetic = self.read_restart(config, coordinate_geodetic, velocity_polar)

    print('--Iteration: ',iteration)

    
    # Reconstruction
    coordinate_cartesian = []
    coordinate_polar     = []
    velocity_cartesian   = []
    #time_elapsed         = []
    for n in range(0,iteration+1):
      # Initial (restart) coordinate and velocity are given by those in geodetic coordinte
      # Those values are converted to in cartesian coordinate
      cartesian_coord_tmp = coordinate_system.convert_geodetic_cartesian(config, coordinate_geodetic[n])
      coordinate_cartesian.append( np.array( cartesian_coord_tmp ) )

      # Set parameters in polar coordinate from cartesian coordinate
      # polar_coord: [radius, 極座標における緯度(beta), 極座標における経度(alpha)]
      polar_coord_tmp = coordinate_system.set_angle_polar(config, cartesian_coord_tmp)
      coordinate_polar.append( np.array( polar_coord_tmp) )

      # Velocity
      longitude = polar_coord_tmp[2]
      latitude  = polar_coord_tmp[1]
      cartesian_veloc_tmp = coordinate_system.convert_polar_carteasian(config, velocity_polar[n] ,longitude, latitude)
      velocity_cartesian.append( np.array( cartesian_veloc_tmp) )

      # Time 
      #time_elapsed.append(n*timestep)

    coordinate_dict = {'geodetic': coordinate_geodetic, 'cartesian':coordinate_cartesian, 'polar':coordinate_polar}
    velocity_dict   = {'cartesian': velocity_cartesian, 'polar':velocity_polar}
    

    # Trajectory properties
    density_trajectory     = []
    temperature_trajectory = []
    knudsen_trajectory     = []
    #for n in range(0,iteration+1):

    trajectory_dict = {'density':  density_trajectory, 'temperature': temperature_trajectory, 'knudsen': knudsen_trajectory}

    return iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict


  def output_restart(self, config, iteration, time_elapsed, coordinate, velocity):

    dir_restart      = config['restart_process']['directory_output']
    file_restart     = config['restart_process']['file_restart']
    flag_time_series = config['restart_process']['flag_time_series']   # --True: stored individually as time series, False: stored by overwriting
    digid_step       = config['restart_process']['digid_step']

    if flag_time_series :
      addfile = '_'+str(iteration).zfill(digid_step)
      filename_tmp = dir_restart + '/' + self.split_file(file_restart,addfile,'.')
    else :
      filename_tmp = dir_restart + '/' + file_restart

    print('Writing restart data...:',filename_tmp)
    
    # File open  
    file = open(filename_tmp, "w")
    file.write('# Restart data (ECEF, cartesian system)' + self.newline_code)
    file.write('# Iteration, Elapsed time' + self.newline_code)
    file.write('# '+ str(iteration) + self.blank_code + str(time_elapsed) + self.newline_code)
    for n in range(0,iteration):
      str_coord = str(coordinate[n][0])+ self.blank_code +str(coordinate[n][1]) + self.blank_code + str(coordinate[n][2]) 
      vel_coord = str(velocity[n][0])  + self.blank_code +str(velocity[n][1])   + self.blank_code + str(velocity[n][2]) 
      file.write( str_coord + self.blank_code + vel_coord + self.newline_code)
    file.close()

    return


  def read_restart(self, config, coordinate, velocity):

    dir_restart      = config['restart_process']['directory_output']
    file_restart     = config['restart_process']['file_restart']
    flag_time_series = config['restart_process']['flag_time_series']   # --True: stored individually as time series, False: stored by overwriting
    digid_step       = config['restart_process']['digid_step']
    restart_step     = config['restart_process']['restart_step']

    if flag_time_series :
      addfile = '_s'+str(restart_step).zfill(digid_step)
      filename_tmp = dir_restart + '/' + self.split_file(file_restart,addfile,'.')
    else :
      filename_tmp = dir_restart + '/' + file_restart
    
    # Open file
    with open(filename_tmp) as f:
      lines = f.readlines()
    # リストとして取得 
    lines_strip = [line.strip() for line in lines]
    # Iteration 
    words = lines_strip[2].split()
    iteration    = int(words[1])
    time_elapsed = float(words[2])

    # Reading data
    for n in range(3,len(lines_strip )):
      words = lines_strip[n].split()
      coordinate.append( float(words[0]), float(words[1]), float(words[2]) )
      velocity.append( float(words[3]), float(words[4]), float(words[5]) )
    
    f.close()

    return iteration, time_elapsed, coordinate, velocity


  def output_tecplot(self, config, iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict):

    if config['post_process']['tecplot']['flag_output'] :

      # Position and velocity
      coordinate_cart = coordinate_dict['cartesian']
      coordinate_geod = coordinate_dict['geodetic']
      velocity_cart   = velocity_dict['cartesian']
      velocity_pola   = velocity_dict['polar']

      # Trajectory properties
      density_traj     = trajectory_dict['density']
      temperature_traj = trajectory_dict['temperature']
      knudsen_traj     = trajectory_dict['knudsen']

      # Config,
      filename_tmp = config['post_process']['directory_output'] + '/' + config['post_process']['tecplot']['filename_output']
      dt = config['time_integration']['timestep_constant']
      frequency_output = config['post_process']['tecplot']['frequency_output']

      # Output
      print('Writing Tecplot file... ', filename_tmp)
      file = open(filename_tmp, "w")
      file.write('# Tecplot data: Tacode' + self.newline_code)
      file.write('Variables = Time[s],X[km],Y[km],Z[km],Long[deg.],Lati[deg.],Alti[km],Upl[m/s],Vpl[m/s],Wpl[m/s],VelplAbs[m/s],Dens[kg/m3],Temp[K],Kn' + self.newline_code)
      file.write('zone t=time i= '+str(iteration/frequency_output)+' f=point' + self.newline_code )
      for n in range(0,iteration):
        if n%frequency_output == 0 :
          time_tmp = float(n)*dt
          velo_abs = np.linalg.norm(velocity_pola)
          str_time = str( time_tmp ) + self.blank_code
          str_coord_cart = ''
          str_coord_geod = ''
          str_veloc_pola = ''
          for m in range(0,3):
            str_coord_cart = str_coord_cart + str(coordinate_cart[n][m]*self.m2km)   + self.blank_code 
            str_coord_geod = str_coord_geod + str(coordinate_geod[n][m]*self.unit_convert_geoditic[m]) + self.blank_code 
            str_veloc_pola = str_veloc_pola + str(velocity_pola[n][m]) + self.blank_code 
          #str_coord_cart = str(coordinate_cart[n][0]*self.m2km)   + self.blank_code +str(coordinate_cart[n][1]*self.m2km)    + self.blank_code + str(coordinate_cart[n][2]*self.m2km) 
          #str_coord_geod = str(coordinate_geod[n][0]*self.rad2deg)+ self.blank_code +str(coordinate_geod[n][1]*self.rad2deg) + self.blank_code + str(coordinate_geod[n][2]*self.m2km) 
          #str_veloc_pola = str(velocity_pola[n][0])               + self.blank_code +str(velocity_pola[n][1])                + self.blank_code + str(velocity_pola[n][2]) + str(velo_abs)
          str_veloc_pola = str_veloc_pola + str(np.linalg.norm(velocity_pola[n])) + self.blank_code
          str_traj       = str(density_traj[n]) + self.blank_code + str(temperature_traj[n]) + self.blank_code + str(knudsen_traj[n]) 
          #
          file.write( str_time  + str_coord_cart + str_coord_geod + str_veloc_pola + str_traj + self.newline_code)
      file.close()

    return



  def routine_postprocess(self, config, iteration, meshnode_dict, meshelem_dict, metrics_dict, gas_property_dict, var_primitiv):

    # Tecplot (not implemented yet)
    #if ( config['post_process']['flag_output_tecplot'] ):
    #  self.output_tecplot(config, dimension_dict, grid_list, geom_dict, iteration, var_primitiv, var_primitiv_bd, var_gradient, var_limiter)

    # VTK
    if config['post_process']['flag_output_vtk'] :
      # -- Setting variables
      scalar_dict, vector_dict = self.prepare_postprocess(config, gas_property_dict, metrics_dict, var_primitiv)

      # -- Output
      self.write_gmsh2vtk(config, iteration, meshnode_dict, meshelem_dict, scalar_dict=scalar_dict, vector_dict=vector_dict)

    return


  def prepare_postprocess(self, config, gas_property_dict, metrics_dict, var_primitiv):

    specfic_heat_ratio = gas_property_dict['specfic_heat_ratio']

    density      = var_primitiv[0,:]
    temperature  = var_primitiv[4,:]
    pressure     = var_primitiv[5,:]
    velocity     = [ var_primitiv[1,:],var_primitiv[2,:],var_primitiv[3,:] ]
    velocity_mag = np.sqrt( velocity[0]**2 + velocity[1]**2 + velocity[2]**2 )
    speedofsound = self.get_speedofsound(specfic_heat_ratio, density, pressure)
    machnumber   = velocity_mag/speedofsound

    #volume = metrics_dict['volume_cell']
    #cellcenter_list = metrics_dict['coord_cellcenter']
    #cellcenter  = [ cellcenter_list[0],cellcenter_list[1],cellcenter_list[2] ]

    scalar_dict = { 'density': density, 'temperature': temperature, 'pressure': pressure, 'machnumber': machnumber }
    #scalar_dict = { 'density': density, 'temperature': temperature, 'pressure': pressure, 'machnumber': machnumber, 'volume':volume }
    vector_dict = { 'velocity':velocity }
    #vector_dict = { 'velocity':velocity, 'cellcenter':cellcenter }

    return scalar_dict, vector_dict

