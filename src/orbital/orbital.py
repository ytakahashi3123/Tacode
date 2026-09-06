#!/usr/bin/env python3

import numpy as np
import os as os
import sys as sys
from general.general import general
import attitude.attitude as attitude
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
  
  KEY_ATTITUDE_QUATERNION = 'quaternion'
  KEY_ATTITUDE_OMEGA      = 'angular_velocity'

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


  def get_directory_path(self, path_specify, default_path, manual_path):
    script_directory = os.path.dirname(os.path.realpath(__file__))
    if path_specify == 'auto' or path_specify == 'default':
      directory_path = script_directory + default_path
    elif path_specify == 'manual':
      directory_path = manual_path
    else :
      directory_path = script_directory + default_path
    return directory_path


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
      sys.exit(1)

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


  def initial_settings_attitude(self, config, coordinate_dict, velocity_dict):
    #
    # 姿勢（6 自由度）の初期条件。attitude.flag_attitude が False（既定）なら
    # None を返し、呼び出し側は従来どおりの質点 3 自由度計算になる。
    #
    #  初期姿勢     : ローカル水平系（地心 NED）基準の 3-2-1 オイラー角 [ヨー, ピッチ, ロール], deg
    #  初期角速度   : 機体軸成分の [p, q, r], deg/s。ECEF に対する角速度として解釈する
    #                （内部状態は慣性系に対する角速度なので、ここで変換する）
    #
    section = attitude.get_setting(config, 'attitude', None)
    if not bool( attitude.get_setting(section, 'flag_attitude', False) ) :
      return None

    print('Setting initial conditions of the attitude (6-DOF)')

    if not config['computational_setup']['flag_initial'] :
      print('Restart of the attitude computation is not implemented in this version.')
      print('Program stopped')
      sys.exit(1)

    section_initial = attitude.get_setting(config, 'initial_settings', None)
    euler_init = np.array( attitude.get_setting(section_initial, 'attitude', [0.0, 0.0, 0.0]), dtype=float )*self.deg2rad
    omega_init = np.array( attitude.get_setting(section_initial, 'angular_velocity', [0.0, 0.0, 0.0]), dtype=float )*self.deg2rad

    # 初期位置の地心経度・緯度でローカル水平系を決める
    coordinate_polar = coordinate_dict[self.KEY_COORD_POLAR][0]
    latitude  = coordinate_polar[1]
    longitude = coordinate_polar[2]

    quaternion_init = attitude.get_quaternion_from_euler(euler_init[0], euler_init[1], euler_init[2], longitude, latitude)

    # ECEF 基準の角速度を慣性系基準に直す
    matrix_be           = attitude.quaternion_to_matrix(quaternion_init)
    omega_inertial_init = np.array(omega_init) \
                        + attitude.get_earth_rate_body(config['planet']['rotation_rate'], matrix_be)

    print('--Euler angle (yaw, pitch, roll, deg.):', euler_init*self.rad2deg)
    print('--Quaternion (ECEF to body):', quaternion_init)

    attitude_dict = {self.KEY_ATTITUDE_QUATERNION: [quaternion_init],
                     self.KEY_ATTITUDE_OMEGA     : [omega_inertial_init]}

    return attitude_dict


  def output_restart(self, config, iteration, time_elapsed, coordinate, velocity, attitude_dict=None):

    dir_restart      = config['restart_process']['directory_output']
    file_restart     = config['restart_process']['file_restart']
    flag_time_series = config['restart_process']['flag_time_series']   # --True: stored individually as time series, False: stored by overwriting
    digid_step       = config['restart_process']['digid_step']
    # 出力間隔。従来は読まれていなかったが、姿勢計算のように時間刻みが細かいと
    # 全ステップ書き出すとファイルが巨大になるので有効にした。既定の 1 では従来と同じ
    frequency_output = config['restart_process'].get('frequency_output', 1)
    if frequency_output < 1 :
      frequency_output = 1

    if flag_time_series :
      addfile = '_'+str(iteration).zfill(digid_step)
      filename_tmp = dir_restart + '/' + self.split_file(file_restart,addfile,'.')
    else :
      filename_tmp = dir_restart + '/' + file_restart

    print('Writing restart data...:',filename_tmp)

    # 6 自由度計算のときだけ姿勢（クォータニオンと慣性系基準の角速度）を続けて書く
    flag_attitude = attitude_dict is not None
    if flag_attitude :
      quaternion_list = attitude_dict[self.KEY_ATTITUDE_QUATERNION]
      omega_list      = attitude_dict[self.KEY_ATTITUDE_OMEGA]

    # File open  
    file = open(filename_tmp, "w")
    file.write('# Restart data (ECEF, cartesian system)' + self.newline_code)
    if flag_attitude :
      file.write('# X, Y, Z, U, V, W, q0, q1, q2, q3, P, Q, R' + self.newline_code)
      file.write('# --Quaternion: ECEF to body. Angular velocity: body axes, relative to the inertial frame' + self.newline_code)
    if frequency_output != 1 :
      file.write('# Output frequency: ' + str(frequency_output) + self.newline_code)
    file.write('# Iteration, Elapsed time' + self.newline_code)
    file.write('# '+ str(iteration) + self.blank_code + str(time_elapsed) + self.newline_code)
    for n in range(0,iteration+1):
      # 最終ステップは間隔によらず必ず書く
      if n%frequency_output != 0 and n != iteration :
        continue
      str_coord = str(coordinate[n][0])+ self.blank_code +str(coordinate[n][1]) + self.blank_code + str(coordinate[n][2]) 
      vel_coord = str(velocity[n][0])  + self.blank_code +str(velocity[n][1])   + self.blank_code + str(velocity[n][2]) 
      str_attitude = ''
      if flag_attitude :
        for m in range(0,4):
          str_attitude = str_attitude + self.blank_code + str(quaternion_list[n][m])
        for m in range(0,3):
          str_attitude = str_attitude + self.blank_code + str(omega_list[n][m])
      file.write( str_coord + self.blank_code + vel_coord + str_attitude + self.newline_code)
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


  def get_attitude_output(self, config, coordinate_cartesian, velocity_cartesian, quaternion, omega_inertial, rotation_rate_planet):
    #
    # Tecplot に書き出す姿勢まわりの量を組み立てる。
    #  クォータニオン（ECEF -> 機体）、ローカル水平基準のオイラー角、
    #  ECEF に対する角速度（機体軸）、空力角。
    # オイラー角と空力角は保存しておらず、ここで毎回作り直す。
    # （そうすることで、位置・速度・姿勢の配列長のずれが入り込む余地が無くなる）
    #
    coordinate_polar = coordinate_system.set_angle_polar(config, coordinate_cartesian)
    latitude  = coordinate_polar[1]
    longitude = coordinate_polar[2]

    yaw, pitch, roll = attitude.get_euler_angle(quaternion, longitude, latitude)

    matrix_be      = attitude.quaternion_to_matrix(quaternion)
    omega_relative = attitude.get_omega_relative(omega_inertial, rotation_rate_planet, matrix_be)
    velocity_body  = np.dot(matrix_be, np.array(velocity_cartesian))

    alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(velocity_body)

    value_output = list(quaternion) \
                 + [yaw*self.rad2deg, pitch*self.rad2deg, roll*self.rad2deg] \
                 + list(np.array(omega_relative)*self.rad2deg) \
                 + [alpha*self.rad2deg, beta*self.rad2deg, alpha_total*self.rad2deg]

    return self.blank_code.join([str(value) for value in value_output])


  def output_tecplot(self, config, iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict, attitude_dict=None):

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

      # Attitude (6-DOF)
      flag_attitude = attitude_dict is not None
      if flag_attitude :
        quaternion_list      = attitude_dict[self.KEY_ATTITUDE_QUATERNION]
        omega_list           = attitude_dict[self.KEY_ATTITUDE_OMEGA]
        rotation_rate_planet = config['planet']['rotation_rate']

      # Config,
      filename_tmp = config['post_process']['directory_output'] + '/' + config['post_process']['tecplot']['filename_output']
      dt = config['time_integration']['timestep_constant']
      frequency_output = config['post_process']['tecplot']['frequency_output']

      # Output
      print('Writing Tecplot file... ', filename_tmp)
      file = open(filename_tmp, "w")
      file.write('# Tecplot data: Tacode' + self.newline_code)
      variables_tmp = 'Variables = Time[s],X[km],Y[km],Z[km],Long[deg.],Lati[deg.],Alti[km],Upl[m/s],Vpl[m/s],Wpl[m/s],VelplAbs[m/s],Dens[kg/m3],Temp[K],Kn'
      if flag_attitude :
        variables_tmp = variables_tmp + ',q0,q1,q2,q3,Yaw[deg.],Pitch[deg.],Roll[deg.],P[deg/s],Q[deg/s],R[deg/s],AoA[deg.],Sideslip[deg.],AoAtotal[deg.]'
      file.write(variables_tmp + self.newline_code)
      # 実際に出力する点数（最終ステップを含む）。zone ヘッダの i= と実点数は一致していなければならない
      num_output = iteration//frequency_output + 1
      file.write('zone t=time i= '+str(num_output)+' f=point' + self.newline_code )
      for n in range(0,iteration+1):
        if n%frequency_output == 0 :
          time_tmp = float(n)*dt
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
          str_attitude = ''
          if flag_attitude :
            str_attitude = self.blank_code + self.get_attitude_output(config, coordinate_cart[n], velocity_cart[n],
                                                                      quaternion_list[n], omega_list[n], rotation_rate_planet)
          file.write( str_time  + str_coord_cart + str_coord_geod + str_veloc_pola + str_traj + str_attitude + self.newline_code)
      file.close()

    return
