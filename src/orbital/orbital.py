#!/usr/bin/env python3

import numpy as np
import os as os
import sys as sys
from general.general import general, vector_norm
import attitude.attitude as attitude
import coordinate_system.coordinate_system as coordinate_system
import epoch.epoch as epoch_module
import wind.wind as wind_module

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

  # read_restart が返す再開点の状態
  KEY_RESTART_ITERATION  = 'iteration'
  KEY_RESTART_TIME       = 'time_elapsed'
  KEY_RESTART_COORDINATE = 'coordinate'
  KEY_RESTART_VELOCITY   = 'velocity'
  KEY_RESTART_QUATERNION = 'quaternion'
  KEY_RESTART_OMEGA      = 'angular_velocity'
  KEY_RESTART_EPOCH      = 'epoch'
  KEY_RESTART_FILE       = 'filename'

  # リスタートファイルのヘッダ。行数は姿勢・エポック・出力間隔の有無で変わるので、
  # 行番号ではなく内容で探す（6 自由度の restart.dat では反復回数の行が 5 行目に来る）
  MARKER_RESTART_EPOCH     = '# Epoch (UTC):'
  MARKER_RESTART_ITERATION = '# Iteration, Elapsed time'

  # データ行の列数。3 自由度は X,Y,Z,U,V,W、6 自由度はそれに q0..q3,P,Q,R が続く
  NUM_COLUMN_RESTART_3DOF = 6
  NUM_COLUMN_RESTART_6DOF = 13

  # クォータニオンの単位長からのずれの許容値（書き出しは丸めていないので本来は 1e-16 の桁）
  TOLERANCE_QUATERNION_NORM = 1.e-6

  newline_code='\n'
  blank_code=' '

  def __init__(self):
    print("Calling class: orbital")

    # リスタートで読んだ再開点の状態（initial_settings が入れ、
    # initial_settings_attitude が姿勢の分を取り出す）
    self.restart_state = None

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


  def initial_settings(self, config, epoch_dict=None):

    print('Setting initial conditions')

    flag_initial = config['computational_setup']['flag_initial']

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

      # 初期条件は測地系・地心ローカル系で与えられるので、ECEF 直交系に直す
      coordinate_cartesian = []
      coordinate_polar     = []
      velocity_cartesian   = []
      for n in range(0,iteration+1):
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

    else :
      print('--From restart file')

      # 再開点の状態。姿勢の分は initial_settings_attitude が受け取る
      # （返り値の数を増やさないため、解析結果をインスタンスに置く）
      restart_state      = self.read_restart(config, epoch_dict)
      self.restart_state = restart_state

      # 使うのは最終行＝再開点だけで、履歴は復元しない。
      # 出力の配列は再開点から始まるので反復回数は 0 に戻し、経過時間だけを継ぐ
      # （時刻列は再開時刻から連続し、風のテーブルも同じ時刻で引かれる）
      print('--Resuming at iteration {:d} of the previous run, elapsed time {:g} s'.format(
            restart_state[self.KEY_RESTART_ITERATION], restart_state[self.KEY_RESTART_TIME]))

      iteration    = 0
      time_elapsed = restart_state[self.KEY_RESTART_TIME]

      # リスタートファイルが持っているのは ECEF 直交系そのものなので、測地系を経由しない。
      # 往復（直交 -> 測地 -> 直交）を通すと最後の桁が動き、続きの計算が
      # 元の計算の続きとビット一致しなくなる
      coordinate_cartesian = [ restart_state[self.KEY_RESTART_COORDINATE] ]
      velocity_cartesian   = [ restart_state[self.KEY_RESTART_VELOCITY] ]

      polar_coord_tmp     = coordinate_system.set_angle_polar(config, coordinate_cartesian[0])
      coordinate_polar    = [ np.array( polar_coord_tmp ) ]
      coordinate_geodetic = [ np.array( coordinate_system.convert_cartesian_geodetic(config, coordinate_cartesian[0]) ) ]
      velocity_polar      = [ np.array( coordinate_system.convert_carteasian_polar(config, velocity_cartesian[0],
                                                                                   polar_coord_tmp[2], polar_coord_tmp[1]) ) ]

      # 姿勢を持つファイルを 3 自由度で読んだときは、捨てていることを知らせる
      section_attitude = attitude.get_setting(config, 'attitude', None)
      if restart_state[self.KEY_RESTART_QUATERNION] is not None and \
         not bool( attitude.get_setting(section_attitude, 'flag_attitude', False) ) :
        print('Warning: the restart file holds attitude data but attitude.flag_attitude is False.')
        print('--The run continues as a 3-DOF computation and the attitude in the file is dropped.')

    print('--Iteration: ',iteration)

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
      # リスタート。クォータニオンと角速度はファイルの値をそのまま使う。
      # ファイルの角速度は「機体軸成分・慣性系基準」で、内部状態と同じ規約なので変換しない
      # （config で与える初期角速度が ECEF 基準なのとは違う。output_restart の
      #   ヘッダにもそう書いてある）
      restart_state = self.restart_state
      if restart_state is None or restart_state[self.KEY_RESTART_QUATERNION] is None :
        print('The restart file has no attitude data, but attitude.flag_attitude is True.')
        if restart_state is not None :
          print('--File:', restart_state[self.KEY_RESTART_FILE])
        print('--A 6-DOF run must restart from a restart file written by a 6-DOF run.')
        print('Program stopped.')
        sys.exit(1)

      quaternion_init     = restart_state[self.KEY_RESTART_QUATERNION]
      omega_inertial_init = restart_state[self.KEY_RESTART_OMEGA]

      print('--Quaternion (ECEF to body):', quaternion_init)
      print('--Angular velocity (body axes, inertial reference, deg/s):', omega_inertial_init*self.rad2deg)

      return {self.KEY_ATTITUDE_QUATERNION: [quaternion_init],
              self.KEY_ATTITUDE_OMEGA     : [omega_inertial_init]}

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


  def output_restart(self, config, iteration, time_elapsed, coordinate, velocity, attitude_dict=None, epoch_dict=None):

    # 出力間隔。従来は読まれていなかったが、姿勢計算のように時間刻みが細かいと
    # 全ステップ書き出すとファイルが巨大になるので有効にした。既定の 1 では従来と同じ
    frequency_output = config['restart_process'].get('frequency_output', 1)
    if frequency_output < 1 :
      frequency_output = 1

    filename_tmp = self.get_restart_filename(config, iteration)

    print('Writing restart data...:',filename_tmp)

    # 6 自由度計算のときだけ姿勢（クォータニオンと慣性系基準の角速度）を続けて書く
    flag_attitude = attitude_dict is not None
    if flag_attitude :
      quaternion_list = attitude_dict[self.KEY_ATTITUDE_QUATERNION]
      omega_list      = attitude_dict[self.KEY_ATTITUDE_OMEGA]

    # File open  
    file = open(filename_tmp, "w")
    file.write('# Restart data (ECEF, cartesian system)' + self.newline_code)
    # エポックを与えたときだけ、経過秒がどの UTC を起点にしているかを添える
    if epoch_dict is not None :
      file.write('# Epoch (UTC): ' + epoch_module.get_string(epoch_dict, 0.0) + self.newline_code)
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


  def get_restart_filename(self, config, step):
    #
    # リスタートファイルの名前。読む側と書く側で規則がずれていると
    # 「書いたファイルが読めない」ので 1 か所に置く
    # （かつて書き出しは restart_0500.dat、読み取りは restart_s0500.dat を探していた）。
    #
    #   flag_time_series: True  -> 1 ステップ 1 ファイル（step を名前に埋める）
    #                     False -> 同じ名前に上書き
    #
    dir_restart  = config['restart_process']['directory_output']
    file_restart = config['restart_process']['file_restart']

    if not config['restart_process']['flag_time_series'] :
      return dir_restart + '/' + file_restart

    addfile = '_' + str(step).zfill(config['restart_process']['digid_step'])

    return dir_restart + '/' + self.split_file(file_restart, addfile, '.')


  def read_restart(self, config, epoch_dict=None):
    #
    # リスタートファイルを読み、再開点（最終行）の状態を返す。
    #
    # ファイルは初期計算からの履歴を全部持っているが、使うのは最終行だけにしてある。
    # 大気量（密度・温度・Kn）は書かれていないので履歴を復元しても軌道量の配列が
    # 揃わず、restart_process.frequency_output で間引いたファイルでは行と反復回数も
    # 対応しないため。再開後の出力は再開点から始まり、時刻は記録された経過時間から続く。
    #
    # ヘッダの行数は姿勢・エポック・出力間隔の有無で変わるので、行番号ではなく
    # 内容で探す。3 自由度か 6 自由度かはデータ行の列数で決める。
    #
    filename_tmp = self.get_restart_filename(config, config['restart_process']['restart_step'])

    print('Reading restart data...:', filename_tmp)

    if not os.path.exists(filename_tmp) :
      print('The restart file is not found:', filename_tmp)
      print('--Check restart_process: directory_output, file_restart, flag_time_series, restart_step.')
      print('Program stopped.')
      sys.exit(1)

    with open(filename_tmp) as f:
      lines_strip = [line.strip() for line in f.readlines()]

    iteration    = None
    time_elapsed = None
    epoch_string = None
    data_rows    = []

    for index in range(0, len(lines_strip)):
      line = lines_strip[index]

      if line.startswith(self.MARKER_RESTART_EPOCH) :
        epoch_string = line[len(self.MARKER_RESTART_EPOCH):].strip()
        continue

      if line.startswith(self.MARKER_RESTART_ITERATION) :
        # 次の行が "# <反復回数> <経過時間>"
        words = lines_strip[index+1].lstrip('#').split() if index+1 < len(lines_strip) else []
        if len(words) < 2 :
          print('The restart file does not carry the iteration and the elapsed time:', filename_tmp)
          print('Program stopped.')
          sys.exit(1)
        try:
          iteration    = int(words[0])
          time_elapsed = float(words[1])
        except ValueError:
          print('The iteration and the elapsed time in the restart file are not numbers:', ' '.join(words))
          print('--File:', filename_tmp)
          print('Program stopped.')
          sys.exit(1)
        continue

      if line.startswith('#') or len(line) == 0 :
        continue

      try:
        data_rows.append( [float(word) for word in line.split()] )
      except ValueError:
        continue

    if iteration is None :
      print('The restart file has no "' + self.MARKER_RESTART_ITERATION + '" line:', filename_tmp)
      print('Program stopped.')
      sys.exit(1)

    if len(data_rows) == 0 :
      print('The restart file has no numerical data:', filename_tmp)
      print('Program stopped.')
      sys.exit(1)

    # 再開点は最終行（出力間隔で間引いても最終ステップは必ず書かれている）
    state = data_rows[-1]
    if len(state) != self.NUM_COLUMN_RESTART_3DOF and len(state) != self.NUM_COLUMN_RESTART_6DOF :
      print('The last line of the restart file has {:d} columns; expected {:d} (3-DOF) or {:d} (6-DOF).'.format(
            len(state), self.NUM_COLUMN_RESTART_3DOF, self.NUM_COLUMN_RESTART_6DOF))
      print('--File:', filename_tmp)
      print('Program stopped.')
      sys.exit(1)

    restart_state = {self.KEY_RESTART_ITERATION : iteration,
                     self.KEY_RESTART_TIME      : time_elapsed,
                     self.KEY_RESTART_COORDINATE: np.array(state[0:3]),
                     self.KEY_RESTART_VELOCITY  : np.array(state[3:6]),
                     self.KEY_RESTART_QUATERNION: None,
                     self.KEY_RESTART_OMEGA     : None,
                     self.KEY_RESTART_EPOCH     : epoch_string,
                     self.KEY_RESTART_FILE      : filename_tmp}

    if len(state) == self.NUM_COLUMN_RESTART_6DOF :
      quaternion = np.array(state[6:10])
      norm_tmp   = vector_norm(quaternion)
      if abs(norm_tmp - 1.0) > self.TOLERANCE_QUATERNION_NORM :
        # 正規化して黙って進めない。単位長から外れているのは書き写しの誤りか、
        # 発散した計算の残骸で、どちらも「続きを計算してよい状態」ではない
        print('The quaternion in the restart file is not of unit length: {:.6e}'.format(norm_tmp))
        print('--File:', filename_tmp)
        print('Program stopped.')
        sys.exit(1)
      restart_state[self.KEY_RESTART_QUATERNION] = quaternion
      restart_state[self.KEY_RESTART_OMEGA]      = np.array(state[10:13])

    self.check_restart_epoch(restart_state, epoch_dict)

    return restart_state


  def check_restart_epoch(self, restart_state, epoch_dict):
    #
    # 経過時間は「エポックからの秒」なので、エポックが変わると同じ経過時間が
    # 別の絶対時刻を指す（風のテーブルの時間内挿がそこで変わる）。
    # 食い違いは黙って進めず、止める。
    #
    epoch_string = restart_state[self.KEY_RESTART_EPOCH]
    epoch_config = epoch_module.get_string(epoch_dict, 0.0) if epoch_dict is not None else None

    if epoch_string is None and epoch_config is None :
      return

    if epoch_string is not None and epoch_config is not None :
      if epoch_string != epoch_config :
        print('The epoch differs between the restart file and the control file.')
        print('--Restart file:', epoch_string)
        print('--Control file:', epoch_config)
        print('--The elapsed time in the restart file is counted from its own epoch.')
        print('Program stopped.')
        sys.exit(1)
      return

    if epoch_string is None :
      print('Warning: the control file sets an epoch but the restart file carries none.')
      print('--The elapsed time of the restart file is taken as seconds from', epoch_config)
    else :
      print('Warning: the restart file carries an epoch (' + epoch_string + ') but the control file sets none.')
      print('--The run continues with the elapsed time alone.')

    return


  def get_attitude_output(self, config, coordinate_cartesian, velocity_cartesian, quaternion, omega_inertial, rotation_rate_planet, velocity_air=None):
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
    # 空力角はソルバーと同じ対気速度から作る。風が無ければ ECEF 速度そのもの
    velocity_aero  = velocity_cartesian if velocity_air is None else velocity_air
    velocity_body  = np.dot(matrix_be, np.array(velocity_aero))

    alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(velocity_body)

    value_output = list(quaternion) \
                 + [yaw*self.rad2deg, pitch*self.rad2deg, roll*self.rad2deg] \
                 + list(np.array(omega_relative)*self.rad2deg) \
                 + [alpha*self.rad2deg, beta*self.rad2deg, alpha_total*self.rad2deg]

    return self.blank_code.join([str(value) for value in value_output])


  def output_tecplot(self, config, iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict, attitude_dict=None, epoch_dict=None, wind_dict=None, time_elapsed_initial=0.0):

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

      # Wind
      # --風は位置（と、いずれ時刻）だけの関数なので保存しておらず、ここで引き直す。
      #   オイラー角や迎角と同じ扱いで、配列長のずれが入り込む余地を無くすため
      flag_wind = wind_dict is not None

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
      # エポックを与えたときだけ、Time[s] の起点を添える。列は増やさない
      # （Tecplot の point 形式に文字列の列を混ぜられないため。各行の UTC は
      #   このエポックに Time[s] を足せば得られる）
      if epoch_dict is not None :
        file.write('# Epoch (UTC): ' + epoch_module.get_string(epoch_dict, 0.0) + self.newline_code)
        file.write('# --UTC of each row is this epoch plus Time[s]' + self.newline_code)
      variables_tmp = 'Variables = Time[s],X[km],Y[km],Z[km],Long[deg.],Lati[deg.],Alti[km],Upl[m/s],Vpl[m/s],Wpl[m/s],VelplAbs[m/s],Dens[kg/m3],Temp[K],Kn'
      # 風を入れたときだけ、風そのものと対気速度の大きさを添える。
      # 対地速度（Upl/Vpl/Wpl）はそのまま残す。両方見えないと風の効きが読めないため
      if flag_wind :
        variables_tmp = variables_tmp + ',WindE[m/s],WindN[m/s],WindU[m/s],VelairAbs[m/s]'
      if flag_attitude :
        variables_tmp = variables_tmp + ',q0,q1,q2,q3,Yaw[deg.],Pitch[deg.],Roll[deg.],P[deg/s],Q[deg/s],R[deg/s],AoA[deg.],Sideslip[deg.],AoAtotal[deg.]'
      file.write(variables_tmp + self.newline_code)
      # 実際に出力する点数（最終ステップを含む）。zone ヘッダの i= と実点数は一致していなければならない
      num_output = iteration//frequency_output + 1
      file.write('zone t=time i= '+str(num_output)+' f=point' + self.newline_code )
      for n in range(0,iteration+1):
        if n%frequency_output == 0 :
          # 時刻はリスタートの再開時刻から連続させる（初期計算では 0）。
          # 累積加算にしないのはソルバーと同じ理由で、丸めを溜めないため
          time_tmp = time_elapsed_initial + float(n)*dt
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
          str_veloc_pola = str_veloc_pola + str(vector_norm(velocity_pola[n])) + self.blank_code
          str_traj       = str(density_traj[n]) + self.blank_code + str(temperature_traj[n]) + self.blank_code + str(knudsen_traj[n]) 
          #
          str_wind     = ''
          velocity_air = None
          if flag_wind :
            wind_local   = wind_module.get_wind_local(coordinate_geod[n], wind_dict, time_tmp)
            velocity_air = wind_module.get_relative_velocity(config, coordinate_cart[n], coordinate_geod[n],
                                                             velocity_cart[n], wind_dict, time_tmp)
            for m in range(0,3):
              str_wind = str_wind + self.blank_code + str(wind_local[m])
            str_wind = str_wind + self.blank_code + str(vector_norm(velocity_air))
          #
          str_attitude = ''
          if flag_attitude :
            str_attitude = self.blank_code + self.get_attitude_output(config, coordinate_cart[n], velocity_cart[n],
                                                                      quaternion_list[n], omega_list[n], rotation_rate_planet,
                                                                      velocity_air)
          file.write( str_time  + str_coord_cart + str_coord_geod + str_veloc_pola + str_traj + str_wind + str_attitude + self.newline_code)
      file.close()

    return
