#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/05/31

import numpy as np
import sys as sys
import atmosphere.atmosphere as atmosphere
import attitude.attitude as attitude
import coordinate_system.coordinate_system as coordinate_system
import force_term.force_term as force_term
import moment_term.moment_term as moment_term
import satellite.satellite as satellite
import wind.wind as wind
from orbital.orbital import orbital
from general.general import get_setting, vector_norm

# Constants
one_sixth = 1.0/6.0
fact_rk   = [0.5, 0.5, 1.0, 0.0]
fact_up   = [1.0, 2.0, 2.0, 1.0]
# 各段が評価する時刻（現在時刻からの dt 倍）。fact_rk が段の終わりに次の仮想状態を
# 作る係数なので、段 m が見ているのは t + fact_time_rk[m]*dt。
# 力が時刻に依らなかったこれまでは要らなかったが、風は時刻に依る
fact_time_rk = [0.0, 0.5, 0.5, 1.0]

# 姿勢の時間刻みの検査を行う間隔, s
# --検査は「まだ警告していない間」ずっと走る。毎ステップ引くと係数表を 2 回引き直すので、
#   6 自由度計算の 9% をそれだけで使う（実測: 12.3 s -> 11.3 s、出力は md5 一致）。
#   周期の見積もりを動かすのは動圧で、動圧が変わる時間尺度は秒の桁なので、
#   それより細かく見直す意味は無い。
INTERVAL_CHECK_TIMESTEP = 1.0



def solve_equation_motion(config, iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict, atmosphere_dict, aerodynamic_dict, attitude_dict=None, wind_dict=None):
  
  print( 'Start calculation of equation of motion...' )

  # Mass-point properties
  mass_satellite      = config['satellite']['mass']
  area_satellite      = config['satellite']['characteristic_area']
  length_satellite    = config['satellite']['characteristic_length']
  density_factor      = config['initial_settings']['density_factor'][0] #Added by Tomoki Sakai 2023/2/3

  # Atmosphere model
  kind_atmosphere_model = config['atmosphere']['kind_atmosphere_model']
  if kind_atmosphere_model == 'constant' :
    # 定数モデルではループ内で更新されないので、ループで使う変数に直接与える
    density     = config['atmosphere']['density']
    temperature = config['atmosphere']['temperature']
    try:
      knudsen   = config['atmosphere']['knudsen']
    except KeyError:
      print('"knudsen" is not given in the atmosphere section of config.')
      print('--It is required for kind_atmosphere_model: constant.')
      print('Program stopped.')
      sys.exit(1)
  elif kind_atmosphere_model == 'fileread' :
    # 内挿・外挿に必要なものは atmosphere_dict に入っているのでそのまま渡す
    pass
  else :
    print('kind_atmosphere_model in config is incorrect.')
    print('Program stopped.')
    sys.exit(1)

  # Aerodynamic model
  kind_aerodynamic_model = config['satellite']['kind_aerodynamic_model']
  if kind_aerodynamic_model == 'constant' :
    # 定数モデルではループ内で更新されないので、ループで使う変数に直接与える
    cdmean = config['satellite']['drag_coefficient']
  elif kind_aerodynamic_model == 'fileread' :
    knudsen_aerodynamic  = aerodynamic_dict['Knudsen_number']
    cdmean_aerodynamic   = aerodynamic_dict['CD_mean']
  else :
    print('kind_aerodynamic_model in config is incorrect.')
    print('Program stopped.')
    sys.exit(1)

  # Force setting
  acceleration = force_term.acceleration_initialsettings(config)

  # Wind setting
  # --wind_dict が None なら大気は地球と共回転し、ECEF 速度がそのまま対気速度になる（従来の仮定）
  flag_wind = wind_dict is not None
  if flag_wind :
    print('--Taking the wind into account (the aerodynamics uses the air-relative velocity)')

  # Attitude (6-DOF) setting
  # --attitude_dict が None なら従来どおりの質点 3 自由度計算
  flag_attitude = attitude_dict is not None
  if flag_attitude :
    print('--Solving the attitude (6-DOF) as well')
    moment, attitude_property = moment_term.moment_initialsettings(config)

    quaternion_list = attitude_dict[attitude.KEY_QUATERNION]
    omega_list      = attitude_dict[attitude.KEY_ANGULAR_VELOCITY]

    rotation_rate_planet = config['planet']['rotation_rate']

    # 迎角依存の係数表が無いと復元モーメントが立たない
    if kind_aerodynamic_model == 'fileread' and len(aerodynamic_dict[satellite.KEY_AOA]) <= 1 :
      print('--Caution: the aerodynamic table has a single angle of attack.')
      print('  The aerodynamic coefficients are then independent of the attitude,')
      print('  so no restoring (static) moment appears.')

    # kind_aerodynamic_model: constant のときの静安定微係数 (1/rad)
    stability_derivative = float(attitude.get_setting(attitude.get_setting(config, 'attitude', None),
                                                      'static_stability_derivative', 0.0))

    # 姿勢振動に対して時間刻みが粗すぎないかの確認（1 度だけ警告する）
    flag_warned_timestep = False
    time_checked_timestep = time_elapsed - INTERVAL_CHECK_TIMESTEP
  else :
    moment               = None
    attitude_property    = None
    quaternion_list      = None
    omega_list           = None
    rotation_rate_planet = None
    stability_derivative = 0.0
    flag_warned_timestep = True
    time_checked_timestep = None

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
  # 経過時間は累積加算せず初期値からのオフセットで求めるため、開始点を保存する
  time_elapsed_initial = time_elapsed
  iteration_initial    = iteration

  # 発散を検知する上限（既定で有効。健全な計算では一度も引っかからない）
  flag_check_divergence, radius_maximum, velocity_maximum \
    = set_divergence_limit(config, coordinate_cart[iteration])

  # 条件を <= にすると time_elapsed == time_elapsed_maximum でも回り、1 ステップ余分に進む
  while time_elapsed < config['computational_setup']['time_elapsed_maximum'] :

    coord_tmp = coordinate_cart[iteration]
    veloc_tmp = velocity_cart[iteration]

    v_res = 0.0
    r_res = 0.0

    if flag_attitude :
      quat_tmp  = quaternion_list[iteration]
      omega_tmp = omega_list[iteration]

    if kind_time_scheme == 'explicit_euler' :
    # Euler explicit 

      # Cartesian --> Geodetic system
      coord_geod = coordinate_system.convert_cartesian_geodetic(config, coord_tmp)

      # Recalculate atmosphere status
      altitude_tmp = coord_geod[2] * orbital.m2km
      if kind_atmosphere_model == 'fileread' :
        density, temperature, knudsen = atmosphere.get_atmosphere_property(altitude_tmp, atmosphere_dict)

      # Aerodynamic coefficient
      if kind_aerodynamic_model == 'fileread' :
        cdmean = satellite.get_aerodynamic_coefficient(knudsen, knudsen_aerodynamic, cdmean_aerodynamic)

      # 対気速度。風が無効なら None のままで、以後は従来どおり ECEF 速度が使われる
      velocity_air = None
      if flag_wind :
        velocity_air = wind.get_relative_velocity(config, coord_tmp, coord_geod, veloc_tmp, wind_dict,
                                                  time_elapsed)

      # 出力用の大気量
      # --オイラー陽解法では現在位置での評価値がそのまま現在位置に対応する
      density_output     = density
      temperature_output = temperature
      knudsen_output     = knudsen

      # Aerodynamics depending on the attitude
      if flag_attitude :
        acceleration_aerodynamic, moment_total, omega_relative \
          = get_aerodynamic_state(config, attitude_property, aerodynamic_dict, kind_aerodynamic_model, \
                                  coord_tmp, veloc_tmp, quat_tmp, omega_tmp, \
                                  mass_satellite, area_satellite, length_satellite, \
                                  density, density_factor, knudsen, cdmean, stability_derivative, \
                                  rotation_rate_planet, moment, velocity_air)
      else :
        acceleration_aerodynamic = None

      # Calculate acceleration
      # --揚力のバンク角は時刻に依り得る（表で与えるとき）。段の時刻で引く
      angle_bank = force_term.get_bank_angle(config, time_elapsed)
      acceleration = force_term.acceleration_routine(config, coord_tmp, veloc_tmp,   \
                                       mass_satellite, area_satellite, \
                                       cdmean, density_factor, density, acceleration, acceleration_aerodynamic, velocity_air,
                                       angle_bank)
      #"density factor" added by Tomoki Sakai 2023/2/3
      acceleration_total = acceleration[0,:]
 
      # Update solution and Calculate residual
      coord_tmp, veloc_tmp, \
      r_res, v_res =  solve_eulerexplicit(delta_time, acceleration_total, \
                                          coord_tmp, veloc_tmp, \
                                          r_res, v_res)

      if flag_attitude :
        quaternion_rate = attitude.quaternion_derivative(quat_tmp, omega_relative)
        omega_rate      = moment_term.solve_angular_acceleration(attitude_property, omega_tmp, moment_total)
        quat_tmp        = attitude.quaternion_normalize( quat_tmp + quaternion_rate*delta_time )
        omega_tmp       = omega_tmp + omega_rate*delta_time


    elif kind_time_scheme == 'runge_kutta' :
    # 4th stage Runge-Kutta method
    
      # Reinitialization
      r_virtual = coord_tmp
      v_virtual = veloc_tmp
      r_virprev = r_virtual
      v_virprev = v_virtual

      if flag_attitude :
        q_virtual = quat_tmp
        w_virtual = omega_tmp
        q_virprev = q_virtual
        w_virprev = w_virtual

      # R.K. 1st - 4th stages
      for m in range(0, 4):
        # Cartesian --> Geodetic system
        coord_geod = coordinate_system.convert_cartesian_geodetic(config, r_virtual)

        # Recalculate atmosphere status
        altitude_tmp = coord_geod[2] * orbital.m2km
        if kind_atmosphere_model == 'fileread' :
          density, temperature, knudsen = atmosphere.get_atmosphere_property(altitude_tmp, atmosphere_dict)

        # Aerodynamic coefficient
        if kind_aerodynamic_model == 'fileread' :
          cdmean = satellite.get_aerodynamic_coefficient(knudsen, knudsen_aerodynamic, cdmean_aerodynamic)

        # 対気速度。各段の仮想位置で引く（大気量と同じ場所）。
        # 風が無効なら None のままで、以後は従来どおり ECEF 速度が使われる
        velocity_air = None
        if flag_wind :
          velocity_air = wind.get_relative_velocity(config, r_virtual, coord_geod, v_virtual, wind_dict,
                                                    time_elapsed + fact_time_rk[m]*delta_time)

        # 出力用の大気量
        # --1 段目は現在位置 r_n での評価値なので、これを出力に使う。
        # --4 段目の値は r_{n+1} 近傍での評価値であり、r_n の行に並べると 1 ステップずれる。
        if m == 0 :
          density_output     = density
          temperature_output = temperature
          knudsen_output     = knudsen

        # External acceleration (but "Delta V" is given here)
        # --Not yet
        # --Convert cartesian coordinate
        #set_angle_longlat_cartesian(r_virtual, longitude_tmp, latitude_tmp)
        #exchanege_axis_xyz2zxy(externalforce_tmp)
        #convert_coordinate_ry( latitude_tmp ,externalforce_tmp)
        #convert_coordinate_rz(-longitude_tmp,externalforce_tmp)

        # Aerodynamics depending on the attitude
        # --各段の仮想状態で評価する（姿勢と軌道は双方向に結合している）
        if flag_attitude :
          acceleration_aerodynamic, moment_total, omega_relative \
            = get_aerodynamic_state(config, attitude_property, aerodynamic_dict, kind_aerodynamic_model, \
                                    r_virtual, v_virtual, q_virtual, w_virtual, \
                                    mass_satellite, area_satellite, length_satellite, \
                                    density, density_factor, knudsen, cdmean, stability_derivative, \
                                    rotation_rate_planet, moment, velocity_air)
        else :
          acceleration_aerodynamic = None

        # Calculate acceleration
        # --揚力のバンク角も各段の時刻で引く（風と同じ扱い）
        angle_bank = force_term.get_bank_angle(config, time_elapsed + fact_time_rk[m]*delta_time)
        acceleration = force_term.acceleration_routine(config, r_virtual, v_virtual, \
                                         mass_satellite, area_satellite, \
                                         cdmean, density_factor, density, acceleration, acceleration_aerodynamic, velocity_air,
                                         angle_bank)
        #"density factor" added by Tomoki Sakai 2023/2/3
        acceleration_total = acceleration[0,:]
 
        # Update solution and Calculate residual
        coord_tmp, veloc_tmp, \
        r_virtual, v_virtual, \
        r_res, v_res = solve_rungekutta(m, delta_time, acceleration_total, \
                                        coord_tmp, veloc_tmp, \
                                        r_virtual, v_virtual, \
                                        r_virprev, v_virprev, \
                                        r_res, v_res)

        if flag_attitude :
          quaternion_rate = attitude.quaternion_derivative(q_virtual, omega_relative)
          omega_rate      = moment_term.solve_angular_acceleration(attitude_property, w_virtual, moment_total)
          quat_tmp, omega_tmp, \
          q_virtual, w_virtual = solve_rungekutta_attitude(m, delta_time, \
                                                           quaternion_rate, omega_rate, \
                                                           quat_tmp, omega_tmp, \
                                                           q_virprev, w_virprev)

    else :
      print('kind_time_scheme in config is incorrect.')
      print('Program stopped.')
      sys.exit(1)

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

    density_traj.append( density_output )
    temperature_traj.append( temperature_output )
    knudsen_traj.append( knudsen_output )

    if flag_attitude :
      quaternion_list.append( quat_tmp )
      omega_list.append( omega_tmp )

      # 姿勢振動の周期に対して時間刻みが粗いと発散するので、最初に条件を割ったところで警告する
      # （INTERVAL_CHECK_TIMESTEP ごと。警告が 1 秒遅れても困らない）
      if not flag_warned_timestep and \
         time_elapsed - time_checked_timestep >= INTERVAL_CHECK_TIMESTEP :
        time_checked_timestep = time_elapsed
        # 動圧は対気速度で作る（風があると振動周期の見積もりが変わる）
        veloc_aero_tmp = wind.get_relative_velocity(config, coord_tmp, coord_geodetic, veloc_tmp, wind_dict,
                                                    time_elapsed + delta_time)
        flag_warned_timestep = check_timestep_attitude(delta_time, attitude_property, aerodynamic_dict, \
                                                       kind_aerodynamic_model, stability_derivative, \
                                                       density_factor*density, vector_norm(veloc_aero_tmp), \
                                                       knudsen, area_satellite, length_satellite)

    iteration    = iteration + 1
    # 累積加算だと丸め誤差が乗るので、初期値からのオフセットで求める
    time_elapsed = time_elapsed_initial + float(iteration - iteration_initial)*delta_time

    coord_geodetic_display = np.multiply(coord_geodetic, orbital.unit_convert_geoditic)
    veloc_polar_mag = vector_norm(veloc_polar)
    #print("Time elapsed:", time_elapsed, 'Longitude:', coord_geodetic_display[0], 'Latitude:', coord_geodetic_display[1], 'Altitude:', coord_geodetic_display[2], 'Velocity (Mag.)',veloc_polar_mag )
    print('{:.1f}'.format(time_elapsed)+', '+'{:.3f}'.format(coord_geodetic_display[0])+', '+'{:.3f}'.format(coord_geodetic_display[1])+', '+'{:.3f}'.format(coord_geodetic_display[2])+', '+'{:.3f}'.format(veloc_polar_mag) )

    if flag_check_divergence :
      # 高度が負になること以外の打ち切り条件が無いので、発散した状態はそのまま
      # 出力に書かれていた。数値が並んでいるので一見それらしく、後処理も読む
      check_divergence(coord_tmp, veloc_tmp, radius_maximum, velocity_maximum, time_elapsed,
                       quat_tmp if flag_attitude else None,
                       omega_tmp if flag_attitude else None)

    if coord_geodetic_display[2] <= 0.0 :
      break

  if kind_atmosphere_model == 'fileread' :
    altitude_top  = atmosphere_dict[atmosphere.KEY_Height][-1]
    altitude_traj = np.array([coord[2] for coord in coordinate_geod])*orbital.m2km
    num_above     = int(np.count_nonzero(altitude_traj > altitude_top))
    if num_above > 0 :
      if atmosphere_dict[atmosphere.KEY_SCALE_HEIGHT] is None :
        treatment = 'clamped to the value at the top of the table'
      else :
        treatment = 'extrapolated with a scale height of {:.2f} km'.format(
                    atmosphere_dict[atmosphere.KEY_SCALE_HEIGHT])
      print('Note: {:d} of {:d} steps are above the atmosphere table ({:.1f} km).'.format(
            num_above, len(altitude_traj), altitude_top))
      print('--Maximum altitude: {:.1f} km. Atmospheric properties there are {:s}.'.format(
            altitude_traj.max(), treatment))

  # ループ終了後の最終位置における大気量を追加し、位置・速度の配列（iteration+1 個）と長さを揃える。
  # これを行わないと出力ルーチンで最終ステップが欠落する。
  coord_geod   = coordinate_system.convert_cartesian_geodetic(config, coordinate_cart[iteration])
  altitude_tmp = coord_geod[2] * orbital.m2km
  if kind_atmosphere_model == 'fileread' :
    density, temperature, knudsen = atmosphere.get_atmosphere_property(altitude_tmp, atmosphere_dict)
  density_traj.append( density )
  temperature_traj.append( temperature )
  knudsen_traj.append( knudsen )

  return iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict


def solve_eulerexplicit(dt, acceleration, coord_tmp, veloc_tmp, r_res, v_res):

  # Update solution and Calculate residual
  dv        = ( acceleration )*dt
  dr        = ( veloc_tmp )*dt
  veloc_tmp = veloc_tmp + dv
  coord_tmp = coord_tmp + dr

  # Discrepancies from the previous step
  v_res = v_res + vector_norm(dv)**2
  r_res = r_res + vector_norm(dr)**2

  return coord_tmp, veloc_tmp, r_res, v_res


def solve_rungekutta(m, dt, acceleration, coord_tmp, veloc_tmp, r_virtual, v_virtual, r_virprev, v_virprev, r_res, v_res):

  # Update solution and Calculate residual
  kv_rk = ( acceleration )*dt
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
  v_res = v_res + vector_norm(dv)**2
  r_res = r_res + vector_norm(dr)**2
 
  return coord_tmp, veloc_tmp, r_virtual, v_virtual, r_res, v_res

def get_aerodynamic_state(config, attitude_property, aerodynamic_dict, kind_aerodynamic_model,
                          coordinate, velocity, quaternion, omega_inertial,
                          mass_satellite, area_satellite, length_satellite,
                          density, density_factor, knudsen, cdmean, stability_derivative,
                          rotation_rate_planet, moment, velocity_air=None):
  #
  # 姿勢に依存する空力（力・モーメント）を求める。
  #
  # 戻り値
  #   acceleration_aerodynamic: 空力による加速度（ECEF 成分, m/s2）。force_term に渡す
  #   moment_total     : 重心まわりのモーメントの合計（機体軸成分, N m）
  #   omega_relative   : ECEF に対する角速度（機体軸成分, rad/s）
  #
  # 対気速度は velocity_air で与える。None なら ECEF 速度をそのまま使う
  # （大気は地球と共回転しているという、風を入れないときの従来の仮定）。
  # 角速度（omega_relative）は風では変わらない。風は並進の量だから。
  #
  quaternion = attitude.quaternion_normalize(quaternion)
  matrix_be  = attitude.quaternion_to_matrix(quaternion)

  velocity_aero = velocity if velocity_air is None else velocity_air
  velocity_body = np.dot(matrix_be, velocity_aero)
  velocity_mag  = vector_norm(velocity_body)

  alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(velocity_body)

  dynamic_pressure = 0.50*density_factor*density*velocity_mag**2

  # 係数表の面（機体軸 x-z 面）における係数
  if kind_aerodynamic_model == 'fileread' :
    coefficient_force_table, coefficient_moment_table \
      = satellite.get_aerodynamic_coefficient_attitude(knudsen, alpha_total*orbital.rad2deg, aerodynamic_dict)
  else :
    # 定数モデル: 軸力は CD 一定、ピッチングモーメントは静安定微係数の線形モデル
    coefficient_force_table  = np.array([cdmean, 0.0, 0.0])
    coefficient_moment_table = np.array([0.0, stability_derivative*alpha_total, 0.0])

  # 実際の横流れ面へ回す
  matrix_roll = attitude.matrix_aerodynamic_roll(phi_aero)
  coefficient_force_body  = np.dot(matrix_roll, coefficient_force_table)
  coefficient_moment_body = np.dot(matrix_roll, coefficient_moment_table)

  # 力の符号: CFx > 0（= 迎角 0 での CD）が -x 方向（後ろ向き）の力になるようにとる。
  # これで迎角 0 では 3 自由度の抗力と一致する。
  force_aerodynamic_body = -dynamic_pressure*area_satellite*coefficient_force_body
  acceleration_aerodynamic      = np.dot(matrix_be.T, force_aerodynamic_body)/mass_satellite

  omega_relative = attitude.get_omega_relative(omega_inertial, rotation_rate_planet, matrix_be)

  moment = moment_term.moment_routine(config, attitude_property, coordinate, matrix_be, omega_relative,
                                      coefficient_moment_body, force_aerodynamic_body,
                                      dynamic_pressure, area_satellite, length_satellite, velocity_mag,
                                      moment)

  return acceleration_aerodynamic, moment[0,:], omega_relative


def solve_rungekutta_attitude(m, dt, quaternion_rate, omega_rate, quat_tmp, omega_tmp, q_virprev, w_virprev):

  # 並進側の solve_rungekutta と同じ重みで姿勢を進める
  kq_rk = quaternion_rate*dt
  kw_rk = omega_rate*dt

  # Update solution
  quat_tmp  = quat_tmp  + fact_up[m]*one_sixth*kq_rk
  omega_tmp = omega_tmp + fact_up[m]*one_sixth*kw_rk

  # Temporary variables
  # --各段のクォータニオンは単位長に戻す（変換行列を直交に保つため）
  q_virtual = attitude.quaternion_normalize( q_virprev + fact_rk[m]*kq_rk )
  w_virtual = w_virprev + fact_rk[m]*kw_rk

  # 4 段目まで積んだ最終値も単位長に戻す
  if m == 3 :
    quat_tmp = attitude.quaternion_normalize(quat_tmp)

  return quat_tmp, omega_tmp, q_virtual, w_virtual


def set_divergence_limit(config, coordinate_initial):
  #
  # 発散を検知する上限を初期状態から作る。
  #
  # computational_setup:
  #   flag_check_divergence   : 検査そのもののon/off（既定True）
  #   factor_radius_maximum   : 地心距離の上限。初期の地心距離の何倍か（既定10倍）
  #   factor_velocity_maximum : 速度の上限。初期位置での脱出速度の何倍か（既定10倍）
  #
  # 初期状態を基準にするのは、上限を絶対値で決めると軌道の高さによって
  # 意味が変わってしまうため。既定の10倍は、まともな計算では届かない一方で
  # 発散した計算はすぐに超える（実測では地心距離 6.6e9 m、速度 1.4e6 m/s）。
  #
  section = config['computational_setup']

  flag_check      = bool( get_setting(section, 'flag_check_divergence', True) )
  factor_radius   = float( get_setting(section, 'factor_radius_maximum', 10.0) )
  factor_velocity = float( get_setting(section, 'factor_velocity_maximum', 10.0) )

  radius_initial = vector_norm(coordinate_initial)

  # 初期位置での脱出速度 sqrt(2 GM/r)
  parameter_gravitational = config['planet']['gravitational_constant']*config['planet']['mass']
  velocity_escape         = np.sqrt( 2.0*parameter_gravitational/radius_initial )

  return flag_check, factor_radius*radius_initial, factor_velocity*velocity_escape


def check_divergence(coordinate, velocity, radius_maximum, velocity_maximum, time_elapsed,
                     quaternion=None, angular_velocity=None):
  #
  # 状態が非物理になっていないかを1ステップに1度だけ確かめ、駄目なら止める。
  #
  # 止めるのが要点である。発散した計算でも出力ファイル（Tecplot・restart・KML）は
  # 最後まで書かれ、終了コードは0だったので、粗すぎる時間刻みの結果が
  # 「計算できた」ものとして残ってしまう。
  #
  state = np.concatenate( (np.asarray(coordinate).ravel(), np.asarray(velocity).ravel()) )
  if quaternion is not None :
    state = np.concatenate( (state, np.asarray(quaternion).ravel(), np.asarray(angular_velocity).ravel()) )

  message = None
  if not np.all( np.isfinite(state) ) :
    message = 'the state is not a finite number (NaN or Inf)'
  else :
    radius = vector_norm(coordinate)
    if radius > radius_maximum :
      message = 'the geocentric distance {:.6g} m is beyond the limit {:.6g} m'.format(radius, radius_maximum)
    else :
      velocity_magnitude = vector_norm(velocity)
      if velocity_magnitude > velocity_maximum :
        message = 'the velocity {:.6g} m/s is beyond the limit {:.6g} m/s'.format(velocity_magnitude, velocity_maximum)

  if message is None :
    return

  print('The solution has diverged at {:.3f} s: {:s}.'.format(time_elapsed, message))
  print('--A time step which is too coarse for the orbit is the usual cause;')
  print('--halve time_integration.timestep_constant and run it again.')
  print('--computational_setup.factor_radius_maximum / factor_velocity_maximum move the limits,')
  print('--and flag_check_divergence: False switches the check off.')
  print('Program stopped.')
  sys.exit(1)


def check_timestep_attitude(delta_time, attitude_property, aerodynamic_dict, kind_aerodynamic_model,
                            stability_derivative, density, velocity_mag, knudsen,
                            area_satellite, length_satellite):
  #
  # 空力復元モーメントによるピッチ振動の周期を見積もり、時間刻みが粗すぎれば警告する。
  # 戻り値は「警告済みかどうか」。
  #
  #   omega_n^2 = |q S L Cm_alpha| / Iyy,  T = 2 pi / omega_n
  #
  angle_reference = 10.0
  if kind_aerodynamic_model == 'fileread' :
    if len(aerodynamic_dict[satellite.KEY_AOA]) <= 1 :
      return True
    coefficient_force_zero, coefficient_moment_zero = satellite.get_aerodynamic_coefficient_attitude(knudsen, 0.0, aerodynamic_dict)
    coefficient_force_ref,  coefficient_moment_ref  = satellite.get_aerodynamic_coefficient_attitude(knudsen, angle_reference, aerodynamic_dict)
    slope_moment = (coefficient_moment_ref[1] - coefficient_moment_zero[1])/(angle_reference*orbital.deg2rad)
  else :
    slope_moment = stability_derivative

  if slope_moment >= 0.0 :
    # 静的に不安定（あるいはモーメントなし）。振動周期は定義できない
    return False

  dynamic_pressure = 0.50*density*velocity_mag**2
  inertia_pitch    = attitude_property[moment_term.KEY_INERTIA][1,1]
  stiffness        = np.abs(dynamic_pressure*area_satellite*length_satellite*slope_moment)
  if stiffness <= 0.0 :
    return False

  period = 2.0*np.pi*np.sqrt(inertia_pitch/stiffness)
  if delta_time > period/20.0 :
    print('Caution: the timestep is coarse for the attitude oscillation.')
    print('--Estimated pitch period: {:.3f} s, timestep: {:.3f} s.'.format(period, delta_time))
    print('--Use a timestep below about {:.3f} s (period/20) to resolve it.'.format(period/20.0))
    return True

  return False
