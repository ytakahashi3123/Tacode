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
from orbital.orbital import orbital

# Constants
one_sixth = 1.0/6.0
fact_rk   = [0.5, 0.5, 1.0, 0.0]
fact_up   = [1.0, 2.0, 2.0, 1.0]



def solve_equation_motion(config, iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict, atmosphere_dict, aerodynamic_dict, attitude_dict=None):
  
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
    altitude_aerodynamic = aerodynamic_dict['Altitude']
    interpolator_aerodynamic = aerodynamic_dict[satellite.KEY_INTERP]
  else :
    print('kind_aerodynamic_model in config is incorrect.')
    print('Program stopped.')
    sys.exit(1)

  # Force setting
  force = force_term.force_initialsettings(config)

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
  else :
    moment               = None
    attitude_property    = None
    quaternion_list      = None
    omega_list           = None
    rotation_rate_planet = None
    stability_derivative = 0.0
    flag_warned_timestep = True

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
        cdmean = satellite.get_aerodynamic_coefficient(knudsen, knudsen_aerodynamic, cdmean_aerodynamic, interpolator_aerodynamic)

      # 出力用の大気量
      # --オイラー陽解法では現在位置での評価値がそのまま現在位置に対応する
      density_output     = density
      temperature_output = temperature
      knudsen_output     = knudsen

      # Aerodynamics depending on the attitude
      if flag_attitude :
        force_aerodynamic, moment_total, omega_relative \
          = get_aerodynamic_state(config, attitude_property, aerodynamic_dict, kind_aerodynamic_model, \
                                  coord_tmp, veloc_tmp, quat_tmp, omega_tmp, \
                                  mass_satellite, area_satellite, length_satellite, \
                                  density, density_factor, knudsen, cdmean, stability_derivative, \
                                  rotation_rate_planet, moment)
      else :
        force_aerodynamic = None

      # Calculate force
      force = force_term.force_routine(config, coord_tmp, veloc_tmp,   \
                                       mass_satellite, area_satellite, \
                                       cdmean, density_factor, density, force, force_aerodynamic)
      #"density factor" added by Tomoki Sakai 2023/2/3
      force_total = force[0,:]
 
      # Update solution and Calculate residual
      coord_tmp, veloc_tmp, \
      r_res, v_res =  solve_eulerexplicit(delta_time, mass_satellite, force_total, \
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
          cdmean = satellite.get_aerodynamic_coefficient(knudsen, knudsen_aerodynamic, cdmean_aerodynamic, interpolator_aerodynamic)

        # 出力用の大気量
        # --1 段目は現在位置 r_n での評価値なので、これを出力に使う。
        # --4 段目の値は r_{n+1} 近傍での評価値であり、r_n の行に並べると 1 ステップずれる。
        if m == 0 :
          density_output     = density
          temperature_output = temperature
          knudsen_output     = knudsen

        # External force (but "Delta V" is given here)
        # --Not yet
        # --Convert cartesian coordinate
        #set_angle_longlat_cartesian(r_virtual, longitude_tmp, latitude_tmp)
        #exchanege_axis_xyz2zxy(externalforce_tmp)
        #convert_coordinate_ry( latitude_tmp ,externalforce_tmp)
        #convert_coordinate_rz(-longitude_tmp,externalforce_tmp)

        # Aerodynamics depending on the attitude
        # --各段の仮想状態で評価する（姿勢と軌道は双方向に結合している）
        if flag_attitude :
          force_aerodynamic, moment_total, omega_relative \
            = get_aerodynamic_state(config, attitude_property, aerodynamic_dict, kind_aerodynamic_model, \
                                    r_virtual, v_virtual, q_virtual, w_virtual, \
                                    mass_satellite, area_satellite, length_satellite, \
                                    density, density_factor, knudsen, cdmean, stability_derivative, \
                                    rotation_rate_planet, moment)
        else :
          force_aerodynamic = None

        # Calculate force
        force = force_term.force_routine(config, r_virtual, v_virtual, \
                                         mass_satellite, area_satellite, \
                                         cdmean, density_factor, density, force, force_aerodynamic)
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
      if not flag_warned_timestep :
        flag_warned_timestep = check_timestep_attitude(delta_time, attitude_property, aerodynamic_dict, \
                                                       kind_aerodynamic_model, stability_derivative, \
                                                       density_factor*density, np.linalg.norm(veloc_tmp), \
                                                       knudsen, area_satellite, length_satellite)

    iteration    = iteration + 1
    # 累積加算だと丸め誤差が乗るので、初期値からのオフセットで求める
    time_elapsed = time_elapsed_initial + float(iteration - iteration_initial)*delta_time

    coord_geodetic_display = np.multiply(coord_geodetic, orbital.unit_convert_geoditic)
    veloc_polar_mag = np.linalg.norm(veloc_polar)
    #print("Time elapsed:", time_elapsed, 'Longitude:', coord_geodetic_display[0], 'Latitude:', coord_geodetic_display[1], 'Altitude:', coord_geodetic_display[2], 'Velocity (Mag.)',veloc_polar_mag )
    print('{:.1f}'.format(time_elapsed)+', '+'{:.3f}'.format(coord_geodetic_display[0])+', '+'{:.3f}'.format(coord_geodetic_display[1])+', '+'{:.3f}'.format(coord_geodetic_display[2])+', '+'{:.3f}'.format(veloc_polar_mag) )

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

def get_aerodynamic_state(config, attitude_property, aerodynamic_dict, kind_aerodynamic_model,
                          coordinate, velocity, quaternion, omega_inertial,
                          mass_satellite, area_satellite, length_satellite,
                          density, density_factor, knudsen, cdmean, stability_derivative,
                          rotation_rate_planet, moment):
  #
  # 姿勢に依存する空力（力・モーメント）を求める。
  #
  # 戻り値
  #   force_aerodynamic: 空力による加速度（ECEF 成分, m/s2）。force_term に渡す
  #   moment_total     : 重心まわりのモーメントの合計（機体軸成分, N m）
  #   omega_relative   : ECEF に対する角速度（機体軸成分, rad/s）
  #
  # ECEF 速度をそのまま対気速度として使う（大気は地球と共回転しているという
  # 3 自由度計算と同じ仮定）。
  #
  quaternion = attitude.quaternion_normalize(quaternion)
  matrix_be  = attitude.quaternion_to_matrix(quaternion)

  velocity_body = np.dot(matrix_be, np.array(velocity))
  velocity_mag  = np.linalg.norm(velocity_body)

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
  force_aerodynamic_body = -dynamic_pressure*area_satellite*np.array(coefficient_force_body)
  force_aerodynamic      = np.dot(matrix_be.T, force_aerodynamic_body)/mass_satellite

  omega_relative = attitude.get_omega_relative(omega_inertial, rotation_rate_planet, matrix_be)

  moment = moment_term.moment_routine(config, attitude_property, coordinate, matrix_be, omega_relative,
                                      coefficient_moment_body, force_aerodynamic_body,
                                      dynamic_pressure, area_satellite, length_satellite, velocity_mag,
                                      moment)

  return force_aerodynamic, moment[0,:], omega_relative


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
