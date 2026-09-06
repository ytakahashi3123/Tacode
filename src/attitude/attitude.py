#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/06
#
# 姿勢（6 自由度）まわりの座標変換とキネマティクス。
#
# 座標系とクォータニオンの定義（ここを取り違えると符号バグになる）:
#
#  ECEF (e)  : 地球中心・地球固定の直交系。並進運動を解いている系。
#  NED  (n)  : 現在位置の地心ローカル水平系 [北, 東, 下]。
#              初期姿勢の入力とオイラー角の出力の基準。
#              地心（幾何緯度ではない）なのは、既存の初期速度 [東, 北, 上] と揃えるため。
#  Body (b)  : 機体固定系 [前方, 右, 下]（航空機の標準的な機体軸）。
#
#  クォータニオン q = [q0, q1, q2, q3]（スカラー先頭）は ECEF -> Body の
#  座標変換（受動回転）を表す。すなわち v_b = C(q) v_e。
#  オイラー角は NED -> Body の 3-2-1（ヨー -> ピッチ -> ロール）。
#
#  角速度は 2 種類あるので混同しないこと:
#   omega_inertial: 慣性系に対する角速度（機体軸成分）。Euler の運動方程式はこちら。
#   omega_relative: ECEF に対する角速度（機体軸成分）。クォータニオンの
#                   キネマティクスと空力減衰はこちら。
#   omega_relative = omega_inertial - C(q) * omega_earth(ECEF)

import numpy as np
import sys as sys

# Dict key
KEY_QUATERNION       = 'quaternion'
KEY_ANGULAR_VELOCITY = 'angular_velocity'


def get_setting(section, key, default):
  # config の姿勢まわりの項目はすべて省略可能にしてあるので、既定値を返せるようにする
  if section is None :
    return default
  try:
    value = section[key]
  except (KeyError, TypeError):
    return default
  if value is None :
    return default
  return value


def quaternion_normalize(quaternion):
  # 数値積分で単位長からずれるので毎ステップ正規化する
  norm = np.linalg.norm(quaternion)
  if norm <= 0.0 :
    print('Quaternion norm became zero.')
    print('Program stopped.')
    sys.exit(1)
  return np.array(quaternion)/norm


def quaternion_to_matrix(quaternion):
  # v_b = C v_e となる変換行列（受動回転）
  q0 = quaternion[0]
  q1 = quaternion[1]
  q2 = quaternion[2]
  q3 = quaternion[3]

  matrix = np.array([
    [q0**2 + q1**2 - q2**2 - q3**2, 2.0*(q1*q2 + q0*q3)          , 2.0*(q1*q3 - q0*q2)          ],
    [2.0*(q1*q2 - q0*q3)          , q0**2 - q1**2 + q2**2 - q3**2, 2.0*(q2*q3 + q0*q1)          ],
    [2.0*(q1*q3 + q0*q2)          , 2.0*(q2*q3 - q0*q1)          , q0**2 - q1**2 - q2**2 + q3**2]
  ])

  return matrix


def matrix_to_quaternion(matrix):
  # Shepperd の方法。対角成分の最大値で分岐して 0 除算を避ける
  trace = matrix[0,0] + matrix[1,1] + matrix[2,2]

  candidate = [trace, matrix[0,0], matrix[1,1], matrix[2,2]]
  index_max = int(np.argmax(candidate))

  if index_max == 0 :
    q0 = 0.5*np.sqrt(1.0 + trace)
    fact = 0.25/q0
    q1 = fact*(matrix[1,2] - matrix[2,1])
    q2 = fact*(matrix[2,0] - matrix[0,2])
    q3 = fact*(matrix[0,1] - matrix[1,0])
  elif index_max == 1 :
    q1 = 0.5*np.sqrt(1.0 + matrix[0,0] - matrix[1,1] - matrix[2,2])
    fact = 0.25/q1
    q0 = fact*(matrix[1,2] - matrix[2,1])
    q2 = fact*(matrix[0,1] + matrix[1,0])
    q3 = fact*(matrix[2,0] + matrix[0,2])
  elif index_max == 2 :
    q2 = 0.5*np.sqrt(1.0 - matrix[0,0] + matrix[1,1] - matrix[2,2])
    fact = 0.25/q2
    q0 = fact*(matrix[2,0] - matrix[0,2])
    q1 = fact*(matrix[0,1] + matrix[1,0])
    q3 = fact*(matrix[1,2] + matrix[2,1])
  else :
    q3 = 0.5*np.sqrt(1.0 - matrix[0,0] - matrix[1,1] + matrix[2,2])
    fact = 0.25/q3
    q0 = fact*(matrix[0,1] - matrix[1,0])
    q1 = fact*(matrix[2,0] + matrix[0,2])
    q2 = fact*(matrix[1,2] + matrix[2,1])

  quaternion = np.array([q0, q1, q2, q3])
  # スカラー部が正になる側を選び、表現の二価性（q と -q）を固定する
  if quaternion[0] < 0.0 :
    quaternion = -quaternion

  return quaternion_normalize(quaternion)


def quaternion_derivative(quaternion, omega_relative):
  # dq/dt = 0.5 * Omega(omega) * q
  # omega_relative は ECEF に対する角速度（機体軸成分）
  p = omega_relative[0]
  q = omega_relative[1]
  r = omega_relative[2]

  omega_matrix = np.array([
    [0.0,  -p,  -q,  -r],
    [  p, 0.0,   r,  -q],
    [  q,  -r, 0.0,   p],
    [  r,   q,  -p, 0.0]
  ])

  return 0.5*np.dot(omega_matrix, quaternion)


def euler_to_matrix(yaw, pitch, roll):
  # 3-2-1 オイラー角（rad）から v_b = C v_n となる変換行列を作る
  sin_psi = np.sin(yaw)
  cos_psi = np.cos(yaw)
  sin_the = np.sin(pitch)
  cos_the = np.cos(pitch)
  sin_phi = np.sin(roll)
  cos_phi = np.cos(roll)

  matrix = np.array([
    [ cos_the*cos_psi,
      cos_the*sin_psi,
     -sin_the ],
    [ sin_phi*sin_the*cos_psi - cos_phi*sin_psi,
      sin_phi*sin_the*sin_psi + cos_phi*cos_psi,
      sin_phi*cos_the ],
    [ cos_phi*sin_the*cos_psi + sin_phi*sin_psi,
      cos_phi*sin_the*sin_psi - sin_phi*cos_psi,
      cos_phi*cos_the ]
  ])

  return matrix


def matrix_to_euler(matrix):
  # v_b = C v_n から 3-2-1 オイラー角（rad）を取り出す
  # ピッチ +-90 度でジンバルロックするので、そこではロールを 0 に倒す
  sin_the = -matrix[0,2]
  sin_the = min(1.0, max(-1.0, sin_the))
  pitch   = np.arcsin(sin_the)

  if np.abs(sin_the) > 1.0 - 1.e-10 :
    yaw  = np.arctan2(-matrix[1,0], matrix[1,1])
    roll = 0.0
  else :
    yaw  = np.arctan2(matrix[0,1], matrix[0,0])
    roll = np.arctan2(matrix[1,2], matrix[2,2])

  return yaw, pitch, roll


def matrix_ecef_to_ned(longitude, latitude):
  # 地心経度・緯度（rad）における ECEF -> NED [北, 東, 下] の変換行列
  sin_lon = np.sin(longitude)
  cos_lon = np.cos(longitude)
  sin_lat = np.sin(latitude)
  cos_lat = np.cos(latitude)

  vector_north = [-sin_lat*cos_lon, -sin_lat*sin_lon,  cos_lat]
  vector_east  = [        -sin_lon,          cos_lon,      0.0]
  vector_down  = [-cos_lat*cos_lon, -cos_lat*sin_lon, -sin_lat]

  return np.array([vector_north, vector_east, vector_down])


def matrix_ecef_to_body(quaternion):
  return quaternion_to_matrix(quaternion)


def get_euler_angle(quaternion, longitude, latitude):
  # ローカル水平系（NED）を基準にしたオイラー角（rad）。出力用
  matrix_be = quaternion_to_matrix(quaternion)
  matrix_ne = matrix_ecef_to_ned(longitude, latitude)
  matrix_bn = np.dot(matrix_be, matrix_ne.T)

  return matrix_to_euler(matrix_bn)


def get_quaternion_from_euler(yaw, pitch, roll, longitude, latitude):
  # ローカル水平系基準のオイラー角（rad）から ECEF 基準のクォータニオンを作る
  matrix_bn = euler_to_matrix(yaw, pitch, roll)
  matrix_ne = matrix_ecef_to_ned(longitude, latitude)
  matrix_be = np.dot(matrix_bn, matrix_ne)

  return matrix_to_quaternion(matrix_be)


def get_earth_rate_body(rotation_rate_planet, matrix_be):
  # 地球自転角速度（ECEF では [0, 0, omega_e]）を機体軸成分にする
  return np.dot(matrix_be, np.array([0.0, 0.0, rotation_rate_planet]))


def get_omega_relative(omega_inertial, rotation_rate_planet, matrix_be):
  # ECEF に対する角速度（機体軸成分）
  return np.array(omega_inertial) - get_earth_rate_body(rotation_rate_planet, matrix_be)


def get_aerodynamic_angle(velocity_body):
  #
  # 機体軸成分の対気速度から空力角を得る。
  #
  #  alpha       : 迎角 atan2(w, u)
  #  beta        : 横滑り角 asin(v/V)
  #  alpha_total : 全迎角 arccos(u/V)。0 以上で、軸対称体の係数表の引数になる
  #  phi_aero    : 空力ロール角 atan2(v, w)。係数表（x-z 面）を実際の
  #                横流れ面に回すのに使う
  #
  velocity_mag = np.linalg.norm(velocity_body)
  if velocity_mag <= 0.0 :
    return 0.0, 0.0, 0.0, 0.0

  vel_u = velocity_body[0]
  vel_v = velocity_body[1]
  vel_w = velocity_body[2]

  alpha = np.arctan2(vel_w, vel_u)
  beta  = np.arcsin( min(1.0, max(-1.0, vel_v/velocity_mag)) )

  cos_alpha_total = min(1.0, max(-1.0, vel_u/velocity_mag))
  alpha_total     = np.arccos(cos_alpha_total)

  velocity_cross = np.sqrt(vel_v**2 + vel_w**2)
  if velocity_cross <= 0.0 :
    phi_aero = 0.0
  else :
    phi_aero = np.arctan2(vel_v, vel_w)

  return alpha, beta, alpha_total, phi_aero


def matrix_aerodynamic_roll(phi_aero):
  #
  # 係数表は横流れが +z 方向（x-z 面）にある状態で作られている。
  # 実際の横流れ方向 [0, sin(phi), cos(phi)] に係数ベクトルを回す行列。
  # 力・モーメントの係数ベクトルの両方に同じものを掛ける（固有回転なので
  # 擬ベクトルであるモーメントも同じ変換で良い）。
  #
  sin_phi = np.sin(phi_aero)
  cos_phi = np.cos(phi_aero)

  return np.array([
    [1.0,      0.0,     0.0],
    [0.0,  cos_phi, sin_phi],
    [0.0, -sin_phi, cos_phi]
  ])
