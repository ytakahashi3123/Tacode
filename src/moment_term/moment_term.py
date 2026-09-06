#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/06
#
# 機体重心まわりのモーメント（トルク）項。force_term.py と同じく
# 「項ごとに配列へ積んでから合計する」構造にしてある。
#
#   moment[0,:] : 合計
#   moment[1,:] : 空力（静的）
#   moment[2,:] : 空力（動的減衰）
#   moment[3,:] : 重力傾斜
#
# いずれも機体軸成分 [N m]。force_term.py の "force" が実際には加速度なのに対し、
# こちらは角加速度ではなくモーメントそのものであることに注意。

import numpy as np
import sys as sys
import attitude.attitude as attitude

# Dict key
KEY_INERTIA      = 'inertia_tensor'
KEY_INERTIA_INV  = 'inertia_tensor_inverse'
KEY_CG           = 'center_of_gravity'
KEY_DAMPING      = 'damping_coefficient'
KEY_FLAG_AERO    = 'flag_moment_aerodynamic'
KEY_FLAG_DAMPING = 'flag_moment_damping'
KEY_FLAG_GRAVITY = 'flag_moment_gravity_gradient'


def moment_initialsettings(config):

  moment = np.zeros(4*3).reshape(4,3)

  section = attitude.get_setting(config, 'attitude', None)

  # 慣性テンソル（機体軸・重心まわり, kg m2）
  inertia_setting = attitude.get_setting(section, 'inertia_tensor', None)
  if inertia_setting is None :
    print('"inertia_tensor" is not given in the attitude section of config.')
    print('--It is required for the 6-DOF (attitude) computation.')
    print('Program stopped.')
    sys.exit(1)

  moment_inertia_xx = float(attitude.get_setting(inertia_setting, 'Ixx', 0.0))
  moment_inertia_yy = float(attitude.get_setting(inertia_setting, 'Iyy', 0.0))
  moment_inertia_zz = float(attitude.get_setting(inertia_setting, 'Izz', 0.0))
  # 慣性乗積は「慣性テンソルの非対角成分」として与える（-Ixy ではない）
  moment_inertia_xy = float(attitude.get_setting(inertia_setting, 'Ixy', 0.0))
  moment_inertia_yz = float(attitude.get_setting(inertia_setting, 'Iyz', 0.0))
  moment_inertia_zx = float(attitude.get_setting(inertia_setting, 'Izx', 0.0))

  inertia = np.array([
    [moment_inertia_xx, moment_inertia_xy, moment_inertia_zx],
    [moment_inertia_xy, moment_inertia_yy, moment_inertia_yz],
    [moment_inertia_zx, moment_inertia_yz, moment_inertia_zz]
  ])

  if np.linalg.det(inertia) <= 0.0 :
    print('The inertia tensor in config is not positive definite.')
    print('--Check Ixx, Iyy, Izz and the products of inertia.')
    print('Program stopped.')
    sys.exit(1)

  # 空力係数の基準点から重心へのベクトル（機体軸, m）
  center_of_gravity = np.array(attitude.get_setting(section, 'center_of_gravity', [0.0, 0.0, 0.0]), dtype=float)

  # 動的減衰係数（機体軸まわり）
  damping_setting = attitude.get_setting(section, 'damping_coefficient', None)
  damping = np.array([float(attitude.get_setting(damping_setting, 'Clp', 0.0)),
                      float(attitude.get_setting(damping_setting, 'Cmq', 0.0)),
                      float(attitude.get_setting(damping_setting, 'Cnr', 0.0))])

  property_dict = {
    KEY_INERTIA     : inertia,
    KEY_INERTIA_INV : np.linalg.inv(inertia),
    KEY_CG          : center_of_gravity,
    KEY_DAMPING     : damping,
    KEY_FLAG_AERO   : bool(attitude.get_setting(section, 'flag_moment_aerodynamic', True)),
    KEY_FLAG_DAMPING: bool(attitude.get_setting(section, 'flag_moment_damping', True)),
    KEY_FLAG_GRAVITY: bool(attitude.get_setting(section, 'flag_moment_gravity_gradient', True))
  }

  return moment, property_dict


def moment_routine(config, property_dict, coordinate, matrix_be, omega_relative,
                   coefficient_moment_body, force_aerodynamic_body,
                   dynamic_pressure, area_satellite, length_satellite, velocity_mag,
                   moment):

  moment[:,:] = 0.0

  # 空力モーメント（静的）
  # --係数は空力データベースの基準点まわりなので、重心へ移す項を加える
  if property_dict[KEY_FLAG_AERO] :
    fact_moment = dynamic_pressure*area_satellite*length_satellite
    moment[1,:] = fact_moment*np.array(coefficient_moment_body) \
                + np.cross(-property_dict[KEY_CG], force_aerodynamic_body)

  # 空力モーメント（動的減衰）
  # --M = q S L * (L/2V) * [Clp*p, Cmq*q, Cnr*r]。角速度は大気（ECEF）に対するもの
  if property_dict[KEY_FLAG_DAMPING] and velocity_mag > 0.0 :
    fact_damping = dynamic_pressure*area_satellite*length_satellite*length_satellite/(2.0*velocity_mag)
    moment[2,:] = fact_damping*np.multiply(property_dict[KEY_DAMPING], omega_relative)

  # 重力傾斜トルク
  # --M = 3 GM/R^3 * (u x I u)、u は重心から地心へ向かう単位ベクトル（機体軸成分）
  # --u は 2 回現れるので符号（地心向き／天頂向き）はどちらでも同じ
  if property_dict[KEY_FLAG_GRAVITY] :
    gravity_const_planet = config['planet']['gravitational_constant']
    mass_planet          = config['planet']['mass']

    radius_pmass = np.linalg.norm(coordinate)
    unit_body    = np.dot(matrix_be, np.array(coordinate)/radius_pmass)
    fact_gravity = 3.0*gravity_const_planet*mass_planet/radius_pmass**3
    moment[3,:]  = fact_gravity*np.cross(unit_body, np.dot(property_dict[KEY_INERTIA], unit_body))

  # Total
  moment[0,:] = moment[1,:] + moment[2,:] + moment[3,:]

  return moment


def solve_angular_acceleration(property_dict, omega_inertial, moment_total):
  #
  # Euler の運動方程式
  #   I domega/dt + omega x (I omega) = M
  # omega は慣性系に対する角速度の機体軸成分。ECEF に対する角速度ではないので注意。
  #
  angular_momentum = np.dot(property_dict[KEY_INERTIA], omega_inertial)
  gyroscopic       = np.cross(omega_inertial, angular_momentum)

  return np.dot(property_dict[KEY_INERTIA_INV], np.array(moment_total) - gyroscopic)
