#!/usr/bin/env python3
#
# 6 自由度の初期姿勢（config の initial_settings.attitude に入れる 3-2-1 オイラー角）を、
# **トリム迎角とバンク角**から作る。
#
# 機体軸まわりにロールさせるだけでは「速度まわりのバンク」にならない。機体 x 軸は
# 速度から迎角のぶん傾いているので、機体軸まわりに回すと横流れ面（速度と機体 x 軸が
# 張る面）が重心オフセットの面からずれ、ずれた分だけロールモーメントが立って
# 機体が回り続ける。ここでは**速度ベクトルまわりに**トリム面ごと回す:
#
#   L_hat = cos(bank) n_up + sin(bank) n_right     揚力を向けたい方向
#   x_b   = cos(alpha) v_hat - sin(alpha) L_hat    機体 x 軸（伝熱面の外向き）
#   z_b   = -( sin(alpha) v_hat + cos(alpha) L_hat )
#   y_b   = z_b x x_b
#
# n_up は速度に直交する鉛直上向き側、n_right = v_hat x n_up は進行方向の右。
# これで機体軸成分の対気速度は (cos alpha, 0, -sin alpha)、すなわち横滑り 0・
# 迎角 -alpha になり、横流れ面が機体 x-z 面と一致する（重心オフセットの面と揃う）。
#
# 使い方:
#   python3 attitude_for_bank.py --alpha 24.4 --bank 45
#   python3 attitude_for_bank.py --alpha 24.4 --bank 0 --file-config config_6dof.yml

import argparse
import os
import sys

import numpy as np

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src'))

from general.general import general       # noqa: E402
import attitude.attitude as attitude      # noqa: E402
import coordinate_system.coordinate_system as coordinate_system   # noqa: E402
from orbital.orbital import orbital       # noqa: E402


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--file-config', type=str, default='config.yml',
                      help='Configuration to take the initial position and velocity from')
  parser.add_argument('--alpha', type=float, required=True,
                      help='Trim angle of attack [deg.], the angle between the body x axis'
                           ' (out of the heat shield) and the velocity')
  parser.add_argument('--bank', type=float, default=0.0,
                      help='Bank angle [deg.] of the lift about the velocity: 0 is lift up,'
                           ' 90 is lift to the right of the flight direction, 180 is lift down')
  args = parser.parse_args()

  os.chdir(DIRECTORY_SCRIPT)
  config = general().read_config_yaml(args.file_config)

  orb = orbital()
  iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = orb.initial_settings(config)
  coordinate = np.array(coordinate_dict['cartesian'][0])
  velocity   = np.array(velocity_dict['cartesian'][0])

  angle_attack = args.alpha*np.pi/180.0
  angle_bank   = args.bank*np.pi/180.0

  direction_velocity = velocity/np.linalg.norm(velocity)
  direction_up       = coordinate/np.linalg.norm(coordinate)

  # 速度に直交する鉛直上向き側と、進行方向の右
  component_up = direction_up - np.dot(direction_up, direction_velocity)*direction_velocity
  component_up = component_up/np.linalg.norm(component_up)
  component_right = np.cross(direction_velocity, component_up)

  direction_lift = np.cos(angle_bank)*component_up + np.sin(angle_bank)*component_right

  axis_x = np.cos(angle_attack)*direction_velocity - np.sin(angle_attack)*direction_lift
  axis_z = -( np.sin(angle_attack)*direction_velocity + np.cos(angle_attack)*direction_lift )
  axis_y = np.cross(axis_z, axis_x)

  # ECEF -> 機体 の変換行列は、行が機体軸の ECEF 成分
  matrix_be = np.array([axis_x, axis_y, axis_z])

  # ローカル水平系（地心 NED）に対する 3-2-1 オイラー角。
  # 基準の経度・緯度は solver と同じく初期位置の**地心**極座標から採る
  coordinate_polar = coordinate_dict[orb.KEY_COORD_POLAR][0]
  matrix_ne = attitude.matrix_ecef_to_ned(coordinate_polar[2], coordinate_polar[1])
  matrix_bn = np.dot(matrix_be, matrix_ne.T)
  yaw, pitch, roll = attitude.matrix_to_euler(matrix_bn)

  velocity_body = np.dot(matrix_be, velocity)
  alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(velocity_body)

  print('Initial attitude for alpha = {:.3f} deg. and bank = {:.3f} deg.'.format(args.alpha, args.bank))
  print('--Check: angle of attack {:+.3f} deg., sideslip {:+.3f} deg., total {:.3f} deg., aerodynamic roll {:.3f} deg.'.format(
        alpha*orbital.rad2deg, beta*orbital.rad2deg, alpha_total*orbital.rad2deg, phi_aero*orbital.rad2deg))
  print('')
  print('  attitude:')
  print('    - {:.4f}'.format(yaw*orbital.rad2deg))
  print('    - {:.4f}'.format(pitch*orbital.rad2deg))
  print('    - {:.4f}'.format(roll*orbital.rad2deg))

  return


if __name__ == '__main__':
  main()
