#!/usr/bin/env python3
#
# 6 自由度の初期姿勢（config の initial_settings.attitude に入れる 3-2-1 オイラー角）と
# 初期角速度を作る。
#
# **Fire の機体にはトリム迎角が無い。**軸対称で重心も軸上にあるので、釣り合う姿勢は
# 迎角 0 のただ一つで、飛行中に見えた迎角はすべて擾乱の名残の振動である。したがって
# ここで与えるのは「トリム姿勢」ではなく**振動の初期振幅**で、報告書（NASA TN D-4183）が
#
#   分離時に機体 x 軸と速度ベクトルの間が 1.52 度、そのあとの coning 半角が 0.8 度、
#   大気に入るまで全迎角は 0.72〜2.32 度の間を往復した
#
# と書いているとおりの 1.52 度を既定にしてある。軸対称なので**どの面に倒すかは任意**で、
# ここでは鉛直面（速度に直交する上向き側）に倒す。
#
# もう一つの初期値が**スピン**で、こちらが姿勢運動の形を決める。分離直後のロールレートは
# テレメトリの信号強度から 171 rpm = 17.91 rad/s と測られている（TN D-3569）。
# 軸対称のスピン機は、静安定による振動と自転が結合して 2 つの周波数に分かれるので、
# スピンを入れないと周期の比較にならない。
#
# config の initial_settings.angular_velocity は **ECEF 基準**の [p, q, r] を deg/s で
# 与える（Tacode の規約）。地球の自転 7.3e-5 rad/s はスピンの 4 万分の 1 なので、
# 慣性基準で測られた 171 rpm をそのまま入れてよい。
#
# 使い方:
#   python3 initial_attitude.py
#   python3 initial_attitude.py --alpha 2.32 --spin 171

import argparse
import os
import sys

import numpy as np

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src'))

from general.general import general       # noqa: E402
import attitude.attitude as attitude      # noqa: E402
from orbital.orbital import orbital       # noqa: E402

# NASA TN D-4183: 分離時の機体軸と速度ベクトルのなす角
ANGLE_ATTACK_DEFAULT = 1.52

# NASA TN D-3569: 分離直後のロールレート、rpm
SPIN_DEFAULT = 171.0


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--file-config', type=str, default='config.yml',
                      help='Configuration to take the initial position and velocity from')
  parser.add_argument('--alpha', type=float, default=ANGLE_ATTACK_DEFAULT,
                      help='Angle between the body x axis and the velocity [deg.]')
  parser.add_argument('--spin', type=float, default=SPIN_DEFAULT,
                      help='Roll rate about the body x axis [rpm]')
  args = parser.parse_args()

  os.chdir(DIRECTORY_SCRIPT)
  config = general().read_config_yaml(args.file_config)

  orb = orbital()
  dummy, dummy2, coordinate_dict, velocity_dict, dummy3 = orb.initial_settings(config)
  coordinate = np.array(coordinate_dict['cartesian'][0])
  velocity   = np.array(velocity_dict['cartesian'][0])

  angle_attack = args.alpha*np.pi/180.0

  direction_velocity = velocity/np.linalg.norm(velocity)
  direction_up       = coordinate/np.linalg.norm(coordinate)

  # 速度に直交する鉛直上向き側。軸対称なので倒す面はどこでもよい
  component_up = direction_up - np.dot(direction_up, direction_velocity)*direction_velocity
  component_up = component_up/np.linalg.norm(component_up)

  axis_x = np.cos(angle_attack)*direction_velocity - np.sin(angle_attack)*component_up
  axis_z = -(np.sin(angle_attack)*direction_velocity + np.cos(angle_attack)*component_up)
  axis_y = np.cross(axis_z, axis_x)

  # ECEF -> 機体 の変換行列は、行が機体軸の ECEF 成分
  matrix_be = np.array([axis_x, axis_y, axis_z])

  coordinate_polar = coordinate_dict[orb.KEY_COORD_POLAR][0]
  matrix_ne = attitude.matrix_ecef_to_ned(coordinate_polar[2], coordinate_polar[1])
  yaw, pitch, roll = attitude.matrix_to_euler(np.dot(matrix_be, matrix_ne.T))

  velocity_body = np.dot(matrix_be, velocity)
  alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(velocity_body)

  spin = args.spin*2.0*np.pi/60.0

  print('Initial attitude for a total angle of attack of {:.3f} deg.'.format(args.alpha))
  print('--Check: angle of attack {:+.3f} deg., sideslip {:+.3f} deg., total {:.3f} deg.'.format(
        alpha*orbital.rad2deg, beta*orbital.rad2deg, alpha_total*orbital.rad2deg))
  print('--Spin : {:.2f} rpm = {:.3f} rad/s = {:.1f} deg/s'.format(
        args.spin, spin, spin*orbital.rad2deg))
  print('')
  print('  attitude:')
  for angle in (yaw, pitch, roll):
    print('    - {:.4f}'.format(angle*orbital.rad2deg))
  print('  angular_velocity:')
  print('    - {:.2f}'.format(spin*orbital.rad2deg))
  print('    - 0.0')
  print('    - 0.0')

  return


if __name__ == '__main__':
  main()
