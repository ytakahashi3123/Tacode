#!/usr/bin/env python3
#
# Mars Pathfinder の突入状態（Spencer et al. 1998, AAS 98-146 表 1、軌道決定の推定値）を
# config.yml が要る形に直す。
#
#   Epoch 1997-07-04 16:51:50.482 UTC   Radial distance 3522.2 km
#   Areocentric latitude 22.6303 deg.   Longitude 338.1691 deg.
#   Inertial flight path angle -14.0610 deg.   Inertial velocity 7.2642 km/s
#   Inertial flight path azimuth 253.1479 deg.
#
# 表 1 は**慣性系**の速さ・経路角・方位角を与える。Tacode の initial_settings.velocity は
# 火星に対する（火星固定系の）速度を地心ローカル系 [東, 北, 上] で与えるので、自転ぶん
# omega x r を引く。惑星は球（config.yml の ellipticity 0）なので、測地と地心の区別は無い。
#
# **なぜ PDS の突入状態を使わないか**: edl_ddr.tab に添えられた突入状態（r = 3597.2 km、
# 緯度 23.0、経度 343.67）は、同じファイルの位置の履歴の上に乗っていない（速さ・経路角・
# 方位角はそのままで、位置だけを合わせると北へ 54 km・西へ 49 km 動く。convert_pds.py が
# 表示する）。表 1 の状態は位置の履歴と平行に、全区間で一定の距離（6.9 km）を保って進む。
# 速さの履歴は両者で 0.02 % 以内に一致する。
#
# 使い方:
#   python3 entry_state.py

import numpy as np

RADIUS_PLANET = 3396.0e3          # m。config.yml の planet.radius
ROTATION_RATE = 7.0882181e-5      # rad/s。config.yml の planet.rotation_rate

RADIUS    = 3522.2e3              # m
LATITUDE  = 22.6303               # deg. N
LONGITUDE = 338.1691              # deg. E
SPEED_INERTIAL   = 7264.2         # m/s
PATH_INERTIAL    = -14.0610       # deg.
AZIMUTH_INERTIAL = 253.1479       # deg.


def main():
  latitude, longitude = np.radians(LATITUDE), np.radians(LONGITUDE)
  east  = np.array([-np.sin(longitude), np.cos(longitude), 0.0])
  north = np.array([-np.sin(latitude)*np.cos(longitude), -np.sin(latitude)*np.sin(longitude), np.cos(latitude)])
  up    = np.array([ np.cos(latitude)*np.cos(longitude),  np.cos(latitude)*np.sin(longitude), np.sin(latitude)])

  path, azimuth = np.radians(PATH_INERTIAL), np.radians(AZIMUTH_INERTIAL)
  velocity_inertial = SPEED_INERTIAL*(np.cos(path)*(np.sin(azimuth)*east + np.cos(azimuth)*north) + np.sin(path)*up)
  velocity = velocity_inertial - np.cross([0.0, 0.0, ROTATION_RATE], RADIUS*up)
  local = [float(velocity @ east), float(velocity @ north), float(velocity @ up)]
  speed = float(np.linalg.norm(velocity))

  longitude_east = (LONGITUDE + 180.0) % 360.0 - 180.0
  print('Relative (Mars-fixed) speed {:.2f} m/s, flight-path angle {:.4f} deg., azimuth {:.4f} deg.'.format(
        speed, np.degrees(np.arcsin(local[2]/speed)), np.degrees(np.arctan2(local[0], local[1])) % 360.0))
  print('')
  print('initial_settings:')
  print('  coordinate:')
  print('    - {:.4f}'.format(longitude_east))
  print('    - {:.4f}'.format(LATITUDE))
  print('    - {:.1f}'.format((RADIUS - RADIUS_PLANET)*1.e-3))
  print('  velocity:')
  for value in local:
    print('    - {:.2f}'.format(value))


if __name__ == '__main__':
  main()
