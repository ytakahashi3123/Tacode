#!/usr/bin/env python3
#
# Project Fire flight II の突入時の状態を config.yml が要る形に直す。
#
# 表 V（reference/fire2_trajectory.dat）が与えるのは**地球相対**の速度・経路角・
# 方位角で、位置は測地座標。Tacode の initial_settings.velocity は **ECEF 速度**を
# **地心**ローカル系 [東, 北, 上] で与えるので、変換は 2 段しかない:
#
#   1. 経路角・方位角から**測地**ローカル系の速度 [東, 北, 上] を作る
#   2. 地心ローカル系 [東, 北, 上] に直す
#
# **Apollo 10 と違って自転を引く処理は要らない。**あちらの表 II は慣性速度だが、
# 表 V は「Earth relative velocity」と「Earth relative flight-path angle」で、
# 既に ECEF の量になっている。ここを取り違えると 400 m/s ほどずれる。
#
# 測地と地心では鉛直が最大 0.19 度違い、ここでは鉛直速度が 30 m/s ほど変わる。
# 報告書は経路角を「relative to the local geodetic horizon」と明記しているので
# 測地を採る。--horizon geocentric で他方も試せる。
#
# 突入界面（t = 1617.75 s、高度 121 920 m）ではなく **t = 1618.25 s から始める**。
# 走査の傷みで 1617.75 s の緯度・経度が読めないためで、0.5 秒ぶんの違いしかない
# （高度 1.4 km、速度 1 m/s）。比較は飛行経過時刻で突き合わせるので影響しない。
#
# 使い方:
#   python3 entry_state.py
#   python3 entry_state.py --time 1620.25
#   python3 entry_state.py --horizon geocentric

import argparse
import os
import sys

import numpy as np

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src'))

from general.general import general                              # noqa: E402
import coordinate_system.coordinate_system as coordinate_system  # noqa: E402

FILE_TRAJECTORY = os.path.join(DIRECTORY_SCRIPT, 'reference/fire2_trajectory.dat')

# 既定の開始時刻。表 V で全項目が読めている最初のレコード
TIME_START = 1618.25

# 打ち上げ時刻（NASA TN D-3569 本文）。飛行経過時刻を UTC に直すのに使う
EPOCH_LIFTOFF = '1965-05-22T21:54:59.703Z'

COLUMN = ['time', 'latitude', 'longitude', 'altitude', 'velocity', 'flightpath',
          'heading', 'dynamic_pressure', 'pressure', 'density', 'temperature',
          'mach', 'acceleration', 'reynolds']


def read_record(file_trajectory, time):
  # 表 V の 1 レコードを読む
  with open(file_trajectory) as stream:
    for line in stream:
      if line.startswith('#') :
        continue
      value = [float(item) for item in line.split()]
      if abs(value[0] - time) < 1.0e-6 :
        return dict(zip(COLUMN, value))
  print('No record at t = {:.2f} s in {:s}.'.format(time, file_trajectory))
  print('Program stopped.')
  sys.exit(1)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--file-config', type=str, default='config.yml',
                      help='Configuration to take the planet constants from')
  parser.add_argument('--time', type=float, default=TIME_START,
                      help='Elapsed flight time of the record to convert, s')
  parser.add_argument('--horizon', type=str, default='geodetic',
                      choices=['geodetic', 'geocentric'],
                      help='Horizon the flight-path angle and the azimuth refer to')
  args = parser.parse_args()

  os.chdir(DIRECTORY_SCRIPT)
  record = read_record(FILE_TRAJECTORY, args.time)
  for name in ('latitude', 'longitude', 'altitude', 'velocity', 'flightpath', 'heading'):
    if record[name] != record[name] :
      print('The record at t = {:.2f} s has no {:s} (the scan lost it).'.format(
            args.time, name))
      print('--Pick another time with --time.')
      print('Program stopped.')
      sys.exit(1)

  config        = general().read_config_yaml(args.file_config)
  longitude     = record['longitude']*np.pi/180.0
  latitude      = record['latitude']*np.pi/180.0
  altitude      = record['altitude']
  coordinate    = np.array(coordinate_system.convert_geodetic_cartesian(
                    config, [longitude, latitude, altitude]))

  # 地心の経度・緯度（ローカル系の基準）
  polar                = coordinate_system.set_angle_polar(config, coordinate)
  longitude_geocentric = polar[2]
  latitude_geocentric  = polar[1]

  # 経路角・方位角を測る地平線
  if args.horizon == 'geodetic' :
    latitude_horizon = latitude
  else:
    latitude_horizon = latitude_geocentric

  # その地平線の [東, 北, 上] の ECEF 成分
  direction_east  = np.array([-np.sin(longitude_geocentric), np.cos(longitude_geocentric), 0.0])
  direction_up    = np.array([np.cos(latitude_horizon)*np.cos(longitude_geocentric),
                              np.cos(latitude_horizon)*np.sin(longitude_geocentric),
                              np.sin(latitude_horizon)])
  direction_north = np.cross(direction_up, direction_east)

  speed         = record['velocity']
  angle_path    = record['flightpath']*np.pi/180.0
  angle_azimuth = record['heading']*np.pi/180.0
  velocity_horizontal = speed*np.cos(angle_path)
  velocity_horizon    = np.array([velocity_horizontal*np.sin(angle_azimuth),
                                  velocity_horizontal*np.cos(angle_azimuth),
                                  speed*np.sin(angle_path)])

  velocity_relative = velocity_horizon[0]*direction_east \
                    + velocity_horizon[1]*direction_north \
                    + velocity_horizon[2]*direction_up

  # 地心ローカル系 [東, 北, 上] へ
  velocity_local = np.array(coordinate_system.convert_carteasian_polar(
                     config, velocity_relative, longitude_geocentric, latitude_geocentric))

  angle_path_geocentric    = np.arcsin(velocity_local[2]/speed)*180.0/np.pi
  angle_azimuth_geocentric = np.arctan2(velocity_local[0], velocity_local[1])*180.0/np.pi

  print('Project Fire flight II, NASA TN D-3569 table V at t = {:.2f} s'.format(args.time))
  print('--Horizon of the tabulated angles : ' + args.horizon)
  print('--Position (geodetic)             : {:.4f} deg. E, {:.4f} deg. N, {:.4f} km'.format(
        record['longitude'], record['latitude'], altitude*1.0e-3))
  print('--Geocentric latitude             : {:.4f} deg.'.format(latitude_geocentric*180.0/np.pi))
  print('--Relative (ECEF) speed           : {:.2f} m/s'.format(speed))
  print('--Relative flight-path angle      : {:.3f} deg. (tabulated, {:s} horizon)'.format(
        record['flightpath'], args.horizon))
  print('--                                  {:.3f} deg. (geocentric horizon)'.format(
        angle_path_geocentric))
  print('--Relative azimuth                : {:.3f} deg. (tabulated)'.format(record['heading']))
  print('--                                  {:.3f} deg. (geocentric horizon)'.format(
        angle_azimuth_geocentric))
  print('--Lift-off                        : ' + EPOCH_LIFTOFF)
  print('')
  print('  coordinate:')
  print('    - {:.4f}'.format(record['longitude']))
  print('    - {:.4f}'.format(record['latitude']))
  print('    - {:.4f}'.format(altitude*1.0e-3))
  print('  velocity:')
  for component in velocity_local:
    print('    - {:.2f}'.format(component))

  return


if __name__ == '__main__':
  main()
