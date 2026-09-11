#!/usr/bin/env python3
#
# Apollo 10 の突入界面の状態ベクトル（NASA TN D-6725 表 II）を config.yml が要る形に直す。
#
# 表 II が与えるのは**慣性系**の速度・経路角・方位角で、位置は測地座標:
#
#   Time 191:48:52.16 (GET)   Velocity 36 309.257 ft/sec   Flight-path angle -6.616 deg.
#   Azimuth 71.928 deg.       Longitude 174.244 deg. E     Latitude 23.652 deg. S
#   Altitude 406 441.29 ft
#
# Tacode の initial_settings.velocity は **ECEF（大気に対する）速度**を
# **地心**ローカル系 [東, 北, 上] で与える。変換は 3 段:
#
#   1. 経路角・方位角から**測地**ローカル系の慣性速度 [東, 北, 上] を作る
#      （経路角は地平線基準。測地と地心では鉛直が最大 0.19 度違うので、どちらを採るかで
#        鉛直成分が 30 m/s ほど変わる。表に断りが無いので測地を採り、
#        --horizon geocentric で他方も試せるようにしてある）
#   2. ECEF 成分に直し、地球の自転ぶん omega x r を引いて ECEF 速度にする
#   3. 地心ローカル系 [東, 北, 上] に直す
#
# 使い方:
#   python3 entry_state.py
#   python3 entry_state.py --horizon geocentric

import argparse
import os
import sys

import numpy as np

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src'))

from general.general import general                              # noqa: E402
import coordinate_system.coordinate_system as coordinate_system  # noqa: E402

FOOT_TO_METRE = 0.3048

# NASA TN D-6725 表 II（best-estimated trajectory）
VELOCITY_INERTIAL   = 36309.257*FOOT_TO_METRE     # m/s
ANGLE_PATH_INERTIAL = -6.616                      # deg.
ANGLE_AZIMUTH       = 71.928                      # deg.
LONGITUDE           = 174.244                     # deg. East
LATITUDE            = -23.652                     # deg. North（表は 23.652 deg. S）
ALTITUDE            = 406441.29*FOOT_TO_METRE     # m


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--file-config', type=str, default='config.yml',
                      help='Configuration to take the planet constants from')
  parser.add_argument('--horizon', type=str, default='geodetic', choices=['geodetic', 'geocentric'],
                      help='Horizon the flight-path angle and the azimuth of table II refer to')
  args = parser.parse_args()

  os.chdir(DIRECTORY_SCRIPT)
  config = general().read_config_yaml(args.file_config)
  rotation_rate = config['planet']['rotation_rate']

  longitude = LONGITUDE*np.pi/180.0
  latitude  = LATITUDE*np.pi/180.0
  coordinate = np.array(coordinate_system.convert_geodetic_cartesian(config, [longitude, latitude, ALTITUDE]))

  # 地心の経度・緯度（ローカル系の基準）
  polar = coordinate_system.set_angle_polar(config, coordinate)
  longitude_geocentric = polar[2]
  latitude_geocentric  = polar[1]

  # 経路角・方位角を測る地平線
  if args.horizon == 'geodetic' :
    latitude_horizon = latitude
  else :
    latitude_horizon = latitude_geocentric

  # その地平線の [東, 北, 上] の ECEF 成分
  direction_east  = np.array([-np.sin(longitude_geocentric), np.cos(longitude_geocentric), 0.0])
  direction_up    = np.array([np.cos(latitude_horizon)*np.cos(longitude_geocentric),
                              np.cos(latitude_horizon)*np.sin(longitude_geocentric),
                              np.sin(latitude_horizon)])
  direction_north = np.cross(direction_up, direction_east)

  angle_path    = ANGLE_PATH_INERTIAL*np.pi/180.0
  angle_azimuth = ANGLE_AZIMUTH*np.pi/180.0
  velocity_horizontal = VELOCITY_INERTIAL*np.cos(angle_path)
  velocity_inertial_local = np.array([velocity_horizontal*np.sin(angle_azimuth),
                                      velocity_horizontal*np.cos(angle_azimuth),
                                      VELOCITY_INERTIAL*np.sin(angle_path)])

  velocity_inertial = velocity_inertial_local[0]*direction_east \
                    + velocity_inertial_local[1]*direction_north \
                    + velocity_inertial_local[2]*direction_up

  # 自転を引いて ECEF 速度に
  velocity_rotation = np.cross(np.array([0.0, 0.0, rotation_rate]), coordinate)
  velocity_relative = velocity_inertial - velocity_rotation

  # 地心ローカル系 [東, 北, 上] へ
  velocity_local = coordinate_system.convert_carteasian_polar(
                     config, velocity_relative, longitude_geocentric, latitude_geocentric)
  velocity_local = np.array(velocity_local)

  speed_relative = np.linalg.norm(velocity_relative)
  angle_path_relative = np.arcsin(velocity_local[2]/speed_relative)*180.0/np.pi
  angle_azimuth_relative = np.arctan2(velocity_local[0], velocity_local[1])*180.0/np.pi

  print('Apollo 10 entry interface (NASA TN D-6725, table II)')
  print('--Horizon of the tabulated angles : ' + args.horizon)
  print('--Position (geodetic)             : {:.4f} deg. E, {:.4f} deg. N, {:.4f} km'.format(
        LONGITUDE, LATITUDE, ALTITUDE*1.0e-3))
  print('--Geocentric latitude             : {:.4f} deg.'.format(latitude_geocentric*180.0/np.pi))
  print('--Inertial speed                  : {:.2f} m/s'.format(VELOCITY_INERTIAL))
  print('--Earth rotation at the point     : {:.2f} m/s'.format(np.linalg.norm(velocity_rotation)))
  print('--Relative (ECEF) speed           : {:.2f} m/s'.format(speed_relative))
  print('--Relative flight-path angle      : {:.3f} deg. (geocentric horizon)'.format(angle_path_relative))
  print('--Relative azimuth                : {:.3f} deg.'.format(angle_azimuth_relative))
  print('')
  print('  coordinate:')
  print('    - {:.4f}'.format(LONGITUDE))
  print('    - {:.4f}'.format(LATITUDE))
  print('    - {:.4f}'.format(ALTITUDE*1.0e-3))
  print('  velocity:')
  for component in velocity_local:
    print('    - {:.2f}'.format(component))

  return


if __name__ == '__main__':
  main()
