# coding:utf-8
#!/usr/bin/env python3

# Tacode: Trajectory analysis code
# Version 2.5.1

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2026/09/14


import sys as sys
import numpy as np
from orbital.orbital import orbital
import atmosphere.atmosphere as atmosphere
import epoch.epoch as epoch
import satellite.satellite as satellite
import solver.solver as solver
import wind.wind as wind
import output_gpsdata.output_gpsdata as output_gpsdata


def main():

  # 設定ファイルの読み込み
  file_control_default = orbital.file_control_default
  arg                  = orbital.argument(file_control_default)
  file_control         = arg.file
  config               = orbital.read_config_yaml(file_control)

  # Set the epoch (absolute time)
  # --None が返れば従来どおり経過秒だけで、出力にも時刻は載らない
  epoch_dict = epoch.initial_settings_epoch(config)

  # Set amosphere parameter
  atmosphere_dict = atmosphere.initial_settings_atmosphere(config)

  # Set satellite parameters
  aerodynamic_dict = satellite.initial_settings_satellite(config)

  # Set the wind
  # --None が返れば大気は地球と共回転し、ECEF 速度がそのまま対気速度になる（従来の仮定）
  wind_dict = wind.initial_settings_wind(config, epoch_dict)

  # Make directories for output data
  orbital.make_directory_output(config)

  # Initial setting
  # --epoch_dict はリスタートのときだけ使う（ファイルのエポックとの食い違いを見る）
  iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = orbital.initial_settings(config, epoch_dict)

  # 出力の時刻列の起点。リスタートでは 0 ではなく、再開時刻から続ける
  time_elapsed_initial = time_elapsed

  # Initial setting of the attitude
  # --None が返れば従来どおりの質点 3 自由度計算
  attitude_dict = orbital.initial_settings_attitude(config, coordinate_dict, velocity_dict)

  # Main routine
  # --attitude_dict は他の辞書と同じくソルバー内で更新される
  iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = solver.solve_equation_motion(config, iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict, atmosphere_dict, aerodynamic_dict, attitude_dict, wind_dict)

  # Output : Geodetic data
  output_gpsdata.output_routine(config, iteration, coordinate_dict, velocity_dict, time_elapsed_initial)

  # Output : Tecplot
  orbital.output_tecplot(config, iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict, attitude_dict, epoch_dict, wind_dict, time_elapsed_initial)

  # Output restart
  orbital.output_restart(config, iteration, time_elapsed, coordinate_dict['cartesian'], velocity_dict['cartesian'], attitude_dict, epoch_dict)

  return


if __name__ == '__main__':

  print('Initializing Tacode')

  # Calling classes
  orbital  = orbital()

  # main 
  main()

  print('Finalizing Tacode')

  sys.exit(0)
