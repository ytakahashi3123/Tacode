# coding:utf-8
#!/usr/bin/env python3

# Tacode: Trajectory analysis code
# Version 2.3.1

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2026/09/09


import sys as sys
import numpy as np
from orbital.orbital import orbital
import atmosphere.atmosphere as atmosphere
import satellite.satellite as satellite
import solver.solver as solver
import output_gpsdata.output_gpsdata as output_gpsdata


def main():

  # 設定ファイルの読み込み
  file_control_default = orbital.file_control_default
  arg                  = orbital.argument(file_control_default)
  file_control         = arg.file
  config               = orbital.read_config_yaml(file_control)

  # Set amosphere parameter
  atmosphere_dict = atmosphere.initial_settings_atmosphere(config)

  # Set satellite parameters
  aerodynamic_dict = satellite.initial_settings_satellite(config)

  # Make directories for output data
  orbital.make_directory_output(config)

  # Initial setting
  iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = orbital.initial_settings(config)

  # Initial setting of the attitude
  # --None が返れば従来どおりの質点 3 自由度計算
  attitude_dict = orbital.initial_settings_attitude(config, coordinate_dict, velocity_dict)

  # Main routine
  # --attitude_dict は他の辞書と同じくソルバー内で更新される
  iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = solver.solve_equation_motion(config, iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict, atmosphere_dict, aerodynamic_dict, attitude_dict)

  # Output : Geodetic data
  output_gpsdata.output_routine(config, iteration, coordinate_dict, velocity_dict)

  # Output : Tecplot
  orbital.output_tecplot(config, iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict, attitude_dict)

  # Output restart
  orbital.output_restart(config, iteration, time_elapsed, coordinate_dict['cartesian'], velocity_dict['cartesian'], attitude_dict)

  return


if __name__ == '__main__':

  print('Initializing Tacode')

  # Calling classes
  orbital  = orbital()

  # main 
  main()

  print('Finalizing Tacode')

  sys.exit(0)
