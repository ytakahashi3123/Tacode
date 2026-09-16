#!/usr/bin/env python3

import sys as sys
import numpy as np
import os as os
import shutil as shutil
import random as random
from orbital.orbital import orbital
from montecarlo.montecarlo import montecarlo

def main():

  # orbital と montecarlo はクラス。インスタンスは局所に持つ（クラス名を
  # 潰すと、クラス属性として参照している側と読み分けが要る）
  orb = orbital()
  mc  = montecarlo()

  # 設定ファイルの読み込み
  file_control_default = orb.file_control_default
  arg                  = orb.argument(file_control_default)
  file_control         = arg.file
  config               = orb.read_config_yaml(file_control)

  # Initial setting
  mc.initial_settings(config)

  # Monte-Carlo simulation
  mc.montecarlo_routine(config)

  return


if __name__ == '__main__':

  print('Initializing Tacode-MonteCarlo')

  # Main
  main()

  print('Finalizing Tacode-MonteCarlo')

  sys.exit(0)
