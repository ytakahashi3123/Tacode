#!/usr/bin/env python3

import numpy as np
import os as os
import shutil as shutil
import random as random
from orbital.orbital import orbital
from montecarlo.montecarlo import montecarlo

def main():

  # 設定ファイルの読み込み
  file_control_default = orbital.file_control_default
  arg                  = orbital.argument(file_control_default)
  file_control         = arg.file
  config               = orbital.read_config_yaml(file_control)

  # Initial setting
  montecarlo.initial_settings(config)

  for n in range(0,10):
    montecarlo.f_tacode(config)

  return


if __name__ == '__main__':

  print('Initializing Tacode-MonteCarlo')

  # Call classes
  orbital = orbital()

  montecarlo = montecarlo()

  # Main
  main()

  print('Finalizing Tacode-MonteCarlo')

  exit()