#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/03/31

import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import numpy as np


class general:

  def __init__(self):
    print("Calling class: general")

# FUnctions
  def argument(self, filename_default):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', action='store', type=str, default=filename_default)
    args = parser.parse_args()
    return args


  def read_config_yaml(self, file_control):

    import yaml as yaml
    import sys as sys
    #import pprint as pprint

    print("Reading control file...:", file_control)

    try:
      with open(file_control) as file:
        config = yaml.safe_load(file)
#        pprint.pprint(config)
    except Exception as e:
      print('Exception occurred while loading YAML...', file=sys.stderr)
      print(e, file=sys.stderr)
      sys.exit(1)

    return config


  def make_directory(self, dir_path):
  
    import os as os
    import shutil as shutil

    if not os.path.exists(dir_path):
      os.mkdir(dir_path)
    
    return


  def make_directory_rm(self, dir_path):
  
    import os as os
    import shutil as shutil

    if not os.path.exists(dir_path):
      os.mkdir(dir_path)
    else:
      shutil.rmtree(dir_path)
      os.mkdir(dir_path)

    return
    

  def check_file_exist(self, dir_path):
  
    import os as os

    if os.path.exists(dir_path):
      flag_file_exist = True
    else: 
      flag_file_exist = False

    return flag_file_exist


  def split_file(self, filename,addfile,splitchar):
    """
    特定の文字の前に'_***'を加える.
    特定文字列が２つ以上ある場合は未対応
    """
  #  import re
  #  splitchar_tmp   ='['+splitchar+']'
  #  filename_split  = re.split(splitchar_tmp, filename)
  #  filename_result = filename_split[0]+addfile+splitchar+filename_split[1]
    splitchar_tmp   = splitchar
    filename_split  = filename.rsplit(splitchar_tmp, 1)
    filename_result = filename_split[0]+addfile+splitchar+filename_split[1]

    return filename_result

    
  def getNearestValue(self, list, num):
      # copied from https://qiita.com/icchi_h/items/fc0df3abb02b51f81657
      """
      概要: リストからある値に最も近い値を返却する関数
      @param list: データ配列
      @param num: 対象値
      @return 対象値に最も近い値
      """
      # リスト要素と対象値の差分を計算し最小値のインデックスを取得
      idx = np.abs(np.asarray(list) - num).argmin()
      return list[idx]
  
  
  def getNearestIndex(self, list, num):
      # リスト要素と対象値の差分を計算し最小値のインデックスを取得し返す
      idx = np.abs(np.asarray(list) - num).argmin()
      return idx