#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/03/31

import os as os
import sys as sys
import numpy as np


# データベースの既定の置き場。このファイルの位置から解決するので、
# カレントディレクトリがどこでもリポジトリ直下の database/ を指す
DIRECTORY_DATABASE = os.path.normpath(
  os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../database'))

# directory_path_specify に書ける値
LIST_PATH_SPECIFY = ['default', 'auto', 'manual']


def get_setting(section, key, default):
  # config の省略可能な項目を既定値つきで読む。
  # 姿勢・エポック・風のセクションはいずれも省略可能にしてあるので、
  # セクションごと無い場合も含めて既定値を返せるようにする。
  if section is None :
    return default
  try:
    value = section[key]
  except (KeyError, TypeError):
    return default
  if value is None :
    return default
  return value


def get_database_directory(section, name_section, name_database, key_directory, directory_manual_default=None):
  #
  # データベース（大気・空力・風）のディレクトリを解決する。**4 か所にあった同じ分岐を
  # ここ 1 つにまとめてある。**
  #
  #   default / auto : リポジトリ直下の database/<name_database>（スクリプト位置基準）
  #   manual         : section[key_directory]（カレントディレクトリ基準の相対パス）
  #   それ以外       : 止める
  #
  # **綴り違いを既定へ落とさない。** 以前は 'Manual' のような値が黙って default 扱いに
  # なり、ケースが自分の database/ に置いたテーブルではなくマスターを読んでいた
  # （しかも分岐ごとに既定のパスが違い、片方は綴りを間違えたディレクトリを指していた）。
  #
  path_specify = get_setting(section, 'directory_path_specify', 'default')

  if path_specify not in LIST_PATH_SPECIFY :
    print('directory_path_specify in the section "' + name_section + '" is incorrect:', path_specify)
    print('--Give one of:', ', '.join(LIST_PATH_SPECIFY))
    print('Program stopped.')
    sys.exit(1)

  if path_specify == 'manual' :
    directory_path = get_setting(section, key_directory, directory_manual_default)
    if directory_path is None :
      print('"' + key_directory + '" is missing in the section "' + name_section + '".')
      print('--It is required by directory_path_specify: manual (a path relative to the')
      print('--current directory), or use directory_path_specify: default instead.')
      print('Program stopped.')
      sys.exit(1)
    return directory_path

  return os.path.join(DIRECTORY_DATABASE, name_database)


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

    # 親ディレクトリごと作る。os.mkdir だと directory_output: nested/output_result で
    # FileNotFoundError の生のトレースバックになり、他のエラー経路と流儀が揃わない
    os.makedirs(dir_path, exist_ok=True)

    return


  def make_directory_rm(self, dir_path):

    import shutil as shutil

    if os.path.exists(dir_path):
      shutil.rmtree(dir_path)
    os.makedirs(dir_path)

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