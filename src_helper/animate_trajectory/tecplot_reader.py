#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/06
#
# Tacode が書いた Tecplot 形式の出力を読む。
# 計算はしない。src/ のソルバーには触らず、出力ファイルだけを見る。

import numpy as np
import os
import sys

# 姿勢（6 自由度）計算のときだけ現れる列
COLUMN_ATTITUDE = ['q0', 'q1', 'q2', 'q3',
                   'Yaw', 'Pitch', 'Roll',
                   'P', 'Q', 'R',
                   'AoA', 'Sideslip', 'AoAtotal']


def strip_unit(name):
  # 'Alti[km]' -> 'Alti'、'Long[deg.]' -> 'Long'
  name = name.strip()
  if '[' in name :
    name = name.split('[')[0]
  return name.strip()


def read_tecplot(filename):
  #
  # Tacode の Tecplot point 形式を読み、列名をキーにした辞書で返す。
  #
  #   data['Time'], data['Alti'], data['q0'], ...
  #
  # 3 自由度の出力（姿勢の列が無いもの）もそのまま読める。
  #
  if not os.path.exists(filename) :
    print('File not found:', filename)
    sys.exit(1)

  name_list = None
  row_list  = []

  with open(filename) as f:
    for line in f:
      stripped = line.strip()
      if len(stripped) == 0 :
        continue
      if stripped.lower().startswith('variables') :
        body = stripped.split('=', 1)[1]
        name_list = [strip_unit(word) for word in body.split(',')]
        continue
      if stripped.startswith('#') or stripped.lower().startswith('zone') :
        continue

      word_list = stripped.split()
      try:
        row_list.append([float(word) for word in word_list])
      except ValueError:
        continue

  if name_list is None :
    print('No "Variables" line was found in', filename)
    sys.exit(1)
  if len(row_list) == 0 :
    print('No data row was found in', filename)
    sys.exit(1)

  table = np.array(row_list)
  if table.shape[1] != len(name_list) :
    print('The number of columns ({:d}) does not match the Variables line ({:d}).'.format(
          table.shape[1], len(name_list)))
    sys.exit(1)

  data = {}
  for index, name in enumerate(name_list):
    data[name] = table[:, index]

  return data


def has_attitude(data):
  # 6 自由度の出力かどうか
  return all(name in data for name in COLUMN_ATTITUDE)


def quaternion_array(data):
  return np.stack([data['q0'], data['q1'], data['q2'], data['q3']], axis=1)


def position_array(data):
  # ECEF 直交座標 [km]
  return np.stack([data['X'], data['Y'], data['Z']], axis=1)


def velocity_geocentric_array(data):
  # 地心ローカル系 [東, 北, 上] の速度 [m/s]
  return np.stack([data['Upl'], data['Vpl'], data['Wpl']], axis=1)
