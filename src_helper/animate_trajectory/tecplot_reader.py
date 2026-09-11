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


def read_tecplot_zone(filename):
  #
  # Tacode の Tecplot point 形式を読み、**ゾーンごとに**列名をキーにした辞書を返す。
  #
  #   zone_list[0]['Time'], zone_list[0]['Alti'], ...
  #
  # 3 自由度の出力（姿勢の列が無いもの）もそのまま読める。
  # モンテカルロの postprocess がまとめたファイルは 1 ケース = 1 ゾーンなので、
  # ここで分けておかないと別ケースの行が 1 本の軌跡として繋がる。
  #
  if not os.path.exists(filename) :
    print('File not found:', filename)
    sys.exit(1)

  name_list = None
  zone_row  = []
  row_list  = None

  with open(filename) as f:
    for line in f:
      stripped = line.strip()
      if len(stripped) == 0 :
        continue
      if stripped.lower().startswith('variables') :
        body = stripped.split('=', 1)[1]
        name_list = [strip_unit(word) for word in body.split(',')]
        continue
      if stripped.lower().startswith('zone') :
        row_list = []
        zone_row.append(row_list)
        continue
      if stripped.startswith('#') :
        continue

      word_list = stripped.split()
      try:
        value_list = [float(word) for word in word_list]
      except ValueError:
        continue
      # zone 行が無いファイル（あるいは行より前のデータ）も 1 つのゾーンとして扱う
      if row_list is None :
        row_list = []
        zone_row.append(row_list)
      row_list.append(value_list)

  if name_list is None :
    print('No "Variables" line was found in', filename)
    sys.exit(1)

  zone_row = [rows for rows in zone_row if len(rows) > 0]
  if len(zone_row) == 0 :
    print('No data row was found in', filename)
    sys.exit(1)

  zone_list = []
  for rows in zone_row:
    table = np.array(rows)
    if table.shape[1] != len(name_list) :
      print('The number of columns ({:d}) does not match the Variables line ({:d}).'.format(
            table.shape[1], len(name_list)))
      sys.exit(1)
    zone_list.append({name: table[:, index] for index, name in enumerate(name_list)})

  return zone_list


def read_tecplot(filename):
  #
  # ゾーンが 1 つだけの Tecplot ファイルを読む（Tacode の 1 ケースの出力）。
  #
  # **複数ゾーンのファイルは受け取らない。** 以前は zone 行を読み飛ばしていたので、
  # モンテカルロのまとめ出力を渡すと別々のケースが 1 本の軌跡として繋がり、
  # 終端点が「最後のケースの最後の点」になっていた（黙って間違う）。
  #
  zone_list = read_tecplot_zone(filename)

  if len(zone_list) > 1 :
    print('The file holds {:d} zones, and this reader takes one:'.format(len(zone_list)), filename)
    print('--A gathered Monte-Carlo result has one zone per case. Read the case files')
    print('--themselves, or use read_tecplot_zone() to take the zones apart.')
    print('Program stopped.')
    sys.exit(1)

  return zone_list[0]


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
