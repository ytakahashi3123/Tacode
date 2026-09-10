#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/11
#
# 3 次元の図に入れる海岸線。
#
# 同梱のテキスト（coastline_110m.txt。Natural Earth 1:110m、パブリックドメイン）を
# 読むだけで、**numpy 以外に何も要らない**。cartopy / shapely / pyproj を持ち込むと
# 計算環境に GEOS/PROJ まで背負わせることになるので、テキスト 1 枚に落としてある
# （作り直しは src_helper/general/generate_coastline.py。ネットワークが要るのはそちら）。
#
# 使う側は 3 つ:
#
#   read_coastline()      -- 折れ線（経度・緯度, deg）の並びを読む。1 度読めば足りる
#   select_range()        -- 窓（経度・緯度の範囲）に入る部分だけを切り出す
#   get_segment_ecef()    -- ECEF の 2 点セグメント（M, 2, 3）にする。
#                            Line3DCollection は**点数の揃ったセグメント**しか受けない
#
# 経度は [-180, 180) で入っているが、窓の範囲は 180 度をまたぐことがある
# （make_planet が中心 ± 余白で作るため）。select_range は ±360 度ずらして合わせる。

import numpy as np
import os as os

# 同梱のテキスト
PATH_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'coastline_110m.txt')

# 地球の赤道半径 [km]（描画用。animate_trajectory / montecarlo_animation と同じ値。
# 一致は test_helper_visualization.py が検査する）
RADIUS_PLANET = 6378.137


def read_coastline(filename=None):
  #
  # 折れ線の並びを返す。各要素は (N, 2) の配列で、列は [経度, 緯度]（deg）。
  #
  if filename is None :
    filename = PATH_DEFAULT

  segment_list = []
  point_list   = []
  with open(filename) as f:
    for line in f:
      text = line.strip()
      if text.startswith('#') :
        continue
      if text == '' :
        # 空行が折れ線の区切り
        if len(point_list) > 1 :
          segment_list.append(np.array(point_list))
        point_list = []
        continue
      word = text.split()
      point_list.append([float(word[0]), float(word[1])])

  if len(point_list) > 1 :
    segment_list.append(np.array(point_list))

  return segment_list


def select_range(segment_list, range_longitude, range_latitude):
  #
  # 窓に入る部分だけを切り出す。窓から出るところで折れ線を切る
  # （切らずに端の点を結ぶと、大陸を横断する直線が 1 本引かれてしまう）。
  #
  # 隣の点まで 1 つずつ延ばして残すので、窓の縁で線が途切れて見えることはない。
  #
  selected = []
  for segment in segment_list:
    longitude = fold_longitude(segment[:, 0], range_longitude)
    inside    = ( (longitude >= range_longitude[0]) & (longitude <= range_longitude[1])
                & (segment[:, 1] >= range_latitude[0]) & (segment[:, 1] <= range_latitude[1]) )
    if not np.any(inside) :
      continue

    # 窓に入る点が続いているところを 1 つの断片にし、**断片ごとに**両隣まで
    # 1 点ずつ延ばす（縁で線が途切れて見えないように）。
    #
    # 「入る点の両隣を残す」という書き方だと、窓の外を 2 点で通り過ぎる部分が
    # 前後とつながってしまい、窓を横切る直線が 1 本引かれる。
    point = np.column_stack([longitude, segment[:, 1]])
    index = np.nonzero(inside)[0]
    start = 0
    for i in range(1, len(index)+1) :
      if i == len(index) or index[i] != index[i-1] + 1 :
        first = max(index[start] - 1, 0)
        last  = min(index[i-1] + 1, len(point) - 1)
        if last > first :
          selected.append(point[first:last+1])
        start = i

  return selected


def fold_longitude(longitude, range_longitude):
  #
  # 窓の範囲に合わせて経度を ±360 度ずらす（範囲は 180 度をまたぐことがある）。
  #
  centre = 0.5*(range_longitude[0] + range_longitude[1])

  return longitude + 360.0*np.round((centre - longitude)/360.0)


def get_view_direction(elevation, azimuth):
  #
  # カメラの向き（単位ベクトル、ECEF）。matplotlib の view_init と同じ規約:
  # 方位角は x 軸から y 軸へ、仰角は xy 平面から z 軸へ。
  #
  elevation = np.radians(elevation)
  azimuth   = np.radians(azimuth)

  return np.array([np.cos(elevation)*np.cos(azimuth),
                   np.cos(elevation)*np.sin(azimuth),
                   np.sin(elevation)])


def select_visible(segment, direction_view):
  #
  # 地球の裏側にあるセグメントを落とす。
  #
  # matplotlib の 3 次元は面で線を隠してくれない（不透明な地表を描いても、
  # zorder が上の線はその上に出る）。裏側の海岸線がそのまま透けて、
  # 大陸が鏡像で重なって見えるので、法線がカメラを向いていないものを外す。
  #
  if len(segment) == 0 :
    return segment

  middle = 0.5*(segment[:, 0, :] + segment[:, 1, :])
  normal = middle/np.linalg.norm(middle, axis=1)[:, np.newaxis]

  return segment[np.dot(normal, direction_view) > 0.0]


def get_segment_ecef(segment_list, radius=RADIUS_PLANET):
  #
  # 折れ線を ECEF の 2 点セグメント (M, 2, 3) にする。
  #
  # Line3DCollection は点数の揃ったセグメントしか受けないので、2 点ずつに割る
  # （不揃いのまま渡すと ValueError）。
  #
  point_list = []
  for segment in segment_list:
    longitude = np.radians(segment[:, 0])
    latitude  = np.radians(segment[:, 1])
    point     = np.column_stack([radius*np.cos(latitude)*np.cos(longitude),
                                 radius*np.cos(latitude)*np.sin(longitude),
                                 radius*np.sin(latitude)])
    if len(point) > 1 :
      point_list.append(np.stack([point[:-1], point[1:]], axis=1))

  if len(point_list) == 0 :
    return np.zeros((0, 2, 3))

  return np.concatenate(point_list, axis=0)
