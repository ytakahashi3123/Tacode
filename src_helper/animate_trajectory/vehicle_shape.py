#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/06
#
# 姿勢を目で見るための機体形状。
#
# 座標は Tacode の機体軸 [前方 x, 右 y, 下 z] で作る。大きさは最大寸法が
# おおよそ 1 になるよう正規化してあり、描画側で好きな倍率に拡大する。
#
# 空力計算には一切使わない。見た目のためだけのもの。

import numpy as np

# 形状の名前
KIND_CAPSULE   = 'capsule'
KIND_SATELLITE = 'satellite'
KIND_AIRCRAFT  = 'aircraft'

LIST_KIND = [KIND_CAPSULE, KIND_SATELLITE, KIND_AIRCRAFT]


def make_box(center, half_size, color):
  # 直方体の 6 面
  cx, cy, cz = center
  hx, hy, hz = half_size

  corner = np.array([[cx-hx, cy-hy, cz-hz], [cx+hx, cy-hy, cz-hz],
                     [cx+hx, cy+hy, cz-hz], [cx-hx, cy+hy, cz-hz],
                     [cx-hx, cy-hy, cz+hz], [cx+hx, cy-hy, cz+hz],
                     [cx+hx, cy+hy, cz+hz], [cx-hx, cy+hy, cz+hz]])

  index_face = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4],
                [2, 3, 7, 6], [1, 2, 6, 5], [0, 3, 7, 4]]

  return [(corner[index], color) for index in index_face]


def make_revolution(generator, color, num_azimuth=24, color_marker=None, angle_marker=np.pi):
  #
  # 母線 generator = [(x, r), ...] を x 軸まわりに回した面。
  # r = 0 の端は自動的に三角形になる。
  #
  # color_marker を与えると、方位角 angle_marker のまわりの数枚だけ色を変える。
  # 回転体は形だけではロールが分からないので、その基準線になる。
  #
  angle = np.linspace(0.0, 2.0*np.pi, num_azimuth, endpoint=False)
  sin_a = np.sin(angle)
  cos_a = np.cos(angle)

  face = []
  for index in range(0, len(generator)-1):
    position_1, radius_1 = generator[index]
    position_2, radius_2 = generator[index+1]
    for m in range(0, num_azimuth):
      n = (m+1) % num_azimuth
      point = [[position_1, radius_1*sin_a[m], radius_1*cos_a[m]],
               [position_1, radius_1*sin_a[n], radius_1*cos_a[n]],
               [position_2, radius_2*sin_a[n], radius_2*cos_a[n]],
               [position_2, radius_2*sin_a[m], radius_2*cos_a[m]]]
      color_tmp = color
      if color_marker is not None :
        difference = np.abs((angle[m] - angle_marker + np.pi) % (2.0*np.pi) - np.pi)
        if difference < 1.5*(2.0*np.pi/float(num_azimuth)) :
          color_tmp = color_marker
      # 半径 0 の側で潰れた点を落とす
      unique = [point[0]]
      for candidate in point[1:]:
        if not np.allclose(candidate, unique[-1]) :
          unique.append(candidate)
      if len(unique) >= 3 :
        face.append((np.array(unique), color_tmp))

  return face


def make_plate(corner, color):
  # 平板（四隅を与える）
  return [(np.array(corner), color)]


def shape_capsule():
  #
  # 球円錐カプセル（database/aerodynamic の球円錐テーブルに合わせた見た目）。
  # ノーズが +x（前方）。
  #
  color_shell  = '#d95f02'
  color_base   = '#7f4010'
  color_marker = '#f7f2e0'   # ロールの基準になる明色のストライプ（機体 -z 側 = 上）

  radius_base = 0.5
  radius_nose = 0.25
  angle_cone  = np.deg2rad(45.0)

  angle_tangent    = 0.5*np.pi - angle_cone
  radius_tangent   = radius_nose*np.sin(angle_tangent)
  position_tangent = (radius_base - radius_tangent)/np.tan(angle_cone)
  position_center  = position_tangent - radius_nose*np.cos(angle_tangent)

  generator = []
  num_nose = 8
  for index in range(0, num_nose+1):
    angle_tmp = angle_tangent*float(index)/float(num_nose)
    generator.append((position_center + radius_nose*np.cos(angle_tmp),
                      radius_nose*np.sin(angle_tmp)))
  generator.append((0.0, radius_base))
  # 母線は先端から底面へ。全体を重心（およそ 0.2）まわりに置き直す
  generator = [(x - 0.2, r) for x, r in generator]

  face = make_revolution(generator, color_shell, color_marker=color_marker)
  face += make_revolution([(-0.2, radius_base), (-0.2, 0.0)], color_base)

  return face


def shape_satellite():
  #
  # 箱型バス＋太陽電池パドル。パドルは機体 y 軸（右）に伸びる。
  #
  color_bus   = '#c9c9c9'
  color_panel = '#1f4e9c'
  color_dish  = '#e6b800'

  face = make_box((0.0, 0.0, 0.0), (0.30, 0.22, 0.22), color_bus)

  for sign in (1.0, -1.0):
    face += make_plate([[-0.30, sign*0.24, -0.01], [0.30, sign*0.24, -0.01],
                        [0.30, sign*1.00, -0.01], [-0.30, sign*1.00, -0.01]], color_panel)
    face += make_plate([[-0.30, sign*0.24, 0.01], [0.30, sign*0.24, 0.01],
                        [0.30, sign*1.00, 0.01], [-0.30, sign*1.00, 0.01]], color_panel)

  # 前方のアンテナ（機首方向が分かるように）
  face += make_revolution([(0.30, 0.03), (0.55, 0.03)], color_dish)
  face += make_revolution([(0.55, 0.02), (0.70, 0.18)], color_dish)

  # 下面（+z）のセンサー。上下が分かるように
  face += make_box((0.0, 0.0, 0.26), (0.10, 0.10, 0.06), '#404040')

  return face


def shape_aircraft():
  #
  # 主翼・水平尾翼・垂直尾翼を持つ機体。上下・左右・前後がすべて見分けられる。
  #
  color_body = '#e8e8e8'
  color_wing = '#2166ac'
  color_tail = '#b2182b'

  generator = [(-0.50, 0.02), (-0.40, 0.09), (0.10, 0.10), (0.40, 0.07), (0.55, 0.0)]
  face = make_revolution(generator, color_body)

  # 主翼（後退角つき）
  for sign in (1.0, -1.0):
    face += make_plate([[0.10, 0.0, 0.0], [-0.12, 0.0, 0.0],
                        [-0.22, sign*0.55, 0.0], [-0.08, sign*0.55, 0.0]], color_wing)

  # 水平尾翼
  for sign in (1.0, -1.0):
    face += make_plate([[-0.34, 0.0, 0.0], [-0.48, 0.0, 0.0],
                        [-0.50, sign*0.22, 0.0], [-0.40, sign*0.22, 0.0]], color_tail)

  # 垂直尾翼（-z が上）
  face += make_plate([[-0.32, 0.0, 0.0], [-0.50, 0.0, 0.0],
                      [-0.50, 0.0, -0.28], [-0.36, 0.0, -0.28]], color_tail)

  return face


def get_shape(kind):
  if kind == KIND_CAPSULE :
    return shape_capsule()
  if kind == KIND_SATELLITE :
    return shape_satellite()
  if kind == KIND_AIRCRAFT :
    return shape_aircraft()

  print('Unknown shape:', kind)
  print('--Available:', ', '.join(LIST_KIND))
  raise SystemExit(1)


def get_vertices_and_colors(kind):
  # 描画側で扱いやすい形（頂点配列のリストと色のリスト）にして返す
  face = get_shape(kind)
  return [np.array(vertex) for vertex, color in face], [color for vertex, color in face]
