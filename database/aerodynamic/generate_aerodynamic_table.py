#!/usr/bin/env python3
#
# 迎角依存の空力係数テーブルを作る（6 自由度計算のサンプルデータ）。
#
# 計測値でも DSMC でもなく、球円錐（sphere-cone）形状に対する解析モデルの出力である。
#
#   連続流側 : 修正ニュートン流理論      Cp = Cp_max sin^2(delta), Ctau = 0
#   自由分子流側: 極超音速極限の完全運動量交換  Cp = 2 sin^2(delta), Ctau = 2 sin(delta) cos(delta)
#   遷移域   : Knudsen 数による sin^2 のブリッジ関数で内挿
#              f = sin^2( pi/2 * (log10(Kn) + 3)/4 )   (1e-3 <= Kn <= 10)
#
# delta は局所の流れの傾斜角（面が流れに向いていれば正）。凸形状なので、
# 影の判定は sin(delta) > 0 だけで足りる。
#
# 機体軸は [前方, 右, 下]。出力の符号は Tacode の規約に合わせてある。
#
#   力      : F_body = -q S [CFx, CFy, CFz]   （迎角 0 で CFx = CD）
#   モーメント: M_body =  q S L [CMx, CMy, CMz]（重心まわり）
#
# 使い方:
#   python3 generate_aerodynamic_table.py -o aerodynamic_spherecone_aoa.txt
#

import argparse
import numpy as np
import os
import sys

# 形状（既定値は tutorial/work_reentry に合わせた直径 1 m の球円錐）
RADIUS_BASE      = 0.5      # 底面半径, m
RADIUS_NOSE      = 0.25     # ノーズ半径, m
ANGLE_CONE_HALF  = 45.0     # 円錐半頂角, deg
POSITION_CG      = 0.20     # 重心位置（底面からの軸方向距離）, m

# 基準量（config.yml の satellite と揃えること）
AREA_REFERENCE   = np.pi*RADIUS_BASE**2   # 底面積, m2
LENGTH_REFERENCE = RADIUS_BASE            # 基準長, m

# 修正ニュートン流の最大圧力係数（比熱比 1.4、極超音速）
CP_MAXIMUM = 1.839

# ブリッジ関数の範囲
KNUDSEN_CONTINUUM   = 1.e-3
KNUDSEN_FREEMOLECULE = 1.e+1

# 分割数
NUM_PANEL_MERIDIONAL = 200
NUM_PANEL_AZIMUTHAL  = 180


def make_panel():
  #
  # 回転体の表面をパネルに分割する。
  # 戻り値: 面心座標 (n,3), 外向き単位法線 (n,3), 面積 (n)
  #
  angle_cone = ANGLE_CONE_HALF*np.pi/180.0

  # 球と円錐の接点（軸からの極角 phi_t = 90 - 円錐半頂角）
  angle_tangent = 0.5*np.pi - angle_cone
  radius_tangent   = RADIUS_NOSE*np.sin(angle_tangent)
  position_tangent = (RADIUS_BASE - radius_tangent)/np.tan(angle_cone)
  position_center  = position_tangent - RADIUS_NOSE*np.cos(angle_tangent)

  if position_center + RADIUS_NOSE <= 0.0 :
    print('The geometry is inconsistent. Check the nose radius and the cone angle.')
    sys.exit(1)

  # 母線: (x, r) の列。ノーズ先端から底面の縁まで
  generator = []
  num_sphere = int(NUM_PANEL_MERIDIONAL*0.4)
  num_cone   = NUM_PANEL_MERIDIONAL - num_sphere
  for index in range(0, num_sphere+1):
    angle_tmp = angle_tangent*float(index)/float(num_sphere)
    generator.append([position_center + RADIUS_NOSE*np.cos(angle_tmp), RADIUS_NOSE*np.sin(angle_tmp)])
  for index in range(1, num_cone+1):
    fact = float(index)/float(num_cone)
    generator.append([position_tangent*(1.0-fact), radius_tangent + (RADIUS_BASE-radius_tangent)*fact])
  generator = np.array(generator)

  angle_azimuth = np.linspace(0.0, 2.0*np.pi, NUM_PANEL_AZIMUTHAL, endpoint=False)
  delta_azimuth = 2.0*np.pi/float(NUM_PANEL_AZIMUTHAL)

  position = []
  normal   = []
  area     = []

  # 側面
  for index in range(0, len(generator)-1):
    position_1, radius_1 = generator[index]
    position_2, radius_2 = generator[index+1]
    length_slant = np.sqrt((position_2-position_1)**2 + (radius_2-radius_1)**2)
    position_mid = 0.5*(position_1+position_2)
    radius_mid   = 0.5*(radius_1+radius_2)
    # 母線に垂直な外向き法線（(x, r) 平面）
    normal_axis   =  (radius_2 - radius_1)/length_slant
    normal_radial = -(position_2 - position_1)/length_slant
    for angle_tmp in angle_azimuth:
      sin_psi = np.sin(angle_tmp)
      cos_psi = np.cos(angle_tmp)
      position.append([position_mid, radius_mid*sin_psi, radius_mid*cos_psi])
      normal.append([normal_axis, normal_radial*sin_psi, normal_radial*cos_psi])
      area.append(radius_mid*delta_azimuth*length_slant)

  # 底面（x = 0 の円板、外向き法線は -x）
  num_radial = 20
  for index in range(0, num_radial):
    radius_1 = RADIUS_BASE*float(index)/float(num_radial)
    radius_2 = RADIUS_BASE*float(index+1)/float(num_radial)
    radius_mid = 0.5*(radius_1+radius_2)
    area_tmp   = 0.5*(radius_2**2 - radius_1**2)*delta_azimuth
    for angle_tmp in angle_azimuth:
      position.append([0.0, radius_mid*np.sin(angle_tmp), radius_mid*np.cos(angle_tmp)])
      normal.append([-1.0, 0.0, 0.0])
      area.append(area_tmp)

  return np.array(position), np.array(normal), np.array(area)


def bridge_function(knudsen):
  # 連続流(0) と自由分子流(1) をつなぐ sin^2 ブリッジ関数
  if knudsen <= KNUDSEN_CONTINUUM :
    return 0.0
  if knudsen >= KNUDSEN_FREEMOLECULE :
    return 1.0
  fact = (np.log10(knudsen) - np.log10(KNUDSEN_CONTINUUM)) \
       / (np.log10(KNUDSEN_FREEMOLECULE) - np.log10(KNUDSEN_CONTINUUM))
  return np.sin(0.5*np.pi*fact)**2


def get_coefficient(position, normal, area, angle_attack, knudsen):
  #
  # 迎角 angle_attack (deg)、Knudsen 数 knudsen における係数を返す。
  # 速度は機体軸 x-z 面内にとる: V = V(cos(alpha), 0, sin(alpha))
  #
  angle_tmp = angle_attack*np.pi/180.0
  velocity_unit = np.array([np.cos(angle_tmp), 0.0, np.sin(angle_tmp)])

  fact_bridge = bridge_function(knudsen)

  sin_delta = np.dot(normal, velocity_unit)
  flag_wind = sin_delta > 0.0

  sin_delta = np.where(flag_wind, sin_delta, 0.0)
  cos_delta = np.sqrt(np.maximum(0.0, 1.0 - sin_delta**2))

  coefficient_pressure = (1.0-fact_bridge)*CP_MAXIMUM*sin_delta**2 + fact_bridge*2.0*sin_delta**2
  coefficient_shear    = fact_bridge*2.0*sin_delta*cos_delta

  # 面に沿う流れの方向（流れは -V 方向に進む）
  vector_tangent = -velocity_unit[np.newaxis,:] + sin_delta[:,np.newaxis]*normal
  norm_tangent   = np.linalg.norm(vector_tangent, axis=1)
  vector_tangent = np.where(norm_tangent[:,np.newaxis] > 1.e-12,
                            vector_tangent/np.maximum(norm_tangent, 1.e-12)[:,np.newaxis], 0.0)

  # 動圧で規格化した力（q = 1 としたもの）
  force_panel = area[:,np.newaxis]*( -coefficient_pressure[:,np.newaxis]*normal
                                    + coefficient_shear[:,np.newaxis]*vector_tangent )

  force_total = np.sum(force_panel, axis=0)

  arm = position - np.array([POSITION_CG, 0.0, 0.0])[np.newaxis,:]
  moment_total = np.sum(np.cross(arm, force_panel), axis=0)

  # Tacode の符号規約に合わせる（F = -q S CF, M = q S L CM）
  coefficient_force  = -force_total/AREA_REFERENCE
  coefficient_moment =  moment_total/(AREA_REFERENCE*LENGTH_REFERENCE)

  return coefficient_force, coefficient_moment


def get_altitude(knudsen_list):
  # Kn に対応する高度を大気テーブルから引く（出力の Altitude 列は参考値）
  script_directory = os.path.dirname(os.path.realpath(__file__))
  sys.path.insert(0, os.path.join(script_directory, '..', '..', 'src'))
  try:
    import atmosphere.atmosphere as atmosphere
  except ImportError:
    return [0.0 for knudsen in knudsen_list]

  config = {'satellite'  : {'characteristic_length': LENGTH_REFERENCE},
            'atmosphere' : {'directory_path_specify': 'default',
                            'filename_atmosphere'   : 'atmospheremodel.txt',
                            'kind_extrapolation'    : 'exponential'}}
  atmosphere_dict = atmosphere.initial_settings_atmosphere(config)

  altitude_table = atmosphere_dict[atmosphere.KEY_Height]
  knudsen_table  = atmosphere_dict[atmosphere.KEY_KN]

  return [float(np.interp(knudsen, knudsen_table, altitude_table)) for knudsen in knudsen_list]


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-o', '--output', type=str, default='aerodynamic_spherecone_aoa.txt')
  args = parser.parse_args()

  position, normal, area = make_panel()

  print('Number of panels:', len(area))
  print('Wetted area (m2):', np.sum(area))

  knudsen_list = [1.e-4, 1.e-3, 1.e-2, 1.e-1, 1.e+0, 1.e+1, 1.e+2, 1.e+3, 1.e+4, 1.e+5]
  angle_list   = list(np.arange(0.0, 181.0, 10.0))

  altitude_list = get_altitude(knudsen_list)

  with open(args.output, 'w') as file:
    file.write('Sphere-cone aerodynamic data (modified Newtonian + free-molecular, analytic sample)\n')
    file.write('variables = Kn, CFx, CFy, CFz, CMx, CMy, CMz, SDV_CFx, SDV_CFy, SDV_CFz, SDV_CMx, SDV_CMy, SDV_CMz, Altitude \n')
    for angle_attack in angle_list:
      file.write('AOA {:g}\n'.format(angle_attack))
      for index, knudsen in enumerate(knudsen_list):
        coefficient_force, coefficient_moment = get_coefficient(position, normal, area, angle_attack, knudsen)
        values = list(coefficient_force) + list(coefficient_moment) + [0.0]*6 + [altitude_list[index]]
        file.write('{:.18e}'.format(knudsen) + '\t'
                   + '\t'.join(['{:.18e}'.format(value) for value in values]) + '\n')

  print('Written:', args.output)

  # 静安定の確認（迎角 0 近傍の Cm 傾き）
  for knudsen in [1.e-3, 1.e+0, 1.e+4]:
    _, moment_zero = get_coefficient(position, normal, area, 0.0, knudsen)
    _, moment_ref  = get_coefficient(position, normal, area, 10.0, knudsen)
    slope = (moment_ref[1] - moment_zero[1])/(10.0*np.pi/180.0)
    force_zero, _ = get_coefficient(position, normal, area, 0.0, knudsen)
    print('Kn = {:9.1e}: CD(0 deg) = {:8.4f}, Cm_alpha = {:9.4f} /rad ({:s})'.format(
          knudsen, force_zero[0], slope, 'stable' if slope < 0.0 else 'UNSTABLE'))


if __name__ == '__main__':
  main()
