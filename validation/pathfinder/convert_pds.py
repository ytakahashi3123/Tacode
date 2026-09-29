#!/usr/bin/env python3
#
# Mars Pathfinder の突入（1997-07-04）の PDS データを、この検証ケースが使う形に直す。
#
# 入力（リポジトリの外、../../../references/20260929_Mars_EDL/pathfinder_pds/）:
#   edl_ddr.tab   ASI/MET の再構成（MPFL-M-ASIMET-4-DDR-EDL-V1.0）。0.25 s ごとの
#                 時刻・高度・緯度・経度・密度・圧力・温度・平均分子量と不確かさ
#   r_sacc_s.tab  科学用加速度計の 32 Hz の較正済みデータ（RAM に記録したもの）
#
# 出力（このディレクトリの下。どれも git に入れる）:
#   reference/pathfinder_ddr.dat            edl_ddr.tab を SI と Tacode の高度に直したもの
#   reference/pathfinder_acceleration.dat   軸方向の加速度 m/s2（レンジ切り替えの過渡を除いたもの）
#   reference/pathfinder_reconstruction.dat 加速度から組み直した速さと、逆算した CD・Mach・Kn
#   database/atmosphere/atmospheremodel_pathfinder_ddr.txt  飛行で測った大気（Tacode の形式）
#   database/aerodynamic/aerodynamic_pathfinder.txt         飛行から逆算した CD(Mach, Kn)
#
# **何が実測で、何がそうでないか**（README の「参照データの層」と同じ）:
#   A 実測       加速度（32 Hz）、軌道決定による突入状態
#   B 再構成     edl_ddr.tab の軌道と大気。A と空力データベース C から NASA が作ったもの
#   C 空力       Braun et al. (1995) の事前データベース。**数表は公開されていない**。
#                B の密度は rho = 2 m a / (CD A V^2) で作られているので、A と B から
#                CD = 2 m a / (rho A V^2) として**再構成が使った CD を逆にたどれる**
#
# 逆算に要る速さ V は edl_ddr.tab に無い（位置は緯度・経度が 0.001 度刻みで、微分すると
# ±数 % 暴れる）。そこで、文書の突入状態から実測加速度を積分して組み直す。これは
# NASA の再構成と同じ手順（edlddrds.htm の Processing 節）を独立に書いたもので、
# 組み直した軌道は edl_ddr.tab の位置と水平 0.2 km 以内で合う（半径は一定の 1.16 km。下記）。
#
# 使い方:
#   python3 convert_pds.py
#   python3 convert_pds.py --pds <edl_ddr.tab と r_sacc_s.tab のあるディレクトリ>

import argparse
import os
import sys

import numpy as np

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
DIRECTORY_PDS    = os.path.join(DIRECTORY_SCRIPT, '../../../references/20260929_Mars_EDL/pathfinder_pds')

# --- 機体（edlddrds.htm の Processing 節）
MASS = 585.3          # kg
AREA = 5.526          # m2

# --- 突入状態（edlddrds.htm の Processing 節。**火星固定（回転）系**）
#   1997-07-04 16:51:12.28 UTC、edl_ddr.tab の時刻 0 がここ（r = 3597.2 km を通過した時刻）
ENTRY_RADIUS    = 3597.2e3        # m
ENTRY_LATITUDE  = 23.0            # deg. N（地心）
ENTRY_LONGITUDE = 343.67          # deg. E
ENTRY_SPEED     = 7444.7          # m/s
ENTRY_PATH      = -16.85          # deg.（地平線より下）
ENTRY_AZIMUTH   = 255.41          # deg.（北から時計回り）

# --- 再構成が使った重力（同じ節）。J2 は正規化係数 8.759e-4 を非正規化して使う
GRAVITATIONAL_PARAMETER = 42828.3748574e9     # m3/s2（JGMRO_120F。config.yml と同じ）
RADIUS_GRAVITY          = 3396.0e3            # m
J2                      = 8.759e-4*np.sqrt(5.0)
ROTATION_RATE           = 7.0882181e-5        # rad/s（IAU 2015。config.yml と同じ）

# --- 時刻と高度の基準
#   edl_ddr.tab の最初の点は SCLK 1246726312.75 で時刻 34.156 s（edlddrds.htm）
SCLK_TIME_ZERO       = 1246726312.75 - 34.156
RADIUS_LANDING_SITE  = 3389.72                # km。edl_ddr.tab の高度の基準（Folkner et al. 1997）
RADIUS_PLANET        = 3396.0                 # km。Tacode の planet.radius（高度はこの球から測る）
GRAVITY_ACCELEROMETER = 9.795433              # m/s2。r_sacc_s.lbl の「1 g」

# --- 加速度計のレンジ切り替え
#   ゲインの列が切り替わった直後、生のカウントはまだ旧レンジの値から落ち着く途中で、
#   較正値が数十倍に跳ねる（82.28 s の 0.8 g -> 40 g で 36 g が 0.1 s 続く）。
#   落ち着くのに 0.3〜0.5 s かかるので 16 点（0.5 s）を捨てて前後から線形に埋める。
#   捨てずに積分すると偽の減速 35 m/s が速さに乗り、CD が 2〜4 % 大きく出る
NUMBER_MASK_SWITCH = 16

# --- 逆算した CD を採る範囲
#   加速度がこれより小さいと量子化とオフセットで CD が読めない（高度 77 km より上）
ACCELERATION_MINIMUM = 0.5        # m/s2
#   Spencer et al. (1998) 表 2 の迫撃砲の点火（突入界面 r = 3522.2 km から 169.6 s）。
#   その時刻は edl_ddr.tab の 38.202 s（下の TIME_ENTRY_INTERFACE）。点火の反動を避ける
TIME_ENTRY_INTERFACE = 50.482 - 12.28         # s。Spencer 表 1 のエポック 16:51:50.482 UTC
TIME_MORTAR          = TIME_ENTRY_INTERFACE + 169.6

# --- 音速（config*.yml の atmosphere 節と同じ値にすること）
#   比熱比は MSL の MEADS の再構成（Karlgaard et al. 2013、図 8）が使った CO2 大気の値、
#   平均分子量は edl_ddr.tab の 100 km 以下の値（Viking の質量分析計、Owen et al. 1977）
SPECIFIC_HEAT_RATIO = 1.335
MOLECULAR_WEIGHT    = 43.49        # g/mol
CONSTANT_GAS        = 8.314462618  # J/(mol K)
MASS_ATOMIC         = 1.66053906660e-27   # kg
CONSTANT_BOLTZMANN  = 1.380649e-23        # J/K

# --- Kn（config*.yml の satellite.characteristic_length と同じ値にすること）
LENGTH_REFERENCE = 1.325           # m。機体の底面半径（チュートリアルの 70 度球円錐と同じ）

# --- 大気の組成
#   100 km 以下は Viking 着陸機の測定（Owen et al. 1977, JGR 82, 4635、体積比）。
#   edl_ddr.tab の平均分子量はそれより上で下がる（光解離と拡散分離）。その差を
#   原子状酸素に割り当てる: x_O = (M_low - M)/(M_low - 16)。Kn にしか効かず、Kn が
#   効くのは CD を自由分子流側へ橋渡しする 77 km より上だけ
COMPOSITION_LOWER = [('CO2', 0.9532, 44.0095), ('N2', 0.027, 28.0134), ('Ar', 0.016, 39.948),
                     ('O2', 0.0013, 31.9988), ('CO', 0.0007, 28.0101)]
WEIGHT_OXYGEN     = 15.9994
# 剛体球の直径 m（src/atmosphere/atmosphere.py の DICT_DIAMETER_MOLECULAR と同じ値）
DIAMETER = {'CO2': 4.59e-10, 'N2': 3.75e-10, 'Ar': 3.64e-10, 'O2': 3.54e-10, 'CO': 3.76e-10, 'O': 3.04e-10}

# --- 空力表の格子
#   Mach の節点（飛行が通った 1.6〜44 の範囲）。各節点の値は、その ±3 % の点に
#   直線を当てはめた節点での値
MACH_NODE = [1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
             11.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 43.0]
WINDOW_MACH = 0.03
#   希薄化が CD を押し上げ始める Kn。これより小さい Kn の点で Mach 依存を読み、
#   これより大きい Kn の点で希薄化の増分を読む（飛行では CD は Kn 5e-3 で最小）
KNUDSEN_RAREFIED = 5.e-3
KNUDSEN_FLIGHT   = [1.e-2, 2.e-2, 3.e-2, 5.e-2]
WINDOW_KNUDSEN   = 0.1            # log の幅
#   自由分子流の CD と、そこへ着く Kn。値はチュートリアルの 70 度球円錐の表
#   （database/aerodynamic/generate_aerodynamic_table.py、完全な運動量交換）。
#   飛行からは読めない範囲なので、生成器と同じ sin^2 の橋渡しで飛行の最後の点からつなぐ
CD_FREE_MOLECULE      = 2.0
KNUDSEN_FREE_MOLECULE = 10.0
KNUDSEN_NODE_LOW      = [1.e-7, 1.e-4, KNUDSEN_RAREFIED]
KNUDSEN_NODE_BRIDGE   = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, KNUDSEN_FREE_MOLECULE, 1.e+2, 1.e+4]


def read_pds(directory):
  ddr = np.loadtxt(os.path.join(directory, 'edl_ddr.tab'), delimiter=',')
  acc = np.loadtxt(os.path.join(directory, 'r_sacc_s.tab'), delimiter=',')
  return ddr, acc


def clean_acceleration(acc):
  #
  # 軸方向（z 軸、機体の対称軸にほぼ沿う）の加速度 m/s2。レンジ切り替えの直後を捨てる
  #
  time  = acc[:,0] - SCLK_TIME_ZERO
  value = acc[:,4]*GRAVITY_ACCELEROMETER
  gain  = acc[:,10]
  switch = np.where(np.diff(gain) != 0.0)[0] + 1
  mask = np.zeros(len(time), dtype=bool)
  for index in switch:
    mask[index:index + NUMBER_MASK_SWITCH] = True
  value_clean = value.copy()
  value_clean[mask] = np.interp(time[mask], time[~mask], value[~mask])
  return time, value_clean, mask, time[switch]


def basis_local(latitude, longitude):
  east  = np.array([-np.sin(longitude), np.cos(longitude), 0.0])
  north = np.array([-np.sin(latitude)*np.cos(longitude), -np.sin(latitude)*np.sin(longitude), np.cos(latitude)])
  up    = np.array([ np.cos(latitude)*np.cos(longitude),  np.cos(latitude)*np.sin(longitude), np.sin(latitude)])
  return east, north, up


def gravity(position):
  radius = np.linalg.norm(position)
  sine2  = (position[2]/radius)**2
  factor = 1.5*J2*GRAVITATIONAL_PARAMETER*RADIUS_GRAVITY**2/radius**5
  return -GRAVITATIONAL_PARAMETER/radius**3*position \
         - factor*position*np.array([1.0 - 5.0*sine2, 1.0 - 5.0*sine2, 3.0 - 5.0*sine2])


def reconstruct(time_acc, acceleration, time_end, step=1.0/64.0):
  #
  # 突入状態から、重力（J2 まで）・コリオリ力・遠心力と、実測の軸方向加速度
  # （対気速度の逆向き。迎角はほぼ 0、Spencer et al. 1998 図 13）で積分する。
  # 大気は火星と共回転（風は無視。再構成も同じ仮定、edlddrds.htm）
  #
  east, north, up = basis_local(np.radians(ENTRY_LATITUDE), np.radians(ENTRY_LONGITUDE))
  path, azimuth = np.radians(ENTRY_PATH), np.radians(ENTRY_AZIMUTH)
  position = ENTRY_RADIUS*up
  velocity = ENTRY_SPEED*(np.cos(path)*(np.sin(azimuth)*east + np.cos(azimuth)*north) + np.sin(path)*up)
  omega = np.array([0.0, 0.0, ROTATION_RATE])

  def derivative(time, state):
    position, velocity = state[:3], state[3:]
    drag = -np.interp(time, time_acc, acceleration)*velocity/np.linalg.norm(velocity)
    return np.concatenate([velocity, gravity(position) + drag
                                     - 2.0*np.cross(omega, velocity) - np.cross(omega, np.cross(omega, position))])

  state = np.concatenate([position, velocity])
  time = 0.0
  history = []
  while time <= time_end:
    history.append([time, *state])
    k1 = derivative(time, state)
    k2 = derivative(time + 0.5*step, state + 0.5*step*k1)
    k3 = derivative(time + 0.5*step, state + 0.5*step*k2)
    k4 = derivative(time + step, state + step*k3)
    state = state + step/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
    time = time + step
  return np.array(history)


def composition(weight_mean):
  #
  # 平均分子量から数密度の割合 {種: x}。100 km 以下の組成に原子状酸素を混ぜて合わせる
  #
  weight_lower = sum(x*m for _, x, m in COMPOSITION_LOWER)/sum(x for _, x, _ in COMPOSITION_LOWER)
  fraction_oxygen = min(max((weight_lower - weight_mean)/(weight_lower - WEIGHT_OXYGEN), 0.0), 1.0)
  total_lower = sum(x for _, x, _ in COMPOSITION_LOWER)
  result = {name: (1.0 - fraction_oxygen)*x/total_lower for name, x, _ in COMPOSITION_LOWER}
  result['O'] = fraction_oxygen
  return result


def knudsen_number(density, weight_mean):
  #
  # src/atmosphere/atmosphere.py の set_knudsen_number と同じ式: Kn = 1/(sqrt(2) pi sum(n d^2) L)
  #
  number = density/(weight_mean*MASS_ATOMIC)
  total = np.zeros(len(density))
  for index in range(len(density)):
    total[index] = sum(number[index]*x*DIAMETER[name]**2 for name, x in composition(weight_mean[index]).items())
  return 1.0/(np.sqrt(2.0)*np.pi*total*LENGTH_REFERENCE)


def mach_number(speed, temperature):
  return speed/np.sqrt(SPECIFIC_HEAT_RATIO*CONSTANT_GAS/(MOLECULAR_WEIGHT*1.e-3)*temperature)


def fit_at(x, y, x0):
  # x0 のまわりの点に直線を当てはめ、x0 での値を返す
  if len(x) < 3:
    return None
  coefficient = np.polyfit(x - x0, y, 1)
  return coefficient[1]


def build_table(mach, knudsen, cd):
  #
  # CD(Mach, Kn) = CD_c(Mach) + dCD(Kn) の分離形で格子を作る。
  #
  # 飛行の中で、希薄化が効く区間（Kn > 5e-3、高度 62 km より上）は Mach 41〜48 の
  # 極超音速だけで、Mach が CD を動かす区間（Mach 38 以下）は Kn 1e-3 未満の連続流だけ。
  # 2 つの効果が別々の区間に現れるので、それぞれを別の区間から読める:
  #   CD_c(Mach) … Kn < 5e-3 の点
  #   dCD(Kn)    … Kn > 5e-3 の点で CD - CD_c(その点の Mach)。CD_c は Mach 43 より上を端の値で止める
  # 飛行が通っていない角（低い Mach で大きい Kn）は、この分離形を延ばしただけの値で、
  # データではない。
  #
  continuum = knudsen < KNUDSEN_RAREFIED
  node_mach, value_mach = [], []
  for node in MACH_NODE:
    select = continuum & (np.abs(np.log(mach/node)) < WINDOW_MACH)
    value = fit_at(mach[select], cd[select], node)
    if value is not None:
      node_mach.append(node)
      value_mach.append(value)
  node_mach, value_mach = np.array(node_mach), np.array(value_mach)

  def cd_continuum(value):
    return np.interp(value, node_mach, value_mach)

  increment = cd - cd_continuum(mach)
  node_knudsen = list(KNUDSEN_NODE_LOW)
  value_knudsen = [0.0]*len(KNUDSEN_NODE_LOW)
  for node in KNUDSEN_FLIGHT:
    select = (~continuum) & (np.abs(np.log(knudsen/node)) < WINDOW_KNUDSEN)
    if select.sum() < 3:
      continue
    node_knudsen.append(node)
    value_knudsen.append(float(np.mean(increment[select])))

  # 飛行の最後の Kn の点から自由分子流へ、生成器と同じ sin^2 で橋渡しする
  knudsen_last, increment_last = node_knudsen[-1], value_knudsen[-1]
  increment_free = CD_FREE_MOLECULE - value_mach[-1]
  for node in KNUDSEN_NODE_BRIDGE:
    if node >= KNUDSEN_FREE_MOLECULE:
      fraction = 1.0
    else:
      fraction = np.sin(0.5*np.pi*np.log(node/knudsen_last)/np.log(KNUDSEN_FREE_MOLECULE/knudsen_last))**2
    node_knudsen.append(node)
    value_knudsen.append(increment_last + (increment_free - increment_last)*fraction)

  return node_mach, value_mach, np.array(node_knudsen), np.array(value_knudsen)


def write_aerodynamic_table(path, node_mach, value_mach, node_knudsen, value_knudsen):
  with open(path, 'w') as f:
    f.write('Mars Pathfinder drag coefficient traced back from the flight (validation/pathfinder/convert_pds.py)\n')
    f.write('# Reference length: {:g} m\n'.format(LENGTH_REFERENCE))
    f.write('# CD = 2 m a/(rho A V^2) with a the measured axial acceleration (PDS r_sacc_s.tab), rho the\n')
    f.write('# reconstructed density (PDS edl_ddr.tab) and V the speed integrated from the entry state.\n')
    f.write('# This is the drag coefficient the NASA reconstruction used (Braun et al. 1995 database),\n')
    f.write('# along the path the vehicle flew. Separable form CD(Mach, Kn) = CD_c(Mach) + dCD(Kn):\n')
    f.write('# CD_c from the continuum part (Kn < {:g}), dCD from the rarefied part up to Kn {:g},\n'.format(
            KNUDSEN_RAREFIED, KNUDSEN_FLIGHT[-1]))
    f.write('# then a sin^2 bridge to the free-molecular {:g} at Kn {:g} (no flight data there).\n'.format(
            CD_FREE_MOLECULE, KNUDSEN_FREE_MOLECULE))
    f.write('# Mach = V/sqrt(gamma R T/M) with gamma {:g}, M {:g} g/mol.\n'.format(SPECIFIC_HEAT_RATIO, MOLECULAR_WEIGHT))
    f.write('variables = Kn, CFx, CFy, CFz, CMx, CMy, CMz, SDV_CFx, SDV_CFy, SDV_CFz, SDV_CMx, SDV_CMy, SDV_CMz, Altitude\n')
    for mach, cd_continuum in zip(node_mach, value_mach):
      f.write('MACH {:g}\n'.format(mach))
      for knudsen, increment in zip(node_knudsen, value_knudsen):
        f.write('{:.6e} {:.6f} 0 0 0 0 0 0 0 0 0 0 0 0\n'.format(knudsen, cd_continuum + increment))


def write_atmosphere_table(path, altitude, density, temperature, weight_mean, offset_radius):
  #
  # 飛行で測った大気を Tacode の CCMC 形式で書く。高度は Tacode の planet.radius の球から。
  # 数密度は rho/(M m_u) を composition() の割合で分けたもの（Kn にだけ効く）
  #
  order = np.argsort(altitude)
  altitude, density, temperature, weight_mean = altitude[order], density[order], temperature[order], weight_mean[order]
  keep = np.concatenate([[True], np.diff(altitude) > 0.0])
  names = ['CO2', 'N2', 'Ar', 'CO', 'O', 'O2']
  with open(path, 'w') as f:
    f.write('Mars Pathfinder entry 1997-07-04: the atmosphere measured in flight\n')
    f.write('Source: PDS MPFL-M-ASIMET-4-DDR-EDL-V1.0, edl_ddr.tab (Magalhaes, Schofield and Seiff 1999)\n')
    f.write('Written by validation/pathfinder/convert_pds.py\n')
    f.write('Height is above the sphere of radius {:.1f} km: the radius of the landing site {:.2f} km plus the\n'.format(
            RADIUS_PLANET, RADIUS_LANDING_SITE))
    f.write('altitude of edl_ddr.tab, less {:.3f} km (edl_ddr.tab puts the vehicle that much higher than the\n'.format(offset_radius))
    f.write('documented entry state and the measured acceleration do, the same at every point)\n')
    f.write('Number densities: rho/(M m_u) split by the Viking lander composition below 100 km (Owen et al. 1977),\n')
    f.write('with atomic oxygen making up the drop of the mean molecular weight M above it\n')
    f.write('\n')
    f.write('Selected parameters are:\n')
    f.write('1 Height, km\n')
    for index, name in enumerate(names):
      f.write('{} {}, cm-3\n'.format(index + 2, name))
    f.write('{} Mass_density, g/cm-3\n'.format(len(names) + 2))
    f.write('{} Temperature_neutral, K\n'.format(len(names) + 3))
    f.write('\n')
    f.write(''.join('{:>11d}'.format(index + 1) for index in range(len(names) + 3)) + '\n')
    for index in np.where(keep)[0]:
      number = density[index]/(weight_mean[index]*MASS_ATOMIC)*1.e-6
      fraction = composition(weight_mean[index])
      value = [altitude[index]] + [number*fraction[name] for name in names] + [density[index]*1.e-3, temperature[index]]
      f.write(''.join('{:11.4e}'.format(v) if i > 0 else '{:11.4f}'.format(v) for i, v in enumerate(value)) + '\n')


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--pds', type=str, default=DIRECTORY_PDS)
  args = parser.parse_args()

  ddr, acc = read_pds(args.pds)
  time_ddr = ddr[:,1]
  time_acc, acceleration, mask, time_switch = clean_acceleration(acc)
  print('Gain switches at t =', ', '.join('{:.3f}'.format(t) for t in time_switch), 's;',
        mask.sum(), 'samples replaced')

  # 軌道を組み直し、edl_ddr.tab の位置と比べる
  history = reconstruct(time_acc, acceleration, time_ddr[-1] + 0.1)
  radius_rec = np.linalg.norm(history[:,1:4], axis=1)
  speed_rec  = np.linalg.norm(history[:,4:7], axis=1)
  radius_ddr = (RADIUS_LANDING_SITE + ddr[:,2])*1.e3
  offset = radius_ddr - np.interp(time_ddr, history[:,0], radius_rec)
  offset_radius = float(np.mean(offset))*1.e-3
  print('edl_ddr.tab radius minus the reconstruction: mean {:.3f} km, spread {:.3f} km (min {:.3f}, max {:.3f})'.format(
        offset_radius, np.std(offset)*1.e-3, offset.min()*1.e-3, offset.max()*1.e-3))
  latitude_rec  = np.degrees(np.arcsin(history[:,3]/radius_rec))
  longitude_rec = np.degrees(np.arctan2(history[:,2], history[:,1])) % 360.0
  north = np.radians(np.interp(time_ddr, history[:,0], latitude_rec) - ddr[:,3])*radius_ddr
  east  = np.radians(np.interp(time_ddr, history[:,0], longitude_rec) - ddr[:,4])*radius_ddr*np.cos(np.radians(ddr[:,3]))
  print('edl_ddr.tab position minus the reconstruction: north {:.1f} km, east {:.1f} km (spread {:.2f}, {:.2f} km)'.format(
        -np.mean(north)*1.e-3, -np.mean(east)*1.e-3, np.std(north)*1.e-3, np.std(east)*1.e-3))

  # 逆算
  speed = np.interp(time_ddr, history[:,0], speed_rec)
  accel = np.interp(time_ddr, time_acc, acceleration)
  density, temperature, weight = ddr[:,5], ddr[:,7], ddr[:,8]
  cd = 2.0*MASS*accel/(density*AREA*speed**2)
  mach = mach_number(speed, temperature)
  knudsen = knudsen_number(density, weight)
  valid = (accel > ACCELERATION_MINIMUM) & (time_ddr < TIME_MORTAR)
  print('CD traced back over t = {:.2f} .. {:.2f} s, Mach {:.2f} .. {:.2f}, Kn {:.2e} .. {:.2e}'.format(
        time_ddr[valid][0], time_ddr[valid][-1], mach[valid].min(), mach[valid].max(),
        knudsen[valid].min(), knudsen[valid].max()))

  node_mach, value_mach, node_knudsen, value_knudsen = build_table(mach[valid], knudsen[valid], cd[valid])
  print('CD_c(Mach):', ', '.join('{:g}:{:.4f}'.format(m, v) for m, v in zip(node_mach, value_mach)))
  print('dCD(Kn)   :', ', '.join('{:.0e}:{:+.4f}'.format(k, v) for k, v in zip(node_knudsen, value_knudsen)))

  # 表がたどった経路の CD をどれだけ再現するか
  import_path = os.path.join(DIRECTORY_SCRIPT, '../../src')
  sys.path.append(import_path)
  import satellite.satellite as satellite   # noqa: E402
  table = {satellite.KEY_MACH: node_mach, satellite.KEY_KN: node_knudsen,
           satellite.KEY_CD_MACH: value_mach[:,None] + value_knudsen[None,:]}
  cd_table = np.array([satellite.get_aerodynamic_coefficient_mach(k, m, table)
                       for k, m in zip(knudsen[valid], mach[valid])])
  residual = cd_table/cd[valid] - 1.0
  print('Table against the traced-back CD: rms {:.2f} %, max {:.2f} %'.format(
        100.0*np.sqrt(np.mean(residual**2)), 100.0*np.max(np.abs(residual))))

  # 書き出し
  altitude_sphere = RADIUS_LANDING_SITE + ddr[:,2] - offset_radius - RADIUS_PLANET
  os.makedirs(os.path.join(DIRECTORY_SCRIPT, 'reference'), exist_ok=True)
  os.makedirs(os.path.join(DIRECTORY_SCRIPT, 'database/atmosphere'), exist_ok=True)
  os.makedirs(os.path.join(DIRECTORY_SCRIPT, 'database/aerodynamic'), exist_ok=True)

  with open(os.path.join(DIRECTORY_SCRIPT, 'reference/pathfinder_ddr.dat'), 'w') as f:
    f.write('# Mars Pathfinder, PDS MPFL-M-ASIMET-4-DDR-EDL-V1.0 edl_ddr.tab in SI\n')
    f.write('# time: s from r = 3597.2 km (1997-07-04 16:51:12.28 UTC); radius: 3389.72 km + the altitude of the file\n')
    f.write('# altitude_tacode: radius - {:.3f} km - 3396.0 km, i.e. above the sphere of config.yml in the frame\n'.format(offset_radius))
    f.write('#   of the entry state (see convert_pds.py)\n')
    f.write('# time[s] radius[km] altitude_tacode[km] latitude[deg] longitude[deg] density[kg/m3] density_plus[kg/m3]'
            ' density_minus[kg/m3] pressure[Pa] temperature[K] temperature_plus[K] temperature_minus[K] molecular_weight\n')
    for index in range(len(time_ddr)):
      row = ddr[index]
      f.write('{:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.4e} {:.4e} {:.4e} {:.4e} {:.3f} {:.4f} {:.4f} {:.3f}\n'.format(
              row[1], RADIUS_LANDING_SITE + row[2], altitude_sphere[index], row[3], row[4],
              row[5], row[9], row[10], row[6]*100.0, row[7], row[13], row[14], row[8]))

  with open(os.path.join(DIRECTORY_SCRIPT, 'reference/pathfinder_acceleration.dat'), 'w') as f:
    f.write('# Mars Pathfinder, PDS MPFL-M-ASIMET-2/3 r_sacc_s.tab: science z-axis accelerometer, 32 Hz, in m/s2\n')
    f.write('# (1 g = {} m/s2 of the label). Time as in pathfinder_ddr.dat. replaced = 1 for the {} samples after\n'.format(
            GRAVITY_ACCELEROMETER, NUMBER_MASK_SWITCH))
    f.write('# each range switch, filled by linear interpolation\n')
    f.write('# time[s] acceleration_axial[m/s2] replaced\n')
    select = (time_acc > -1.0) & (time_acc < time_ddr[-1] + 1.0)
    for t, a, m in zip(time_acc[select], acceleration[select], mask[select]):
      f.write('{:.5f} {:.6e} {:d}\n'.format(t, a, int(m)))

  with open(os.path.join(DIRECTORY_SCRIPT, 'reference/pathfinder_reconstruction.dat'), 'w') as f:
    f.write('# Mars Pathfinder: speed integrated from the documented entry state with the measured acceleration,\n')
    f.write('# and the drag coefficient traced back from it, CD = 2 m a/(rho A V^2) (convert_pds.py)\n')
    f.write('# valid = 1 where the acceleration is above {:g} m/s2 and before the mortar fire ({:.1f} s)\n'.format(
            ACCELERATION_MINIMUM, TIME_MORTAR))
    f.write('# time[s] speed[m/s] acceleration[m/s2] CD Mach Kn valid\n')
    for index in range(len(time_ddr)):
      f.write('{:.3f} {:.3f} {:.5e} {:.5f} {:.4f} {:.4e} {:d}\n'.format(
              time_ddr[index], speed[index], accel[index], cd[index], mach[index], knudsen[index], int(valid[index])))

  write_atmosphere_table(os.path.join(DIRECTORY_SCRIPT, 'database/atmosphere/atmospheremodel_pathfinder_ddr.txt'),
                         altitude_sphere, density, temperature, weight, offset_radius)
  write_aerodynamic_table(os.path.join(DIRECTORY_SCRIPT, 'database/aerodynamic/aerodynamic_pathfinder.txt'),
                          node_mach, value_mach, node_knudsen, value_knudsen)
  print('Written: reference/*.dat, database/atmosphere/atmospheremodel_pathfinder_ddr.txt,'
        ' database/aerodynamic/aerodynamic_pathfinder.txt')


if __name__ == '__main__':
  main()
