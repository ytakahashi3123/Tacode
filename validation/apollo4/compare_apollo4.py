#!/usr/bin/env python3
#
# Apollo 4（AS-501）の飛行データと Tacode を突き合わせる。
#
# 3 つを別々に見る（混ぜると、どこがどれだけ効いているか分からなくなる）:
#
#   (1) 軌道      Tacode を突入状態から自由飛行させ、高度と対地速度の履歴を飛行値と比べる。
#                 抗力のみ（config.yml）、揚力つき（config_lift.yml）、6 自由度
#                 （config_6dof.yml、重心オフセットでトリムさせる）の 3 本を並べる。
#                 **抗力だけではスキップアウトが出ない**ので、比較になるのは最初の
#                 ディップまで。揚力を入れると全域が比較になる（ただしバンク角は仮定）
#   (2) 大気       飛行値の高度をそのまま入力に、Tacode の大気テーブル（NRLMSISE-00）の
#                 密度を飛行値（よどみ圧から導かれた自由流密度）と比べる。
#                 **軌道計算を通さない**ので、揚力の有無に影響されない
#   (3) 加熱       よどみ点加熱の相関式（Detra-Kemp-Riddell と Sutton-Graves）を
#                 「飛行値の密度・速度」と「テーブルの密度・飛行値の速度」で評価し、
#                 飛行値と比べる。相関式は Tacode の一部ではないので、これは
#                 「Tacode の密度・速度を加熱評価に使えるか」を見る検算
#
# 加熱の比較で気をつけること:
#   - 飛行値の加熱は **S/R = 0.732**（放射計 CA3363K の位置、迎角 24.4 度でのよどみ域）の
#     冷壁値で、よどみ点値ではない
#   - 参照 [1] は迎角 24.4 度の Apollo 司令船の正面衝撃波条件が **半径 10 ft
#     （3.048 m）の球**に相当するとしている。機体の伝熱面の曲率半径 4.694 m とは
#     別物で、加熱率は 1/sqrt(Rn) で効く。既定は 4.694 m、--nose-radius で変えられる
#
# 使い方:
#   cd validation/apollo4
#   ./run_tacode.sh                               # 抗力のみ
#   ./run_tacode.sh -file config_lift.yml         # 揚力つき
#   ./run_tacode.sh -file config_6dof.yml         # 6 自由度（重心オフセットでトリム）
#   python3 compare_apollo4.py                    # 在るものを全部並べて図と要約
#   python3 compare_apollo4.py --no-figure        # 数字だけ
#   python3 compare_apollo4.py --nose-radius 3.048
#   python3 compare_apollo4.py --case output_result/tecplot.dat "drag only"

import argparse
import os
import sys

import numpy as np


# Tacode 本体と後処理の部品を使う（規約を二重に持たない）
DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
DIRECTORY_SRC    = os.path.join(DIRECTORY_SCRIPT, '../../src')
DIRECTORY_HELPER = os.path.join(DIRECTORY_SCRIPT, '../../src_helper')
sys.path.append(DIRECTORY_SRC)
sys.path.append(os.path.join(DIRECTORY_HELPER, 'animate_trajectory'))

import tecplot_reader as tecplot_reader   # noqa: E402
from general.general import general       # noqa: E402
import atmosphere.atmosphere as atmosphere_module   # noqa: E402


# 相関式の定数
DENSITY_SEA_LEVEL   = 1.225        # kg/m3
VELOCITY_CIRCULAR   = 7924.8       # m/s, Detra-Kemp-Riddell の基準速度
COEFFICIENT_DKR     = 1.1035e8     # W/m2, Rn [m], rho/rho_sl, V/V_c
COEFFICIENT_SUTTON  = 1.7415e-4    # W/m2, rho [kg/m3], Rn [m], V [m/s]

# 飛行値の時刻（打ち上げからの秒）と計算の経過時間の対応
TIME_ENTRY_DEFAULT  = 29968.0

# 既定で並べるケース（在るものだけを使う）
CASE_DEFAULT = (
  ('output_result/tecplot.dat',      '3-DOF, drag only'),
  ('output_result_lift/tecplot.dat', '3-DOF, measured vertical lift'),
  ('output_result_6dof/tecplot.dat', '6-DOF, trimmed capsule (constant bank 55 deg.)'),
)

# 図の色（ケースの並び順）
COLOUR_CASE = ('tab:red', 'tab:purple', 'tab:orange', 'tab:brown')


def heating_dkr(density, velocity, nose_radius):
  # Detra-Kemp-Riddell のよどみ点対流加熱, W/m2
  return COEFFICIENT_DKR/np.sqrt(nose_radius) \
         * np.sqrt(density/DENSITY_SEA_LEVEL) * (velocity/VELOCITY_CIRCULAR)**3.15


def heating_sutton_graves(density, velocity, nose_radius):
  # Sutton-Graves のよどみ点対流加熱, W/m2
  return COEFFICIENT_SUTTON*np.sqrt(density/nose_radius)*velocity**3


def statistics_relative(value_model, value_reference, mask_valid=None):
  # 相対差の統計。参照値が 0 の点（表に値が無い行）は落とす
  mask = value_reference > 0.0
  if mask_valid is not None :
    mask = np.logical_and(mask, mask_valid)
  if not np.any(mask) :
    return None
  ratio = value_model[mask]/value_reference[mask]
  return {
    'number': int(np.sum(mask)),
    'mean':   float(np.mean(ratio)),
    'min':    float(np.min(ratio)),
    'max':    float(np.max(ratio)),
    'rms':    float(np.sqrt(np.mean((ratio - 1.0)**2))),
  }


def print_statistics(title, statistics):
  if statistics is None :
    print('{:<34s} (no point to compare)'.format(title))
    return
  print('{:<34s} mean {:6.3f}, min {:6.3f}, max {:6.3f}, rms of (ratio-1) {:6.3f}  [{:d} points]'.format(
        title, statistics['mean'], statistics['min'], statistics['max'],
        statistics['rms'], statistics['number']))


def compare_trajectory(data_case, data_flight, time_entry, label):
  # 飛行値の時刻に計算値を内挿する。計算が終わった後（地表に達した後）は比べない
  time_flight   = data_flight['Time'] - time_entry
  time_case     = data_case['Time']
  mask          = time_flight <= time_case[-1]

  altitude_case = np.interp(time_flight[mask], time_case, data_case['Alti'])
  velocity_case = np.interp(time_flight[mask], time_case, data_case['VelplAbs'])

  difference_altitude = altitude_case - data_flight['Alti'][mask]
  difference_velocity = velocity_case - data_flight['VelrelAbs'][mask]

  # 突入直後（最大加熱までの区間）はどうか
  time_peak = data_flight['Time'][np.argmax(data_flight['Qconv'])] - time_entry
  window    = time_flight[mask] <= time_peak

  print('')
  print('--' + label)
  print('  points compared            : {:d} of {:d} (the run ends at {:.1f} s)'.format(
        int(np.sum(mask)), len(time_flight), time_case[-1]))
  print('  altitude difference        : max {:+.1f} km, rms {:.1f} km'.format(
        difference_altitude[np.argmax(np.abs(difference_altitude))], np.sqrt(np.mean(difference_altitude**2))))
  print('  velocity difference        : max {:+.0f} m/s, rms {:.0f} m/s'.format(
        difference_velocity[np.argmax(np.abs(difference_velocity))], np.sqrt(np.mean(difference_velocity**2))))
  print('  up to the peak heating at {:.0f} s:'.format(time_peak))
  print('    altitude difference      : max {:+.1f} km, rms {:.1f} km'.format(
        difference_altitude[window][np.argmax(np.abs(difference_altitude[window]))],
        np.sqrt(np.mean(difference_altitude[window]**2))))
  print('    velocity difference      : max {:+.0f} m/s, rms {:.0f} m/s'.format(
        difference_velocity[window][np.argmax(np.abs(difference_velocity[window]))],
        np.sqrt(np.mean(difference_velocity[window]**2))))

  # 最初のディップ（スキップアウトの有無がここに出る）
  index_dip = int(np.argmin(data_case['Alti']))
  print('  lowest point of the run    : {:.2f} km at {:.1f} s'.format(
        data_case['Alti'][index_dip], time_case[index_dip]))

  return


def density_from_table(config, altitude_flight):
  # Tacode の大気テーブルを Tacode 自身の読み取り・内挿で引く
  # （ここで別の実装を書くと、検証しているものが変わってしまう）
  atmosphere_dict = atmosphere_module.initial_settings_atmosphere(config)
  density_list     = []
  temperature_list = []
  for altitude in altitude_flight:
    # 高度は km で渡す（テーブルの Height がそのまま km、solver も m2km を掛けて渡す）
    density, temperature, knudsen = atmosphere_module.get_atmosphere_property(altitude, atmosphere_dict)
    density_list.append(density)
    temperature_list.append(temperature)
  return np.array(density_list), np.array(temperature_list)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--case-directory', type=str, default=DIRECTORY_SCRIPT,
                      help='Directory of the case (config.yml and output_result are read from here)')
  parser.add_argument('--file-config', type=str, default='config.yml')
  parser.add_argument('--case', action='append', nargs=2, metavar=('PATH', 'LABEL'), default=None,
                      help='A Tecplot output and the label to draw it with; may be given more than once')
  parser.add_argument('--file-reference', type=str, default='reference/apollo4_flight.dat')
  parser.add_argument('--time-entry', type=float, default=TIME_ENTRY_DEFAULT,
                      help='Flight time [s from lift-off] of the initial state of the run')
  parser.add_argument('--nose-radius', type=float, default=4.694,
                      help='Nose radius [m] for the heating correlations')
  parser.add_argument('--output-directory', type=str, default='output_comparison')
  parser.add_argument('--no-figure', action='store_true',
                      help='Write the numbers only (matplotlib is then not needed)')
  args = parser.parse_args()

  # config の database は "manual"（カレントディレクトリ基準）なので、ケースへ移る
  os.chdir(args.case_directory)

  for filename in [args.file_config, args.file_reference]:
    if not os.path.exists(filename) :
      print('File not found:', filename)
      print('Program stopped.')
      sys.exit(1)

  # 並べるケース。指定が無ければ在るものを既定の一覧から採る
  case_list = args.case if args.case is not None else \
              [ [path, label] for path, label in CASE_DEFAULT if os.path.exists(path) ]
  if len(case_list) == 0 :
    print('No result to compare with.')
    print('--Run the case first: ./run_tacode.sh (and ./run_tacode.sh -file config_lift.yml)')
    print('Program stopped.')
    sys.exit(1)
  for path, label in case_list:
    if not os.path.exists(path) :
      print('File not found:', path)
      print('Program stopped.')
      sys.exit(1)

  config = general().read_config_yaml(args.file_config)

  data_case_list = [ (tecplot_reader.read_tecplot(path), label) for path, label in case_list ]
  data_flight    = tecplot_reader.read_tecplot(args.file_reference)

  print('Apollo 4 (AS-501) entry: flight data against Tacode')
  for path, label in case_list:
    print('--Case      :', os.path.join(args.case_directory, path), '(' + label + ')')
  print('--Reference :', os.path.join(args.case_directory, args.file_reference))
  print('--Nose radius for the correlations: {:.3f} m'.format(args.nose_radius))

  print('')
  print('(1) Trajectory against the flight data')
  for data_case, label in data_case_list:
    compare_trajectory(data_case, data_flight, args.time_entry, label)

  # (2) 大気: 飛行値の高度でテーブルを引く
  altitude_flight = data_flight['Alti']
  density_table, temperature_table = density_from_table(config, altitude_flight)
  density_flight = data_flight['Dens']

  print('')
  print('(2) Atmosphere: the table of the case at the altitudes of the flight data')
  print_statistics('--Density, table / flight   :', statistics_relative(density_table, density_flight))
  print('--The flight density is derived from the measured stagnation pressure; the three')
  print('  points before the first pressure reading carry no value and are left out.')

  # (3) 加熱: 相関式を飛行値の密度とテーブルの密度で評価する
  velocity_flight = data_flight['VelrelAbs']
  heating_flight  = data_flight['Qconv']

  # 表の 3 行（突入直後）にはよどみ圧の計測が無く、密度が 0 で入っている
  mask_density_flight = density_flight > 0.0

  heating_dkr_flight    = heating_dkr(density_flight, velocity_flight, args.nose_radius)
  heating_dkr_table     = heating_dkr(density_table,  velocity_flight, args.nose_radius)
  heating_sutton_flight = heating_sutton_graves(density_flight, velocity_flight, args.nose_radius)

  print('')
  print('(3) Stagnation-point convective heating correlations against the flight value')
  print('    (the flight value is the cold-wall heating at S/R = 0.732, not at the stagnation point)')
  print_statistics('--DKR with the flight density:',
                   statistics_relative(heating_dkr_flight, heating_flight, mask_density_flight))
  print_statistics('--DKR with the table density :',
                   statistics_relative(heating_dkr_table, heating_flight, mask_density_flight))
  print_statistics('--Sutton-Graves, flight dens.:',
                   statistics_relative(heating_sutton_flight, heating_flight, mask_density_flight))

  index_peak = int(np.argmax(heating_flight))
  print('--At the peak of the flight heating ({:.0f} s, {:.1f} km, {:.0f} m/s):'.format(
        data_flight['Time'][index_peak], altitude_flight[index_peak], velocity_flight[index_peak]))
  print('  flight {:.3e} W/m2, DKR {:.3e} W/m2, Sutton-Graves {:.3e} W/m2'.format(
        heating_flight[index_peak], heating_dkr_flight[index_peak], heating_sutton_flight[index_peak]))

  if not args.no_figure :
    make_figure(args, data_case_list, data_flight,
                density_table, heating_dkr_flight, heating_dkr_table, heating_sutton_flight)

  return


def make_figure(args, data_case_list, data_flight,
                density_table, heating_dkr_flight, heating_dkr_table, heating_sutton_flight):
  try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
  except ImportError:
    print('')
    print('Caution: matplotlib is not installed, so no figure is written.')
    return

  os.makedirs(args.output_directory, exist_ok=True)

  colour_flight = 'black'
  colour_case   = COLOUR_CASE[0]

  def finalise(axis, filename, legend_location='best'):
    axis.grid(linewidth=0.5, alpha=0.7)
    axis.legend(loc=legend_location, frameon=False)
    plt.tight_layout()
    path = os.path.join(args.output_directory, filename)
    plt.savefig(path, dpi=200)
    plt.close()
    print('--Figure:', path)

  print('')
  print('Figures')

  # 高度履歴
  figure, axis = plt.subplots(figsize=(6, 4))
  axis.plot(data_flight['Time'], data_flight['Alti'], color=colour_flight, marker='o',
            markerfacecolor='white', markersize=4, linewidth=1.5, label='Apollo 4 flight')
  for index, (data_case, label) in enumerate(data_case_list):
    axis.plot(data_case['Time'] + args.time_entry, data_case['Alti'],
              color=COLOUR_CASE[index % len(COLOUR_CASE)], linewidth=2.0, label='Tacode, ' + label)
  axis.set_xlabel('Time from lift-off, s')
  axis.set_ylabel('Altitude, km')
  axis.set_xlim(data_flight['Time'][0], data_flight['Time'][-1])
  axis.set_ylim(0.0, 130.0)
  finalise(axis, 'profile_time-altitude.png')

  # 速度履歴
  figure, axis = plt.subplots(figsize=(6, 4))
  axis.plot(data_flight['Time'], data_flight['VelrelAbs'], color=colour_flight, marker='o',
            markerfacecolor='white', markersize=4, linewidth=1.5, label='Apollo 4 flight')
  for index, (data_case, label) in enumerate(data_case_list):
    axis.plot(data_case['Time'] + args.time_entry, data_case['VelplAbs'],
              color=COLOUR_CASE[index % len(COLOUR_CASE)], linewidth=2.0, label='Tacode, ' + label)
  axis.set_xlabel('Time from lift-off, s')
  axis.set_ylabel('Relative velocity, m/s')
  axis.set_xlim(data_flight['Time'][0], data_flight['Time'][-1])
  finalise(axis, 'profile_time-velocity.png')

  # 速度-高度
  figure, axis = plt.subplots(figsize=(6, 4))
  axis.plot(data_flight['VelrelAbs'], data_flight['Alti'], color=colour_flight, marker='o',
            markerfacecolor='white', markersize=4, linewidth=1.5, label='Apollo 4 flight')
  for index, (data_case, label) in enumerate(data_case_list):
    axis.plot(data_case['VelplAbs'], data_case['Alti'],
              color=COLOUR_CASE[index % len(COLOUR_CASE)], linewidth=2.0, label='Tacode, ' + label)
  axis.set_xlabel('Relative velocity, m/s')
  axis.set_ylabel('Altitude, km')
  axis.set_ylim(0.0, 130.0)
  finalise(axis, 'profile_velocity-altitude.png')

  # 密度（飛行値 vs テーブル）
  mask_density = data_flight['Dens'] > 0.0
  figure, axis = plt.subplots(figsize=(6, 4))
  axis.semilogx(data_flight['Dens'][mask_density], data_flight['Alti'][mask_density],
                color=colour_flight, marker='o', markerfacecolor='white', markersize=4,
                linewidth=1.5, label='Apollo 4 flight (from the stagnation pressure)')
  axis.semilogx(density_table[mask_density], data_flight['Alti'][mask_density],
                color='tab:blue', linewidth=2.0, label='NRLMSISE-00 table of the case')
  axis.set_xlabel(r'Free-stream density, kg/m$^3$')
  axis.set_ylabel('Altitude, km')
  finalise(axis, 'profile_density-altitude.png')

  # 対流加熱
  figure, axis = plt.subplots(figsize=(6, 4))
  axis.plot(data_flight['Time'], data_flight['Qconv'], color=colour_flight, marker='o',
            markerfacecolor='white', markersize=4, linewidth=1.5,
            label='Apollo 4 flight, S/R = 0.732')
  axis.plot(data_flight['Time'], heating_dkr_flight, color='tab:blue', linestyle='--',
            linewidth=2.0, label='DKR, flight density')
  axis.plot(data_flight['Time'], heating_dkr_table, color=colour_case, linestyle='-.',
            linewidth=2.0, label='DKR, table density')
  axis.plot(data_flight['Time'], heating_sutton_flight, color='tab:green', linestyle=':',
            linewidth=2.0, label='Sutton-Graves, flight density')
  axis.set_xlabel('Time from lift-off, s')
  axis.set_ylabel(r'Heat flux, W/m$^2$')
  axis.set_xlim(data_flight['Time'][0], data_flight['Time'][-1])
  finalise(axis, 'profile_time-heating.png')

  # 放射加熱（飛行計測のみ。相関式は持っていない）
  mask_radiation = data_flight['Qrad'] > 0.0
  if np.any(mask_radiation) :
    figure, axis = plt.subplots(figsize=(6, 4))
    axis.plot(data_flight['Time'][mask_radiation], data_flight['Qrad'][mask_radiation],
              color=colour_flight, marker='o', markerfacecolor='white', markersize=4,
              linewidth=1.5, label='Apollo 4 radiometer CA3363K')
    axis.plot(data_flight['Time'], data_flight['Qconv'], color='tab:blue', linestyle='--',
              linewidth=2.0, label='Convective, S/R = 0.732')
    axis.set_xlabel('Time from lift-off, s')
    axis.set_ylabel(r'Heat flux, W/m$^2$')
    axis.set_xlim(data_flight['Time'][0], data_flight['Time'][-1])
    finalise(axis, 'profile_time-radiation.png')

  return


if __name__ == '__main__':
  main()
