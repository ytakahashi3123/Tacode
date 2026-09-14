#!/usr/bin/env python3
#
# Project Fire flight II と Tacode を突き合わせる。
#
# **比べている相手が 2 種類あることに注意。**
#
#   reference/fire2_trajectory.dat  NASA TN D-3569 表 V。レーダで拘束して実測大気を
#                                   飛ばした **NASA 自身の 3 自由度計算**であって
#                                   実測ではない。これと合うことは「同じ入力から
#                                   同じ軌道が出る」ことの確認（ソルバーの検証）
#   database/atmosphere/...         Tacode が飛ぶ NRLMSISE-00 のテーブル。表 V の
#                                   密度（アセンション島のゾンデ実測）と比べる
#
# 機体は**完全に無制御**で軸対称なので、揚力の向きもバンク角も仮定していない。
# 自由度として残っているのは一つも無い（質量も直径も報告書の実測値）。
#
# 使い方:
#   cd validation/fire2 && ./run_tacode.sh
#   python3 compare_fire2.py
#   python3 compare_fire2.py --no-figure

import argparse
import os
import sys

import numpy as np

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src'))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src_helper/animate_trajectory'))

import tecplot_reader as tecplot_reader      # noqa: E402
from general.general import general          # noqa: E402

CASE_DEFAULT = (('output_result/tecplot.dat', '3-DOF, drag only'),)

COLOUR_CASE = ('tab:red', 'tab:purple', 'tab:orange', 'tab:brown')

FILE_TRAJECTORY = 'reference/fire2_trajectory.dat'
DIRECTORY_OUTPUT = 'output_comparison'

COLUMN_TRAJECTORY = ['time', 'latitude', 'longitude', 'altitude', 'velocity',
                     'flightpath', 'heading', 'dynamic_pressure', 'pressure',
                     'density', 'temperature', 'mach', 'acceleration', 'reynolds']

# 表 V の時刻は**飛行経過時刻**、Tacode の時刻は計算開始からの経過時間。
# config.yml の開始時刻を足して揃える
TIME_START_DEFAULT = 1618.25

# 全区間のほかに、この区間だけを別に出す。**極超音速のうちだけ**を切り出すため:
# t = 1700 s より後は Mach 3 を割っており、Tacode の係数表にも表 V の一定の
# 弾道係数にも Mach 依存が無いので、そこは空力の欠落を見ていることになる
WINDOW_HYPERSONIC = (1618.25, 1700.0)


def read_trajectory(file_trajectory):
  # 表 V。読めなかった値は NaN なので、量ごとに生きている点だけを使う
  data = {name: [] for name in COLUMN_TRAJECTORY}
  with open(file_trajectory) as stream:
    for line in stream:
      if line.startswith('#') :
        continue
      value = [float(item) for item in line.split()]
      for name, item in zip(COLUMN_TRAJECTORY, value):
        data[name].append(item)
  return {name: np.array(data[name]) for name in data}


def read_atmosphere_table(file_table):
  # CCMC 形式のテーブルから高度 [m] と質量密度 [kg/m3] を取る
  altitude, density = [], []
  flag_data = False
  with open(file_table) as stream:
    for line in stream:
      item = line.split()
      if len(item) == 7 and item[0] == '1' and item[1] == '2' :
        flag_data = True
        continue
      if not flag_data or len(item) != 7 :
        continue
      try:
        value = [float(entry) for entry in item]
      except ValueError:
        continue
      altitude.append(value[0]*1.0e3)
      density.append(value[4]*1.0e3)      # g/cm3 -> kg/m3
  return np.array(altitude), np.array(density)


def difference(time_reference, value_reference, time_case, value_case, window=None):
  #
  # 参照の点ごとに、計算値との差を返す。参照が NaN の点と、計算の範囲外の点は外す
  #
  mask = np.isfinite(value_reference) \
         & (time_reference >= time_case[0]) & (time_reference <= time_case[-1])
  if window is not None :
    mask = mask & (time_reference >= window[0]) & (time_reference <= window[1])
  if not mask.any() :
    return np.array([]), np.array([])
  interpolated = np.interp(time_reference[mask], time_case, value_case)
  return time_reference[mask], interpolated - value_reference[mask]


def report(name, unit, time, error):
  if len(error) == 0 :
    print('  {:<16s} no overlapping points'.format(name))
    return
  index = np.argmax(np.abs(error))
  print('  {:<16s} {:4d} points, max {:+9.2f} {:s} at t = {:.2f} s, rms {:8.2f} {:s}'.format(
        name, len(error), error[index], unit, time[index], np.sqrt(np.mean(error**2)), unit))


def compare_atmosphere(data_reference, file_table):
  #
  # Tacode が飛ぶ大気テーブルと、表 V の密度（ゾンデ実測）の比。
  # **軌道を通していない**ので、この比は計算結果に依らない
  #
  altitude_table, density_table = read_atmosphere_table(file_table)
  mask = np.isfinite(data_reference['density']) & np.isfinite(data_reference['altitude']) \
         & (data_reference['altitude'] <= altitude_table[-1])
  altitude = data_reference['altitude'][mask]
  ratio    = data_reference['density'][mask]/np.interp(altitude, altitude_table, density_table)

  print('Atmosphere table against the density measured in flight (table V)')
  print('  {:d} points, mean ratio flight/table {:.4f}, rms deviation {:.1f} %'.format(
        len(ratio), ratio.mean(), 100.0*np.std(ratio/ratio.mean())))
  for lower, upper in ((0.0, 20.0), (20.0, 40.0), (40.0, 60.0),
                       (60.0, 80.0), (80.0, 100.0), (100.0, 125.0)):
    select = (altitude >= lower*1.0e3) & (altitude < upper*1.0e3)
    if select.any() :
      print('    {:3.0f} to {:3.0f} km: {:4d} points, mean ratio {:.3f}'.format(
            lower, upper, int(select.sum()), ratio[select].mean()))
  return altitude, ratio


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--case', type=str, nargs='+', default=None,
                      help='Tecplot outputs to compare, in place of the default')
  parser.add_argument('--label', type=str, nargs='+', default=None,
                      help='Legend entry for each --case')
  parser.add_argument('--time-start', type=float, default=TIME_START_DEFAULT,
                      help='Elapsed flight time the run starts at, s')
  parser.add_argument('--file-config', type=str, default='config.yml')
  parser.add_argument('--no-figure', action='store_true')
  args = parser.parse_args()

  os.chdir(DIRECTORY_SCRIPT)
  config = general().read_config_yaml(args.file_config)

  if args.case is None :
    case_list = list(CASE_DEFAULT)
  else:
    label_list = args.label if args.label is not None else args.case
    case_list  = list(zip(args.case, label_list))

  data_reference = read_trajectory(FILE_TRAJECTORY)

  data_case_list = []
  for file_case, label in case_list:
    if not os.path.exists(file_case) :
      print('No such output: ' + file_case)
      print('--Run ./run_tacode.sh first.')
      print('Program stopped.')
      sys.exit(1)
    data = tecplot_reader.read_tecplot(file_case)
    data['Time'] = data['Time'] + args.time_start
    data_case_list.append((data, label))

  print('Tacode against NASA TN D-3569 table V')
  print('--Table V is NASA\'s own 3-DOF simulation, anchored to the radar and flown')
  print('  through the measured atmosphere. Agreement with it is a check of the')
  print('  solver, not of the physics of the flight.')
  for data, label in data_case_list:
    for window, title in ((None, 'whole entry'),
                          (WINDOW_HYPERSONIC, 'hypersonic, to t = {:.0f} s'.format(
                             WINDOW_HYPERSONIC[1]))):
      print('{:s} ({:s}):'.format(label, title))
      for name, unit, scale in (('altitude', 'm', 1.0e3), ('velocity', 'm/s', 1.0),
                                ('flightpath', 'deg.', 1.0)):
        if name == 'altitude' :
          value_case = data['Alti']*scale
        elif name == 'velocity' :
          value_case = data['VelplAbs']
        else:
          value_case = np.degrees(np.arcsin(data['Wpl']/data['VelplAbs']))
        time, error = difference(data_reference['time'], data_reference[name],
                                 data['Time'], value_case, window)
        report(name, unit, time, error)

  file_table = os.path.join(config['atmosphere']['directory_atmosphere'],
                            config['atmosphere']['filename_atmosphere'])
  altitude_ratio, ratio = compare_atmosphere(data_reference, file_table)

  if not args.no_figure :
    make_figure(data_case_list, data_reference, altitude_ratio, ratio)

  return


def make_figure(data_case_list, data_reference, altitude_ratio, ratio):
  try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
  except ImportError:
    print('Caution: matplotlib is not installed, so no figure is written.')
    return

  if not os.path.isdir(DIRECTORY_OUTPUT) :
    os.makedirs(DIRECTORY_OUTPUT)

  # 高度と速度の履歴
  for name, unit, scale, column in (('altitude', 'km', 1.0e-3, 'Alti'),
                                    ('velocity', 'm/s', 1.0, 'VelplAbs')):
    figure, axis = plt.subplots(2, 1, figsize=(7.0, 6.0), sharex=True,
                                gridspec_kw={'height_ratios': [2, 1]})
    reference = data_reference[name]*(scale if name == 'altitude' else 1.0)
    axis[0].plot(data_reference['time'], reference, 'k.', markersize=2,
                 label='NASA TN D-3569 table V')
    for index, (data, label) in enumerate(data_case_list):
      value = data[column]*(1.0 if name == 'altitude' else scale)
      axis[0].plot(data['Time'], value, color=COLOUR_CASE[index % len(COLOUR_CASE)],
                   linewidth=1.2, label='Tacode, ' + label)
      time, error = difference(data_reference['time'], data_reference[name],
                               data['Time'],
                               data[column]*(1.0e3 if name == 'altitude' else 1.0))
      axis[1].plot(time, error*(scale if name == 'altitude' else 1.0),
                   color=COLOUR_CASE[index % len(COLOUR_CASE)], linewidth=1.0)
    axis[0].set_ylabel(name.capitalize() + ' [' + unit + ']')
    axis[0].legend(fontsize=8)
    axis[0].grid(alpha=0.3)
    axis[1].axhline(0.0, color='k', linewidth=0.8)
    axis[1].set_ylabel('Tacode - table V [' + unit + ']')
    axis[1].set_xlabel('Elapsed flight time [s]')
    axis[1].grid(alpha=0.3)
    figure.tight_layout()
    path = os.path.join(DIRECTORY_OUTPUT, 'profile_time-' + name + '.png')
    figure.savefig(path, dpi=150)
    plt.close(figure)
    print('Wrote ' + path)

  # 大気テーブルと飛行実測の密度比
  figure, axis = plt.subplots(figsize=(5.0, 5.0))
  axis.plot(ratio, altitude_ratio*1.0e-3, 'k.', markersize=3)
  axis.axvline(1.0, color='tab:red', linewidth=1.0)
  axis.set_xlabel('Density, flight (table V) / NRLMSISE-00 table')
  axis.set_ylabel('Altitude [km]')
  axis.set_xlim(0.5, 1.5)
  axis.grid(alpha=0.3)
  figure.tight_layout()
  path = os.path.join(DIRECTORY_OUTPUT, 'profile_density-ratio.png')
  figure.savefig(path, dpi=150)
  plt.close(figure)
  print('Wrote ' + path)


if __name__ == '__main__':
  main()
