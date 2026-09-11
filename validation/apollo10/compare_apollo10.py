#!/usr/bin/env python3
#
# Apollo 10（AS-505）の飛行データと Tacode を突き合わせる。
#
# Apollo 4 のケース（validation/apollo4）と違い、**揚力の向きを仮定していない**。
# NASA TN D-6725 は誘導が実際に打ったロール角の履歴（図 12）を載せており、それを
# そのまま config の bank_angle_table に入れている。したがってここで見ているのは
#
#   「突入界面の状態ベクトルと、実機のバンク角の履歴を与えたとき、
#     Tacode は実機の高度・速度の履歴を再現するか」
#
# であって、合わせ込む自由度は残っていない（弾道係数の元になる質量だけが仮定）。
#
# 速度の扱い:
#   図 13 の速度は**慣性系**（突入界面で 36 309 ft/s）。Tacode が出すのは ECEF 速度
#   なので、比較のときに自転ぶんを足して慣性速度にする:
#
#     |v_in|^2 = |v_ecef|^2 + 2 v_rot u_east + v_rot^2,  v_rot = omega sqrt(X^2+Y^2)
#
#   （v_rot は回転系そのものの東向き速度、u_east は ECEF 速度の東成分）
#
# 使い方:
#   cd validation/apollo10 && ./run_tacode.sh
#   python3 compare_apollo10.py
#   python3 compare_apollo10.py --no-figure

import argparse
import os
import sys

import numpy as np

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src'))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src_helper/animate_trajectory'))

import tecplot_reader as tecplot_reader      # noqa: E402
from general.general import general          # noqa: E402

CASE_DEFAULT = (
  ('output_result/tecplot.dat', '3-DOF, measured bank angle'),
)

COLOUR_CASE = ('tab:red', 'tab:purple', 'tab:orange', 'tab:brown')


def velocity_inertial(data, rotation_rate):
  # ECEF 速度を慣性速度の大きさに直す
  radius_equatorial = np.sqrt(data['X']**2 + data['Y']**2)*1.0e3
  velocity_rotation = rotation_rate*radius_equatorial
  return np.sqrt(data['VelplAbs']**2 + 2.0*velocity_rotation*data['Upl'] + velocity_rotation**2)


def compare_trajectory(data_case, data_flight, label, rotation_rate):
  time_flight = data_flight['Time']
  time_case   = data_case['Time']
  mask        = time_flight <= time_case[-1]

  altitude_case = np.interp(time_flight[mask], time_case, data_case['Alti'])
  velocity_case = np.interp(time_flight[mask], time_case, velocity_inertial(data_case, rotation_rate))

  difference_altitude = altitude_case - data_flight['Alti'][mask]
  difference_velocity = velocity_case - data_flight['VelinAbs'][mask]

  print('')
  print('--' + label)
  print('  points compared            : {:d} of {:d} (the run ends at {:.0f} s)'.format(
        int(np.sum(mask)), len(time_flight), time_case[-1]))
  print('  altitude difference        : max {:+.1f} km, rms {:.1f} km'.format(
        difference_altitude[np.argmax(np.abs(difference_altitude))], np.sqrt(np.mean(difference_altitude**2))))
  print('  velocity difference        : max {:+.0f} m/s, rms {:.0f} m/s'.format(
        difference_velocity[np.argmax(np.abs(difference_velocity))], np.sqrt(np.mean(difference_velocity**2))))

  # 最初のディップ（第 1 ピーク加熱・加速度の場所）
  index_case   = int(np.argmin(data_case['Alti'][data_case['Time'] < 250.0]))
  index_flight = int(np.argmin(data_flight['Alti'][data_flight['Time'] < 250.0]))
  print('  first minimum of altitude  : {:.1f} km at {:.0f} s (flight: {:.1f} km at {:.0f} s)'.format(
        data_case['Alti'][index_case], data_case['Time'][index_case],
        data_flight['Alti'][index_flight], data_flight['Time'][index_flight]))

  return


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--case-directory', type=str, default=DIRECTORY_SCRIPT)
  parser.add_argument('--file-config', type=str, default='config.yml')
  parser.add_argument('--file-reference', type=str, default='reference/apollo10_flight.dat')
  parser.add_argument('--file-roll', type=str, default='reference/apollo10_roll.dat')
  parser.add_argument('--case', action='append', nargs=2, metavar=('PATH', 'LABEL'), default=None,
                      help='A Tecplot output and the label to draw it with; may be given more than once')
  parser.add_argument('--output-directory', type=str, default='output_comparison')
  parser.add_argument('--no-figure', action='store_true')
  args = parser.parse_args()

  os.chdir(args.case_directory)

  case_list = args.case if args.case is not None else \
              [ [path, label] for path, label in CASE_DEFAULT if os.path.exists(path) ]
  if len(case_list) == 0 :
    print('No result to compare with.')
    print('--Run the case first: ./run_tacode.sh')
    print('Program stopped.')
    sys.exit(1)

  config = general().read_config_yaml(args.file_config)
  rotation_rate = config['planet']['rotation_rate']

  data_flight = tecplot_reader.read_tecplot(args.file_reference)
  data_roll   = tecplot_reader.read_tecplot(args.file_roll)
  data_case_list = [ (tecplot_reader.read_tecplot(path), label) for path, label in case_list ]

  print('Apollo 10 (AS-505) entry: flight data against Tacode')
  print('--Reference :', os.path.join(args.case_directory, args.file_reference))
  print('--Bank angle: the roll angle measured in flight (an input, not a result)')

  print('')
  print('(1) Trajectory against the flight data')
  for data_case, label in data_case_list:
    compare_trajectory(data_case, data_flight, label, rotation_rate)

  if not args.no_figure :
    make_figure(args, data_case_list, data_flight, data_roll, rotation_rate)

  return


def make_figure(args, data_case_list, data_flight, data_roll, rotation_rate):
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

  def finalise(axis, filename):
    axis.grid(linewidth=0.5, alpha=0.7)
    axis.legend(frameon=False)
    plt.tight_layout()
    path = os.path.join(args.output_directory, filename)
    plt.savefig(path, dpi=200)
    plt.close()
    print('--Figure:', path)

  print('')
  print('Figures')

  # 高度
  figure, axis = plt.subplots(figsize=(6, 4))
  axis.plot(data_flight['Time'], data_flight['Alti'], color=colour_flight, linewidth=2.0,
            label='Apollo 10 flight')
  for index, (data_case, label) in enumerate(data_case_list):
    axis.plot(data_case['Time'], data_case['Alti'], color=COLOUR_CASE[index % len(COLOUR_CASE)],
              linestyle='--', linewidth=2.0, label='Tacode, ' + label)
  axis.set_xlabel('Time from entry interface, s')
  axis.set_ylabel('Altitude, km')
  axis.set_xlim(0.0, data_flight['Time'][-1])
  axis.set_ylim(0.0, 130.0)
  finalise(axis, 'profile_time-altitude.png')

  # 慣性速度
  figure, axis = plt.subplots(figsize=(6, 4))
  axis.plot(data_flight['Time'], data_flight['VelinAbs'], color=colour_flight, linewidth=2.0,
            label='Apollo 10 flight')
  for index, (data_case, label) in enumerate(data_case_list):
    axis.plot(data_case['Time'], velocity_inertial(data_case, rotation_rate),
              color=COLOUR_CASE[index % len(COLOUR_CASE)], linestyle='--', linewidth=2.0,
              label='Tacode, ' + label)
  axis.set_xlabel('Time from entry interface, s')
  axis.set_ylabel('Inertial velocity, m/s')
  axis.set_xlim(0.0, data_flight['Time'][-1])
  finalise(axis, 'profile_time-velocity.png')

  # 速度-高度
  figure, axis = plt.subplots(figsize=(6, 4))
  axis.plot(data_flight['VelinAbs'], data_flight['Alti'], color=colour_flight, linewidth=2.0,
            label='Apollo 10 flight')
  for index, (data_case, label) in enumerate(data_case_list):
    axis.plot(velocity_inertial(data_case, rotation_rate), data_case['Alti'],
              color=COLOUR_CASE[index % len(COLOUR_CASE)], linestyle='--', linewidth=2.0,
              label='Tacode, ' + label)
  axis.set_xlabel('Inertial velocity, m/s')
  axis.set_ylabel('Altitude, km')
  axis.set_ylim(0.0, 130.0)
  finalise(axis, 'profile_velocity-altitude.png')

  # 入力したバンク角
  figure, axis = plt.subplots(figsize=(6, 3.2))
  axis.plot(data_roll['Time'], data_roll['Bank'], color='tab:blue', linewidth=1.8,
            label='Roll angle measured in flight (input)')
  axis.set_xlabel('Time from entry interface, s')
  axis.set_ylabel('Bank angle, deg.')
  axis.set_xlim(0.0, data_flight['Time'][-1])
  axis.set_yticks([-180, -90, 0, 90, 180])
  finalise(axis, 'profile_time-bank.png')

  return


if __name__ == '__main__':
  main()
