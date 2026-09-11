#!/usr/bin/env python3
#
# Apollo 10 の飛行データを NASA TN D-6725 の図から読み取る。
#
#   図 12: ロール角（実機の姿勢。**これが揚力の向き＝バンク角**）と荷重倍数
#   図 13: 高度と（慣性系の）速度
#   表 II: 突入界面の状態ベクトル（本文に数値があるので読み取りは不要）
#
# **これがこの検証ケースの肝**である。Apollo 4（TM X-58091）には軌道の表があっても
# バンク角が無く、定数バンクを仮定するしかなかった。Apollo 10 は誘導が実際に打った
# ロール角の履歴が載っているので、**揚力の向きを仮定せずに軌道を突き合わせられる**。
#
# 読み取りの手順（どれも figure の画素から決めており、目分量はしていない）:
#
#   1. pdftoppm で 18 ページ目を 600 dpi の PNG にする
#   2. 軸線・目盛りの画素位置から座標変換を作る（下の CALIBRATION_*）
#   3. 列ごとに黒画素の塊（クラスタ）を取り、直前までの点に直線を当てて予測し、
#      予測に最も近いクラスタを採る
#   4. 図 13 は高度と速度が t = 3.3 min あたりで交差する。速度を先に追跡し、
#      高度は**速度の曲線から 13 px 以内のクラスタを除いて**前後両方向から追跡して繋ぐ
#      （交差の前後で上下が入れ替わるので、片方向だけだと乗り移る）
#
# 出力（このディレクトリの reference/）:
#
#   apollo10_flight.dat   時刻 [s]・高度 [km]・慣性速度 [m/s]
#   apollo10_roll.dat     時刻 [s]・バンク角 [deg.]
#
# 併せて config.yml に貼る bank_angle_table の YAML も標準出力に書く。
#
# 必要なもの: poppler の pdftoppm、Pillow、numpy。**参照 PDF はリポジトリの外**に
# あるので、読み取り済みの reference/*.dat をリポジトリに置いてある。
#
# 使い方:
#   python3 digitize_figure.py                    # reference/ を作り直す
#   python3 digitize_figure.py --check figure.png # 読み取りを図に重ねた確認図を書く

import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np


DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
DIRECTORY_REFERENCE_DEFAULT = os.path.join(DIRECTORY_SCRIPT, '../../../references/20260423_Apollo_radiation')
FILE_REPORT = 'Apollo_FlightPathAngle/19720013191.pdf'
PAGE_FIGURE = 26          # pdftoppm のページ番号（本文 18 ページ）
RESOLUTION  = 600         # dpi

# 600 dpi で描いたページから図を切り出す矩形 (left, upper, right, lower)
CROP_FIGURE_12 = (500, 2350, 2500, 4300)
CROP_FIGURE_13 = (2500, 2350, 4700, 4300)

# 図 12（ロール角）の較正: 画素 -> データ
#   縦軸 0 deg. の行、1 度あたりの画素、時刻 0 の列、1 分あたりの画素
CALIBRATION_ROLL = {'row_zero': 544.0, 'pixel_per_degree': 1.95,
                    'column_zero': 350.0, 'pixel_per_minute': 193.3}

# 図 13（高度・速度）の較正
CALIBRATION_TRAJECTORY = {'column_zero': 720.0, 'pixel_per_minute': 153.5,
                          'row_altitude_zero': 1777.0, 'pixel_per_altitude': 154.1/40000.0,
                          'row_velocity_zero': 1772.0, 'pixel_per_velocity': 155.7/4000.0}

FOOT_TO_METRE = 0.3048


def render_page(file_pdf, directory_work):
  # PDF の 1 ページを PNG にする
  prefix = os.path.join(directory_work, 'page')
  command = ['pdftoppm', '-f', str(PAGE_FIGURE), '-l', str(PAGE_FIGURE),
             '-r', str(RESOLUTION), '-png', file_pdf, prefix]
  try:
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
  except (OSError, subprocess.CalledProcessError) as instance:
    print('Cannot render the report with pdftoppm:', instance)
    print('--Install poppler-utils, or keep the reference files already in reference/.')
    print('Program stopped.')
    sys.exit(1)

  for name in sorted(os.listdir(directory_work)):
    if name.startswith('page') and name.endswith('.png') :
      return os.path.join(directory_work, name)

  print('pdftoppm wrote no page.')
  print('Program stopped.')
  sys.exit(1)


def load_image(path, crop):
  try:
    from PIL import Image
  except ImportError:
    print('Pillow is not installed; it is needed to read the figures.')
    print('--Install it with: pip install pillow')
    print('Program stopped.')
    sys.exit(1)
  return np.array(Image.open(path).convert('L').crop(crop))


def cluster_column(mask, column):
  # 1 列の黒画素を、途切れ 10 px 以内でまとめた塊の中心にする
  row_list = np.where(mask[:, column])[0]
  if len(row_list) == 0 :
    return []
  group = []
  for row in row_list:
    if group and row - group[-1][-1] <= 10 :
      group[-1].append(row)
    else:
      group.append([row])
  return [float(np.mean(item)) for item in group]


def track_curve(mask, column_start, row_start, column_end, step=3,
                window=8, tolerance=45.0, forbidden=None):
  #
  # 曲線を 1 本たどる。直前 window 点に直線を当てて次の列の位置を予測し、
  # 予測に最も近いクラスタを採る。予測から tolerance 以上離れていれば見送る
  # （破線の隙間や、ほぼ垂直な区間で他の曲線に乗り移らないようにするため）。
  # forbidden(column) が返す行の近傍は候補から外す（交差する別の曲線を避ける）。
  #
  direction = 1 if column_end > column_start else -1
  column_list = [column_start]
  row_list    = [row_start]
  column = column_start
  while (column + step*direction - column_end)*direction < 0 :
    column += step*direction
    candidate = cluster_column(mask, column)
    if forbidden is not None :
      row_forbidden = forbidden(column)
      candidate = [row for row in candidate if abs(row - row_forbidden) > 13.0]
    if len(candidate) == 0 :
      continue
    number = min(window, len(column_list))
    if number >= 4 :
      coefficient = np.polyfit(column_list[-number:], row_list[-number:], 1)
      prediction  = np.polyval(coefficient, column)
    else :
      prediction = row_list[-1]
    row = min(candidate, key=lambda value: abs(value - prediction))
    if abs(row - prediction) > tolerance :
      continue
    column_list.append(column)
    row_list.append(row)

  return np.array(column_list, dtype=float), np.array(row_list, dtype=float)


def read_roll(image):
  # 図 12 の上のパネル。ロール角（実機）と指令はほぼ重なっており、離れるのは
  # 1.5-2.4 min の反転のあいだだけ。そこでは指令が +180 に張り付き、実機は
  # -180 側を回る。**バンク角としては同じ姿勢**（cos も sin も一致する）なので、
  # 連続な枝をたどればよい。
  calibration = CALIBRATION_ROLL
  mask = image < 128
  panel = np.zeros_like(mask)
  panel[180:905, 356:1925] = True
  work = mask & panel

  column_list = []
  row_list    = []
  previous    = 0.0
  for column in range(356, 1925, 5):
    candidate = cluster_column(work, column)
    # 凡例（右上）を落とす
    candidate = [row for row in candidate
                 if not (column > 950 and (calibration['row_zero'] - row)/calibration['pixel_per_degree'] > 115.0)]
    if len(candidate) == 0 or len(candidate) >= 4 :
      continue                      # ほぼ垂直な区間は飛ばし、前後から補間させる
    value = [(calibration['row_zero'] - row)/calibration['pixel_per_degree'] for row in candidate]
    angle = min(value, key=lambda item: abs(item - previous))
    # +-180 のラップ（同じ姿勢）
    if abs(angle - previous) > 150.0 and abs(abs(angle) - 180.0) < 40.0 and abs(abs(previous) - 180.0) < 40.0 :
      angle = np.sign(previous)*abs(angle)
    column_list.append(column)
    row_list.append(angle)
    previous = angle

  time = (np.array(column_list) - calibration['column_zero'])/calibration['pixel_per_minute']*60.0
  return time, np.array(row_list)


def read_trajectory(image):
  # 図 13。高度と速度は t = 3.3 min あたりで交差する
  calibration = CALIBRATION_TRAJECTORY
  mask = image < 128
  panel = np.zeros_like(mask)
  panel[190:1745, 745:1800] = True
  work = mask & panel
  work[300:560, 1150:1800] = False     # 凡例
  work[190:300, 745:960]   = False     # 「40 X 10^3」のラベル

  # 速度: 突入界面では 36 000 ft/s で平ら
  column_velocity, row_velocity = track_curve(work, 760, 371.0, 1798)

  def row_of_velocity(column):
    return np.interp(column, column_velocity, row_velocity)

  # 高度: 速度の曲線を避けて、前後両方向から
  column_forward,  row_forward  = track_curve(work, 750, 438.0, 1240, forbidden=row_of_velocity)
  column_backward, row_backward = track_curve(work, 1780, 1440.0, 1320, forbidden=row_of_velocity)
  column_altitude = np.concatenate([column_forward, column_backward[::-1]])
  row_altitude    = np.concatenate([row_forward, row_backward[::-1]])

  def time_of(column):
    return (column - calibration['column_zero'])/calibration['pixel_per_minute']*60.0

  altitude = (calibration['row_altitude_zero'] - row_altitude)/calibration['pixel_per_altitude']*FOOT_TO_METRE*1.0e-3
  velocity = (calibration['row_velocity_zero'] - row_velocity)/calibration['pixel_per_velocity']*FOOT_TO_METRE

  return (time_of(column_altitude), altitude), (time_of(column_velocity), velocity)


def write_reference(directory, time_altitude, altitude, time_velocity, velocity, time_roll, roll):
  # 高度と速度は時刻が違うので、高度の時刻に速度を内挿してひとつの表にする
  velocity_on_altitude = np.interp(time_altitude, time_velocity, velocity)

  path = os.path.join(directory, 'apollo10_flight.dat')
  with open(path, 'w') as f:
    f.write('# Apollo 10 (AS-505) entry, 1969-05-26, flight data in SI units\n')
    f.write('# --Source: NASA TN D-6725 (Graves and Harpold, 1972), figure 13,\n')
    f.write('#   digitised by validation/apollo10/digitize_figure.py.\n')
    f.write('# --Time is seconds from the entry interface (400 000 ft).\n')
    f.write('# --The velocity is the INERTIAL speed, as plotted in the report.\n')
    f.write('Variables = Time[s],Alti[km],VelinAbs[m/s]\n')
    f.write('zone t="Apollo 10 flight" i= {:d} f=point\n'.format(len(time_altitude)))
    for index in range(0, len(time_altitude)):
      f.write('{:15.7e}{:15.7e}{:15.7e}\n'.format(
              time_altitude[index], altitude[index], velocity_on_altitude[index]))
  print('Reference file written:', path)

  path = os.path.join(directory, 'apollo10_roll.dat')
  with open(path, 'w') as f:
    f.write('# Apollo 10 (AS-505) entry, roll (bank) angle measured in flight\n')
    f.write('# --Source: NASA TN D-6725, figure 12 (the "Roll angle" trace),\n')
    f.write('#   digitised by validation/apollo10/digitize_figure.py.\n')
    f.write('# --Time is seconds from the entry interface. 0 deg. is lift up.\n')
    f.write('Variables = Time[s],Bank[deg.]\n')
    f.write('zone t="Apollo 10 roll" i= {:d} f=point\n'.format(len(time_roll)))
    for index in range(0, len(time_roll)):
      f.write('{:15.7e}{:15.7e}\n'.format(time_roll[index], roll[index]))
  print('Reference file written:', path)


def print_bank_table(time_roll, roll, interval):
  # config.yml の satellite.bank_angle_table に貼る形で書き出す
  time_node = np.arange(0.0, time_roll[-1] + 0.5*interval, interval)
  angle     = np.interp(time_node, time_roll, roll)
  print('')
  print('  # Bank angle [deg.] against the time from entry [s]')
  print('  bank_angle_table:')
  for index in range(0, len(time_node)):
    print('    - [{:6.1f}, {:8.2f}]'.format(time_node[index], angle[index]))


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--reference-directory', type=str, default=DIRECTORY_REFERENCE_DEFAULT,
                      help='Directory holding the scanned report (outside the repository)')
  parser.add_argument('--output-directory', type=str,
                      default=os.path.join(DIRECTORY_SCRIPT, 'reference'))
  parser.add_argument('--bank-interval', type=float, default=5.0,
                      help='Step [s] of the bank table printed for config.yml')
  parser.add_argument('--check', type=str, default=None,
                      help='Write a figure with the digitised curves drawn over the scan')
  args = parser.parse_args()

  file_pdf = os.path.join(args.reference_directory, FILE_REPORT)
  if not os.path.exists(file_pdf) :
    print('The report was not found:', file_pdf)
    print('--Give the directory which holds it with --reference-directory.')
    print('--The digitised result is already in reference/, so this script is only')
    print('--needed to reproduce it.')
    print('Program stopped.')
    sys.exit(1)

  with tempfile.TemporaryDirectory() as directory_work:
    path_page = render_page(file_pdf, directory_work)
    image_roll       = load_image(path_page, CROP_FIGURE_12)
    image_trajectory = load_image(path_page, CROP_FIGURE_13)

  time_roll, roll = read_roll(image_roll)
  (time_altitude, altitude), (time_velocity, velocity) = read_trajectory(image_trajectory)

  print('Roll angle : {:d} points, {:.0f} to {:.0f} s, {:.0f} to {:.0f} deg.'.format(
        len(time_roll), time_roll[0], time_roll[-1], roll.min(), roll.max()))
  print('Altitude   : {:d} points, {:.0f} to {:.0f} s, {:.1f} to {:.1f} km'.format(
        len(time_altitude), time_altitude[0], time_altitude[-1], altitude.min(), altitude.max()))
  print('Velocity   : {:d} points, {:.0f} to {:.0f} s, {:.0f} to {:.0f} m/s'.format(
        len(time_velocity), time_velocity[0], time_velocity[-1], velocity.min(), velocity.max()))

  write_reference(args.output_directory, time_altitude, altitude, time_velocity, velocity,
                  time_roll, roll)
  print_bank_table(time_roll, roll, args.bank_interval)

  if args.check is not None :
    draw_check(args.check, image_roll, image_trajectory, time_roll, roll,
               time_altitude, altitude, time_velocity, velocity)

  return


def draw_check(filename, image_roll, image_trajectory, time_roll, roll,
               time_altitude, altitude, time_velocity, velocity):
  try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
  except ImportError:
    print('matplotlib is not installed, so no check figure is written.')
    return

  cal_r = CALIBRATION_ROLL
  cal_t = CALIBRATION_TRAJECTORY
  figure, axis = plt.subplots(1, 2, figsize=(16, 6))

  axis[0].imshow(image_roll[150:950, 300:1980], cmap='gray', extent=[300, 1980, 950, 150])
  axis[0].plot(cal_r['column_zero'] + time_roll/60.0*cal_r['pixel_per_minute'],
               cal_r['row_zero'] - roll*cal_r['pixel_per_degree'], 'r-', linewidth=1.0)
  axis[0].set_title('figure 12: roll angle')

  axis[1].imshow(image_trajectory[180:1800, 650:1900], cmap='gray', extent=[650, 1900, 1800, 180])
  axis[1].plot(cal_t['column_zero'] + time_altitude/60.0*cal_t['pixel_per_minute'],
               cal_t['row_altitude_zero'] - altitude*1.0e3/FOOT_TO_METRE*cal_t['pixel_per_altitude'],
               'r-', linewidth=1.2, label='altitude')
  axis[1].plot(cal_t['column_zero'] + time_velocity/60.0*cal_t['pixel_per_minute'],
               cal_t['row_velocity_zero'] - velocity/FOOT_TO_METRE*cal_t['pixel_per_velocity'],
               'b-', linewidth=1.2, label='velocity')
  axis[1].legend(loc='lower left')
  axis[1].set_title('figure 13: altitude and velocity')

  plt.tight_layout()
  plt.savefig(filename, dpi=130)
  plt.close()
  print('Check figure written:', filename)

  return


if __name__ == '__main__':
  main()
