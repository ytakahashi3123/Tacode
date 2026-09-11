#!/usr/bin/env python3
#
# Apollo 4 が実際に飛んだ**揚力の向き**を NASA TN D-5399 の図から読み取る。
#
#   E. R. Hillje, "Entry Aerodynamics at Lunar Return Conditions Obtained From the
#   Flight of Apollo 4 (AS-501)", NASA TN D-5399, 1969.
#
#   図 6(b)  "Bank-angle time history"（本文 34 ページ）
#            誘導計算機の制御履歴としてのバンク角
#   図 16(a) "Vertical lift-to-drag ratio (L/D)v"（本文 59 ページ）
#            **飛行データから求めた鉛直方向の揚抗比**。加速度計と再構成軌道から出ており、
#            軌道を決めるのはこの量そのもの
#
# **ケースを駆動するのは図 16(a) のほう**である。図 6(b) と図 16(a) は
# (L/D)v = (L/D)_R cos(bank) の関係で概ね一致する（差の rms は 0.03）が、ロール反転の
# 直後、機体が最も深いところに居て動圧が最大のあたりでは、バンク角の読み取りが数秒
# ずれるだけで (L/D)v が倍ほど変わる。そこは軌道が最も敏感な区間でもあるので、
# **鉛直揚力そのものを読んだほうが素性がよい**:
#
#   bank(t) = arccos( (L/D)v(t) / (L/D)_R ),   (L/D)_R = 0.368（同報告書の飛行由来値）
#
# 横方向（左右）の符号は落ちるが、軌道の上下動には効かない。
#
# 読み取りの手順（どちらの図も同じ）:
#
#   1. pdftoppm でページを 300/600 dpi にし、**90 度回して**縦置きを直す
#   2. 軸ラベルの文字の塊を画素で測って座標変換を作る（目分量では 20 秒ずれた）
#   3. データは丸記号。罫線が密なので、**局所の黒画素密度**で記号と罫線を分ける
#   4. 塊の重心を点にし、1 秒ごとの中央値にまとめる
#
# 出力（reference/）:
#
#   apollo4_bank.dat          図 6(b) のバンク角（記録として）
#   apollo4_liftvertical.dat  図 16(a) の鉛直揚抗比（**ケースを駆動する量**）
#
# 併せて config_lift.yml に貼る bank_angle_table の YAML を標準出力に書く。
#
import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np


DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
DIRECTORY_REFERENCE_DEFAULT = os.path.join(DIRECTORY_SCRIPT, '../../../references/20260423_Apollo_radiation')
FILE_REPORT = 'Apollo_EntryAerodynamics/19690029435.pdf'
PAGE_BANK   = 42          # pdftoppm のページ番号（本文 34 ページ、図 6(b)）
PAGE_LIFT   = 67          # 同（本文 59 ページ、図 16(a)）
RESOLUTION_BANK = 300     # dpi
RESOLUTION_LIFT = 600     # dpi（記号が小さいので細かく）

# 報告書が飛行データから求めた合力の揚抗比（トリム後から第 1 ピーク g までの平均）
LIFT_DRAG_RESULTANT = 0.368
# 同じ報告書が挙げる最大値（M = 6 で 0.41）。鉛直成分はこれを超え得ない
LIFT_DRAG_MAXIMUM = 0.41

# 90 度回したあとの画素 -> データ
#   時刻: GET 30 040 s が column_reference、1 秒あたり pixel_per_second
#   角度: 0 度が row_zero、1 度あたり pixel_per_degree
# 較正は目分量ではなく、軸ラベルの文字の塊を画素で測って決めた:
#   時刻  "29 960" の中心 x = 753、"30 040" は 1030、"30 120" は 1307.5  -> 3.4656 px/s
#   角度  160 -> y 586、0 -> 1151、-160 -> 1718.5                        -> 3.540 px/deg
CALIBRATION_BANK = {'column_reference': 1030.0, 'time_reference': 30040.0, 'pixel_per_second': 3.4656,
                    'row_zero': 1151.0, 'pixel_per_value': 3.540}
WINDOW_BANK = {'row': (625, 1755), 'column': (780, 2875)}

# 図 16(a)（600 dpi）: "29 960" の中心 x = 1253、80 秒あたり 449.2 px -> 5.615 px/s
#                      (L/D)v = 0 は y = 2506、1 単位あたり 2352 px
CALIBRATION_LIFT = {'column_reference': 1253.0, 'time_reference': 29960.0, 'pixel_per_second': 5.615,
                    'row_zero': 2506.0, 'pixel_per_value': 2352.0}
WINDOW_LIFT = {'row': (1340, 3672), 'column': (1490, 5400)}

# 突入界面の時刻（GET, s）。TN D-5399 表 I
TIME_ENTRY_INTERFACE = 29968.54


def render_page(file_pdf, directory_work, page, resolution):
  prefix = os.path.join(directory_work, 'page{:d}'.format(page))
  command = ['pdftoppm', '-f', str(page), '-l', str(page),
             '-r', str(resolution), '-png', file_pdf, prefix]
  try:
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
  except (OSError, subprocess.CalledProcessError) as instance:
    print('Cannot render the report with pdftoppm:', instance)
    print('--Install poppler-utils, or keep the reference file already in reference/.')
    print('Program stopped.')
    sys.exit(1)
  for name in sorted(os.listdir(directory_work)):
    if name.startswith('page{:d}'.format(page)) and name.endswith('.png') :
      return os.path.join(directory_work, name)
  print('pdftoppm wrote no page.')
  print('Program stopped.')
  sys.exit(1)


def read_symbol(path_image, window, calibration, density_size, density_threshold, size_minimum):
  #
  # 丸記号で描かれた曲線を読む。罫線は細く、記号は太い輪なので、
  # 局所の黒画素密度で分けられる。戻り値は (時刻 [s], 値) と回した画像。
  #
  try:
    from PIL import Image
  except ImportError:
    print('Pillow is not installed; it is needed to read the figures.')
    print('Program stopped.')
    sys.exit(1)
  try:
    from scipy import ndimage
  except ImportError:
    print('SciPy is not installed; it is needed to read the figures.')
    print('Program stopped.')
    sys.exit(1)

  # 図は縦置きに刷られているので回す
  image = np.array(Image.open(path_image).rotate(-90, expand=True).convert('L'))
  row_0, row_1       = window['row']
  column_0, column_1 = window['column']
  region_dark = (image[row_0:row_1, column_0:column_1] < 150).astype(float)

  density = ndimage.uniform_filter(region_dark, size=density_size)
  blob    = density > density_threshold

  label, number = ndimage.label(blob, structure=np.ones((3,3)))
  size = ndimage.sum(blob, label, range(1, number+1))
  slice_list = ndimage.find_objects(label)

  point = []
  for index, item in enumerate(slice_list):
    if size[index] < size_minimum :
      continue
    region = (label[item] == index+1)
    height = item[0].stop - item[0].start
    width  = item[1].stop - item[1].start
    if height <= 3*density_size and width <= 3*density_size :
      centre_row, centre_column = ndimage.center_of_mass(region)
      point.append((item[1].start + centre_column, item[0].start + centre_row))
    else :
      # 記号が連なった塊は列ごとに中心を採る
      for column in range(0, region.shape[1]):
        row_list = np.where(region[:, column])[0]
        if len(row_list) == 0 :
          continue
        point.append((item[1].start + column, item[0].start + row_list.mean()))

  point = np.array(sorted(point))
  time = calibration['time_reference'] \
       + (point[:,0] + column_0 - calibration['column_reference'])/calibration['pixel_per_second']
  value = (calibration['row_zero'] - (point[:,1] + row_0))/calibration['pixel_per_value']

  return time, value, image


def reduce_per_second(time, value):
  # 1 点につき記号の画素が何本も出るので、1 秒ごとの中央値にまとめる
  time_round = np.round(time).astype(int)
  time_node  = np.unique(time_round)
  value_node = np.array([np.median(value[time_round == item]) for item in time_node])
  return time_node.astype(float), value_node


def read_bank(path_image):
  # 図 6(b)。ロール反転のあいだは記号が疎で、機体は +-180 を通って回っている
  time, bank, image = read_symbol(path_image, WINDOW_BANK, CALIBRATION_BANK, 13, 0.28, 20)
  time_node, bank_node = reduce_per_second(time, bank)

  # 孤立した外れ値を落とす。前後 2 点のどれからも 80 度以上離れている点は、
  # 網掛けや枠を拾ったもの（機体が回っている区間の点は前後と 40 度程度でつながる）
  keep = np.ones(len(bank_node), dtype=bool)
  for index in range(0, len(bank_node)):
    neighbour = [bank_node[other] for other in range(max(0, index-2), min(len(bank_node), index+3))
                 if other != index]
    if len(neighbour) == 0 :
      continue
    distance = [abs((bank_node[index] - item + 180.0) % 360.0 - 180.0) for item in neighbour]
    if min(distance) > 80.0 :
      keep[index] = False
  time_node = time_node[keep]
  bank_node = bank_node[keep]

  # 最初のまばらな記号（大気がまだ効いていない区間で、読み取りが枠と紛れる）を落とし、
  # **点が 1-2 秒間隔で続く区間から**採る
  index_start = 0
  for index in range(0, len(time_node)-4):
    if np.all(np.diff(time_node[index:index+5]) <= 2.0) :
      index_start = index
      break
  time_node = time_node[index_start:]
  bank_node = bank_node[index_start:]
  while len(bank_node) > 1 and abs(bank_node[0] - bank_node[1]) > 60.0 :
    time_node = time_node[1:]
    bank_node = bank_node[1:]

  return time_node, bank_node, image


def read_lift_vertical(path_image):
  # 図 16(a)。図の上下の網掛けと、左上の凡例の枠を落とす。
  # **合力の揚抗比 0.41（報告書が挙げる最大値）を超える鉛直成分は物理的に有り得ない**ので、
  # そこで切るのが素直（凡例の枠はちょうど +0.43 に乗っていて、これに引っかかっていた）。
  time, value, image = read_symbol(path_image, WINDOW_LIFT, CALIBRATION_LIFT, 21, 0.34, 60)
  keep = (value > -LIFT_DRAG_MAXIMUM) & (value < LIFT_DRAG_MAXIMUM)
  time_node, value_node = reduce_per_second(time[keep], value[keep])

  return time_node, value_node, image


def write_reference(filename, time, value, title, column, note):
  with open(filename, 'w') as f:
    f.write('# Apollo 4 (AS-501) entry, ' + title + '\n')
    f.write('# --Source: NASA TN D-5399 (Hillje, 1969), digitised by\n')
    f.write('#   validation/apollo4/digitize_bank.py\n')
    for line in note:
      f.write('# --' + line + '\n')
    f.write('Variables = Time[s],' + column + '\n')
    f.write('zone t="Apollo 4" i= {:d} f=point\n'.format(len(time)))
    for index in range(0, len(time)):
      f.write('{:15.7e}{:15.7e}\n'.format(time[index], value[index]))
  print('Reference file written:', filename)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--reference-directory', type=str, default=DIRECTORY_REFERENCE_DEFAULT)
  parser.add_argument('--output-directory', type=str,
                      default=os.path.join(DIRECTORY_SCRIPT, 'reference'))
  parser.add_argument('--bank-interval', type=float, default=5.0,
                      help='Step [s] of the bank table printed for config_lift.yml')
  parser.add_argument('--lift-drag', type=float, default=LIFT_DRAG_RESULTANT,
                      help='Resultant lift-to-drag ratio the bank angle is computed with')
  parser.add_argument('--check', type=str, default=None,
                      help='Write a figure with the digitised points over the scans')
  args = parser.parse_args()

  file_pdf = os.path.join(args.reference_directory, FILE_REPORT)
  if not os.path.exists(file_pdf) :
    print('The report was not found:', file_pdf)
    print('--NASA TN D-5399 is public; it can be downloaded from')
    print('--https://ntrs.nasa.gov/citations/19690029435')
    print('--The digitised result is already in reference/, so this script is only')
    print('--needed to reproduce it.')
    print('Program stopped.')
    sys.exit(1)

  with tempfile.TemporaryDirectory() as directory_work:
    path_bank = render_page(file_pdf, directory_work, PAGE_BANK, RESOLUTION_BANK)
    time_bank_get, bank, image_bank = read_bank(path_bank)
    path_lift = render_page(file_pdf, directory_work, PAGE_LIFT, RESOLUTION_LIFT)
    time_lift_get, lift_vertical, image_lift = read_lift_vertical(path_lift)

  time_bank = time_bank_get - TIME_ENTRY_INTERFACE
  mask = time_bank >= 0.0
  time_bank, bank = time_bank[mask], bank[mask]

  time_lift = time_lift_get - TIME_ENTRY_INTERFACE
  mask = time_lift >= 0.0
  time_lift, lift_vertical = time_lift[mask], lift_vertical[mask]

  print('Bank angle (figure 6(b))      : {:d} points, {:.0f} to {:.0f} s, {:.0f} to {:.0f} deg.'.format(
        len(time_bank), time_bank[0], time_bank[-1], bank.min(), bank.max()))
  print('Vertical L/D (figure 16(a))   : {:d} points, {:.0f} to {:.0f} s, {:.2f} to {:.2f}'.format(
        len(time_lift), time_lift[0], time_lift[-1], lift_vertical.min(), lift_vertical.max()))

  # 2 つの図は (L/D)v = (L/D)_R cos(bank) で結ばれる。読み取りが互いに整合するか見る
  bank_on_lift = np.interp(time_lift, time_bank, bank)
  model = args.lift_drag*np.cos(np.radians(bank_on_lift))
  window = (time_lift > 20.0) & (time_lift < 560.0)
  difference = model[window] - lift_vertical[window]
  print('Cross-check, (L/D)_R cos(bank) against figure 16(a): mean {:+.3f}, rms {:.3f}'.format(
        np.mean(difference), np.sqrt(np.mean(difference**2))))

  write_reference(os.path.join(args.output_directory, 'apollo4_bank.dat'),
                  time_bank, bank, 'bank angle from the entry-control history',
                  'Bank[deg.]',
                  ['Figure 6(b). Time is seconds from the entry interface (GET 29 968.54 s).',
                   '0 deg. is lift up, positive towards the right of the flight direction.'])
  write_reference(os.path.join(args.output_directory, 'apollo4_liftvertical.dat'),
                  time_lift, lift_vertical, 'vertical lift-to-drag ratio derived from the flight',
                  'LiftDragVertical',
                  ['Figure 16(a). Time is seconds from the entry interface (GET 29 968.54 s).',
                   'Positive is lift up. This is what drives config_lift.yml.'])

  # config に貼る表。**鉛直揚抗比から**バンク角を作る
  angle = np.degrees(np.arccos(np.clip(lift_vertical/args.lift_drag, -1.0, 1.0)))
  time_node = np.arange(0.0, time_lift[-1] + 0.5*args.bank_interval, args.bank_interval)
  angle_node = np.interp(time_node, time_lift, angle)
  print('')
  print('  # Bank angle [deg.] against the time from the entry interface [s],')
  print('  # from the flight-derived vertical lift-to-drag ratio:')
  print('  #   bank = arccos( (L/D)v / {:.3f} )'.format(args.lift_drag))
  print('  bank_angle_table:')
  for index in range(0, len(time_node)):
    print('    - [{:6.1f}, {:8.2f}]'.format(time_node[index], angle_node[index]))

  if args.check is not None :
    draw_check(args.check, image_bank, time_bank_get, bank, image_lift, time_lift_get, lift_vertical)

  return


def draw_check(filename, image_bank, time_bank, bank, image_lift, time_lift, lift_vertical):
  try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
  except ImportError:
    print('matplotlib is not installed, so no check figure is written.')
    return

  figure, axis = plt.subplots(2, 1, figsize=(14, 11))

  calibration = CALIBRATION_BANK
  axis[0].imshow(image_bank[520:1820, 760:2900], cmap='gray', extent=[760, 2900, 1820, 520])
  axis[0].plot(calibration['column_reference'] + (time_bank - calibration['time_reference'])*calibration['pixel_per_second'],
               calibration['row_zero'] - bank*calibration['pixel_per_value'], 'r.', markersize=2.5)
  axis[0].set_xticks([calibration['column_reference'] + (item - calibration['time_reference'])*calibration['pixel_per_second']
                      for item in range(29960, 30681, 80)])
  axis[0].set_xticklabels(range(29960, 30681, 80), rotation=90, fontsize=7)
  axis[0].set_yticks([calibration['row_zero'] - item*calibration['pixel_per_value'] for item in range(-160, 161, 40)])
  axis[0].set_yticklabels(range(-160, 161, 40), fontsize=7)
  axis[0].set_title('figure 6(b): bank angle, deg.')

  calibration = CALIBRATION_LIFT
  axis[1].imshow(image_lift[1300:3750, 1450:5450], cmap='gray', extent=[1450, 5450, 3750, 1300])
  axis[1].plot(calibration['column_reference'] + (time_lift - calibration['time_reference'])*calibration['pixel_per_second'],
               calibration['row_zero'] - lift_vertical*calibration['pixel_per_value'], 'r.', markersize=2)
  axis[1].set_xticks([calibration['column_reference'] + (item - calibration['time_reference'])*calibration['pixel_per_second']
                      for item in range(29960, 30681, 80)])
  axis[1].set_xticklabels(range(29960, 30681, 80), rotation=90, fontsize=7)
  axis[1].set_yticks([calibration['row_zero'] - item*calibration['pixel_per_value'] for item in (-0.4, -0.2, 0.0, 0.2, 0.4)])
  axis[1].set_yticklabels((-0.4, -0.2, 0.0, 0.2, 0.4), fontsize=7)
  axis[1].set_title('figure 16(a): vertical lift-to-drag ratio')

  plt.tight_layout()
  plt.savefig(filename, dpi=120)
  plt.close()
  print('Check figure written:', filename)

  return


if __name__ == '__main__':
  main()
