#!/usr/bin/env python3
#
# Project Fire flight II の飛行データを NASA TN D-4183 の図から読み取る。
#
#   図 17 TN D-3569 の軌道の動圧（破線）と、**報告書自身の 6 自由度シミュレーションが
#         実測のピッチ・ヨーレートの周波数に合わせるために必要とした動圧**（実線の
#         短い区間）。周波数 omega は omega^2 ∝ -Cma q S d / Iy なので、この実線は
#         **実測の振動周波数を動圧に翻訳したもの**であり、Cma・S・d・Iy が同じなら
#         Tacode の動圧と直接比べられる
#   図 19 全迎角の包絡線（振動の最大と最小）。擾乱ごとに階段状に上がる
#
# 読み取らないと決めた図が 3 つある。**読めなかったのではなく、読んでも分解能が
# 足りないか、公平な標本にならないから**である。
#
#   図 8（ロールレート）  commutate した点のばらつきが大きい。報告書が本文に数値を
#         書いている（分離直後 171 rpm = 17.9 rad/s、t = 1658 s で 13.6 rad/s、
#         1662 s で 6.6 rad/s）ので、そちらを README に引く
#   図 9（機上加速度計）  白抜きの丸が重なっていて、連結成分に分けても 200 個ほどの
#         うち 30 個しか切り離せない。しかも切り離せるのは**記号が疎な区間だけ**なので、
#         30 個を並べると疎な区間に偏った標本になる。同じ物理（q S CD / m）は
#         表 V との突き合わせがもっと細かく見ている
#   TN D-3569 の図 17（レーダ実測）  高度・速度・経路角の 3 本と 4 種類の追尾装置の
#         記号が重なっていて、読み取り分解能は高度で 1.5 km ほど。表 V との一致が
#         0.12 km なので、10 倍粗い物差しでは何も判別できない
#
# 較正はすべて**図そのものの罫線と、PDF の文字層にある目盛りラベル**から取っており、
# 目分量はしていない。罫線の間隔が目盛り間隔と合うことを数で確かめてある。
#
# 必要なもの: poppler の pdftoppm、Pillow、numpy。**参照 PDF はリポジトリの外**に
# あるので、読み取り済みの reference/*.dat をリポジトリに置いてある。
#
# 使い方:
#   python3 digitize_figure.py                    # reference/ を作り直す
#   python3 digitize_figure.py --check check      # 読み取りを図に重ねた確認図

import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np


DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
FILE_REPORT_DEFAULT = os.path.join(
  DIRECTORY_SCRIPT, '../../../references/pdf',
  'NASA-TN-D-4183_Scallion-Lewis_1967_Fire-II-Body-Motions.pdf')
FILE_TRAJECTORY = os.path.join(DIRECTORY_SCRIPT, 'reference/fire2_trajectory.dat')

RESOLUTION = 600          # dpi

STEP_COLUMN    = 4        # 読み取る列の間隔（画素）
SIZE_CLUSTER   = 4        # これより小さい黒画素の塊は染みとみなす

# 図ごとの較正（600 dpi の画素）。
#
#   page      pdftoppm のページ番号
#   x_time    (画素, 飛行経過時刻) の 2 点。横軸の目盛りラベルの中心から
#   y_value   (画素, 値) の 2 点。**横罫線の本数が目盛り間隔と合うことを確かめてある**
#   grid_row     横罫線 (基準行, 間隔, 消す半幅, 列の範囲) の一覧。曲線と紛れるので消す
#   grid_column  縦罫線 (基準列, 間隔, 消す半幅, 行の範囲) の一覧
#   box       読み取る範囲 (row_min, row_max, column_min, column_max)
#   cluster_minimum  これより小さい塊は捨てる
#   grid_reject  (基準行, 間隔, 距離, 太さ)。罫線の位置にある細い塊を捨てる
#   exclude   読み取りから外す矩形の一覧（凡例など）
FIGURE = {
  'dynamic_pressure': {
    'page': 81,
    # 横軸は**目盛りラベルではなく縦罫線**から取った。罫線は 2 秒ごと 16 本あり、
    # ラベルの中心から取ると 0.3 秒ずれる（走査の文字の中心決めの誤差）
    'x_time':  ((1792.0, 1640.0), (5015.5, 1670.0)),
    'y_value': ((3600.0, 0.0), (477.5, 140000.0)),
    # 横罫線は 10 000 N/m2 ごと 15 本。走査が傾いていて左右で 11 画素ずれるので
    # 消す幅を広めに取る（曲線が罫線を横切るのは急な区間だけなので損は小さい）
    'grid_row': [(477.5, 223.04, 14, 1700, 5100)],
    'grid_column': [(1792.0, 214.9, 13, 420, 3640)],
    'grid_reject': None,
    'cluster_minimum': SIZE_CLUSTER,
    'box': (420, 3640, 1800, 5020),
    'exclude': [(950, 1350, 3700, 5100)],        # 凡例
  },
  'angle_attack': {
    'page': 83,
    'x_time':  ((1063.0, 1639.0), (5292.0, 1675.0)),
    'y_value': ((3328.0, 0.0), (869.0, 20.0)),
    # 横罫線は 1 度ごと 21 本だが、**消さない**。曲線が 0.8 度や 19.5 度という
    # 罫線のすぐそばを長く走るので、帯ごと消すと曲線まで持っていかれる。
    # 代わりに太さで分ける: 罫線の塊は 1 列あたり 5〜7 画素、曲線は 11 画素以上ある
    # 横罫線は 1 度ごと 21 本だが、**帯ごと消してはいけない**。曲線が 0.8 度や
    # 19.5 度という罫線のすぐそばを長く走るので、消すと曲線まで持っていかれる。
    # 代わりに「罫線の位置にある細いインク」だけを捨てる（grid_reject）。
    # 罫線の塊は 1 列あたり 5〜7 画素、曲線は 11 画素以上、曲線が罫線に重なれば
    # 16 画素以上になるので、位置と太さの両方で分けられる
    'grid_row': [],
    'grid_reject': (869.0, 122.95, 5, 10),
    # 目盛りの短い縦線が 1 秒ごとに下端から 3.3 度あたりまで伸びていて、
    # 下側の曲線と繋がってしまう。そこだけ消す
    'grid_column': [(1063.0, 117.47, 7, 2840, 3310)],
    'cluster_minimum': 5,
    # 下端は 0 度の枠線の手前で切る（枠線を曲線と取り違えないため）
    'box': (850, 3310, 1070, 5290),
    'exclude': [],
  },
}

# 図 17 で、破線（表 V の動圧）とみなす予測からのずれ（画素）。これより離れた
# インクを実線の区間（実測の周波数）とみなす
TOLERANCE_DASHED = 30.0

# 図 17 の実線は短い区間に分かれている。この幅より狭い塊は染みとして捨てる
WIDTH_SEGMENT_MINIMUM = 60.0

def render_page(file_report, page, directory_work):
  if not os.path.exists(file_report) :
    print('No such report:', file_report)
    print('--The scanned report lives outside the repository. Give its path with')
    print('  --report, or keep the reference/ files already in the case.')
    print('Program stopped.')
    sys.exit(1)

  prefix  = os.path.join(directory_work, 'page{:d}'.format(page))
  command = ['pdftoppm', '-f', str(page), '-l', str(page), '-r', str(RESOLUTION),
             '-png', file_report, prefix]
  try:
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
  except (OSError, subprocess.CalledProcessError) as instance:
    print('Cannot render the report with pdftoppm:', instance)
    print('--Install poppler-utils, or keep the reference files already in reference/.')
    print('Program stopped.')
    sys.exit(1)

  for name in sorted(os.listdir(directory_work)):
    if name.startswith('page{:d}'.format(page)) and name.endswith('.png') :
      return os.path.join(directory_work, name)

  print('pdftoppm wrote no page.')
  print('Program stopped.')
  sys.exit(1)


def load_image(path):
  try:
    from PIL import Image
  except ImportError:
    print('Pillow is not installed; it is needed to read the figures.')
    print('--Install it with: pip install pillow')
    print('Program stopped.')
    sys.exit(1)
  return np.array(Image.open(path).convert('L'))


def linear_map(pair):
  # (画素, 値) の 2 点から画素 -> 値 の直線を作る
  (pixel_a, value_a), (pixel_b, value_b) = pair
  slope = (value_b - value_a)/(pixel_b - pixel_a)
  return lambda pixel: value_a + slope*(pixel - pixel_a)


def inverse_map(pair):
  (pixel_a, value_a), (pixel_b, value_b) = pair
  slope = (pixel_b - pixel_a)/(value_b - value_a)
  return lambda value: pixel_a + slope*(value - value_a)


def build_mask(image, setting):
  # 読み取る範囲を切り出し、横罫線と除外矩形を消す
  row_min, row_max, column_min, column_max = setting['box']
  mask = image < 128

  keep = np.zeros(mask.shape, dtype=bool)
  keep[row_min:row_max, column_min:column_max] = True
  mask = mask & keep

  for row_a, row_b, column_a, column_b in setting['exclude']:
    mask[row_a:row_b, column_a:column_b] = False

  for base, spacing, half, column_a, column_b in setting['grid_row']:
    index = 0
    while True:
      row = int(round(base + index*spacing))
      if row > row_max :
        break
      if row >= row_min :
        mask[max(0, row - half):row + half + 1, column_a:column_b] = False
      index += 1

  for base, spacing, half, row_a, row_b in setting['grid_column']:
    index = 0
    while True:
      column = int(round(base + index*spacing))
      if column > column_max :
        break
      if column >= column_min :
        mask[row_a:row_b, max(0, column - half):column + half + 1] = False
      index += 1

  return mask


def cluster_column(mask, column, minimum=SIZE_CLUSTER):
  # 1 列の黒画素を、途切れ 8 画素以内でまとめた (中心, 画素数) にする
  row_list = np.where(mask[:, column])[0]
  if len(row_list) == 0 :
    return []
  group = []
  for row in row_list:
    if group and row - group[-1][-1] <= 8 :
      group[-1].append(row)
    else:
      group.append([row])
  return [(float(np.mean(item)), len(item)) for item in group
          if len(item) >= minimum]


def reject_gridline(candidate, setting):
  # 罫線の位置にある細い塊を落とす。曲線が罫線に重なった塊は太いので残る
  if setting['grid_reject'] is None :
    return candidate
  base, spacing, distance, thickness = setting['grid_reject']
  keep = []
  for row, size in candidate:
    index = round((row - base)/spacing)
    if abs(row - (base + index*spacing)) < distance and size < thickness :
      continue
    keep.append((row, size))
  return keep


def read_envelope(mask, setting):
  #
  # 各列の**一番上と一番下**の塊を採る。図 19 の 2 本は交わらないので、
  # 曲線を 1 本ずつ追いかけるより頑丈（1648 s の跳びも素通りできる）
  #
  to_time  = linear_map(setting['x_time'])
  to_value = linear_map(setting['y_value'])
  row_min, row_max, column_min, column_max = setting['box']

  time_list, upper_list, lower_list = [], [], []
  for column in range(column_min, column_max, STEP_COLUMN):
    candidate = reject_gridline(
                  cluster_column(mask, column, setting['cluster_minimum']), setting)
    # 2 本とも見えている列だけを採る。片方が罫線や目盛りに呑まれた列を混ぜると
    # 上下の包絡線が入れ替わる
    if len(candidate) < 2 :
      continue
    time_list.append(to_time(column))
    upper_list.append(to_value(min(item[0] for item in candidate)))
    lower_list.append(to_value(max(item[0] for item in candidate)))
  return np.array(time_list), np.array(upper_list), np.array(lower_list)


def read_trajectory_dynamic_pressure():
  # 表 V の動圧。図 17 の破線はこれそのものなので、破線を見分けるのに使う
  time, value = [], []
  with open(FILE_TRAJECTORY) as stream:
    for line in stream:
      if line.startswith('#') :
        continue
      item = [float(entry) for entry in line.split()]
      if np.isfinite(item[7]) :
        time.append(item[0])
        value.append(item[7])
  return np.array(time), np.array(value)


def read_solid_segments(mask, setting):
  #
  # 図 17 の**実線の短い区間**だけを拾う。
  #
  # 破線は表 V の動圧そのものなので、その画素上の位置は reference から**予測できる**。
  # 予測から TOLERANCE_DASHED 以内のインクを破線として除き、残ったものを実線とみなす。
  # 図から破線を読み直して実線と見分けるより確かで、較正の確認にもなる
  #
  to_time    = linear_map(setting['x_time'])
  to_value   = linear_map(setting['y_value'])
  from_value = inverse_map(setting['y_value'])
  row_min, row_max, column_min, column_max = setting['box']

  time_table, value_table = read_trajectory_dynamic_pressure()

  time_dashed, value_dashed = [], []
  point = []
  for column in range(column_min, column_max, STEP_COLUMN):
    time      = to_time(column)
    predicted = from_value(np.interp(time, time_table, value_table))
    for row, dummy in cluster_column(mask, column, setting['cluster_minimum']):
      if abs(row - predicted) <= TOLERANCE_DASHED :
        time_dashed.append(time)
        value_dashed.append(to_value(row))
      else:
        point.append((column, time, to_value(row)))

  # 残ったインクを、列が続いている塊にまとめる
  segment, current = [], []
  for column, time, value in point:
    if current and column - current[-1][0] <= 3*STEP_COLUMN :
      current.append((column, time, value))
    else:
      if current :
        segment.append(current)
      current = [(column, time, value)]
  if current :
    segment.append(current)

  segment = [item for item in segment
             if item[-1][0] - item[0][0] >= WIDTH_SEGMENT_MINIMUM]

  return (np.array(time_dashed), np.array(value_dashed),
          [(np.array([entry[1] for entry in item]),
            np.array([entry[2] for entry in item])) for item in segment])


def write_columns(path, header, column_list):
  directory = os.path.dirname(path)
  if directory and not os.path.isdir(directory) :
    os.makedirs(directory)
  name_list, data_list = zip(*column_list)
  with open(path, 'w') as stream:
    stream.write(header)
    stream.write('#' + ' '.join('{:>15s}'.format(name) for name in name_list) + '\n')
    for index in range(len(data_list[0])):
      stream.write(' ' + ' '.join('{:15.6f}'.format(data[index])
                                  for data in data_list) + '\n')
  print('Wrote ' + os.path.relpath(path))


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--report', type=str, default=FILE_REPORT_DEFAULT,
                      help='Scanned NASA TN D-4183')
  parser.add_argument('--check', type=str, default=None,
                      help='Prefix for figures with the digitised points drawn on')
  args = parser.parse_args()

  os.chdir(DIRECTORY_SCRIPT)
  directory_work = tempfile.mkdtemp()
  try:
    mark = {}

    # 図 19: 全迎角の包絡線
    setting = FIGURE['angle_attack']
    image   = load_image(render_page(args.report, setting['page'], directory_work))
    mask    = build_mask(image, setting)
    time, upper, lower = read_envelope(mask, setting)
    print('Figure 19: {:d} columns read, t = {:.1f} to {:.1f} s, '
          'eta = {:.2f} to {:.2f} deg.'.format(
          len(time), time[0], time[-1], lower.min(), upper.max()))
    write_columns('reference/fire2_angle_attack.dat',
                  '# Project Fire flight II total angle-of-attack envelope\n'
                  '#\n'
                  '# NASA TN D-4183 (Scallion and Lewis, 1967), figure 19,\n'
                  '# digitised by digitize_figure.py.\n'
                  '#\n'
                  '# The maximum and the minimum of the oscillation, derived in the\n'
                  '# report from the rate-gyro records and the Mach 35 static\n'
                  '# aerodynamic data of its figure 4.\n'
                  '#\n'
                  '# time       s, elapsed flight time\n'
                  '# eta_max    deg., upper bound of the oscillation\n'
                  '# eta_min    deg., lower bound\n',
                  [('time', time), ('eta_max', upper), ('eta_min', lower)])
    mark['angle_attack'] = (setting, [(time, upper), (time, lower)])

    # 図 17: 実測の周波数が要求する動圧
    setting = FIGURE['dynamic_pressure']
    image   = load_image(render_page(args.report, setting['page'], directory_work))
    mask    = build_mask(image, setting)
    time_dashed, value_dashed, segment_list = read_solid_segments(mask, setting)
    time_table, value_table = read_trajectory_dynamic_pressure()
    residual = value_dashed - np.interp(time_dashed, time_table, value_table)
    print('Figure 17: the dashed curve is table V to {:.0f} N/m2 rms over {:d} points '
          '(a check of the calibration)'.format(
          np.sqrt(np.mean(residual**2)), len(time_dashed)))
    print('           {:d} solid segments found'.format(len(segment_list)))
    time_solid  = np.concatenate([item[0] for item in segment_list])
    value_solid = np.concatenate([item[1] for item in segment_list])
    order = np.argsort(time_solid)
    time_solid, value_solid = time_solid[order], value_solid[order]
    ratio = value_solid/np.interp(time_solid, time_table, value_table)
    print('           required / table V: {:.3f} to {:.3f}'.format(ratio.min(), ratio.max()))
    write_columns('reference/fire2_dynamic_pressure.dat',
                  '# Project Fire flight II: the dynamic pressure the measured\n'
                  '# pitch and yaw frequencies require\n'
                  '#\n'
                  '# NASA TN D-4183 (Scallion and Lewis, 1967), figure 17, the SOLID\n'
                  '# segments, digitised by digitize_figure.py.\n'
                  '#\n'
                  '# The report matched six windows of its own six-degree-of-freedom\n'
                  '# simulation to the measured rate-gyro traces by varying the\n'
                  '# dynamic pressure at the wind-tunnel value of Cm_alpha. Since the\n'
                  '# oscillation frequency obeys omega^2 = -Cm_alpha q S d / Iy at\n'
                  '# fixed spin and inertia, these values ARE the measured frequency,\n'
                  '# expressed as a dynamic pressure.\n'
                  '#\n'
                  '# time       s, elapsed flight time\n'
                  '# q_required N/m2\n'
                  '# q_tableV   N/m2, the dashed curve of the same figure\n',
                  [('time', time_solid), ('q_required', value_solid),
                   ('q_tableV', np.interp(time_solid, time_table, value_table))])
    mark['dynamic_pressure'] = (setting, [(item[0], item[1]) for item in segment_list])

    if args.check is not None :
      write_check(args.check, args.report, directory_work, mark)
  finally:
    for name in os.listdir(directory_work):
      os.unlink(os.path.join(directory_work, name))
    os.rmdir(directory_work)

  return


def write_check(prefix, file_report, directory_work, mark):
  try:
    from PIL import Image, ImageDraw
  except ImportError:
    print('Pillow is needed for --check.')
    return
  for name in mark:
    setting, series_list = mark[name]
    image = Image.open(render_page(file_report, setting['page'], directory_work)).convert('RGB')
    draw  = ImageDraw.Draw(image)
    from_time  = inverse_map(setting['x_time'])
    from_value = inverse_map(setting['y_value'])
    for time, value in series_list:
      for index in range(len(time)):
        x, y = from_time(time[index]), from_value(value[index])
        draw.ellipse([x - 6, y - 6, x + 6, y + 6], outline=(255, 0, 0), width=3)
    path = prefix + '_' + name + '.png'
    image.save(path)
    print('Wrote ' + path)


if __name__ == '__main__':
  main()
