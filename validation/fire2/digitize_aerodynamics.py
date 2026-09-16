#!/usr/bin/env python3
#
# Project Fire flight II の空力係数を NASA TN D-4183 の図 4 から読み取り、
# Tacode の係数表 database/aerodynamic/aerodynamic_fire2.txt を書く。
#
# 図 4 は 3 枚組で、迎角 0〜80 度に対する
#
#   Cx  軸力係数（図では負に振ってある。alpha = 0 で -1.51）
#   CR  Yb-Zb 面内の合力係数（＝法線力係数）
#   Cm  ピッチングモーメント係数
#
# を与える。**これは解析モデルではなく実測**で、Ames の極超音速自由飛行装置で
# 「basic Apollo body shape」を M = 35 で撃ったもの（D-4183 の本文）。
# D-4183 自身が飛行データから迎角を出すのにこの図を使っているので、
# **Tacode に入れる係数と、実測の迎角を出した係数が同じもの**になる。
#
# Tacode の規約へ:
#
#   F_body = -q S [CFx, CFy, CFz]   ->  CFx = |Cx|（alpha = 0 で CD）, CFz = CR
#   M_body =  q S L [CMx, CMy, CMz] ->  CMy = Cm（静安定なので alpha > 0 で負）
#
# 基準量は S = pi d^2 / 4、L = d（報告書の Cm の定義に合わせる）。
# **モーメントの基準点は図 4 の基準点 x_mc** で、表 I の脚注が
# `x_mc = x_r - r + 0.2875 d`（Fire Station）と与えている。重心はそこから
# ずれているので、config の `attitude.center_of_gravity` で移す（README 参照）。
#
# **Knudsen 依存は入っていない。**図 4 は単一の試験条件（M = 35、連続流）なので、
# 表は全 Kn に同じ値を書く。希薄の効きは README に見積もりを書いてある。
#
# 必要なもの: poppler の pdftoppm、Pillow、numpy。**参照 PDF はリポジトリの外**に
# あるので、読み取り済みの reference/fire2_aerodynamics.dat をリポジトリに置いてある。
#
# 使い方:
#   python3 digitize_aerodynamics.py                     # reference/ と database/ を作り直す
#   python3 digitize_aerodynamics.py --check check.png   # 読み取りを図に重ねた確認図

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
FILE_REFERENCE = os.path.join(DIRECTORY_SCRIPT, 'reference/fire2_aerodynamics.dat')

# 係数表はケースの database/ と**リポジトリのマスター**の両方に書く。
# ケースが自分の使うテーブルを持ち運ぶのがこのリポジトリの方針で、両者が
# バイト一致することを test/test_database_path.py が md5 で検査する
FILE_TABLE_LIST = [os.path.join(DIRECTORY_SCRIPT, 'database/aerodynamic/aerodynamic_fire2.txt'),
                   os.path.join(DIRECTORY_SCRIPT, '../../database/aerodynamic/aerodynamic_fire2.txt')]

PAGE_FIGURE = 23          # pdftoppm のページ番号（本文 21 ページ）
RESOLUTION  = 600         # dpi

# 図 4 の較正（600 dpi の画素）。
#
# **横軸の原点は枠の左端ではない。**縦罫線は 5 度ごとに 16 本あり、左端の 1 本は
# alpha = 0 ではなく **alpha = 5 度**で、alpha = 0 はその 144 画素左にある。
# 罫線の間隔から求めた alpha = 80 度の位置（3698.5）は、PDF の文字層にある
# 目盛りラベル "80" の中心（3699）と 1 画素で一致する。ラベル "0" は 1382 で
# 12 画素ずれるが、これはラベルの中心決めの誤差で、罫線のほうが確か。
#
# ここを枠の左端と取り違えると、係数がまるごと 5 度ずれる（実際に一度やった。
# CD が 1.47 に、微小迎角の Cm の傾きが 2 倍近くに出た）。
COLUMN_ANGLE_ZERO    = 1394.5     # alpha = 0
COLUMN_ANGLE_MAXIMUM = 3698.5     # alpha = 80 deg.
ANGLE_MAXIMUM        = 80.0

# 縦罫線は 5 度ごと、alpha = 5 度から 80 度まで 16 本
GRIDLINE_ANGLE = [5.0*(index + 1) for index in range(16)]

# パネルごとの (枠の上の行, 枠の下の行, 上の値, 下の値)
PANEL = {'CX': (614.0, 1723.0,  0.00, -1.60),
         'CR': (1854.0, 3237.0, 0.50,  0.00),
         'CM': (3368.0, 4749.0, 0.02, -0.08)}

# 縦罫線はこの幅で避ける（画素）
WIDTH_GRIDLINE = 8

# 横罫線は、パネルの幅のこれ以上にわたって黒い行として消す。曲線が平らな区間より
# 長く、罫線より短い値にしてある（Cx の alpha < 10 度はおよそ 12 %）
FRACTION_GRIDLINE = 0.40

# 追跡の刻みと、予測から許す外れ（画素）
STEP_COLUMN   = 6
WINDOW_FIT    = 8
TOLERANCE_ROW = 40.0

# 出力する迎角（度）。表は全迎角で引かれるので 0〜180 度を張る必要がある
ANGLE_OUTPUT = list(range(0, 185, 5))

# Knudsen の格子（既存の表と同じ）
KNUDSEN_LIST = [1.e-4, 1.e-3, 1.e-2, 1.e-1, 1.e0, 1.e1, 1.e2, 1.e3, 1.e4, 1.e5]

# 機体の基準量（NASA TN D-4183 表 I、complete reentry package）
DIAMETER_REFERENCE = 0.672


def render_page(file_report, directory_work):
  # PDF の 1 ページを PNG にする
  if not os.path.exists(file_report) :
    print('No such report:', file_report)
    print('--The scanned report lives outside the repository. Give its path with')
    print('  --report, or keep the reference/ file already in the case.')
    print('Program stopped.')
    sys.exit(1)

  prefix  = os.path.join(directory_work, 'page')
  command = ['pdftoppm', '-f', str(PAGE_FIGURE), '-l', str(PAGE_FIGURE),
             '-r', str(RESOLUTION), '-png', file_report, prefix]
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


def load_image(path):
  try:
    from PIL import Image
  except ImportError:
    print('Pillow is not installed; it is needed to read the figures.')
    print('--Install it with: pip install pillow')
    print('Program stopped.')
    sys.exit(1)
  return np.array(Image.open(path).convert('L'))


def clean_mask(image, row_top, row_bottom):
  #
  # パネルを切り出して黒画素の真偽値にし、罫線を消す。
  # 縦罫線は 5 度ごとに全高を走り、横罫線は幅の大半を走る。どちらも曲線より長い
  #
  mask = image[int(row_top):int(row_bottom) + 1,
               int(COLUMN_ANGLE_ZERO):int(COLUMN_ANGLE_MAXIMUM) + 1] < 128

  width = mask.shape[1]
  for row in range(mask.shape[0]):
    if mask[row, :].sum() > FRACTION_GRIDLINE*width :
      mask[row, :] = False

  for angle in GRIDLINE_ANGLE:
    column = int(round(angle/ANGLE_MAXIMUM*(COLUMN_ANGLE_MAXIMUM - COLUMN_ANGLE_ZERO)))
    mask[:, max(0, column - WIDTH_GRIDLINE):column + WIDTH_GRIDLINE + 1] = False

  return mask


def cluster_column(mask, column):
  # 1 列の黒画素を、途切れ 10 px 以内でまとめた塊の (中心, 画素数) にする
  row_list = np.where(mask[:, column])[0]
  if len(row_list) == 0 :
    return []
  group = []
  for row in row_list:
    if group and row - group[-1][-1] <= 10 :
      group[-1].append(row)
    else:
      group.append([row])
  return [(float(np.mean(item)), len(item)) for item in group]


def track_curve(mask, column_start, row_start):
  #
  # 曲線を左から右へたどる。直前 WINDOW_FIT 点に直線を当てて次の列を予測し、
  # 予測に最も近い塊を採る。TOLERANCE_ROW 以上離れていれば見送る（罫線を消した
  # 隙間や、ほかの線に乗り移らないようにするため）
  #
  column_list = [column_start]
  row_list    = [row_start]
  column      = column_start
  while column + STEP_COLUMN < mask.shape[1]:
    column += STEP_COLUMN
    candidate = cluster_column(mask, column)
    if len(column_list) >= 2 :
      fit = np.polyfit(column_list[-WINDOW_FIT:], row_list[-WINDOW_FIT:], 1)
      predicted = np.polyval(fit, column)
    else:
      predicted = row_list[-1]
    if not candidate :
      continue
    row = min(candidate, key=lambda item: abs(item[0] - predicted))[0]
    if abs(row - predicted) > TOLERANCE_ROW :
      continue
    column_list.append(column)
    row_list.append(row)
  return np.array(column_list, dtype=float), np.array(row_list, dtype=float)


def find_seed(mask):
  #
  # 追跡を始める列。**塊がちょうど 1 つしか無い最初の列**にする。
  # 縦軸のすぐ右には目盛りの数字の滲みが残っていて、そこで始めると数字を
  # 曲線と取り違える。0 度のすぐ右では曲線しか無い列がすぐ見つかる
  #
  for column in range(WIDTH_GRIDLINE + 2, mask.shape[1]):
    candidate = [item for item in cluster_column(mask, column) if item[1] >= 4]
    if len(candidate) == 1 :
      return column, candidate[0][0]
  print('No single curve found anywhere in the panel.')
  print('Program stopped.')
  sys.exit(1)


def digitize_panel(image, name, value_origin):
  row_top, row_bottom, value_top, value_bottom = PANEL[name]
  mask = clean_mask(image, row_top, row_bottom)

  column_start, row_start = find_seed(mask)
  column_list, row_list   = track_curve(mask, column_start, row_start)

  angle = column_list/(COLUMN_ANGLE_MAXIMUM - COLUMN_ANGLE_ZERO)*ANGLE_MAXIMUM
  value = value_top + (value_bottom - value_top)*row_list/(row_bottom - row_top)

  # alpha = 0 の値は**対称性で決まっている**ので、読み取りではなくそちらを使う。
  # 軸対称の機体は迎角 0 で法線力もピッチングモーメントも持たず（value_origin = 0）、
  # 軸力は迎角の偶関数なので傾きが 0 になる（value_origin = None で端の値を延ばす）。
  # 図の左端は縦軸や Cm = 0 の罫線と重なっていて、そこだけ読み取りが効かない
  angle = np.concatenate([[0.0], angle])
  value = np.concatenate([[value[0] if value_origin is None else value_origin], value])

  return angle, value, (column_list, row_list, row_top)


def resample(angle, value, angle_output):
  # 読み取り点を出力の迎角へ線形内挿する。範囲外は端の値でクランプ
  return np.interp(angle_output, angle, value)


def write_reference(file_output, angle, coefficient):
  directory = os.path.dirname(file_output)
  if directory and not os.path.isdir(directory) :
    os.makedirs(directory)
  with open(file_output, 'w') as stream:
    stream.write('# Project Fire flight II static aerodynamic coefficients\n'
                 '#\n'
                 '# NASA TN D-4183 (Scallion and Lewis, 1967), figure 4:\n'
                 '#   "Static longitudinal aerodynamic data used in the\n'
                 '#    angle-of-attack analysis".\n'
                 '#\n'
                 '# Measured in the Ames hypersonic free-flight facility at\n'
                 '# Mach 35 on a basic Apollo body shape. Digitised from the\n'
                 '# scan by digitize_aerodynamics.py.\n'
                 '#\n'
                 '# The moment is referenced to the station x_mc of table I,\n'
                 '# not to the centre of gravity.\n'
                 '#\n'
                 '# alpha  deg.\n'
                 '# CX     axial force, as plotted (negative)\n'
                 '# CR     resultant force in the Yb-Zb plane\n'
                 '# CM     pitching moment, q S d\n'
                 '#{:>11s} {:>15s} {:>15s} {:>15s}\n'.format('alpha', 'CX', 'CR', 'CM'))
    for index in range(len(angle)):
      stream.write(' {:11.2f} {:15.6f} {:15.6f} {:15.6f}\n'.format(
                   angle[index], coefficient['CX'][index],
                   coefficient['CR'][index], coefficient['CM'][index]))


def write_table(file_output, angle, coefficient):
  #
  # Tacode の係数表。**図 4 は 80 度までしか無い**ので、そこから先は
  # 端の値でクランプする（この機体は静安定で、実際の迎角は 20 度を超えない。
  # 表の残りは内挿器が範囲外にならないようにするためだけのもの）
  #
  directory = os.path.dirname(file_output)
  if directory and not os.path.isdir(directory) :
    os.makedirs(directory)

  with open(file_output, 'w') as stream:
    stream.write('Project Fire II aerodynamic data (NASA TN D-4183 figure 4, Mach 35, measured)\n')
    # Kn 軸を作った代表長さ（読み取り側が config の characteristic_length と突き合わせる）
    stream.write('# Reference length: {:g} m\n'.format(DIAMETER_REFERENCE))
    stream.write('variables = Kn, CFx, CFy, CFz, CMx, CMy, CMz, '
                 'SDV_CFx, SDV_CFy, SDV_CFz, SDV_CMx, SDV_CMy, SDV_CMz, Altitude \n')
    for index, angle_value in enumerate(angle):
      stream.write('AOA {:g}\n'.format(angle_value))
      value_x = abs(coefficient['CX'][index])
      value_z = coefficient['CR'][index]
      value_m = coefficient['CM'][index]
      for knudsen in KNUDSEN_LIST:
        stream.write('\t'.join(['{:.18e}'.format(item) for item in
                                [knudsen, value_x, 0.0, value_z,
                                 0.0, value_m, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]) + '\n')


def write_check(file_check, path_image, trace):
  try:
    from PIL import Image, ImageDraw
  except ImportError:
    print('Pillow is needed for --check.')
    return
  image = Image.open(path_image).convert('RGB')
  draw  = ImageDraw.Draw(image)
  for name in trace:
    column_list, row_list, row_top = trace[name]
    for index in range(len(column_list)):
      x = COLUMN_ANGLE_ZERO + column_list[index]
      y = row_top + row_list[index]
      draw.ellipse([x - 5, y - 5, x + 5, y + 5], outline=(255, 0, 0), width=3)
  image.save(file_check)
  print('Wrote ' + file_check)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--report', type=str, default=FILE_REPORT_DEFAULT,
                      help='Scanned NASA TN D-4183')
  parser.add_argument('--check', type=str, default=None,
                      help='Write the digitised points back onto the figure')
  args = parser.parse_args()

  directory_work = tempfile.mkdtemp()
  try:
    path_image = render_page(args.report, directory_work)
    image      = load_image(path_image)

    trace       = {}
    coefficient = {}
    for name, value_origin in (('CX', None), ('CR', 0.0), ('CM', 0.0)):
      angle, value, mark = digitize_panel(image, name, value_origin)
      trace[name]       = mark
      coefficient[name] = resample(angle, value, ANGLE_OUTPUT)
      print('{:s}: {:d} points read, alpha = {:.1f} to {:.1f} deg., '
            'value = {:.4f} to {:.4f}'.format(
            name, len(angle), angle[0], angle[-1], value[0], value[-1]))

    if args.check is not None :
      write_check(args.check, path_image, trace)
  finally:
    for name in os.listdir(directory_work):
      os.unlink(os.path.join(directory_work, name))
    os.rmdir(directory_work)

  angle_output = np.array(ANGLE_OUTPUT, dtype=float)
  write_reference(FILE_REFERENCE, angle_output, coefficient)
  print('Wrote ' + os.path.relpath(FILE_REFERENCE))
  for file_table in FILE_TABLE_LIST:
    write_table(file_table, angle_output, coefficient)
    print('Wrote ' + os.path.relpath(file_table))

  print('--CD at alpha = 0        : {:.4f}'.format(abs(coefficient['CX'][0])))
  slope = (coefficient['CM'][2] - coefficient['CM'][0])/(10.0*np.pi/180.0)
  print('--dCm/dalpha at alpha = 0: {:.4f} per radian'.format(slope))
  print('--Reference length       : {:.3f} m'.format(DIAMETER_REFERENCE))

  return


if __name__ == '__main__':
  main()
