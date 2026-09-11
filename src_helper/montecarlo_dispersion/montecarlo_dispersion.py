#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/10
#
# モンテカルロ（src/tacode-montecarlo.py）が作ったケース群の終端点をまとめる。
# 計算はしない。各ケースの出力（Tecplot）と制御ファイルを読むだけ。
#
#   python3 src_helper/montecarlo_dispersion/montecarlo_dispersion.py work_montecarlo_wind
#
# 設定はカレントディレクトリの config_helper.yml の montecarlo_dispersion セクションに
# 書ける（コマンドラインが優先。--save-config でいまの設定を書き出せる）。
#
# 各ケースの最後の点（再突入なら着地点）を集め、全ケースの平均位置からのずれを
# ローカル水平系（東・北）の km で出す。ばらつかせた入力は、モンテカルロが
# 残す case_template の制御ファイルとの差分から拾うので、風でも密度でも
# 初期速度でも同じ扱いになる（ツール側にモンテカルロの設定は持たない）。
#
# ずれの分解に使う東・北は**地心**のローカル水平系。ソルバーの入出力
# （initial_settings.velocity、風、Tecplot の Upl/Vpl/Wpl）と同じ規約である。

import argparse
import csv
import glob
import re
import os
import sys

import numpy as np
import yaml

# Tecplot の読み取りは animate_trajectory のものを使う（形式を二重に持たないため）
HELPER_DIR = os.path.dirname(os.path.abspath(__file__))
READER_DIR = os.path.normpath(os.path.join(HELPER_DIR, '..', 'animate_trajectory'))
GENERAL_DIR = os.path.normpath(os.path.join(HELPER_DIR, '..', 'general'))
for directory in (READER_DIR, GENERAL_DIR) :
  if directory not in sys.path :
    sys.path.insert(0, directory)

import tecplot_reader as tecplot_reader  # noqa: E402
import helper_config as helper_config    # noqa: E402

# 設定ファイルの中でこのツールが読むセクション
NAME_SECTION = 'montecarlo_dispersion'

# モンテカルロが作るディレクトリの名前（montecarlo.case_dir + 連番、および
# その隣に置かれるテンプレート）
SUFFIX_TEMPLATE = '_template'

# ケースの中の出力（config.yml の post_process.directory_output と tecplot.filename_output）
PATH_TECPLOT_DEFAULT = os.path.join('output_result', 'tecplot.dat')
NAME_CONTROL_DEFAULT = 'config.yml'


def find_case_directory(directory):
  #
  # <directory> の直下から、ケースディレクトリを名前順に集める。
  #
  # ケースの名前は montecarlo.get_case_directory が作る「<case_dir> + 4 桁」なので、
  # **末尾が 4 桁の数字のディレクトリだけ**を採る。「テンプレート以外の全部」だと、
  # 作業ディレクトリの中に置いた結果ディレクトリなどもケースとして数えてしまう。
  #
  path_list = sorted(glob.glob(os.path.join(directory, '*')))
  case_list = []
  for path in path_list :
    if not os.path.isdir(path) :
      continue
    name = os.path.basename(path)
    if name.endswith(SUFFIX_TEMPLATE) :
      continue
    if re.search(r'\d{4}$', name) is None :
      continue
    case_list.append(path)
  return case_list


def find_template_directory(directory):
  path_list = sorted(glob.glob(os.path.join(directory, '*'+SUFFIX_TEMPLATE)))
  for path in path_list :
    if os.path.isdir(path) :
      return path
  return None


def flatten_config(node, prefix=''):
  #
  # 入れ子の辞書・リストを 'wind.velocity[0]' のようなキーの平坦な辞書にする。
  # 制御ファイルどうしの差分を取るため。
  #
  flat = {}
  if isinstance(node, dict) :
    for key, value in node.items() :
      name = str(key) if prefix == '' else prefix+'.'+str(key)
      flat.update(flatten_config(value, name))
  elif isinstance(node, list) :
    for i, value in enumerate(node) :
      flat.update(flatten_config(value, prefix+'['+str(i)+']'))
  else :
    flat[prefix] = node
  return flat


def get_dispersed_variable(case_list, template_path, name_control):
  #
  # テンプレートの制御ファイルと各ケースのそれを比べ、値の違うキーを拾う。
  # 返すのは (キーの一覧, ケースごとの値の辞書)。
  #
  if template_path is None :
    return [], {}

  file_template = os.path.join(template_path, name_control)
  if not os.path.exists(file_template) :
    return [], {}

  with open(file_template) as f:
    flat_template = flatten_config(yaml.safe_load(f))

  name_list  = []
  value_dict = {}
  for case_path in case_list :
    file_case = os.path.join(case_path, name_control)
    if not os.path.exists(file_case) :
      continue
    with open(file_case) as f:
      flat_case = flatten_config(yaml.safe_load(f))
    value_dict[case_path] = {}
    for name, value in flat_case.items() :
      if name not in flat_template or flat_template[name] == value :
        continue
      value_dict[case_path][name] = value
      if name not in name_list :
        name_list.append(name)

  return sorted(name_list), value_dict


def get_terminal_state(filename):
  #
  # Tecplot 出力の最後の行（軌道の終端。再突入なら着地点）を辞書で返す。
  #
  data = tecplot_reader.read_tecplot(filename)
  state = {}
  for name in data :
    state[name] = float(data[name][-1])
  return state


def get_local_horizon(position):
  #
  # 地心のローカル水平系の単位ベクトル [東, 北, 上] を、ECEF の位置から作る。
  #
  radius_xy = np.hypot(position[0], position[1])
  if radius_xy == 0.0 :
    # 極の上。東西が定義できないので x 軸を東にとる（実用上は起きない）
    unit_east = np.array([1.0, 0.0, 0.0])
  else :
    unit_east = np.array([-position[1], position[0], 0.0])/radius_xy
  unit_up    = position/np.linalg.norm(position)
  unit_north = np.cross(unit_up, unit_east)
  return unit_east, unit_north, unit_up


def get_dispersion(state_list, position_reference=None):
  #
  # 終端点のずれを、基準点のローカル水平系で東・北・上に分解する（km）。
  # 基準点は既定では全ケースの平均。--reference を与えたときはそのケースの終端点
  # （たとえば無風のケース）で、そのときは平均のずれ自体が風の効きになる。
  #
  position = np.array([[state['X'], state['Y'], state['Z']] for state in state_list])
  position_mean = np.mean(position, axis=0)

  position_origin = position_mean if position_reference is None else position_reference

  unit_east, unit_north, unit_up = get_local_horizon(position_origin)

  offset = position - position_origin
  east   = offset.dot(unit_east)
  north  = offset.dot(unit_north)
  up     = offset.dot(unit_up)

  return position_origin, east, north, up


def report(case_list, state_list, name_list, value_dict, position_origin, east, north, up,
           name_reference=None):

  horizontal = np.hypot(east, north)

  header = ['Case', 'Time[s]', 'Long[deg.]', 'Lati[deg.]', 'Alti[km]',
            'East[km]', 'North[km]', 'Range[km]'] + name_list
  row_list = []
  for i, case_path in enumerate(case_list) :
    state = state_list[i]
    row = [os.path.basename(case_path),
           '{:.1f}'.format(state['Time']),
           '{:.5f}'.format(state['Long']),
           '{:.5f}'.format(state['Lati']),
           '{:.4f}'.format(state['Alti']),
           '{:+.3f}'.format(east[i]),
           '{:+.3f}'.format(north[i]),
           '{:.3f}'.format(horizontal[i])]
    for name in name_list :
      value = value_dict.get(case_path, {}).get(name, '')
      row.append('{:.6g}'.format(value) if isinstance(value, float) else str(value))
    row_list.append(row)

  width = [max(len(header[j]), max(len(row[j]) for row in row_list)) for j in range(0, len(header))]
  line  = '  '.join(header[j].rjust(width[j]) for j in range(0, len(header)))
  print(line)
  print('-'*len(line))
  for row in row_list :
    print('  '.join(row[j].rjust(width[j]) for j in range(0, len(header))))

  print()
  print('Number of cases: ', len(case_list))
  if name_reference is None :
    print('Origin: the mean of the cases')
  else :
    print('Origin: ', name_reference)
  print('Origin (ECEF): {:.4f}, {:.4f}, {:.4f} km'.format(*position_origin))
  print('Mean offset from the origin:  East {:+.3f} km,  North {:+.3f} km'
        .format(float(np.mean(east)), float(np.mean(north))))
  print('Standard deviation:  East {:.3f} km,  North {:.3f} km,  Up {:.4f} km'
        .format(np.std(east, ddof=1) if len(east) > 1 else 0.0,
                np.std(north, ddof=1) if len(north) > 1 else 0.0,
                np.std(up, ddof=1) if len(up) > 1 else 0.0))
  print('Horizontal distance from the origin: median {:.3f} km,  maximum {:.3f} km'
        .format(float(np.median(horizontal)), float(np.max(horizontal))))
  print('Terminal altitude: {:.4f} to {:.4f} km'
        .format(min(state['Alti'] for state in state_list),
                max(state['Alti'] for state in state_list)))

  return header, row_list


def write_csv(filename, header, row_list):
  with open(filename, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(row_list)
  print('Written: ', filename)

  return


def get_covariance_ellipse(east, north):
  #
  # 散らばりの共分散楕円。固有値分解して (長半径, 短半径, 傾き[deg]) を返す。
  # 円ではなく楕円で描くのは、風のばらつきが等方でないため（東西と南北で
  # 効き方が違う）。1 sigma の楕円は、正規分布なら点の 39% を含む。
  #
  covariance = np.cov(np.vstack([east, north]))
  eigenvalue, eigenvector = np.linalg.eigh(covariance)

  index = int(np.argmax(eigenvalue))
  semi_major = float(np.sqrt(max(eigenvalue[index], 0.0)))
  semi_minor = float(np.sqrt(max(eigenvalue[1-index], 0.0)))
  angle      = float(np.degrees(np.arctan2(eigenvector[1, index], eigenvector[0, index])))

  return semi_major, semi_minor, angle


def shorten_path(path, depth=3):
  #
  # 図に書き込むための短いパス。絶対パスをそのまま置くと枠からはみ出す。
  #
  part = os.path.normpath(path).split(os.sep)
  if len(part) <= depth :
    return os.path.normpath(path)
  return os.path.join('...', *part[-depth:])


def get_mark(text_list, unit_east, unit_north, position_origin):
  #
  # --mark で与えられた「ラベル=Tecplot ファイル」を読み、基準点からのずれ（東, 北）に直す。
  # モンテカルロの外で走らせた計算（たとえば実データの風のケース）を散布図に
  # 重ねるためのもので、統計には入れない。
  #
  mark_list = []
  for text in text_list :
    if '=' in text :
      label, path = text.split('=', 1)
    else :
      label, path = os.path.basename(os.path.dirname(os.path.dirname(text))), text
    if not os.path.exists(path) :
      print('Mark file not found: ', path)
      sys.exit(1)
    state    = get_terminal_state(path)
    offset   = np.array([state['X'], state['Y'], state['Z']])-position_origin
    mark_list.append((label.strip(), float(offset.dot(unit_east)), float(offset.dot(unit_north))))
  return mark_list


def plot_dispersion(filename, case_list, east, north, dpi, text_origin=None, mark_list=None):
  #
  # 着地点の散布図と分散円。matplotlib は Tacode の依存ではないので、ここでだけ import する。
  #
  try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
  except ImportError:
    print('matplotlib is not available; skipping the plot')
    return

  east_mean  = float(np.mean(east))
  north_mean = float(np.mean(north))
  radius     = np.hypot(east-east_mean, north-north_mean)

  figure = plt.figure(figsize=(7.0, 7.0))
  axis   = figure.add_subplot(1, 1, 1)

  # 基準点（--reference なら参照ケースの着地点、そうでなければケースの平均）
  axis.axhline(0.0, color='black', linewidth=0.5, zorder=1)
  axis.axvline(0.0, color='black', linewidth=0.5, zorder=1)
  axis.plot(0.0, 0.0, marker='+', markersize=14, markeredgewidth=1.6, color='black',
            linestyle='none', zorder=5, label='origin')

  # 分散円（k sigma の共分散楕円）と CEP 50%
  if len(east) > 2 :
    semi_major, semi_minor, angle = get_covariance_ellipse(east, north)
    for k in (1, 2, 3) :
      ellipse = Ellipse((east_mean, north_mean), 2.0*k*semi_major, 2.0*k*semi_minor,
                        angle=angle, facecolor='#1f77b4',
                        alpha=0.16 if k == 1 else 0.08, edgecolor='#1f77b4',
                        linestyle='--', linewidth=1.0, zorder=2,
                        label='1, 2, 3 sigma ellipse' if k == 1 else None)
      axis.add_patch(ellipse)

    cep = float(np.median(radius))
    circle = np.linspace(0.0, 2.0*np.pi, 181)
    axis.plot(east_mean+cep*np.cos(circle), north_mean+cep*np.sin(circle),
              color='#d62728', linestyle='-', linewidth=1.4, zorder=4,
              label='CEP 50% ({:.2f} km)'.format(cep))

  axis.scatter(east, north, s=38, color='#1f77b4', edgecolor='black', linewidth=0.5,
               zorder=6, label='cases ({:d})'.format(len(east)))
  axis.plot(east_mean, north_mean, marker='x', markersize=12, markeredgewidth=2.0,
            color='#d62728', linestyle='none', zorder=7, label='mean of the cases')

  # ケースが少ないうちは番号を振る（多いと読めなくなるだけなので止める）
  if len(case_list) <= 25 :
    for i, case_path in enumerate(case_list) :
      label = os.path.basename(case_path).replace('case', '').lstrip('0')
      axis.annotate(label, (east[i], north[i]),
                    textcoords='offset points', xytext=(5, 4), fontsize=7, color='gray')

  # モンテカルロの外の計算（統計には入れない）
  color_list = ['#2ca02c', '#9467bd', '#8c564b', '#e377c2']
  for i, (label, mark_east, mark_north) in enumerate(mark_list or []) :
    axis.plot(mark_east, mark_north, marker='*', markersize=17, linestyle='none',
              color=color_list[i % len(color_list)], markeredgecolor='black',
              markeredgewidth=0.5, zorder=8, label=label)

  text = 'mean offset  {:+.2f}, {:+.2f} km\nsigma  {:.2f} km (E), {:.2f} km (N)'.format(
    east_mean, north_mean,
    float(np.std(east, ddof=1)) if len(east) > 1 else 0.0,
    float(np.std(north, ddof=1)) if len(north) > 1 else 0.0)
  if text_origin is not None :
    text = 'origin: '+text_origin+'\n'+text
  axis.text(0.02, 0.02, text, transform=axis.transAxes, fontsize=9, va='bottom',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.9))

  axis.set_xlabel('East from the origin [km]')
  axis.set_ylabel('North from the origin [km]')
  axis.set_title('Terminal point dispersion')
  axis.set_aspect('equal', adjustable='datalim')
  axis.grid(True, linestyle=':', linewidth=0.5)
  axis.legend(loc='upper right', fontsize=9, framealpha=0.9)

  figure.tight_layout()
  figure.savefig(filename, dpi=dpi)
  plt.close(figure)
  print('Written: ', filename)

  return


def main():

  parser = argparse.ArgumentParser(
    description='Summarise the terminal points of a Monte-Carlo run made by tacode-montecarlo.py')
  parser.add_argument('directory', nargs='?', default=None,
                      help='the Monte-Carlo working directory (montecarlo.work_dir), '
                           'which holds case0001, case0002, ... and case_template')
  parser.add_argument('--tecplot', default=PATH_TECPLOT_DEFAULT,
                      help='path of the Tecplot output inside a case directory '
                           '(default: %(default)s)')
  parser.add_argument('--control', default=NAME_CONTROL_DEFAULT,
                      help='name of the control file inside a case directory '
                           '(default: %(default)s)')
  parser.add_argument('--reference', default=None,
                      help='Tecplot output of a reference case (for instance the same '
                           'entry without wind). The offsets are then measured from its '
                           'terminal point instead of from the mean of the cases')
  parser.add_argument('--mark', action='append', default=None, metavar='LABEL=PATH',
                      help='Tecplot output of another run to mark on the plot, as '
                           'LABEL=PATH. It is drawn as a star and left out of the '
                           'statistics. May be given more than once')
  parser.add_argument('-o', '--output', default=None,
                      help='write the table to this CSV file as well')
  parser.add_argument('--plot', default=None,
                      help='write a scatter plot of the dispersion to this file '
                           '(needs matplotlib)')
  parser.add_argument('--dpi', type=int, default=140, help='resolution of the plot')
  helper_config.add_argument(parser)
  args = helper_config.get_setting(parser, NAME_SECTION)

  if helper_config.save_file(args, NAME_SECTION) :
    return

  helper_config.require(args, 'directory', NAME_SECTION)

  if not os.path.isdir(args.directory) :
    print('Directory not found: ', args.directory)
    sys.exit(1)

  case_list_all = find_case_directory(args.directory)

  case_list  = []
  state_list = []
  for case_path in case_list_all :
    filename = os.path.join(case_path, args.tecplot)
    if not os.path.exists(filename) :
      print('Skipping (no output): ', case_path)
      continue
    case_list.append(case_path)
    state_list.append(get_terminal_state(filename))

  if len(case_list) == 0 :
    print('No case with an output was found in: ', args.directory)
    sys.exit(1)

  template_path = find_template_directory(args.directory)
  name_list, value_dict = get_dispersed_variable(case_list, template_path, args.control)

  position_reference = None
  if args.reference is not None :
    if not os.path.exists(args.reference) :
      print('Reference file not found: ', args.reference)
      sys.exit(1)
    state_reference    = get_terminal_state(args.reference)
    position_reference = np.array([state_reference['X'], state_reference['Y'], state_reference['Z']])

  position_origin, east, north, up = get_dispersion(state_list, position_reference)

  header, row_list = report(case_list, state_list, name_list, value_dict,
                            position_origin, east, north, up, args.reference)

  if args.output is not None :
    write_csv(args.output, header, row_list)

  mark_list = []
  if args.mark is not None :
    unit_east, unit_north, unit_up = get_local_horizon(position_origin)
    mark_list = get_mark(args.mark, unit_east, unit_north, position_origin)
    for label, mark_east, mark_north in mark_list :
      print('Mark: {:s}  East {:+.3f} km,  North {:+.3f} km'.format(label, mark_east, mark_north))

  if args.plot is not None :
    text_origin = 'mean of the cases' if args.reference is None else shorten_path(args.reference)
    plot_dispersion(args.plot, case_list, east, north, args.dpi, text_origin, mark_list)

  return


if __name__ == '__main__':

  main()
