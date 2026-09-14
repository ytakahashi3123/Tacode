#!/usr/bin/env python3
#
# 熱遮蔽の投棄で機体が変わるのを、**区間ごとの計算を繋いで**飛ばす。
#
# Fire の再突入体は飛行中に 3 枚のベリリウム熱量計と 2 枚のフェノール遮蔽を順に捨て、
# 質量が 24 %、ピッチ慣性が 34 %、基準直径が 13 % 減る（NASA TN D-4183 表 I）。
# Tacode は 1 回の計算の途中で機体を変えられないので、表 I の構成ごとに計算を切り、
# **前の区間の最終状態を次の区間の初期条件にして**繋ぐ。
#
# Tacode のリスタート（`flag_initial: False`）は未実装なので、繋ぎはケース側で行う。
# 引き継ぐのは Tecplot 出力の最終行:
#
#   位置   Long, Lati, Alti      -> initial_settings.coordinate
#   速度   Upl, Vpl, Wpl         -> initial_settings.velocity（地心ローカル [東, 北, 上]）
#   姿勢   Yaw, Pitch, Roll      -> initial_settings.attitude（地心 NED に対する 3-2-1）
#   角速度 P, Q, R               -> initial_settings.angular_velocity（**ECEF 基準**, deg/s）
#
# どれも config が要求する規約そのままの量が出力に載っている。Tecplot 出力は倍精度で
# 書かれるので、繋ぎ目で状態が丸められることはない。
#
# **区間の計算は 1 ステップごとに出力させる。**そうしないと最後の出力が区間の終わりより
# 手前になり、その差（最大で出力間隔ぶん）だけ次の区間が早く始まってしまう。0.5 秒
# ずれると高度が 700 m ずれるので、これは無視できない。繋いだあとで元の出力間隔に
# 間引く。
#
# **投棄のときの角力積は入れていない。**報告書は 1640.38 s に 6.440 N・m・s、
# 1648.18 s に 9.49 N・m・s の角力積があったと測っているが、向きが分からない。
# 入れれば迎角の包絡線は合うが、それは合わせ込みであって再現ではない。
# したがってこの計算は**質量・慣性・寸法の階段状の変化だけ**を持つ。
#
# 使い方:
#   python3 run_segments.py                        # config_6dof.yml を区間に割って走らせる
#   python3 run_segments.py --file-config config.yml
#   python3 run_segments.py --python ~/venvs/myenv/bin/python

import argparse
import copy
import os
import subprocess
import sys

import numpy as np
import yaml

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src_helper/animate_trajectory'))

import tecplot_reader as tecplot_reader      # noqa: E402

DIRECTORY_OUTPUT = 'output_segments'
FILE_SOLVER      = '../../src/tacode.py'

# 計算を始める飛行経過時刻（config.yml と揃えること）
TIME_START = 1618.25

#
# NASA TN D-4183 表 I。構成ごとの質量・重心・慣性・基準直径と、係数表の
# モーメント基準点 x_mc。
#
#   time         その構成になる飛行経過時刻
#   mass         kg
#   diameter     m（基準直径。基準面積は pi d^2 / 4）
#   inertia      Ix, Iy, Iz, kg m2
#   offset_cg    x_mc - x_cg, m。重心は基準点より**前**にあるので正
#
# **質量の 86.568 kg は TN D-3569 表 III の値**を採る。TN D-4183 表 I は 86.586 と
# 印字しているが、TN D-3569 表 I の積み上げ（94.73 - 8.16 = 86.57）は 86.568 に
# 丸まり 86.586 には丸まらないので、後者は数字の入れ替わりと見る。差は 0.02 %。
#
# x_mc は表の印字値をそのまま使う。脚注の式 x_mc = x_r - r + 0.2875 d は最初の 3 行では
# 印字値と 0.4 mm 以内で合うが、後の 2 行では 6.5 mm ずれる（報告書側の丸め）。
#
# 3 枚目の熱量計は完全には溶けなかったので、1647.53 s 以降は質量も慣性も変えない
# （表 I の脚注 2）。
#
SEGMENT = [
  {'time': 1617.75, 'mass': 86.568, 'diameter': 0.672,
   'inertia': (3.511, 2.806, 2.874), 'offset_cg': 0.306 - 0.277,
   'name': 'complete reentry package'},
  {'time': 1640.48, 'mass': 83.189, 'diameter': 0.651,
   'inertia': (3.281, 2.644, 2.698), 'offset_cg': 0.307 - 0.282,
   'name': 'less first calorimeter'},
  {'time': 1642.12, 'mass': 76.022, 'diameter': 0.630,
   'inertia': (2.806, 2.305, 2.359), 'offset_cg': 0.313 - 0.293,
   'name': 'less first phenolic layer'},
  {'time': 1646.10, 'mass': 72.166, 'diameter': 0.607,
   'inertia': (2.562, 2.128, 2.183), 'offset_cg': 0.319 - 0.299,
   'name': 'less second calorimeter'},
  {'time': 1647.53, 'mass': 66.179, 'diameter': 0.587,
   'inertia': (2.183, 1.857, 1.925), 'offset_cg': 0.322 - 0.309,
   'name': 'less second phenolic layer'},
]


def read_config(file_config):
  with open(file_config) as stream:
    return yaml.safe_load(stream)


def write_config(config, path):
  directory = os.path.dirname(path)
  if directory and not os.path.isdir(directory) :
    os.makedirs(directory)
  with open(path, 'w') as stream:
    stream.write('# Written by run_segments.py. Do not edit: the comments of the\n'
                 '# configuration it was built from are lost here on purpose, because\n'
                 '# this file is an intermediate product, not a case configuration.\n')
    yaml.safe_dump(config, stream, default_flow_style=False, sort_keys=False)


def final_state(file_tecplot):
  # 直前の区間の最終行を、config が要求する形で返す。経過時間も返す
  data  = tecplot_reader.read_tecplot(file_tecplot)
  index = -1
  state = {'coordinate': [float(data['Long'][index]), float(data['Lati'][index]),
                          float(data['Alti'][index])],
           'velocity': [float(data['Upl'][index]), float(data['Vpl'][index]),
                        float(data['Wpl'][index])]}
  if 'Yaw' in data :
    state['attitude'] = [float(data['Yaw'][index]), float(data['Pitch'][index]),
                         float(data['Roll'][index])]
    state['angular_velocity'] = [float(data['P'][index]), float(data['Q'][index]),
                                 float(data['R'][index])]
  return state, float(data['Time'][index])


def build_segment(config_base, segment, time_begin, time_end, state, directory):
  config = copy.deepcopy(config_base)
  config['computational_setup']['time_elapsed_maximum'] = time_end - time_begin

  diameter = segment['diameter']
  config['satellite']['mass'] = segment['mass']
  config['satellite']['characteristic_area'] = 0.25*np.pi*diameter*diameter
  config['satellite']['characteristic_length'] = diameter

  if state is not None :
    for key in state:
      config['initial_settings'][key] = state[key]

  if config.get('attitude', {}).get('flag_attitude', False) :
    inertia = segment['inertia']
    config['attitude']['inertia_tensor'] = {'Ixx': inertia[0], 'Iyy': inertia[1],
                                            'Izz': inertia[2], 'Ixy': 0.0,
                                            'Iyz': 0.0, 'Izx': 0.0}
    config['attitude']['center_of_gravity'] = [segment['offset_cg'], 0.0, 0.0]

  config['post_process']['directory_output'] = directory
  config['restart_process']['directory_output'] = directory
  config['post_process']['kml']['flag_output'] = False
  config['post_process']['tecplot']['frequency_output'] = 1
  return config


def gather(file_list, time_list, file_output, frequency):
  #
  # 区間ごとの Tecplot 出力を 1 本に繋ぐ。**時刻は最初の区間の開始からの経過時間**に
  # 直し、繋ぎ目の重複（次の区間の 0 秒目）は落とす
  #
  header, variable, row_list = None, None, []
  for path, time_offset in zip(file_list, time_list):
    with open(path) as stream:
      line_list = stream.readlines()
    index = 0
    while index < len(line_list) and not line_list[index].lower().startswith('variables') :
      index += 1
    if header is None :
      header   = line_list[:index]
      variable = line_list[index]
    elif line_list[index] != variable :
      print('The segments were written with different variables:', path)
      print('Program stopped.')
      sys.exit(1)
    count = 0
    for line in line_list[index + 2:]:
      item = line.split()
      if not item :
        continue
      keep  = (count % frequency == 0)
      count += 1
      if keep and count == 1 and row_list :
        continue                      # 繋ぎ目の重複
      if not keep :
        continue
      item[0] = '{:.10g}'.format(float(item[0]) + time_offset)
      row_list.append(' '.join(item) + '\n')

  with open(file_output, 'w') as stream:
    for line in header:
      stream.write(line)
    stream.write(variable)
    stream.write('zone t=time i= {:d} f=point\n'.format(len(row_list)))
    for line in row_list:
      stream.write(line)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--file-config', type=str, default='config_6dof.yml',
                      help='Configuration the segments are built from')
  parser.add_argument('--python', type=str, default=os.environ.get('TACODE_PYTHON', 'python3'),
                      help='Interpreter to run the solver with')
  parser.add_argument('--output', type=str, default=DIRECTORY_OUTPUT)
  args = parser.parse_args()

  os.chdir(DIRECTORY_SCRIPT)
  config_base = read_config(args.file_config)
  time_final  = TIME_START + config_base['computational_setup']['time_elapsed_maximum']

  if not os.path.isdir(args.output) :
    os.makedirs(args.output)

  file_list, time_list = [], []
  state, time_begin_next = None, None
  for index, segment in enumerate(SEGMENT):
    time_begin = max(segment['time'], TIME_START)
    if time_begin_next is not None :
      # 前の区間が実際に終わった時刻。表 I の切り替え時刻と出力の刻みは割り切れない
      time_begin = time_begin_next
    time_end   = SEGMENT[index + 1]['time'] if index + 1 < len(SEGMENT) else time_final
    time_end   = min(time_end, time_final)
    if time_end <= time_begin :
      continue

    directory   = os.path.join(args.output, 'segment{:d}'.format(index))
    file_config = os.path.join(args.output, 'config_segment{:d}.yml'.format(index))
    write_config(build_segment(config_base, segment, time_begin, time_end, state,
                               directory), file_config)

    print('Segment {:d}: t = {:.2f} to {:.2f} s, {:s}, {:.3f} kg, d = {:.3f} m'.format(
          index, time_begin, time_end, segment['name'], segment['mass'],
          segment['diameter']))
    with open(os.path.join(args.output, 'log_segment{:d}'.format(index)), 'w') as stream:
      result = subprocess.run([args.python, FILE_SOLVER, '-file', file_config],
                              stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode != 0 :
      print('The solver failed on segment {:d}; see {:s}'.format(
            index, os.path.join(args.output, 'log_segment{:d}'.format(index))))
      print('Program stopped.')
      sys.exit(1)

    path = os.path.join(directory, config_base['post_process']['tecplot']['filename_output'])
    file_list.append(path)
    time_list.append(time_begin - TIME_START)
    state, elapsed = final_state(path)
    time_begin_next = time_begin + elapsed

  file_output = os.path.join(args.output, 'tecplot.dat')
  gather(file_list, time_list, file_output,
         config_base['post_process']['tecplot']['frequency_output'])
  print('Wrote ' + file_output)

  return


if __name__ == '__main__':
  main()
