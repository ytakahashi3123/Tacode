#!/usr/bin/env python3
#
# レポート（tacode_v2_en.tex / tacode_v2_ja.tex）の図を作る。
#
# **図はすべてリポジトリの中のデータから作る。**チュートリアルの参照出力、
# validation/ の参照データと計算結果、それに時間刻みを振った計算そのものである。
# 外から持ってきた画像は座標系の模式図 4 枚だけで、これは manuscript_ver1.2_en
# から複製して PDF に直したもの（figure/*_coordinate.pdf ほか）。
#
# 計算結果が要る図は、その出力が無ければ**作らずに、作り方を表示して次へ進む**。
#
# 使い方:
#   python3 generate_figure.py                 # 全部
#   python3 generate_figure.py --only fire2    # 名前に fire2 を含む図だけ
#   python3 generate_figure.py --list          # 図の一覧
#   python3 generate_figure.py --skip-run      # 計算を回す図（収束次数）を飛ばす

import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
DIRECTORY_ROOT   = os.path.normpath(os.path.join(DIRECTORY_SCRIPT, '..'))
DIRECTORY_FIGURE = os.path.join(DIRECTORY_SCRIPT, 'figure')

sys.path.append(os.path.join(DIRECTORY_ROOT, 'src_helper/animate_trajectory'))
import tecplot_reader as tecplot_reader      # noqa: E402

# 図の体裁。1 段組の本文幅に 1 枚か 2 枚並べる前提
SIZE_SINGLE = (5.4, 3.6)
SIZE_WIDE   = (7.0, 4.2)
SIZE_TALL   = (5.4, 6.4)
RESOLUTION  = 300

COLOUR_CASE   = 'tab:red'
COLOUR_SECOND = 'tab:blue'

# 突入界面からの経過時刻に直すための、飛行経過時刻の起点
TIME_START_FIRE2 = 1618.25


def path_root(*name):
  return os.path.join(DIRECTORY_ROOT, *name)


def require(*path_list):
  # 必要な入力が揃っているか。揃っていなければ足りないものを返す
  return [path for path in path_list if not os.path.exists(path)]


def announce_missing(name, missing, command):
  print('  {:s}: skipped, {:d} input(s) missing'.format(name, len(missing)))
  for path in missing:
    print('    ' + os.path.relpath(path, DIRECTORY_ROOT))
  print('    produce them with: ' + command)


def save(figure, name):
  if not os.path.isdir(DIRECTORY_FIGURE) :
    os.makedirs(DIRECTORY_FIGURE)
  path = os.path.join(DIRECTORY_FIGURE, name + '.pdf')
  figure.tight_layout()
  figure.savefig(path, dpi=RESOLUTION)
  print('  wrote figure/' + name + '.pdf')


def read_columns(path):
  # `# name name ...` のヘッダを持つ素のテキスト表
  name_list = None
  row_list  = []
  with open(path) as stream:
    for line in stream:
      if line.startswith('#') :
        item = line[1:].split()
        if item and all(part.replace('_', '').isalpha() for part in item) :
          name_list = item
        continue
      row_list.append([float(item) for item in line.split()])
  data = np.array(row_list)
  return {name: data[:, index] for index, name in enumerate(name_list)}


# ---------------------------------------------------------------- schematic

def figure_body_axes(plt):
  #
  # 機体軸と空力角の模式図。左が迎角・横滑り角・全迎角と空力ロール角、
  # 右が 3 自由度の揚力とバンク角
  #
  figure, axis = plt.subplots(1, 2, figsize=(7.4, 3.4))

  # 左: 機体軸から見た対気速度
  left = axis[0]
  left.annotate('', xy=(1.0, 0.0), xytext=(0.0, 0.0),
                arrowprops=dict(arrowstyle='-|>', color='k', lw=1.4))
  left.annotate('', xy=(0.0, -0.85), xytext=(0.0, 0.0),
                arrowprops=dict(arrowstyle='-|>', color='k', lw=1.4))
  left.annotate('', xy=(-0.35, 0.45), xytext=(0.0, 0.0),
                arrowprops=dict(arrowstyle='-|>', color='k', lw=1.4))
  left.text(1.03, 0.0, r'$x_b$ (forward)', va='center', fontsize=9)
  left.text(0.0, -0.93, r'$z_b$ (down)', ha='center', va='top', fontsize=9)
  left.text(-0.40, 0.50, r'$y_b$ (right)', ha='right', fontsize=9)

  angle = np.radians(28.0)
  left.annotate('', xy=(0.86*np.cos(angle), -0.86*np.sin(angle)), xytext=(0.0, 0.0),
                arrowprops=dict(arrowstyle='-|>', color=COLOUR_CASE, lw=1.8))
  left.text(0.90*np.cos(angle) + 0.02, -0.90*np.sin(angle),
            r'$\boldsymbol{v}_{\rm air}$', color=COLOUR_CASE, fontsize=10)
  sweep = np.linspace(0.0, -angle, 40)
  left.plot(0.45*np.cos(sweep), 0.45*np.sin(sweep), color=COLOUR_CASE, lw=0.9)
  left.text(0.50*np.cos(0.5*angle), -0.50*np.sin(0.5*angle) + 0.06,
            r'$\eta$', color=COLOUR_CASE, fontsize=11)

  # 迎角と横滑り角は全迎角を分けたもの
  left.text(0.62, -0.10, r'$\eta^2 = \alpha^2 + \beta^2$', color=COLOUR_CASE,
            fontsize=9)
  left.set_xlim(-0.75, 1.75)
  left.set_ylim(-1.15, 0.75)
  left.set_aspect('equal')
  left.axis('off')
  left.set_title('Body axes and the total angle of attack', fontsize=9)

  # 右: 速度に直交する面内での揚力の向き
  right = axis[1]
  circle = plt.Circle((0.0, 0.0), 1.0, fill=False, color='0.6', lw=0.9)
  right.add_patch(circle)
  right.annotate('', xy=(0.0, 1.0), xytext=(0.0, 0.0),
                 arrowprops=dict(arrowstyle='-|>', color='k', lw=1.2))
  right.annotate('', xy=(1.0, 0.0), xytext=(0.0, 0.0),
                 arrowprops=dict(arrowstyle='-|>', color='k', lw=1.2))
  right.text(0.03, 1.05, 'up (geocentric)', fontsize=9)
  right.text(1.04, 0.0, 'right of flight', va='center', fontsize=9)

  bank = np.radians(50.0)
  right.annotate('', xy=(np.sin(bank), np.cos(bank)), xytext=(0.0, 0.0),
                 arrowprops=dict(arrowstyle='-|>', color=COLOUR_CASE, lw=1.8))
  right.text(1.04*np.sin(bank), 1.06*np.cos(bank), r'$\boldsymbol{F}_{\rm lift}$',
             color=COLOUR_CASE, fontsize=10)
  sweep = np.linspace(0.0, bank, 40)
  right.plot(0.35*np.sin(sweep), 0.35*np.cos(sweep), color=COLOUR_CASE, lw=0.9)
  right.text(0.42*np.sin(0.5*bank), 0.42*np.cos(0.5*bank), r'$\sigma$',
             color=COLOUR_CASE, fontsize=11)
  right.plot(0.0, 0.0, 'ko', markersize=4)
  right.text(0.0, -1.18, r'$\boldsymbol{v}_{\rm air}$ out of the page',
             ha='center', fontsize=9)

  right.set_xlim(-1.30, 2.05)
  right.set_ylim(-1.45, 1.35)
  right.set_aspect('equal')
  right.axis('off')
  right.set_title('Bank angle of the three-degree-of-freedom lift', fontsize=9)

  save(figure, 'body_axes')
  plt.close(figure)
  return True


# ------------------------------------------------------------ order of accuracy

def run_timestep_study(directory_case, file_config, timestep_list, python,
                       time_maximum, time_probe):
  #
  # 時間刻みを振って走らせ、**共通の時刻での状態**を返す。
  #
  # 最終状態で比べてはいけない。時間ループは `time < time_maximum` で止まるので
  # 終了時刻が刻みによって最大 dt ずれ、地表近くの急な軌道ではその差が km に育つ。
  # 1 ステップごとに出力させ、全部の計算に共通の時刻へ内挿して比べる。
  #
  # **ケースのディレクトリで走らせる**（database/ のパスがカレントディレクトリ
  # 基準の manual なので）。
  #
  import yaml
  directory_work = os.path.join(directory_case, 'output_order')
  if not os.path.isdir(directory_work) :
    os.makedirs(directory_work)

  with open(os.path.join(directory_case, file_config)) as stream:
    config_base = yaml.safe_load(stream)

  result = {}
  for timestep in timestep_list:
    tag = 'dt{:g}'.format(timestep)
    config = dict(config_base)
    config['time_integration'] = dict(config_base['time_integration'])
    config['time_integration']['timestep_constant'] = timestep
    config['computational_setup'] = dict(config_base['computational_setup'])
    config['computational_setup']['time_elapsed_maximum'] = time_maximum
    config['post_process'] = dict(config_base['post_process'])
    config['post_process']['directory_output'] = os.path.join('output_order', tag)
    config['post_process']['kml'] = dict(config_base['post_process']['kml'])
    config['post_process']['kml']['flag_output'] = False
    config['post_process']['tecplot'] = dict(config_base['post_process']['tecplot'])
    config['post_process']['tecplot']['frequency_output'] = 1
    config['restart_process'] = dict(config_base['restart_process'])
    config['restart_process']['directory_output'] = os.path.join('output_order', tag)

    name_config = os.path.join('output_order', 'config_' + tag + '.yml')
    with open(os.path.join(directory_case, name_config), 'w') as stream:
      yaml.safe_dump(config, stream, default_flow_style=False, sort_keys=False)

    with open(os.path.join(directory_work, 'log_' + tag), 'w') as stream:
      answer = subprocess.run([python, '../../src/tacode.py', '-file', name_config],
                              cwd=directory_case, stdout=stream,
                              stderr=subprocess.STDOUT)
    if answer.returncode != 0 :
      print('    the solver failed at dt = {:g}'.format(timestep))
      return None
    data = tecplot_reader.read_tecplot(
             os.path.join(directory_work, tag,
                          config['post_process']['tecplot']['filename_output']))
    state = [np.interp(time_probe, data['Time'], data[name])
             for name in ('X', 'Y', 'Z')]
    if 'q0' in data :
      state += [np.interp(time_probe, data['Time'], data[name])
                for name in ('q0', 'q1', 'q2', 'q3')]
    result[timestep] = np.array(state)
  return result


def figure_order_of_accuracy(plt, python, flag_run):
  if not flag_run :
    print('  order_of_accuracy: skipped (--skip-run)')
    return False

  directory = path_root('validation/fire2')
  missing = require(os.path.join(directory, 'config.yml'),
                    os.path.join(directory, 'config_6dof.yml'))
  if missing :
    announce_missing('order_of_accuracy', missing, 'see validation/fire2/README.md')
    return False

  print('  order_of_accuracy: running the solver at several time steps...')
  # **減速のピーク（73 g, t = 1650 s）を含む窓を使う。**最初の数十秒だけだと
  # 空力がほとんど効かず、誤差が丸めと表の内挿の床に当たる
  step_3dof = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625]
  step_6dof = [0.08, 0.04, 0.02, 0.01, 0.005, 0.0025]

  result_3dof = run_timestep_study(directory, 'config.yml', step_3dof, python,
                                   time_maximum=40.0, time_probe=39.0)
  result_6dof = run_timestep_study(directory, 'config_6dof.yml', step_6dof, python,
                                   time_maximum=35.0, time_probe=34.0)
  if result_3dof is None or result_6dof is None :
    return False

  figure, axis = plt.subplots(figsize=SIZE_SINGLE)
  for result, step_list, label, colour, marker in (
      (result_3dof, step_3dof, 'three degrees of freedom, position', COLOUR_CASE, 'o'),
      (result_6dof, step_6dof, 'six degrees of freedom, attitude', COLOUR_SECOND, 's')):
    finest = result[step_list[-1]]
    error, step = [], []
    for timestep in step_list[:-1]:
      if label.startswith('three') :
        difference = np.linalg.norm(result[timestep][0:3] - finest[0:3])*1.0e3
      else:
        difference = np.linalg.norm(result[timestep][3:7] - finest[3:7])
      error.append(difference)
      step.append(timestep)
    axis.loglog(step, error, marker + '-', color=colour, label=label, markersize=4)
    # **細かいほうの 3 点で傾きを測る。**粗いほうの端は振動や減速のピークを
    # 解けていないので漸近域に入っていない（図の左端が折れているのがそれ）
    order = np.polyfit(np.log(step[-3:]), np.log(error[-3:]), 1)[0]
    print('    {:s}: observed order {:.2f}'.format(label, order))
    axis.annotate('order {:.1f}'.format(order),
                  xy=(step[-2], error[-2]), xytext=(1.7*step[-2], 0.18*error[-2]),
                  color=colour, fontsize=8)

  reference = np.array([step_3dof[1], step_3dof[4]])
  axis.loglog(reference, 1.0e-1*(reference/reference[0])**4, 'k--', lw=0.8)
  axis.text(reference[-1], 1.0e-1*(reference[-1]/reference[0])**4,
            r'  slope 4', fontsize=8)
  axis.set_xlabel(r'Time step $\Delta t$ [s]')
  axis.set_ylabel('Difference from the finest step')
  axis.grid(alpha=0.3, which='both')
  axis.legend(fontsize=8)
  save(figure, 'order_of_accuracy')
  plt.close(figure)
  return True


# --------------------------------------------------------------- examples

def figure_example_orbit(plt):
  path = path_root('tutorial/work/output_result/tecplot.dat')
  missing = require(path)
  if missing :
    announce_missing('example_orbit', missing, 'cd tutorial/work && ./run_tacode.sh')
    return False

  data = tecplot_reader.read_tecplot(path)
  figure, axis = plt.subplots(1, 2, figsize=SIZE_WIDE)
  axis[0].plot(data['Time'], data['Alti'], color=COLOUR_CASE, lw=1.0)
  axis[0].set_xlabel('Time [s]')
  axis[0].set_ylabel('Altitude [km]')
  axis[0].grid(alpha=0.3)
  axis[1].plot(data['Long'], data['Lati'], '.', color=COLOUR_CASE, markersize=1.0)
  axis[1].set_xlabel('Longitude [deg.]')
  axis[1].set_ylabel('Latitude [deg.]')
  axis[1].set_xlim(-180.0, 180.0)
  axis[1].set_ylim(-90.0, 90.0)
  axis[1].grid(alpha=0.3)
  save(figure, 'example_orbit')
  plt.close(figure)
  return True


def figure_example_reentry(plt):
  path = path_root('tutorial/work_reentry/output_result/tecplot.dat')
  missing = require(path)
  if missing :
    announce_missing('example_reentry', missing,
                     'cd tutorial/work_reentry && ./run_tacode.sh')
    return False

  data = tecplot_reader.read_tecplot(path)
  figure, axis = plt.subplots(1, 2, figsize=SIZE_WIDE)
  axis[0].plot(data['VelplAbs']*1.0e-3, data['Alti'], color=COLOUR_CASE, lw=1.2)
  axis[0].set_xlabel('Relative velocity [km/s]')
  axis[0].set_ylabel('Altitude [km]')
  axis[0].grid(alpha=0.3)

  pressure = 0.5*data['Dens']*data['VelplAbs']**2
  axis[1].plot(data['Time'], pressure*1.0e-3, color=COLOUR_CASE, lw=1.2)
  axis[1].set_xlabel('Time [s]')
  axis[1].set_ylabel(r'Dynamic pressure [kN/m$^2$]')
  axis[1].grid(alpha=0.3)
  save(figure, 'example_reentry')
  plt.close(figure)
  return True


def figure_example_reentry_6dof(plt):
  path = path_root('tutorial/work_reentry_6dof/output_result/tecplot.dat')
  missing = require(path)
  if missing :
    announce_missing('example_reentry_6dof', missing,
                     'cd tutorial/work_reentry_6dof && ./run_tacode.sh')
    return False

  data = tecplot_reader.read_tecplot(path)
  figure, axis = plt.subplots(2, 1, figsize=SIZE_TALL, sharex=True)
  axis[0].plot(data['Time'], data['AoAtotal'], color=COLOUR_CASE, lw=0.8)
  axis[0].set_ylabel('Total angle of attack [deg.]')
  axis[0].grid(alpha=0.3)
  for name, colour, label in (('P', 'k', r'$p$'), ('Q', COLOUR_CASE, r'$q$'),
                              ('R', COLOUR_SECOND, r'$r$')):
    axis[1].plot(data['Time'], data[name], color=colour, lw=0.7, label=label)
  axis[1].set_xlabel('Time [s]')
  axis[1].set_ylabel('Angular velocity [deg/s]')
  axis[1].legend(fontsize=8, ncol=3)
  axis[1].grid(alpha=0.3)
  save(figure, 'example_reentry_6dof')
  plt.close(figure)
  return True


# ------------------------------------------------------------- validation

def figure_apollo4(plt):
  directory = path_root('validation/apollo4')
  path_case   = os.path.join(directory, 'output_result_lift/tecplot.dat')
  path_flight = os.path.join(directory, 'reference/apollo4_flight.dat')
  missing = require(path_case, path_flight)
  if missing :
    announce_missing('apollo4', missing,
                     'cd validation/apollo4 && ./run_tacode.sh -file config_lift.yml')
    return False

  data   = tecplot_reader.read_tecplot(path_case)
  flight = tecplot_reader.read_tecplot(path_flight)
  time_flight = flight['Time'] - flight['Time'][0]

  figure, axis = plt.subplots(1, 2, figsize=SIZE_WIDE)
  axis[0].plot(time_flight, flight['Alti'], 'k.', markersize=2, label='flight')
  axis[0].plot(data['Time'], data['Alti'], color=COLOUR_CASE, lw=1.2,
               label='Tacode, measured vertical lift')
  axis[0].set_xlabel('Time from the entry interface [s]')
  axis[0].set_ylabel('Altitude [km]')
  axis[0].legend(fontsize=8)
  axis[0].grid(alpha=0.3)

  axis[1].plot(flight['VelrelAbs']*1.0e-3, flight['Alti'], 'k.', markersize=2,
               label='flight')
  axis[1].plot(data['VelplAbs']*1.0e-3, data['Alti'], color=COLOUR_CASE, lw=1.2,
               label='Tacode')
  axis[1].set_xlabel('Relative velocity [km/s]')
  axis[1].set_ylabel('Altitude [km]')
  axis[1].legend(fontsize=8)
  axis[1].grid(alpha=0.3)
  save(figure, 'apollo4_trajectory')
  plt.close(figure)
  return True


def velocity_inertial(data, rotation_rate=7.292115e-5):
  # ECEF 速度を慣性速度の大きさに直す（apollo10 の図 13 は慣性速度で描かれている）
  radius = np.sqrt(data['X']**2 + data['Y']**2)*1.0e3
  speed  = rotation_rate*radius
  return np.sqrt(data['VelplAbs']**2 + 2.0*speed*data['Upl'] + speed**2)


def figure_apollo10(plt):
  directory = path_root('validation/apollo10')
  path_case   = os.path.join(directory, 'output_result/tecplot.dat')
  path_flight = os.path.join(directory, 'reference/apollo10_flight.dat')
  path_roll   = os.path.join(directory, 'reference/apollo10_roll.dat')
  missing = require(path_case, path_flight, path_roll)
  if missing :
    announce_missing('apollo10', missing, 'cd validation/apollo10 && ./run_tacode.sh')
    return False

  data   = tecplot_reader.read_tecplot(path_case)
  flight = tecplot_reader.read_tecplot(path_flight)
  roll   = tecplot_reader.read_tecplot(path_roll)

  figure, axis = plt.subplots(3, 1, figsize=(5.4, 7.4), sharex=True)
  axis[0].plot(flight['Time'], flight['Alti'], 'k.', markersize=2, label='flight')
  axis[0].plot(data['Time'], data['Alti'], color=COLOUR_CASE, lw=1.2, label='Tacode')
  axis[0].set_ylabel('Altitude [km]')
  axis[0].legend(fontsize=8)
  axis[0].grid(alpha=0.3)

  axis[1].plot(flight['Time'], flight['VelinAbs']*1.0e-3, 'k.', markersize=2,
               label='flight')
  axis[1].plot(data['Time'], velocity_inertial(data)*1.0e-3, color=COLOUR_CASE,
               lw=1.2, label='Tacode')
  axis[1].set_ylabel('Inertial velocity [km/s]')
  axis[1].grid(alpha=0.3)

  axis[2].plot(roll['Time'], roll['Bank'], 'k-', lw=0.9)
  axis[2].set_xlabel('Time from the entry interface [s]')
  axis[2].set_ylabel('Bank angle, input [deg.]')
  axis[2].grid(alpha=0.3)
  save(figure, 'apollo10_trajectory')
  plt.close(figure)
  return True


def figure_fire2_aerodynamics(plt):
  path = path_root('validation/fire2/reference/fire2_aerodynamics.dat')
  missing = require(path)
  if missing :
    announce_missing('fire2_aerodynamics', missing,
                     'cd validation/fire2 && python3 digitize_aerodynamics.py')
    return False

  data = read_columns(path)
  select = data['alpha'] <= 80.0
  figure, axis = plt.subplots(1, 3, figsize=(7.4, 2.6))
  for index, (name, label) in enumerate((('CX', r'$C_X$'), ('CR', r'$C_R$'),
                                         ('CM', r'$C_m$'))):
    axis[index].plot(data['alpha'][select], data[name][select], 'o-',
                     color=COLOUR_CASE, markersize=2.5, lw=1.0)
    axis[index].set_xlabel(r'$\alpha$ [deg.]')
    axis[index].set_ylabel(label)
    axis[index].grid(alpha=0.3)
  save(figure, 'fire2_aerodynamics')
  plt.close(figure)
  return True


def read_fire2_trajectory():
  path = path_root('validation/fire2/reference/fire2_trajectory.dat')
  name = ['time', 'latitude', 'longitude', 'altitude', 'velocity', 'flightpath',
          'heading', 'dynamic_pressure', 'pressure', 'density', 'temperature',
          'mach', 'acceleration', 'reynolds']
  row = []
  with open(path) as stream:
    for line in stream:
      if line.startswith('#') :
        continue
      row.append([float(item) for item in line.split()])
  value = np.array(row)
  return {key: value[:, index] for index, key in enumerate(name)}


def figure_fire2_trajectory(plt):
  directory = path_root('validation/fire2')
  path_case    = os.path.join(directory, 'output_result/tecplot.dat')
  path_matched = os.path.join(directory, 'output_matched/tecplot.dat')
  missing = require(path_case, path_matched,
                    os.path.join(directory, 'reference/fire2_trajectory.dat'))
  if missing :
    announce_missing('fire2_trajectory', missing,
                     'cd validation/fire2 && ./run_tacode.sh '
                     '&& ./run_tacode.sh -file config_matched.yml')
    return False

  reference = read_fire2_trajectory()
  figure, axis = plt.subplots(2, 1, figsize=SIZE_TALL, sharex=True,
                              gridspec_kw={'height_ratios': [2, 1]})
  axis[0].plot(reference['time'], reference['altitude']*1.0e-3, 'k.', markersize=1.5,
               label='NASA TN D-3569 table V')

  for path, label, colour in ((path_case, r'Tacode, $C_D$ from figure 4', COLOUR_CASE),
                              (path_matched, r'Tacode, $C_D$ matched to table V',
                               COLOUR_SECOND)):
    data = tecplot_reader.read_tecplot(path)
    time = data['Time'] + TIME_START_FIRE2
    axis[0].plot(time, data['Alti'], color=colour, lw=1.0, label=label)
    select = np.isfinite(reference['altitude']) \
             & (reference['time'] >= time[0]) & (reference['time'] <= time[-1])
    error = np.interp(reference['time'][select], time, data['Alti']*1.0e3) \
            - reference['altitude'][select]
    axis[1].plot(reference['time'][select], error, color=colour, lw=0.9)

  axis[0].set_ylabel('Altitude [km]')
  axis[0].legend(fontsize=8)
  axis[0].grid(alpha=0.3)
  axis[1].axhline(0.0, color='k', lw=0.8)
  axis[1].set_xlabel('Elapsed flight time [s]')
  axis[1].set_ylabel('Tacode - table V [m]')
  axis[1].grid(alpha=0.3)
  save(figure, 'fire2_trajectory')
  plt.close(figure)
  return True


def figure_fire2_density(plt):
  directory = path_root('validation/fire2')
  path_table = os.path.join(directory, 'database/atmosphere/atmospheremodel_fire2.txt')
  missing = require(path_table, os.path.join(directory, 'reference/fire2_trajectory.dat'))
  if missing :
    announce_missing('fire2_density', missing, 'see validation/fire2/README.md')
    return False

  sys.path.append(os.path.join(DIRECTORY_ROOT, 'validation/fire2'))
  import compare_fire2                                   # noqa: E402
  altitude_table, density_table = compare_fire2.read_atmosphere_table(path_table)

  reference = read_fire2_trajectory()
  select = np.isfinite(reference['density']) & np.isfinite(reference['altitude']) \
           & (reference['altitude'] <= altitude_table[-1])
  altitude = reference['altitude'][select]
  ratio = reference['density'][select] \
          / np.interp(altitude, altitude_table, density_table)

  figure, axis = plt.subplots(figsize=(4.2, 4.6))
  axis.plot(ratio, altitude*1.0e-3, 'k.', markersize=2.5)
  axis.axvline(1.0, color=COLOUR_CASE, lw=1.0)
  axis.set_xlabel('Density, flight (table V) / NRLMSISE-00 table')
  axis.set_ylabel('Altitude [km]')
  axis.set_xlim(0.5, 1.5)
  axis.grid(alpha=0.3)
  save(figure, 'fire2_density')
  plt.close(figure)
  return True


def figure_fire2_attitude(plt):
  directory = path_root('validation/fire2')
  path_case = os.path.join(directory, 'output_result_6dof/tecplot.dat')
  path_angle = os.path.join(directory, 'reference/fire2_angle_attack.dat')
  path_pressure = os.path.join(directory, 'reference/fire2_dynamic_pressure.dat')
  missing = require(path_case, path_angle, path_pressure)
  if missing :
    announce_missing('fire2_attitude', missing,
                     'cd validation/fire2 && ./run_tacode.sh -file config_6dof.yml')
    return False

  sys.path.append(os.path.join(DIRECTORY_ROOT, 'validation/fire2'))
  import compare_attitude                                # noqa: E402

  data = tecplot_reader.read_tecplot(path_case)
  data['Time'] = data['Time'] + TIME_START_FIRE2
  envelope  = np.loadtxt(path_angle)
  reference = np.loadtxt(path_pressure)

  inertia_x, inertia_y = 3.511, 2.806
  spin      = np.radians(np.median(data['P']))
  spin_term = inertia_x*spin/(2.0*inertia_y)
  measurement = compare_attitude.measure_frequency(data)
  slope = compare_attitude.fit_relation(measurement, spin, spin_term)[0]

  figure, axis = plt.subplots(2, 1, figsize=SIZE_TALL, sharex=True)
  axis[0].plot(envelope[:, 0], envelope[:, 1], 'k.', markersize=1.5,
               label='measured envelope, figure 19')
  axis[0].plot(envelope[:, 0], envelope[:, 2], 'k.', markersize=1.5)
  axis[0].plot(data['Time'], data['AoAtotal'], color=COLOUR_CASE, lw=0.7,
               label='Tacode, six degrees of freedom')
  axis[0].set_ylabel('Total angle of attack [deg.]')
  axis[0].set_ylim(0.0, 22.0)
  axis[0].legend(fontsize=8)
  axis[0].grid(alpha=0.3)

  omega_case = compare_attitude.frequency_model(
                 compare_attitude.dynamic_pressure(data), spin, spin_term, slope)
  omega_flight = compare_attitude.frequency_model(reference[:, 1], spin, spin_term,
                                                  slope)
  axis[1].plot(data['Time'], omega_case, color=COLOUR_CASE, lw=1.0, label='Tacode')
  axis[1].plot(measurement[:, 0], measurement[:, 1], 'o', color=COLOUR_CASE,
               markersize=3, fillstyle='none', label='Tacode, from the rates')
  axis[1].plot(reference[:, 0], omega_flight, 'k.', markersize=2.5,
               label='flight, figure 17')
  axis[1].set_xlabel('Elapsed flight time [s]')
  axis[1].set_ylabel('Pitch-yaw frequency [rad/s]')
  axis[1].set_xlim(envelope[0, 0], envelope[-1, 0])
  axis[1].legend(fontsize=8)
  axis[1].grid(alpha=0.3)
  save(figure, 'fire2_attitude')
  plt.close(figure)
  return True


# -------------------------------------------------------------------- main

FIGURE = [
  ('body_axes',            lambda plt, python, run: figure_body_axes(plt)),
  ('order_of_accuracy',    lambda plt, python, run: figure_order_of_accuracy(plt, python, run)),
  ('example_orbit',        lambda plt, python, run: figure_example_orbit(plt)),
  ('example_reentry',      lambda plt, python, run: figure_example_reentry(plt)),
  ('example_reentry_6dof', lambda plt, python, run: figure_example_reentry_6dof(plt)),
  ('apollo4_trajectory',   lambda plt, python, run: figure_apollo4(plt)),
  ('apollo10_trajectory',  lambda plt, python, run: figure_apollo10(plt)),
  ('fire2_aerodynamics',   lambda plt, python, run: figure_fire2_aerodynamics(plt)),
  ('fire2_trajectory',     lambda plt, python, run: figure_fire2_trajectory(plt)),
  ('fire2_density',        lambda plt, python, run: figure_fire2_density(plt)),
  ('fire2_attitude',       lambda plt, python, run: figure_fire2_attitude(plt)),
]


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--only', type=str, default=None,
                      help='Only the figures whose name contains this')
  parser.add_argument('--list', action='store_true', help='List the figures')
  parser.add_argument('--skip-run', action='store_true',
                      help='Skip the figures which run the solver')
  parser.add_argument('--python', type=str,
                      default=os.environ.get('TACODE_PYTHON', sys.executable),
                      help='Interpreter to run the solver with')
  args = parser.parse_args()

  if args.list :
    for name, dummy in FIGURE:
      print(name)
    return

  try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
  except ImportError:
    print('matplotlib is not installed; it is needed to draw the report figures.')
    print('Program stopped.')
    sys.exit(1)
  plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9,
                       'savefig.bbox': 'tight'})

  print('Writing the report figures into ' + os.path.relpath(DIRECTORY_FIGURE))
  count = 0
  for name, builder in FIGURE:
    if args.only is not None and args.only not in name :
      continue
    if builder(plt, args.python, not args.skip_run) :
      count += 1
  print('{:d} figure(s) written'.format(count))

  return


if __name__ == '__main__':
  main()
