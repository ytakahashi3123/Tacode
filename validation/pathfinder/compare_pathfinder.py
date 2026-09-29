#!/usr/bin/env python3
#
# Mars Pathfinder（1997-07-04）の飛行データと Tacode を突き合わせる。
#
# 3 つの config は 1 つずつ要素を替えてある（README.md）:
#
#   config.yml        飛行の大気 + 飛行の CD(Mach, Kn)   -> ソルバー
#   config_mcd.yml    MCD の大気 + 飛行の CD(Mach, Kn)   -> 大気モデル
#   config_panel.yml  飛行の大気 + Tacode のパネル法 CD(Kn) -> 空力モデル
#
# 比べる相手:
#   高度      reference/pathfinder_ddr.dat の altitude_tacode（NASA の再構成。突入状態の枠に
#             直してある。convert_pds.py）
#   速さ      reference/pathfinder_reconstruction.dat（実測加速度を突入状態から積分したもの）
#   減速      reference/pathfinder_acceleration.dat（**実測**、32 Hz の軸方向加速度）
#   密度      reference/pathfinder_ddr.dat（NASA の再構成。実測加速度と空力データベースから）
#
# Tacode の減速は 0.5 rho V^2 CD A / m で、CD は計算と同じ表を同じ関数
# （src/satellite）で、出力の Kn と Mach から引き直す。
#
# 時刻: Tacode の 0 は突入界面（Spencer et al. 1998 表 1 のエポック）で、PDS の時刻では
# 38.202 s。比較は迫撃砲の点火（突入界面から 169.6 s）まで。
#
# 使い方:
#   cd validation/pathfinder && ./run_tacode.sh
#   python3 compare_pathfinder.py
#   python3 compare_pathfinder.py --no-figure

import argparse
import os
import sys

import numpy as np
import yaml

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src'))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src_helper/animate_trajectory'))

import tecplot_reader as tecplot_reader      # noqa: E402
import satellite.satellite as satellite      # noqa: E402

TIME_ENTRY_INTERFACE = 50.482 - 12.28        # s（PDS の時刻）
TIME_MORTAR          = 169.6                 # s（突入界面から）

CASE_LIST = (
  ('config.yml',       'flight atmosphere, flight CD(Mach, Kn)'),
  ('config_mcd.yml',   'MCD atmosphere, flight CD(Mach, Kn)'),
  ('config_panel.yml', 'flight atmosphere, panel-method CD(Kn)'),
)

COLOUR_CASE = ('tab:red', 'tab:blue', 'tab:orange')


class quiet:
  # Tacode の読み取り関数の print を黙らせる
  def __enter__(self):
    self.stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
  def __exit__(self, *args):
    sys.stdout.close()
    sys.stdout = self.stdout


def read_config(name):
  with open(os.path.join(DIRECTORY_SCRIPT, name)) as f:
    config = yaml.safe_load(f)
  section = config['satellite']
  section['directory_aerodynamic'] = os.path.join(DIRECTORY_SCRIPT, section['directory_aerodynamic'])
  return config


def drag_deceleration(config, data):
  # Tacode の出力から、計算と同じ CD で 0.5 rho V^2 CD A / m を作る
  with quiet():
    aerodynamic_dict = satellite.initial_settings_satellite(config)
  if aerodynamic_dict[satellite.KEY_MACH] is None:
    cd = np.array([satellite.get_aerodynamic_coefficient(k, aerodynamic_dict[satellite.KEY_KN],
                                                         aerodynamic_dict[satellite.KEY_CD_MEAN]) for k in data['Kn']])
  else:
    cd = np.array([satellite.get_aerodynamic_coefficient_mach(k, m, aerodynamic_dict)
                   for k, m in zip(data['Kn'], data['Mach'])])
  area, mass = config['satellite']['characteristic_area'], config['satellite']['mass']
  return 0.5*data['Dens']*data['VelplAbs']**2*cd*area/mass, cd


def read_reference():
  ddr = np.loadtxt(os.path.join(DIRECTORY_SCRIPT, 'reference/pathfinder_ddr.dat'))
  rec = np.loadtxt(os.path.join(DIRECTORY_SCRIPT, 'reference/pathfinder_reconstruction.dat'))
  acc = np.loadtxt(os.path.join(DIRECTORY_SCRIPT, 'reference/pathfinder_acceleration.dat'))
  return ddr, rec, acc


def peak(time, value):
  index = int(np.argmax(value))
  if 0 < index < len(value) - 1:
    # 3 点の放物線で頂点を取る（サンプル間隔より細かく）
    y0, y1, y2 = value[index-1:index+2]
    shift = 0.5*(y0 - y2)/(y0 - 2.0*y1 + y2)
    step = time[index+1] - time[index]
    return time[index] + shift*step, y1 - 0.25*(y0 - y2)*shift
  return time[index], value[index]


def compare(name, label, ddr, rec, acc):
  config = read_config(name)
  directory = config['post_process']['directory_output']
  path = os.path.join(DIRECTORY_SCRIPT, directory, 'tecplot.dat')
  if not os.path.exists(path):
    print('{}: no output ({}). Run ./run_tacode.sh {}'.format(name, path, name))
    return None
  data = tecplot_reader.read_tecplot(path)
  time = data['Time'] + TIME_ENTRY_INTERFACE
  deceleration, cd = drag_deceleration(config, data)

  window = (ddr[:,0] >= time[0]) & (ddr[:,0] <= TIME_ENTRY_INTERFACE + TIME_MORTAR)
  time_ref = ddr[window,0]
  altitude = np.interp(time_ref, time, data['Alti']) - ddr[window,2]
  speed    = np.interp(time_ref, time, data['VelplAbs']) - rec[window,1]

  window_acc = (acc[:,0] >= time[0]) & (acc[:,0] <= time[-1])
  time_acc = acc[window_acc,0]
  measured = acc[window_acc,1]
  model = np.interp(time_acc, time, deceleration)
  # 減速の比は、減速が大きい区間（1 g 以上）で見る。小さいところは実測の分解能で割り算が暴れる
  strong = measured > 9.8
  ratio = model[strong]/measured[strong]

  time_peak_model, value_peak_model = peak(time, deceleration)
  time_peak_meas,  value_peak_meas  = peak(time_acc, measured)

  # 迫撃砲の点火の時点
  index_end = len(data['Time']) - 1
  time_end = time[index_end]
  altitude_ref_end = np.interp(time_end, ddr[:,0], ddr[:,2])
  speed_ref_end    = np.interp(time_end, rec[:,0], rec[:,1])

  print('')
  print('{} -- {}'.format(name, label))
  print('  Altitude minus the reconstruction : rms {:.3f} km, max {:+.3f} km'.format(
        np.sqrt(np.mean(altitude**2)), altitude[np.argmax(np.abs(altitude))]))
  print('  Speed minus the reconstruction    : rms {:.1f} m/s, max {:+.1f} m/s'.format(
        np.sqrt(np.mean(speed**2)), speed[np.argmax(np.abs(speed))]))
  print('  Deceleration / measured (above 1 g): mean {:.4f}, rms of (ratio - 1) {:.2f} %, range {:.3f} .. {:.3f}'.format(
        np.mean(ratio), 100.0*np.sqrt(np.mean((ratio - 1.0)**2)), ratio.min(), ratio.max()))
  print('  Peak deceleration                 : {:.2f} m/s2 ({:.2f} g0) at {:.2f} s; measured {:.2f} m/s2 at {:.2f} s'.format(
        value_peak_model, value_peak_model/9.80665, time_peak_model - TIME_ENTRY_INTERFACE,
        value_peak_meas, time_peak_meas - TIME_ENTRY_INTERFACE))
  print('  At the mortar fire ({:.1f} s)       : altitude {:.3f} km (reconstruction {:.3f}), speed {:.1f} m/s ({:.1f})'.format(
        time_end - TIME_ENTRY_INTERFACE, data['Alti'][index_end], altitude_ref_end, data['VelplAbs'][index_end], speed_ref_end))

  return {'name': name, 'label': label, 'data': data, 'time': time, 'deceleration': deceleration, 'cd': cd,
          'time_ref': time_ref, 'altitude': altitude, 'speed': speed, 'time_acc': time_acc,
          'ratio_time': time_acc[strong], 'ratio': ratio}


def compare_density(result, ddr):
  # MCD の密度と飛行の密度の比を、高度帯ごとに（同じ高度で比べる）
  #
  # 飛行の密度の高度は 2 通りある。altitude_tacode は突入状態の枠（Tacode が飛ぶ枠）で、
  # edl_ddr.tab の高度そのものより 1.2 km 低い（convert_pds.py）。**どちらを採るかで比が
  # 1 割以上変わる**（スケールハイトが 8〜10 km なので 1.2 km で 13 %）ので、両方を出す
  data = result['data']
  order = np.argsort(data['Alti'])
  altitude_file = ddr[:,1] - 3396.0
  ratio = np.interp(ddr[:,2], data['Alti'][order], data['Dens'][order])/ddr[:,5]
  ratio_file = np.interp(altitude_file, data['Alti'][order], data['Dens'][order])/ddr[:,5]
  print('')
  print('Density of {} against the flight, at the same altitude'.format(result['name']))
  print('  (entry-state frame, as flown by Tacode | the altitude of edl_ddr.tab as it is):')
  for low, high in ((100.0, 130.0), (60.0, 100.0), (40.0, 60.0), (20.0, 40.0), (10.0, 20.0), (1.5, 10.0)):
    select = (ddr[:,2] >= low) & (ddr[:,2] < high)
    if select.any():
      print('  {:5.1f} - {:5.1f} km: ratio mean {:.3f}, range {:.3f} .. {:.3f} | mean {:.3f}'.format(
            low, high, np.mean(ratio[select]), ratio[select].min(), ratio[select].max(), np.mean(ratio_file[select])))
  return ratio


def plot(results, ddr, rec, acc, density_ratio, directory):
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt

  os.makedirs(directory, exist_ok=True)
  t0 = TIME_ENTRY_INTERFACE

  fig, axis = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
  axis[0].plot(ddr[:,0] - t0, ddr[:,2], 'k-', lw=2.5, alpha=0.4, label='NASA reconstruction (PDS)')
  for result, colour in zip(results, COLOUR_CASE):
    axis[0].plot(result['time'] - t0, result['data']['Alti'], color=colour, lw=1, label=result['label'])
    axis[1].plot(result['time_ref'] - t0, result['altitude']*1.e3, color=colour, lw=1)
  axis[0].set_ylabel('Altitude above 3396 km [km]')
  axis[0].legend(fontsize=8)
  axis[1].set_ylabel('Tacode - reconstruction [m]')
  axis[1].set_xlabel('Time from the entry interface [s]')
  for a in axis:
    a.grid(alpha=0.3)
  fig.tight_layout()
  fig.savefig(os.path.join(directory, 'profile_time-altitude.png'), dpi=150)
  plt.close(fig)

  fig, axis = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
  select = (acc[:,0] > t0) & (acc[:,0] < t0 + TIME_MORTAR)
  axis[0].plot(acc[select,0] - t0, acc[select,1], 'k-', lw=2.5, alpha=0.4, label='measured (32 Hz)')
  for result, colour in zip(results, COLOUR_CASE):
    axis[0].plot(result['time'] - t0, result['deceleration'], color=colour, lw=1, label=result['label'])
    axis[1].plot(result['ratio_time'] - t0, result['ratio'], color=colour, lw=1)
  axis[0].set_ylabel('Axial deceleration [m/s$^2$]')
  axis[0].legend(fontsize=8)
  axis[1].axhline(1.0, color='k', lw=0.5)
  axis[1].set_ylabel('Tacode / measured (above 1 g)')
  axis[1].set_xlabel('Time from the entry interface [s]')
  for a in axis:
    a.grid(alpha=0.3)
  fig.tight_layout()
  fig.savefig(os.path.join(directory, 'profile_time-deceleration.png'), dpi=150)
  plt.close(fig)

  fig, axis = plt.subplots(figsize=(7, 4.5))
  valid = rec[:,6] > 0
  axis.plot(rec[valid,4], rec[valid,3], 'k.', ms=2, alpha=0.5, label='traced back from the flight')
  for result, colour in zip(results, COLOUR_CASE):
    if result['name'] == 'config_mcd.yml':
      continue
    mach = result['data']['Mach'] if 'Mach' in result['data'] else rec[valid,4].max()*0 + np.interp(
           result['time'], rec[:,0], rec[:,4])
    axis.plot(mach, result['cd'], color=colour, lw=1, label=result['label'])
  axis.set_xscale('log')
  axis.set_xlabel('Mach number')
  axis.set_ylabel('CD')
  axis.legend(fontsize=8)
  axis.grid(alpha=0.3, which='both')
  fig.tight_layout()
  fig.savefig(os.path.join(directory, 'profile_mach-cd.png'), dpi=150)
  plt.close(fig)

  if density_ratio is not None:
    fig, axis = plt.subplots(figsize=(5, 6))
    axis.plot(density_ratio, ddr[:,2], 'b-', lw=1)
    axis.axvline(1.0, color='k', lw=0.5)
    axis.set_xlabel('MCD density / flight density')
    axis.set_ylabel('Altitude above 3396 km [km]')
    axis.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(directory, 'profile_density-ratio.png'), dpi=150)
    plt.close(fig)

  print('')
  print('Figures written to', directory)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--no-figure', action='store_true')
  parser.add_argument('--output', type=str, default=os.path.join(DIRECTORY_SCRIPT, 'output_comparison'))
  args = parser.parse_args()

  ddr, rec, acc = read_reference()
  results = []
  for name, label in CASE_LIST:
    result = compare(name, label, ddr, rec, acc)
    if result is not None:
      results.append(result)

  density_ratio = None
  for result in results:
    if result['name'] == 'config_mcd.yml':
      density_ratio = compare_density(result, ddr)

  if not args.no_figure and len(results) > 0:
    try:
      import matplotlib  # noqa: F401
    except ImportError:
      print('matplotlib is not installed; no figure')
      return
    plot(results, ddr, rec, acc, density_ratio, args.output)


if __name__ == '__main__':
  main()
