#!/usr/bin/env python3
#
# Tacode の 6 自由度計算と、Project Fire flight II で測られた姿勢運動を突き合わせる。
#
# **比べるのは振動の「周期」と「振幅」で、扱いが違う。**
#
#   周期  比較できる。報告書（NASA TN D-4183）の図 17 は、自分の 6 自由度計算を
#         実測のピッチ・ヨーレートに合わせるために必要とした動圧を描いている。
#         スピンする軸対称体の横方向角速度は 2 つのモードに分かれ、**機体軸で見た**
#         周波数は
#
#           lambda = p - a + sqrt(a^2 + k q),   a = Ix p / (2 Iy),  k = -Cma S d / Iy
#
#         になる（Murphy BRL-1216 の tricyclic 解を機体軸へ戻したもの。動圧 0 では
#         lambda = p - a + a = p、つまり無トルクのときのスピンそのものに落ちる）。
#         **したがってその動圧は実測の周波数そのもの**である。ここでは Tacode の
#         出力から k を同定し（この式が実際に成り立つことの確認も兼ねる）、
#         それで q_required を周波数に直して Tacode の周波数と比べる
#
#   振幅  **比較にならない。**実測の全迎角の包絡線（図 19）は 3 度から 7.7 度、13 度、
#         19.5 度へ階段状に育つが、その段はすべて離散的な擾乱（熱量計の非対称な溶融、
#         熱遮蔽の投棄）が作ったもので、Tacode は擾乱を持たず機体を軸対称としている。
#         **したがって図 19 は「合うはずのないもの」として並べて描くだけ**にし、
#         rms のような成績にはしない
#
# 周波数の測り方: 機体軸の横方向角速度を複素数 xi = q + i r にまとめ、窓ごとに
# FFT を取って正の周波数側の山を拾う（スピンする軸対称体の 2 つのモードは複素平面で
# 正負の周波数に分かれるので、実数信号のスペクトルでは分けられない）。山は放物線
# 内挿で細かく決める。
#
# 使い方:
#   cd validation/fire2 && ./run_tacode.sh -file config_6dof.yml
#   python3 compare_attitude.py
#   python3 compare_attitude.py --case output_segments/tecplot.dat --no-figure

import argparse
import os
import sys

import numpy as np

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src'))
sys.path.append(os.path.join(DIRECTORY_SCRIPT, '../../src_helper/animate_trajectory'))

import tecplot_reader as tecplot_reader      # noqa: E402
from general.general import general          # noqa: E402

CASE_DEFAULT     = 'output_result_6dof/tecplot.dat'
DIRECTORY_OUTPUT = 'output_comparison'
TIME_START       = 1618.25

FILE_ANGLE    = 'reference/fire2_angle_attack.dat'
FILE_PRESSURE = 'reference/fire2_dynamic_pressure.dat'

# 周波数を測る窓の半幅（秒）と刻み。動圧はこの窓の中で 2 割ほど動くので、
# これ以上長くすると周波数が窓の中でぼける
WINDOW_HALF = 1.0
WINDOW_STEP = 0.5

# 窓の中の点がこれより少なければ測らない
POINT_MINIMUM = 40

# k を同定するのに使う窓の条件。**線形な範囲だけを使う**:
# 動圧が小さいうちは FFT が雑音を拾い、全迎角が大きくなると係数表の非線形が効く
PRESSURE_FIT = 5.0e3
ANGLE_FIT    = 5.0


def read_case(file_case):
  data = tecplot_reader.read_tecplot(file_case)
  data['Time'] = data['Time'] + TIME_START
  return data


def dynamic_pressure(data):
  return 0.5*data['Dens']*data['VelplAbs']**2


def peak_frequency(signal, timestep, previous=None):
  #
  # 複素信号の正の周波数側の山を放物線内挿で求める（rad/s）。
  # 山が見つからなければ None。
  #
  # **previous を与えるとその近くの山を選ぶ。**スピンする軸対称体は 2 つのモードを
  # 持ち、動圧が下がると速いほうの山が痩せて遅いほうに主役が移る。一番高い山を
  # 毎回選ぶと、そこで測っている量が入れ替わってしまう
  #
  window    = np.hanning(len(signal))
  spectrum  = np.fft.fftshift(np.fft.fft((signal - signal.mean())*window))
  frequency = np.fft.fftshift(np.fft.fftfreq(len(signal), timestep))*2.0*np.pi
  magnitude = np.abs(spectrum)

  select = frequency > 0.0
  if not select.any() :
    return None
  index_list = np.arange(len(frequency))[select]
  if previous is None :
    index = index_list[np.argmax(magnitude[select])]
  else:
    # 山（両隣より高い点）のうち、前の窓に最も近いもの
    peak = [item for item in index_list
            if 0 < item < len(frequency) - 1
            and magnitude[item] >= magnitude[item - 1]
            and magnitude[item] >= magnitude[item + 1]]
    if not peak :
      return None
    index = min(peak, key=lambda item: abs(frequency[item] - previous))
  if index <= 0 or index >= len(frequency) - 1 :
    return float(frequency[index])

  # 対数振幅に放物線を当てて山の位置を細かく決める
  left, middle, right = (np.log(max(magnitude[index + offset], 1.0e-300))
                         for offset in (-1, 0, 1))
  denominator = left - 2.0*middle + right
  shift = 0.0 if abs(denominator) < 1.0e-30 else 0.5*(left - right)/denominator
  spacing = frequency[1] - frequency[0]
  return float(frequency[index] + shift*spacing)


def measure_frequency(data):
  # 窓ごとの (時刻, 周波数, 動圧)
  time = data['Time']
  step = np.median(np.diff(time))
  xi   = np.radians(data['Q']) + 1j*np.radians(data['R'])
  qbar = dynamic_pressure(data)

  angle    = data['AoAtotal']
  result   = []
  previous = None
  centre   = time[0] + WINDOW_HALF
  while centre + WINDOW_HALF <= time[-1]:
    select = (time >= centre - WINDOW_HALF) & (time <= centre + WINDOW_HALF)
    if select.sum() >= POINT_MINIMUM :
      value = peak_frequency(xi[select], step, previous)
      if value is not None and value > 0.0 :
        previous = value
        result.append((centre, value, float(qbar[select].mean()),
                       float(angle[select].mean())))
    centre += WINDOW_STEP
  return np.array(result)


def frequency_model(qbar, spin, spin_term, slope):
  # lambda = p - a + sqrt(a^2 + k q)
  return spin - spin_term + np.sqrt(np.maximum(spin_term**2 + slope*qbar, 0.0))


def fit_relation(measurement, spin, spin_term):
  #
  # lambda = p - a + sqrt(a^2 + k q) の k を最小二乗で決める。a と p は慣性と
  # スピンから決まっている。**この形が Tacode の出力に実際に乗ることの確認**
  # でもあるので、残差と使った点数も返す
  #
  omega, qbar, angle = measurement[:, 1], measurement[:, 2], measurement[:, 3]
  select = (qbar > PRESSURE_FIT) & (angle < ANGLE_FIT)
  if select.sum() < 4 :
    return None, None, 0
  # (lambda - p + a)^2 - a^2 = k q を最小二乗で
  left  = (omega[select] - spin + spin_term)**2 - spin_term**2
  slope = float(np.dot(left, qbar[select])/np.dot(qbar[select], qbar[select]))
  model = frequency_model(qbar[select], spin, spin_term, slope)
  residual = float(np.sqrt(np.mean((omega[select] - model)**2)))
  return slope, residual, int(select.sum())


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--case', type=str, default=CASE_DEFAULT,
                      help='Tecplot output of a six-degree-of-freedom run')
  parser.add_argument('--file-config', type=str, default='config_6dof.yml')
  parser.add_argument('--no-figure', action='store_true')
  args = parser.parse_args()

  os.chdir(DIRECTORY_SCRIPT)
  if not os.path.exists(args.case) :
    print('No such output: ' + args.case)
    print('--Run ./run_tacode.sh -file config_6dof.yml first.')
    print('Program stopped.')
    sys.exit(1)

  config = general().read_config_yaml(args.file_config)
  inertia = config['attitude']['inertia_tensor']
  data    = read_case(args.case)

  spin      = np.radians(np.median(data['P']))
  spin_term = inertia['Ixx']*spin/(2.0*inertia['Iyy'])
  print('Six-degree-of-freedom run: ' + args.case)
  print('--Spin {:.3f} rad/s, Ix/Iy {:.4f}, spin term Ix p / (2 Iy) = {:.2f} rad/s'.format(
        spin, inertia['Ixx']/inertia['Iyy'], spin_term))

  measurement = measure_frequency(data)
  slope, residual, count = fit_relation(measurement, spin, spin_term)
  if slope is None :
    print('Not enough of the run is inside the atmosphere to measure a frequency.')
    print('Program stopped.')
    sys.exit(1)
  print('--lambda = p - a + sqrt(a^2 + k q) fits {:d} windows with k = {:.4e}'
        ' to {:.2f} rad/s rms'.format(count, slope, residual))
  print('  (k = -Cm_alpha S d / Iy, so this run behaves as Cm_alpha = {:.4f} per radian;'.format(
        -slope*inertia['Iyy']/(config['satellite']['characteristic_area']
                               *config['satellite']['characteristic_length'])))
  print('   the aerodynamic table gives -0.1290 about its own reference point and'
        ' -0.1423 about the centre of gravity)')

  reference = np.loadtxt(FILE_PRESSURE)
  print('')
  print('The oscillation frequency against the measured one (NASA TN D-4183 figure 17)')
  print('  {:>9s} {:>12s} {:>12s} {:>8s} {:>10s} {:>10s} {:>8s}'.format(
        'window', 'q Tacode', 'q required', 'ratio', 'omega calc', 'omega meas', 'ratio'))
  boundary = np.where(np.diff(reference[:, 0]) > 0.5)[0]
  lower = np.concatenate([[0], boundary + 1])
  upper = np.concatenate([boundary, [len(reference) - 1]])
  ratio_list = []
  for first, last in zip(lower, upper):
    time_mean = reference[first:last + 1, 0].mean()
    required  = reference[first:last + 1, 1].mean()
    case      = np.interp(time_mean, data['Time'], dynamic_pressure(data))
    omega_case   = frequency_model(case, spin, spin_term, slope)
    omega_flight = frequency_model(required, spin, spin_term, slope)
    ratio_list.append(omega_flight/omega_case)
    print('  {:9.1f} {:12.0f} {:12.0f} {:8.3f} {:10.2f} {:10.2f} {:8.3f}'.format(
          time_mean, case, required, required/case, omega_case, omega_flight,
          omega_flight/omega_case))
  ratio_list = np.array(ratio_list)
  print('  The measured oscillation is {:.1f} to {:.1f} % faster than the computed one'
        ' ({:.1f} % on average).'.format(
        100.0*(ratio_list.min() - 1.0), 100.0*(ratio_list.max() - 1.0),
        100.0*(ratio_list.mean() - 1.0)))

  envelope = np.loadtxt(FILE_ANGLE)
  select   = (data['Time'] >= envelope[0, 0]) & (data['Time'] <= envelope[-1, 0])
  print('')
  print('The total angle of attack against the measured envelope (figure 19)')
  print('  Tacode over the same interval: {:.2f} to {:.2f} deg.'.format(
        data['AoAtotal'][select].min(), data['AoAtotal'][select].max()))
  print('  Measured envelope            : {:.2f} to {:.2f} deg.'.format(
        envelope[:, 2].min(), envelope[:, 1].max()))
  print('  **These are not comparable.** Every step of the measured envelope was made')
  print('  by a discrete disturbance (asymmetric melting of a calorimeter, ejection of')
  print('  a heat shield) that Tacode has no way to represent; see README.md.')

  if not args.no_figure :
    make_figure(data, measurement, reference, envelope, spin, spin_term, slope)

  return


def make_figure(data, measurement, reference, envelope, spin, spin_term, slope):
  try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
  except ImportError:
    print('Caution: matplotlib is not installed, so no figure is written.')
    return

  if not os.path.isdir(DIRECTORY_OUTPUT) :
    os.makedirs(DIRECTORY_OUTPUT)

  figure, axis = plt.subplots(2, 1, figsize=(7.0, 7.0), sharex=True)

  axis[0].plot(data['Time'], data['AoAtotal'], color='tab:red', linewidth=0.6,
               label='Tacode, 6-DOF')
  axis[0].plot(envelope[:, 0], envelope[:, 1], 'k.', markersize=1.5,
               label='measured envelope (figure 19)')
  axis[0].plot(envelope[:, 0], envelope[:, 2], 'k.', markersize=1.5)
  axis[0].set_ylabel('Total angle of attack [deg.]')
  axis[0].set_ylim(0.0, 22.0)
  axis[0].legend(fontsize=8)
  axis[0].grid(alpha=0.3)
  axis[0].set_title('The envelope cannot be reproduced: its steps are disturbances',
                    fontsize=9)

  omega_case = frequency_model(dynamic_pressure(data), spin, spin_term, slope)
  axis[1].plot(data['Time'], omega_case, color='tab:red', linewidth=1.0,
               label='Tacode, from its own dynamic pressure')
  axis[1].plot(measurement[:, 0], measurement[:, 1], 'o', color='tab:red',
               markersize=3, fillstyle='none',
               label='Tacode, measured from the rates')
  omega_flight = frequency_model(reference[:, 1], spin, spin_term, slope)
  axis[1].plot(reference[:, 0], omega_flight, 'k.', markersize=2,
               label='flight (figure 17)')
  axis[1].set_ylabel('Pitch-yaw frequency [rad/s]')
  axis[1].set_xlabel('Elapsed flight time [s]')
  axis[1].set_xlim(envelope[0, 0], envelope[-1, 0])
  axis[1].legend(fontsize=8)
  axis[1].grid(alpha=0.3)

  figure.tight_layout()
  path = os.path.join(DIRECTORY_OUTPUT, 'profile_attitude.png')
  figure.savefig(path, dpi=150)
  plt.close(figure)
  print('Wrote ' + path)


if __name__ == '__main__':
  main()
