#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/29
#
# Tacode の軌道出力（tecplot.dat）から、軌道に沿った流れの量を作って Tecplot 形式で書く。
# 計算はしない（ソルバーは走らせない）。出力を読むだけ。
#
#   python3 src_helper/flow_properties/flow_properties.py output_result/tecplot.dat
#
# 出す量（1 行が軌道の 1 点）:
#   Time[s] Alti[km] Dens[kg/m3] Temp[K] Pres[Pa] Velair[m/s] Mach Re Kn Qdyn[Pa]
#   HeatFlux[W/m2] HeatLoad[J/m2]
#
# Tacode v1 の src_helper/compute_flowfield_fromTrajectory と同じ量を出す。v1 との違い:
#   - 速さは対気速度（風があれば VelairAbs、無ければ VelplAbs）
#   - Kn は**ソルバーが出した列をそのまま使う**（v1 は直径 4e-10 m の剛体球で作り直していた。
#     ソルバーの Kn は大気表の分子種ごとの直径と satellite.characteristic_length から作る）
#   - Mach 数の比熱比と平均分子量は、既定でケースの config.yml の atmosphere 節
#     （Mach 軸のある空力表でソルバーが使うのと同じ値）から読む。そうすれば出力の
#     Mach 列（あれば）と一致する
#   - 熱流束は W/m2（v1 は Tauber の式の係数を kW/m2 のまま使い、見出しだけ W/m2 と
#     書いていた）。総加熱量の列を見出しにも書く（v1 は値だけ書いて見出しに無かった）
#
# 設定はカレントディレクトリの config_helper.yml の flow_properties セクションに書ける
# （コマンドラインが優先。--save-config でいまの設定を書き出せる）。

import argparse
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
NAME_SECTION = 'flow_properties'

# 一般気体定数 J/(mol K)（CODATA 2018。src/satellite/satellite.py と同じ値）
CONSTANT_GAS_UNIVERSAL = 8.314462618

# 気体ごとの既定値:
#   比熱比・平均分子量 g/mol（config.yml にもコマンドラインにも無いときに使う）、
#   粘性係数の Sutherland 式 mu = mu0 (T/T0)^1.5 (T0 + S)/(T + S) の mu0 Pa s, T0 K, S K、
#   Sutton-Graves の熱流束 q = k sqrt(rho/Rn) V^3 の k kg^0.5/m（Sutton and Graves 1971）
#
#   co2: 比熱比・分子量・Sutherland は v1 の compute_flowfield_fromTrajectory.yaml と同じ
#        （CO2 95.57 %, N2 2.70 %, Ar 1.60 %, O2 0.13 % を自由度 7/5/3/5 で混ぜたもの。
#         分子量 43.49 は Pathfinder の ASI/MET が 100 km 以下で使った値とも一致する）
#   air: 1.4 と 28.9644（国際標準大気）、Sutherland は White, Viscous Fluid Flow の空気の値
GAS = {
  'co2': {'specific_heat_ratio': 1.290723, 'molecular_weight': 43.49,
          'sutherland': [1.37e-5, 293.0, 240.0], 'sutton_graves': 1.9027e-4},
  'air': {'specific_heat_ratio': 1.4, 'molecular_weight': 28.9644,
          'sutherland': [1.716e-5, 273.15, 110.4], 'sutton_graves': 1.7415e-4},
}

# 熱流束 q = C rho^a V^b / sqrt(Rn)（W/m2、rho kg/m3、V m/s、Rn m）
#   sutton-graves: C = k（気体ごと）, a = 0.5, b = 3
#   tauber-mars  : C = 1.35e-4, a = 0.5, b = 3.04（Tauber, Bowles and Yang 1989 の火星の対流加熱。
#                  原典は 1.35e-8 W/cm2。v1 はこれを 1.35e-7 と書いて kW/m2 を出していた）
#   custom       : --heat-flux-law C a b
HEATING = ('sutton-graves', 'tauber-mars', 'custom')
LAW_TAUBER_MARS = [1.35e-4, 0.5, 3.04]

# ノーズ半径を与えなかったときに使う値, m。根拠のある値ではないので、使ったら必ず知らせる
NOSE_RADIUS_DEFAULT = 1.0

VARIABLES = ['Time[s]', 'Alti[km]', 'Dens[kg/m3]', 'Temp[K]', 'Pres[Pa]', 'Velair[m/s]', 'Mach', 'Re', 'Kn',
             'Qdyn[Pa]', 'HeatFlux[W/m2]', 'HeatLoad[J/m2]']


def read_control(filename):
  #
  # ソルバーの制御ファイルから比熱比・平均分子量・代表長さを拾う。無ければ空の辞書
  #
  if filename is None or not os.path.exists(filename) :
    return {}
  with open(filename) as f:
    config = yaml.safe_load(f)
  result = {}
  section = config.get('atmosphere', {}) or {}
  for key in ('specific_heat_ratio', 'molecular_weight'):
    if key in section :
      result[key] = float(section[key])
  section = config.get('satellite', {}) or {}
  if 'characteristic_length' in section :
    result['length'] = float(section['characteristic_length'])
  return result


def resolve_setting(args):
  #
  # 比熱比・平均分子量・Reynolds 数の代表長さを、コマンドライン（と設定ファイル）>
  # ソルバーの制御ファイル > 気体の既定値 の順で決め、どこから取ったかも返す
  #
  control = read_control(args.control)
  gas = GAS[args.gas]
  value, source = {}, {}
  for key in ('specific_heat_ratio', 'molecular_weight'):
    if getattr(args, key) is not None :
      value[key], source[key] = float(getattr(args, key)), 'given'
    elif key in control :
      value[key], source[key] = control[key], args.control
    else :
      value[key], source[key] = gas[key], 'the default for ' + args.gas
  if args.length is not None :
    value['length'], source['length'] = float(args.length), 'given'
  elif 'length' in control :
    value['length'], source['length'] = control['length'], args.control + ' (satellite.characteristic_length)'
  else :
    print('Give the length for the Reynolds number with --length, or a control file with')
    print('satellite.characteristic_length (--control, default: config.yml).')
    sys.exit(1)

  value['sutherland'] = [float(v) for v in (args.sutherland if args.sutherland is not None else gas['sutherland'])]
  if args.heating == 'sutton-graves' :
    value['heat_flux_law'] = [gas['sutton_graves'], 0.5, 3.0]
  elif args.heating == 'tauber-mars' :
    value['heat_flux_law'] = list(LAW_TAUBER_MARS)
  else :
    if args.heat_flux_law is None :
      print('--heating custom needs --heat-flux-law C a b (q = C rho^a V^b / sqrt(Rn), W/m2).')
      sys.exit(1)
    value['heat_flux_law'] = [float(v) for v in args.heat_flux_law]

  if args.nose_radius is None :
    value['nose_radius'], source['nose_radius'] = NOSE_RADIUS_DEFAULT, 'not given'
    print('Caution: the nose radius is not given (--nose-radius, or nose_radius in the settings file).')
    print('--The heat flux and the heat load are computed with Rn = {:g} m, which is only a placeholder;'.format(
          NOSE_RADIUS_DEFAULT))
    print('  they scale as 1/sqrt(Rn).')
  else :
    value['nose_radius'], source['nose_radius'] = float(args.nose_radius), 'given'

  if value['specific_heat_ratio'] <= 1.0 or value['molecular_weight'] <= 0.0 or value['length'] <= 0.0 \
     or value['nose_radius'] <= 0.0 :
    print('The ratio of specific heats must be above 1, and the molecular weight, the length and')
    print('the nose radius positive.')
    sys.exit(1)

  return value, source


def viscosity_sutherland(temperature, mu0, temperature_reference, constant):
  return mu0*(temperature/temperature_reference)**1.5*(temperature_reference + constant)/(temperature + constant)


def compute(data, value, nose_radius):
  #
  # 軌道に沿った流れの量。data は tecplot_reader.read_tecplot が返す辞書（単位を落とした列名）
  #
  time        = data['Time']
  density     = data['Dens']
  temperature = data['Temp']
  speed       = data['VelairAbs'] if 'VelairAbs' in data else data['VelplAbs']

  constant_gas = CONSTANT_GAS_UNIVERSAL/(value['molecular_weight']*1.e-3)
  pressure = density*constant_gas*temperature
  mach = speed/np.sqrt(value['specific_heat_ratio']*constant_gas*temperature)
  viscosity = viscosity_sutherland(temperature, *value['sutherland'])
  reynolds = density*speed*value['length']/viscosity
  pressure_dynamic = 0.5*density*speed**2
  coefficient, exponent_density, exponent_speed = value['heat_flux_law']
  heat_flux = coefficient*density**exponent_density*speed**exponent_speed/np.sqrt(nose_radius)
  heat_load = np.concatenate([[0.0], np.cumsum(0.5*(heat_flux[1:] + heat_flux[:-1])*np.diff(time))])

  return {'Time[s]': time, 'Alti[km]': data['Alti'], 'Dens[kg/m3]': density, 'Temp[K]': temperature,
          'Pres[Pa]': pressure, 'Velair[m/s]': speed, 'Mach': mach, 'Re': reynolds, 'Kn': data['Kn'],
          'Qdyn[Pa]': pressure_dynamic, 'HeatFlux[W/m2]': heat_flux, 'HeatLoad[J/m2]': heat_load}


def write_tecplot(filename, result, header_comment, fmt='.6e'):
  with open(filename, 'w') as f:
    for line in header_comment :
      f.write('# ' + line + '\n')
    # 見出しは Tacode 本体の tecplot.dat と同じ書き方（引用符なし、カンマ区切り）。
    # src_helper の tecplot_reader でも Tecplot でも読める
    f.write('Variables = ' + ','.join(VARIABLES) + '\n')
    f.write('zone t=trajectory i= {:d} f=point\n'.format(len(result['Time[s]'])))
    for row in zip(*[result[name] for name in VARIABLES]):
      f.write(' '.join(format(v, fmt) for v in row) + '\n')


def summarise(result):
  index = {name: int(np.argmax(result[name])) for name in ('HeatFlux[W/m2]', 'Qdyn[Pa]')}
  i = index['HeatFlux[W/m2]']
  print('Peak heat flux   : {:.4g} W/m2 at {:.1f} s, {:.2f} km (Mach {:.2f}, Re {:.3g}, Kn {:.3g})'.format(
        result['HeatFlux[W/m2]'][i], result['Time[s]'][i], result['Alti[km]'][i],
        result['Mach'][i], result['Re'][i], result['Kn'][i]))
  i = index['Qdyn[Pa]']
  print('Peak dynamic pressure: {:.4g} Pa at {:.1f} s, {:.2f} km (Mach {:.2f})'.format(
        result['Qdyn[Pa]'][i], result['Time[s]'][i], result['Alti[km]'][i], result['Mach'][i]))
  print('Heat load        : {:.4g} J/m2 over {:.1f} s'.format(result['HeatLoad[J/m2]'][-1],
        result['Time[s]'][-1] - result['Time[s]'][0]))


def argument():
  parser = argparse.ArgumentParser(
    description='Flow properties along a Tacode trajectory (Mach, Reynolds and Knudsen numbers, '
                'dynamic pressure, stagnation-point heat flux and heat load)')
  parser.add_argument('filename', nargs='?', default=os.path.join('output_result', 'tecplot.dat'),
                      help='Tecplot output of Tacode (default: %(default)s)')
  parser.add_argument('-o', '--output', default='tecplot_output.dat',
                      help='file to write (default: %(default)s)')
  parser.add_argument('--control', default='config.yml',
                      help='control file of the run, for atmosphere.specific_heat_ratio, '
                           'atmosphere.molecular_weight and satellite.characteristic_length '
                           '(default: %(default)s; a missing file is not an error)')
  parser.add_argument('--gas', default='co2', choices=sorted(GAS),
                      help='gas for the defaults below and for the viscosity (default: %(default)s)')
  parser.add_argument('--specific-heat-ratio', type=float, default=None,
                      help='ratio of specific heats for the speed of sound (default: from the control '
                           'file, else the gas)')
  parser.add_argument('--molecular-weight', type=float, default=None,
                      help='mean molecular weight, g/mol (default: from the control file, else the gas)')
  parser.add_argument('--length', type=float, default=None,
                      help='length of the Reynolds number, m (default: satellite.characteristic_length)')
  parser.add_argument('--sutherland', type=float, nargs=3, default=None, metavar=('MU0', 'T0', 'S'),
                      help='Sutherland law of the viscosity, Pa s, K, K (default: that of the gas)')
  parser.add_argument('--nose-radius', type=float, default=None,
                      help='nose radius for the heat flux, m (default: {:g}, with a caution, since '
                           'the heat flux depends on it)'.format(NOSE_RADIUS_DEFAULT))
  parser.add_argument('--heating', default='sutton-graves', choices=HEATING,
                      help='stagnation-point heat flux: sutton-graves (the coefficient of the gas), '
                           'tauber-mars, or custom (default: %(default)s)')
  parser.add_argument('--heat-flux-law', type=float, nargs=3, default=None, metavar=('C', 'A', 'B'),
                      help='with --heating custom: q = C rho^A V^B / sqrt(Rn) in W/m2')
  helper_config.add_argument(parser)
  return parser


def main():
  parser = argument()
  args = helper_config.get_setting(parser, NAME_SECTION)
  if helper_config.save_file(args, NAME_SECTION) :
    return

  if not os.path.exists(args.filename) :
    print('File not found: ', args.filename)
    sys.exit(1)

  value, source = resolve_setting(args)
  data = tecplot_reader.read_tecplot(args.filename)
  result = compute(data, value, value['nose_radius'])

  header = ['Flow properties along the trajectory of ' + os.path.abspath(args.filename),
            'Written by src_helper/flow_properties/flow_properties.py',
            'Speed of sound: specific heat ratio {:g} ({}), molecular weight {:g} g/mol ({})'.format(
              value['specific_heat_ratio'], source['specific_heat_ratio'],
              value['molecular_weight'], source['molecular_weight']),
            'Reynolds number: length {:g} m ({}), Sutherland mu0 {:g} Pa s, T0 {:g} K, S {:g} K ({})'.format(
              value['length'], source['length'], *value['sutherland'], args.gas),
            'Knudsen number: the column of the Tacode output (satellite.characteristic_length)',
            'Heat flux: {} with q = {:g} rho^{:g} V^{:g} / sqrt(Rn), Rn = {:g} m; heat load: its time integral'.format(
              args.heating, *value['heat_flux_law'], value['nose_radius']) +
              ('' if source['nose_radius'] == 'given' else ' (NOT GIVEN: a placeholder)'),
            'Velair is the air-relative speed (VelairAbs when the wind is on, VelplAbs otherwise)']
  for line in header[2:] :
    print(line)
  write_tecplot(args.output, result, header)
  summarise(result)
  print('Written: ', args.output)


if __name__ == '__main__':
  main()
