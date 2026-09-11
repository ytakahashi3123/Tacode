#!/usr/bin/env python3
#
# 大気テーブルを作る。
#
# Tacode のソルバー（src/）はテキストのテーブルしか読まない。大気モデルそのものを
# 呼ぶライブラリはここに閉じ込めてある（計算環境に依存を持ち込まないため。
# database/wind/generate_wind_table.py と同じ方針）。
#
# 出力は **CCMC(VITMO) 形式**（src/atmosphere/atmosphere.py が内容から判定する 2 形式の
# 片方）。1 点のプロファイル、あるいは複数点の平均を書く:
#
#   Selected parameters are:
#   1 Height, km
#   2 O, cm-3
#   ...
#   （空行）
#   （列番号の行）
#   （データ）
#
# 列は Height / O / N2 / O2 / Mass_density / Temperature_neutral / N の 7 つで、
# 同梱の atmospheremodel.txt と同じ並びにしてある。Knudsen 数は N2・O2・N・O の
# 数密度から solver 側が作るので、この 4 種を必ず書く。
#
# データ源:
#
#   msis  NRLMSISE-00（既定）または MSIS 2.1 を pymsis で呼ぶ。**pymsis はこの
#         スクリプトだけの依存**で、requirements.txt には入れない（使うときに
#         `pip install pymsis`）。version 0 が NRLMSISE-00、2.1 が MSIS 2.1
#
# 使い方:
#   # Apollo 4 の突入時刻・突入点の NRLMSISE-00 プロファイル（0-400 km、1 km 刻み）
#   python3 generate_atmosphere_table.py --datetime 1967-11-09T20:19:29Z \
#           --longitude 187.5 --latitude 27.0 --f107 150 --f107a 150 --ap 4 \
#           -o atmospheremodel_apollo4.txt
#
#   # 緯度経度を振って平均したプロファイル（atmospheremodel_700km.txt と同じ作り）
#   python3 generate_atmosphere_table.py --datetime 2024-01-01T00:00:00Z \
#           --longitude 0 315 --longitude-step 45 --latitude -50 50 --latitude-step 25 \
#           --altitude 0 700 -o atmospheremodel_mean.txt
#
# 注意:
#   - **F10.7 と Ap は外から与える**（既定は太陽活動中庸の 150 / 150 / 4）。
#     pymsis.utils.get_f107_ap は外部データの取得が要るので、ここでは呼ばない。
#     高度 100 km 以下の密度は F10.7 にほとんど依らないが、熱圏では効く
#   - O・N・NO は MSIS の下限（72.5 km 付近）より下で NaN が返る。テーブルには
#     0 を書く（同梱の CCMC 版テーブルと同じ扱い）
#   - 高度は**ジオメトリック高度**（MSIS の入力そのまま）。Tacode の測地高度と同じ扱いで読む

import argparse
import datetime
import sys

import numpy as np


# 出力する列（名前, 単位, 変換係数: MSIS の SI から単位へ）
# --MSIS は数密度 m-3、質量密度 kg/m3、温度 K で返す
COLUMN_LIST = [
  ('Height',              'km',      None),
  ('O',                   'cm-3',    1.0e-6),
  ('N2',                  'cm-3',    1.0e-6),
  ('O2',                  'cm-3',    1.0e-6),
  ('Mass_density',        'g/cm-3',  1.0e-3),
  ('Temperature_neutral', 'K',       1.0),
  ('N',                   'cm-3',    1.0e-6),
]

# pymsis.Variable の名前（COLUMN_LIST の並びに対応。Height は別扱い）
VARIABLE_LIST = ['O', 'N2', 'O2', 'MASS_DENSITY', 'TEMPERATURE', 'N']


def parse_datetime(string_time):
  # ISO 8601。末尾の Z は Python 3.11 未満の fromisoformat が読めないので置き換える
  # （src/epoch/epoch.py と同じ扱い）
  text = string_time.strip()
  if text.endswith('Z') or text.endswith('z') :
    text = text[:-1] + '+00:00'
  try:
    time_parsed = datetime.datetime.fromisoformat(text)
  except ValueError:
    print('Cannot read the date and time:', string_time)
    print('--Give it in ISO 8601, for example 1967-11-09T20:19:29Z.')
    print('Program stopped.')
    sys.exit(1)
  if time_parsed.tzinfo is None :
    time_parsed = time_parsed.replace(tzinfo=datetime.timezone.utc)
  return time_parsed


def make_axis(bound, step, name):
  # [下限, 上限] と刻みから節点を作る。1 点だけのときは下限をそのまま返す
  if len(bound) == 1 or bound[0] == bound[-1] :
    return np.array([bound[0]])
  if step <= 0.0 :
    print('The step of the', name, 'axis must be positive:', step)
    print('Program stopped.')
    sys.exit(1)
  num_point = int(round((bound[-1] - bound[0])/step)) + 1
  return bound[0] + step*np.arange(0, num_point)


def run_msis(time_table, longitude, latitude, altitude, f107, f107a, ap, version):
  # pymsis は「使うときだけ」import する（このスクリプトの実行時にしか要らない）
  try:
    import pymsis
  except ImportError:
    print('pymsis is not installed; it is needed to generate a table with --source msis.')
    print('--Install it with: pip install pymsis')
    print('Program stopped.')
    sys.exit(1)

  # aps は [Ap, 3hr x 6] の 7 要素。日平均だけを与える使い方に合わせて同じ値を並べる
  aps = [[ap]*7]
  # pymsis は datetime64 に落とすのでタイムゾーンを持てない。UTC に直して外す
  time_naive = time_table.astimezone(datetime.timezone.utc).replace(tzinfo=None)
  output = pymsis.calculate(time_naive, longitude, latitude, altitude,
                            f107, f107a, aps, version=version)
  # 形は (時刻, 経度, 緯度, 高度, 変数)
  return np.squeeze(output, axis=0)


def gather_profile(output, variable_index):
  # 経度・緯度について平均し、高度のプロファイルにする
  # --MSIS は下限より下で O / N / NO を NaN にする。テーブルには 0 を書く
  profile = {}
  for name in VARIABLE_LIST:
    data = output[:, :, :, variable_index[name]]
    data = np.where(np.isfinite(data), data, 0.0)
    profile[name] = np.mean(data, axis=(0, 1))
  return profile


def write_table(filename, header_line, altitude, profile):
  with open(filename, 'w') as f:
    for line in header_line:
      f.write(line + '\n')
    f.write('\n')
    f.write('   Selected parameters are:\n')
    for index, (name, unit, factor) in enumerate(COLUMN_LIST):
      f.write('{:d} {:s}, {:s}\n'.format(index+1, name, unit))
    f.write('\n')
    f.write(''.join(['{:>11d}'.format(index+1) for index in range(0, len(COLUMN_LIST))]) + '\n')

    for n in range(0, len(altitude)):
      value_list = [altitude[n]]
      for name, unit, factor in COLUMN_LIST[1:]:
        key = 'MASS_DENSITY' if name == 'Mass_density' else \
              'TEMPERATURE'  if name == 'Temperature_neutral' else name
        value_list.append(profile[key][n]*factor)
      f.write('{:11.1f}'.format(value_list[0]))
      f.write(''.join(['{:11.3e}'.format(value) for value in value_list[1:]]) + '\n')

  print('Atmosphere table written:', filename)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--source', type=str, default='msis', choices=['msis'],
                      help='Data source (only the MSIS family is implemented)')
  parser.add_argument('-o', '--output', type=str, default='atmospheremodel.txt')
  parser.add_argument('--datetime', type=str, required=True,
                      help='UTC of the profile in ISO 8601, e.g. 1967-11-09T20:19:29Z')
  parser.add_argument('--longitude', type=float, nargs='+', default=[0.0],
                      help='Longitude [deg. East]; one value, or a lower and an upper bound')
  parser.add_argument('--latitude', type=float, nargs='+', default=[0.0],
                      help='Latitude [deg. North]; one value, or a lower and an upper bound')
  parser.add_argument('--longitude-step', type=float, default=45.0)
  parser.add_argument('--latitude-step', type=float, default=25.0)
  parser.add_argument('--altitude', type=float, nargs=2, default=[0.0, 400.0],
                      help='Altitude range [km]')
  parser.add_argument('--altitude-step', type=float, default=1.0)
  parser.add_argument('--f107', type=float, default=150.0,
                      help='F10.7 of the previous day')
  parser.add_argument('--f107a', type=float, default=150.0,
                      help='81-day averaged F10.7')
  parser.add_argument('--ap', type=float, default=4.0,
                      help='Daily Ap index')
  parser.add_argument('--version', type=float, default=0,
                      help='MSIS version: 0 for NRLMSISE-00 (default), 2.1 for MSIS 2.1')
  args = parser.parse_args()

  time_profile = parse_datetime(args.datetime)

  longitude = make_axis(args.longitude, args.longitude_step, 'longitude')
  latitude  = make_axis(args.latitude,  args.latitude_step,  'latitude')
  altitude  = make_axis(args.altitude,  args.altitude_step,  'altitude')

  print('Calling MSIS (version {:})...'.format(args.version))
  print('--Date and time :', time_profile.isoformat())
  print('--Longitude     : {:} point(s), {:} to {:} deg. East'.format(len(longitude), longitude[0], longitude[-1]))
  print('--Latitude      : {:} point(s), {:} to {:} deg. North'.format(len(latitude), latitude[0], latitude[-1]))
  print('--Altitude      : {:} point(s), {:} to {:} km'.format(len(altitude), altitude[0], altitude[-1]))
  print('--F10.7, F10.7A, Ap: {:}, {:}, {:}'.format(args.f107, args.f107a, args.ap))

  import pymsis
  output = run_msis(time_profile, longitude, latitude, altitude,
                    args.f107, args.f107a, args.ap, args.version)
  variable_index = {name: pymsis.Variable[name].value for name in VARIABLE_LIST}
  profile = gather_profile(output, variable_index)

  name_model = 'NRLMSISE-00' if float(args.version) == 0.0 else 'MSIS {:}'.format(args.version)
  header_line = [
    'Generated by database/atmosphere/generate_atmosphere_table.py (pymsis)',
    '',
    'Model: ' + name_model,
    'Date and time (UTC): ' + time_profile.isoformat(),
    'Longitude [deg. East]: {:} point(s), {:} to {:}'.format(len(longitude), longitude[0], longitude[-1]),
    'Latitude [deg. North]: {:} point(s), {:} to {:}'.format(len(latitude), latitude[0], latitude[-1]),
    'F10.7 = {:}, F10.7A = {:}, Ap = {:}'.format(args.f107, args.f107a, args.ap),
    'The species below the lower bound of the model are written as zero.',
  ]
  write_table(args.output, header_line, altitude, profile)

  return


if __name__ == '__main__':
  main()
