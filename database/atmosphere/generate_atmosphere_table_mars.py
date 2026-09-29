#!/usr/bin/env python3
#
# 火星の大気テーブルを作る（Mars Climate Database の Web 版から）。
#
# Tacode のソルバー（src/）はテキストのテーブルしか読まない。地球の
# generate_atmosphere_table.py が pymsis を閉じ込めているのと同じく、火星の大気モデルを
# 呼ぶ部分はここに閉じ込めてある。**依存は標準ライブラリだけ**だが、**ネットワークが要る**
# （作ったテーブルはリポジトリに置くので、ソルバーの実行には要らない）。
#
# データ源:
#
#   Mars Climate Database (MCD) v6 の Web 版 https://www-mars.lmd.jussieu.fr/mcd_python/
#   LMD/OU/IAA/ESA/CNES。利用条件は https://www-mars.lmd.jussieu.fr/mars/access.html:
#   科学利用は自由（出典を正しく引用し、利用を MCD チームに知らせること）、商用利用は
#   許可が要る、無保証。**Web 版は「適度な利用」が条件**（大量・精密な用途には、
#   登録して配布版を使うよう求めている）。1 本のテーブルで 3 x (高度 35 点ごと) 回
#   問い合わせる（0-200 km・1 km 刻みで 18 回）ので、1 回ごとに間を置く。
#   引用: Millour et al., The Mars Climate Database (version 6.1), EPSC 2022;
#         Forget et al., J. Geophys. Res. 104, 24155 (1999)
#
# **作ったテーブルはリポジトリに入れない**（.gitignore 済み）。MCD のデータを再配布しない
# ためで、使う人がそれぞれこのスクリプトで取得する（tutorial/work_reentry_mars/run_tacode.sh が
# テーブルの無いときに取得のコマンドを示す）。
#
# 出力は **CCMC(VITMO) 形式**（src/atmosphere/atmosphere.py が内容から判定する 2 形式の
# 片方）。列は
#
#   Height / CO2 / N2 / Ar / CO / O / O2 / Mass_density / Temperature_neutral
#
# 数密度は MCD の圧力・温度・体積混合比から n_i = x_i p/(k T) で作る（MCD は数密度を
# 直接は返さない）。Mass_density は MCD の密度そのもの。両者の食い違い
# （書いた種の sum(n_i m_i) と密度の比）はヘッダに書く。Knudsen 数は solver 側が
# これらの数密度から作る。
#
# 高度の基準:
#
#   **半径 3396.0 km の球からの高さ**（MCD の zkey = 5）。Tacode の高度は
#   planet.radius と planet.ellipticity で決まる楕円体からの高さなので、
#   **config は radius: 3396.0e+3 / ellipticity: 0.0 にすること**（そうすれば両者は
#   同じ量になる。緯度も火星の標準である planetocentric になる）。
#   MOLA の標高や MCD の「sea level」は areoid 基準で、これとは km の桁でずれる。
#
# 使い方:
#   # Mars Pathfinder の着陸点・着陸時の季節と地方時（tutorial/work_reentry_mars が使うもの）。
#   # ケースのディレクトリで:
#   python3 ../../database/atmosphere/generate_atmosphere_table_mars.py --latitude 19.13 \
#           --longitude -33.22 --ls 142.7 --local-time 3.0 --altitude 0 200 \
#           -o database/atmosphere/atmospheremodel_mars.txt
#
# 注意:
#   - MCD は 1 回の問い合わせで 1 つの高度範囲を 35 点に等分して返す。ここでは
#     35 点ずつ区切って刻みを --altitude-step にそろえる
#   - 地表より下は MCD が NaN を返す。**範囲の中に NaN があれば止まる**
#     （地表が 3396 km の球より上にある場所では、下端を地表より上に上げること）
#   - 高度 1 点の値は 1 地点・1 時刻の気候値（平均的な年の日変化の中の 1 点）で、
#     日々の変動（MCD の rmsrho）は含まない

import argparse
import math
import os
import sys
import time
import urllib.parse
import urllib.request


URL_MCD = 'https://www-mars.lmd.jussieu.fr/mcd_python/cgi-bin/mcdcgi.py'
URL_TEXT = 'https://www-mars.lmd.jussieu.fr/mcd_python/'

# MCD は 1 回に 4 変数までしか返さない
VARIABLE_GROUP_LIST = [
  ['rho', 't', 'p', 'vmr_co2'],
  ['vmr_n2', 'vmr_ar', 'vmr_co', 'vmr_o'],
  ['vmr_o2', 'vmr_he', 'vmr_h', 'vmr_h2'],
]

# 返ってくるテキストの「Column 2 is ...」の書き出し（変数の取り違えを検査する）
DICT_LABEL = {
  'rho'    : 'Density',
  't'      : 'Temperature',
  'p'      : 'Pressure',
  'vmr_co2': '[CO2] volume mixing ratio',
  'vmr_n2' : '[N2] volume mixing ratio',
  'vmr_ar' : '[Ar] volume mixing ratio',
  'vmr_co' : '[CO] volume mixing ratio',
  'vmr_o'  : '[O] volume mixing ratio',
  'vmr_o2' : '[O2] volume mixing ratio',
  'vmr_he' : '[He] volume mixing ratio',
  'vmr_h'  : '[H] volume mixing ratio',
  'vmr_h2' : '[H2] volume mixing ratio',
}

# テーブルに書く種（名前は src/atmosphere/atmosphere.py の KEY_* と一致させる）と
# MCD の変数、分子量。He / H / H2 は Kn の和に入れないので書かない（200 km までの
# 混合比は下のヘッダに書く「書いていない種の割合」で分かる）
SPECIES_LIST = [
  ('CO2', 'vmr_co2', 44.0095),
  ('N2',  'vmr_n2',  28.0134),
  ('Ar',  'vmr_ar',  39.948),
  ('CO',  'vmr_co',  28.0101),
  ('O',   'vmr_o',   15.9994),
  ('O2',  'vmr_o2',  31.9988),
]

# MCD は 1 回の問い合わせで範囲を 35 点に等分する
NUM_POINT_QUERY = 35

CONSTANT_BOLTZMANN = 1.380649e-23    # J/K
MASS_ATOMIC        = 1.66053907e-27  # kg

# 返ってきた高度が頼んだ格子から外れていないかの許容幅, m
TOLERANCE_ALTITUDE = 1.0

# 問い合わせの間隔, s（Web 版は適度な利用が条件）
INTERVAL_QUERY = 1.0

# 利用条件（https://www-mars.lmd.jussieu.fr/mars/access.html）。テーブルのヘッダに書く。
# 作ったテーブルは MCD の派生データなので、Tacode の MIT ライセンスの外にある
LINE_TERMS = ('Terms of the MCD: scientific use is free provided the origin of the data is quoted and the MCD team '
              'is kept informed; no commercial use without their authorization; no warranty. '
              "This file is not covered by Tacode's MIT license.")


def make_axis(bound, step):
  if step <= 0.0 or bound[1] <= bound[0] :
    print('The altitude range must increase and the step must be positive:', bound, step)
    print('Program stopped.')
    sys.exit(1)
  num_point = int(round((bound[1] - bound[0])/step)) + 1
  return [bound[0] + step*index for index in range(0, num_point)]


def query_mcd(variables, altitude_lower, altitude_upper, args):
  # 1 回の問い合わせ。返り値は {変数: [(高度 m, 値), ...]} と、先頭の見出し行
  parameter = {
    'var1': variables[0], 'var2': variables[1], 'var3': variables[2], 'var4': variables[3],
    'datekeyhtml': '1', 'ls': '{:g}'.format(args.ls), 'localtime': '{:g}'.format(args.local_time),
    'latitude': '{:g}'.format(args.latitude), 'longitude': '{:g}'.format(args.longitude),
    'altitude': '{:.1f} {:.1f}'.format(altitude_lower, altitude_upper),
    'zkey': '5', 'isfixedlt': 'off', 'dust': str(args.dust), 'hrkey': '1',
    'averaging': args.averaging, 'dpi': '80', 'islog': 'off', 'colorm': 'jet', 'proj': 'cyl',
    'iswind': 'off',
  }
  url = URL_MCD + '?' + urllib.parse.urlencode(parameter)
  page = fetch(url)

  # 結果のページはテキストファイル（../txt/<hash>.txt）へのリンクを持つ
  marker = '../txt/'
  index = page.find(marker)
  if index < 0 :
    print('The Mars Climate Database did not return a result.')
    print('--URL:', url)
    print('Program stopped.')
    sys.exit(1)
  index_end = page.find('.txt', index)
  text = fetch(URL_TEXT + page[index+3:index_end+4])

  return parse_text(text, variables)


def fetch(url):
  request = urllib.request.Request(url, headers={'User-Agent': 'Tacode atmosphere table generator'})
  with urllib.request.urlopen(request, timeout=300) as response:
    return response.read().decode('utf-8', errors='replace')


def parse_text(text, variables):
  # テキストは変数ごとのブロックが並ぶ:
  #   ### MCD_v6.2 with climatology average solar scenario.
  #   ### Ls 142.7deg. Latitude 19.13N. Longitude -33.22E. Local time 3.0h
  #   ### Column 1 is altitude above mean Mars radius (m)
  #   ### Column 2 is Density (kg/m3)
  #   ...
  #   <高度> <値>
  block_list = []
  heading = []
  flag_data = False
  for line in text.splitlines():
    stripped = line.strip()
    if stripped.startswith('### Column 2 is') :
      block_list.append({'label': stripped[len('### Column 2 is'):].strip(), 'data': []})
      continue
    if stripped.startswith('#') :
      # 最初のデータ行より前の見出し（版・条件・取得日時）だけを残す。罫線は落とす
      content = stripped.strip('#').strip()
      if not flag_data and content != '' and not content.startswith('-') :
        heading.append(content)
      continue
    words = stripped.split()
    if len(words) != 2 or len(block_list) == 0 :
      continue
    block_list[-1]['data'].append((float(words[0]), float(words[1])))
    flag_data = True

  if len(block_list) != len(variables) :
    print('Unexpected answer from the Mars Climate Database:', len(block_list), 'blocks for', variables)
    print('Program stopped.')
    sys.exit(1)

  result = {}
  for variable, block in zip(variables, block_list):
    if not block['label'].startswith(DICT_LABEL[variable]) :
      print('The Mars Climate Database returned "' + block['label'] + '" where', variable, 'was asked for.')
      print('Program stopped.')
      sys.exit(1)
    result[variable] = block['data']

  return result, heading


def gather_profile(altitude_list, args):
  # 高度 35 点ごと・変数 4 つごとに問い合わせ、変数ごとの高度プロファイルにまとめる
  step = args.altitude_step
  profile = {variable: [] for group in VARIABLE_GROUP_LIST for variable in group}
  heading = None

  for index_start in range(0, len(altitude_list), NUM_POINT_QUERY):
    altitude_chunk = altitude_list[index_start:index_start+NUM_POINT_QUERY]
    # MCD は範囲を常に 35 点に割るので、上端は下端 + 34 刻み（余りは捨てる）
    altitude_lower = altitude_chunk[0]*1.e3
    altitude_upper = (altitude_chunk[0] + step*(NUM_POINT_QUERY-1))*1.e3

    for group in VARIABLE_GROUP_LIST:
      print('--Asking MCD for {:s} at {:.1f}-{:.1f} km'.format(
            ', '.join(group), altitude_lower*1.e-3, altitude_upper*1.e-3))
      result, heading_tmp = query_mcd(group, altitude_lower, altitude_upper, args)
      if heading is None :
        heading = heading_tmp
      for variable in group:
        data = result[variable]
        if len(data) != NUM_POINT_QUERY :
          print('The Mars Climate Database returned', len(data), 'points instead of', NUM_POINT_QUERY)
          print('Program stopped.')
          sys.exit(1)
        for index, altitude in enumerate(altitude_chunk):
          altitude_returned, value = data[index]
          if abs(altitude_returned - altitude*1.e3) > TOLERANCE_ALTITUDE :
            print('The altitude grid of the answer does not match the request:',
                  altitude_returned, 'm against', altitude*1.e3, 'm')
            print('Program stopped.')
            sys.exit(1)
          profile[variable].append(value)
      time.sleep(INTERVAL_QUERY)

  return profile, heading


def check_profile(altitude_list, profile):
  # 地表より下（MCD が NaN を返す）や負の値が混じっていないこと
  for variable, values in profile.items():
    for altitude, value in zip(altitude_list, values):
      if not math.isfinite(value) :
        print('The Mars Climate Database gives no value of', variable, 'at', altitude, 'km.')
        print('--This is below the local surface. Raise the lower bound of --altitude.')
        print('Program stopped.')
        sys.exit(1)
      if value < 0.0 :
        print('A negative value of', variable, 'at', altitude, 'km:', value)
        print('Program stopped.')
        sys.exit(1)


def make_number_density(profile):
  # n_i = x_i p/(k T), m-3
  number_density = {}
  for name, variable, mass_molecular in SPECIES_LIST:
    number_density[name] = [fraction*pressure/(CONSTANT_BOLTZMANN*temperature)
                            for fraction, pressure, temperature
                            in zip(profile[variable], profile['p'], profile['t'])]
  return number_density


def summarise_closure(altitude_list, profile, number_density):
  # 書いた種の sum(n_i m_i) と MCD の密度の比、書いていない種の混合比
  ratio_list = []
  fraction_missing_list = []
  for index in range(0, len(altitude_list)):
    density_species = sum(number_density[name][index]*mass_molecular*MASS_ATOMIC
                          for name, variable, mass_molecular in SPECIES_LIST)
    ratio_list.append(density_species/profile['rho'][index])
    fraction_written = sum(profile[variable][index] for name, variable, mass_molecular in SPECIES_LIST)
    fraction_missing_list.append(1.0 - fraction_written)
  return ratio_list, fraction_missing_list


def write_table(filename, header_line, altitude_list, profile, number_density):
  column_list = [('Height', 'km')] + [(name, 'cm-3') for name, variable, mass in SPECIES_LIST] \
              + [('Mass_density', 'g/cm-3'), ('Temperature_neutral', 'K')]

  directory = os.path.dirname(filename)
  if directory != '' :
    os.makedirs(directory, exist_ok=True)

  with open(filename, 'w') as f:
    for line in header_line:
      f.write(line + '\n')
    f.write('\n')
    f.write('   Selected parameters are:\n')
    for index, (name, unit) in enumerate(column_list):
      f.write('{:d} {:s}, {:s}\n'.format(index+1, name, unit))
    f.write('\n')
    f.write(''.join(['{:>11d}'.format(index+1) for index in range(0, len(column_list))]) + '\n')

    for index, altitude in enumerate(altitude_list):
      value_list = [number_density[name][index]*1.e-6 for name, variable, mass in SPECIES_LIST]
      value_list.append(profile['rho'][index]*1.e-3)
      value_list.append(profile['t'][index])
      f.write('{:11.1f}'.format(altitude))
      f.write(''.join(['{:11.3e}'.format(value) for value in value_list]) + '\n')

  print('Atmosphere table written:', filename)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-o', '--output', type=str, default='atmospheremodel_mars.txt')
  parser.add_argument('--latitude', type=float, required=True,
                      help='Planetocentric latitude [deg. North]')
  parser.add_argument('--longitude', type=float, required=True,
                      help='Longitude [deg. East]')
  parser.add_argument('--ls', type=float, required=True,
                      help='Solar longitude Ls [deg.] (the Martian season)')
  parser.add_argument('--local-time', type=float, required=True,
                      help='Local true solar time at the point [Martian hours]')
  parser.add_argument('--dust', type=int, default=1,
                      help='MCD dust and solar scenario (1: climatology, average solar; default)')
  parser.add_argument('--averaging', type=str, default='off', choices=['off', 'loct'],
                      help='off: the local time given by --local-time (default);'
                           ' loct: the mean over all local times (diurnal mean; --local-time is ignored)')
  parser.add_argument('--altitude', type=float, nargs=2, default=[0.0, 200.0],
                      help='Altitude range above the 3396 km sphere [km]')
  parser.add_argument('--altitude-step', type=float, default=1.0)
  args = parser.parse_args()

  altitude_list = make_axis(args.altitude, args.altitude_step)

  print('Calling the Mars Climate Database (web interface)...')
  print('--Point      : {:g} deg. N, {:g} deg. E'.format(args.latitude, args.longitude))
  print('--Ls, LT     : {:g} deg., {:g} h'.format(args.ls, args.local_time))
  print('--Scenario   :', args.dust)
  print('--Altitude   : {:} point(s), {:} to {:} km'.format(len(altitude_list), altitude_list[0], altitude_list[-1]))

  profile, heading = gather_profile(altitude_list, args)
  check_profile(altitude_list, profile)
  number_density = make_number_density(profile)
  ratio_list, fraction_missing_list = summarise_closure(altitude_list, profile, number_density)

  print('--sum(n_i m_i)/rho      : {:.4f} to {:.4f}'.format(min(ratio_list), max(ratio_list)))
  print('--Species not written   : up to {:.2e} of the molecules'.format(max(fraction_missing_list)))

  header_line = [
    'Generated by database/atmosphere/generate_atmosphere_table_mars.py',
    '',
    'Model: Mars Climate Database, web interface (' + URL_TEXT + ')',
  ] + ['MCD: ' + line for line in heading] + [
    'Latitude (planetocentric) [deg. North]: {:g}'.format(args.latitude),
    'Longitude [deg. East]: {:g}'.format(args.longitude),
    ('Ls [deg.]: {:g}, local true solar time [h]: {:g}, scenario: {:d}'.format(args.ls, args.local_time, args.dust)
     if args.averaging == 'off' else
     'Ls [deg.]: {:g}, diurnal mean over all local times, scenario: {:d}'.format(args.ls, args.dust)),
    'Height is above the sphere of radius 3396.0 km (planet.radius 3396.0e+3, ellipticity 0.0).',
    'Number densities are x_i p/(k T) from the MCD pressure, temperature and volume mixing ratios.',
    'sum(n_i m_i) of the species written / MCD density: {:.4f} to {:.4f}'.format(min(ratio_list), max(ratio_list)),
    'Molecules of the species not written (He, H, H2, H2O, ...): up to {:.2e}'.format(max(fraction_missing_list)),
    'Mars Climate Database (c) LMD/OU/IAA/ESA/CNES; Millour et al., EPSC 2022; Forget et al., JGR 104, 24155 (1999)',
    LINE_TERMS,
  ]
  write_table(args.output, header_line, altitude_list, profile, number_density)

  return


if __name__ == '__main__':
  main()
