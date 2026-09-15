#!/usr/bin/env python3

import numpy as np
import os as os
import sys as sys
import scipy.interpolate

import re

from general.general import get_database_directory

# 単位変換
# ’read_atmosphere_file’でKey errorがあったときは下記の単位変換が正しいかチェックする
unit_convert ={'km':1.0, 'cm-3':1.e6, 'g/cm-3': 1.e3, 'K': 1.0}

# Dict key
KEY_DATA   = 'Number_Data'
KEY_ATM    = 'Number_Atmosphere'
KEY_KN     = 'Knudsen_number'
KEY_INTERP = 'Interpolator'
KEY_SCALE_HEIGHT = 'Scale_height'

KEY_Height = 'Height'
KEY_N2     = 'N2'
KEY_O2     = 'O2'
KEY_N      = 'N'
KEY_O      = 'O'
KEY_Mass_density = 'Mass_density'
KEY_Temperature_neutral = 'Temperature_neutral'

# テーブル上端でスケールハイトをフィットする高度幅, km
RANGE_FIT_SCALE_HEIGHT = 50.0

# 指数外挿の指数の上限（スケールハイトの何個分まで外挿するか）
# --exp(-100) = 3.7e-44 で密度は事実上 0、Kn は 2.7e43 倍。ここから先は
#   真空・自由分子流の極限で値が変わらないので、同じ値で扱う。
#   上限を置かないと exp がアンダーフローして Kn が Inf になり、出力に Inf が載る
EXPONENT_MAXIMUM = 100.0

# テーブル範囲外の扱い
KIND_EXTRAPOLATION_EXPONENTIAL = 'exponential'
KIND_EXTRAPOLATION_CLAMP       = 'clamp'

# 大気モデルファイルの形式
# --ccmc:    CCMC(VITMO) の Web から取得したもの。"Selected parameters are:" ブロックを持つ
# --fortran: NRLMSISE-00 の Fortran 版 (NRLMSISE-00_readctl.FOR) が直接出力したもの
KIND_FILE_CCMC    = 'ccmc'
KIND_FILE_FORTRAN = 'fortran'

MARKER_CCMC    = 'Selected parameters are:'
MARKER_FORTRAN = 'ALTITUDE'

# Fortran 版の見出し行（NRLMSISE-00_readctl.FOR にリテラルで埋め込まれている）から
# 変数名と単位への対応。見出しは "<番号>, <D(n) - >?<説明>(<単位>)" の形。
DICT_FORTRAN_PARAMETER = {
  'ALTITUDE'           : (KEY_Height,              'km'),
  'O NUMBER DENSITY'   : (KEY_O,                   'cm-3'),
  'N2 NUMBER DENSITY'  : (KEY_N2,                  'cm-3'),
  'O2 NUMBER DENSITY'  : (KEY_O2,                  'cm-3'),
  'N NUMBER DENSITY'   : (KEY_N,                   'cm-3'),
  'TOTAL MASS DENSITY' : (KEY_Mass_density,        'g/cm-3'),
  'TEMPERATURE AT ALT' : (KEY_Temperature_neutral, 'K'),
}

LIST_MOLECULAR_KIND = [KEY_N2, KEY_O2, KEY_N, KEY_O]
DICT_DIAMETER_MOLECULAR = {KEY_N2: 3.75e-10, KEY_O2: 3.54e-10, KEY_N: 3.10e-10, KEY_O: 3.04e-10}


def initial_settings_atmosphere(config):

  atmosphere_dict = read_atmosphere_file(config)

  atmosphere_dict = set_knudsen_number(config, atmosphere_dict)

  atmosphere_dict = set_interpolator(atmosphere_dict)

  atmosphere_dict = set_scale_height(config, atmosphere_dict)

  return atmosphere_dict


def set_scale_height(config, atmosphere_dict):
  # テーブル上端より上での外挿に使うスケールハイトを、テーブル自身から求める。
  # 熱圏上部は等温かつ拡散平衡なので密度は高度に対して指数関数で減る:
  #   rho(z) = rho_top * exp( -(z - z_top)/H )
  # H を定数として埋め込まず毎回フィットするので、テーブルを差し替えても追従する。

  kind_extrapolation = config['atmosphere'].get('kind_extrapolation', KIND_EXTRAPOLATION_EXPONENTIAL)

  if kind_extrapolation == KIND_EXTRAPOLATION_CLAMP :
    atmosphere_dict[KEY_SCALE_HEIGHT] = None
    return atmosphere_dict

  if kind_extrapolation != KIND_EXTRAPOLATION_EXPONENTIAL :
    print('kind_extrapolation in config is incorrect:', kind_extrapolation)
    print('--Available: ' + KIND_EXTRAPOLATION_EXPONENTIAL + ', ' + KIND_EXTRAPOLATION_CLAMP)
    print('Program stopped.')
    sys.exit(1)

  altitude_atm = atmosphere_dict[KEY_Height]
  density_atm  = atmosphere_dict[KEY_Mass_density]

  mask = altitude_atm >= altitude_atm[-1] - RANGE_FIT_SCALE_HEIGHT
  if np.count_nonzero(mask) < 2 :
    mask = np.ones(len(altitude_atm), dtype=bool)

  slope = np.polyfit(altitude_atm[mask], np.log(density_atm[mask]), 1)[0]

  if slope >= 0.0 :
    # 上端で密度が減っていないテーブルでは指数外挿が使えない
    print('Warning: the density does not decrease at the top of the atmosphere table.')
    print('--Values above the table are clamped instead of extrapolated.')
    atmosphere_dict[KEY_SCALE_HEIGHT] = None
    return atmosphere_dict

  scale_height = -1.0/slope
  print('Setting scale height for extrapolation...: {:.2f} km above {:.1f} km'.format(
        scale_height, altitude_atm[-1]))

  atmosphere_dict[KEY_SCALE_HEIGHT] = scale_height

  return atmosphere_dict


def set_interpolator(atmosphere_dict):
  # 大気テーブルの補間器を初期化時に 1 度だけ構築する。
  # 毎ステップ構築すると RK4 の段数だけスプライン生成が走り、計算時間を支配してしまう。

  print('Setting atmosphere interpolator...')

  altitude_atm = atmosphere_dict[KEY_Height]

  # scipy.interpolate.interp1d は SciPy 1.10 以降 legacy 扱いなので使わない。
  # interp1d(kind='cubic') は内部で make_interp_spline(k=3) を作っていたので、
  # これを直接呼ぶと係数も評価経路も同じで、値はビット単位で変わらない
  # （CubicSpline は同じ not-a-knot でも評価が PPoly になり 4e-16 ずれる）。
  interpolator = {}
  for key_tmp in [KEY_Mass_density, KEY_Temperature_neutral, KEY_KN]:
    interpolator[key_tmp] = scipy.interpolate.make_interp_spline(altitude_atm, atmosphere_dict[key_tmp], k=3)

  atmosphere_dict[KEY_INTERP] = interpolator

  return atmosphere_dict


def read_atmosphere_file(config):

  directory_path = get_database_directory(config['atmosphere'], 'atmosphere', 'atmosphere', 'directory_atmosphere')

  filename_tmp = directory_path + '/' + config['atmosphere']['filename_atmosphere']
  print('Reading atmosphere model...:', filename_tmp)

  # File open
  with open(filename_tmp) as f:
    lines = f.readlines()
  # リストとして取得 
  lines_strip = [line.strip() for line in lines]

  kind_file = detect_atmosphere_file_kind(lines_strip, filename_tmp)
  print('--File format:', kind_file)

  if kind_file == KIND_FILE_CCMC :
    atmosphere_name, atmosphere_unit, index_data = read_header_ccmc(lines_strip)
  else :
    atmosphere_name, atmosphere_unit, index_data = read_header_fortran(lines_strip)

  num_atmosphere_var = len(atmosphere_name)

  print('Atmosphere name, unit')
  for i in range(0,num_atmosphere_var):
    print(i,atmosphere_name[i],',', atmosphere_unit[i])

  # Check unit convert  
  for i in range(0,num_atmosphere_var):
    try:
      unit_convert[atmosphere_unit[i]]
    except KeyError as instance:
      print(instance)
      print('Key error for unit convert. Check module: atmosphere.py, or atmosphere file')
      print('Program stopped.')
      sys.exit(1)

  # 大気モデルデータの取得
  # --数値が変数の数だけ並んでいる行のみを採る。空行や末尾のゴミ行があっても落ちない
  data_rows = []
  for i in range(index_data, len(lines_strip) ):
    words = lines_strip[i].split()
    if len(words) != num_atmosphere_var :
      continue
    try:
      data_rows.append( [float(word) for word in words] )
    except ValueError:
      continue

  if len(data_rows) == 0 :
    print('There is no numerical data in the atmosphere model file:', filename_tmp)
    print('Program stopped.')
    sys.exit(1)

  num_array_tmp    = len(data_rows)
  atmosphere_model = np.array(data_rows).transpose()
  for j in range(0,num_atmosphere_var):
    atmosphere_model[j,:] = atmosphere_model[j,:]*unit_convert[atmosphere_unit[j]]

  atmosphere_dict = {KEY_DATA:num_array_tmp, KEY_ATM:num_atmosphere_var}
  for i in range(0,num_atmosphere_var):
    atmosphere_dict[ atmosphere_name[i] ] = atmosphere_model[i]

  return atmosphere_dict


def detect_atmosphere_file_kind(lines_strip, filename_tmp):
  # 大気モデルファイルの形式を見出しから判定する

  for line in lines_strip:
    if line == MARKER_CCMC :
      return KIND_FILE_CCMC

  for line in lines_strip:
    if line.startswith('1,') and MARKER_FORTRAN in line :
      return KIND_FILE_FORTRAN

  print('Cannot recognize the format of the atmosphere model file:', filename_tmp)
  print('--Expected either "' + MARKER_CCMC + '" (CCMC/VITMO) or "1, ' + MARKER_FORTRAN + ' ..." (NRLMSISE-00 Fortran).')
  print('Program stopped.')
  sys.exit(1)


def read_header_ccmc(lines_strip):
  # CCMC(VITMO) 形式:
  #   Selected parameters are:
  #   1 Height, km
  #   2 O, cm-3
  #   ...
  #   （空行）
  #   （データ）

  atmosphere_name = []
  atmosphere_unit = []

  index_start = None
  for i in range(0, len(lines_strip) ):
    if lines_strip[i] == MARKER_CCMC :
      index_start = i+1
      break

  if index_start is None :
    print('There is no specific word: "' + MARKER_CCMC + '", Check atmosphere model file.')
    print('Program stopped.')
    sys.exit(1)

  for i in range(index_start, len(lines_strip) ):
    # --ブランク行が来たらループから出る
    if lines_strip[i] == '':
      break
    # --空白で分割
    words = lines_strip[i].split()
    # --コンマ削除
    words = [n.replace(",","") for n in words]
    atmosphere_name.append( words[1] )
    atmosphere_unit.append( words[2] )

  index_data = index_start + len(atmosphere_name) + 2

  return atmosphere_name, atmosphere_unit, index_data


def read_header_fortran(lines_strip):
  # NRLMSISE-00 の Fortran 版が出す形式:
  #   1, ALTITUDE (KM)
  #   2, D(2) - O NUMBER DENSITY(CM-3)
  #   ...
  #   1  2  3  4  5  6  7      <- 列番号の行
  #   （データ）
  #
  # 見出しは NRLMSISE-00_readctl.FOR にリテラルで書かれているので、
  # 説明文から DICT_FORTRAN_PARAMETER で変数名・単位に対応づける。

  atmosphere_name = []
  atmosphere_unit = []
  index_last = None

  for i in range(0, len(lines_strip) ):
    line = lines_strip[i]
    if not re.match(r'^\d+\s*,', line) :
      continue

    # "<番号>," を落とし、"D(n) - " があればその後ろを採る
    description = line.split(',', 1)[1].strip()
    if ' - ' in description :
      description = description.split(' - ', 1)[1].strip()
    # 末尾の "(単位)" を落とす
    description = re.sub(r'\s*\([^()]*\)\s*$', '', description).strip().upper()

    try:
      name_tmp, unit_tmp = DICT_FORTRAN_PARAMETER[description]
    except KeyError:
      print('Unknown parameter in the atmosphere model file:', line)
      print('--Add it to DICT_FORTRAN_PARAMETER in atmosphere.py if it is needed.')
      print('Program stopped.')
      sys.exit(1)

    atmosphere_name.append( name_tmp )
    atmosphere_unit.append( unit_tmp )
    index_last = i

  if len(atmosphere_name) == 0 :
    print('There is no parameter description in the atmosphere model file.')
    print('Program stopped.')
    sys.exit(1)

  # 見出しの次は列番号の行。その次からデータ
  index_data = index_last + 2

  return atmosphere_name, atmosphere_unit, index_data


def set_knudsen_number(config, atmosphere_dict):

  print('Setting Knudsen number...')

  length   = config['satellite']['characteristic_length'] 
  num_data = atmosphere_dict[KEY_DATA]
  num_atm  = atmosphere_dict[KEY_ATM]

  d2_nd_total = np.zeros(num_data).reshape(num_data)
  for m in range(0, len(LIST_MOLECULAR_KIND) ):
    key_tmp      = LIST_MOLECULAR_KIND[m]
    try:
      numb_density = atmosphere_dict[key_tmp]
      diamter      = DICT_DIAMETER_MOLECULAR[key_tmp]
    except KeyError as instance:
      continue
    # Calculate sum( nd*diamter^2 )
    for n in range(0,num_data):
      d2_nd_total[n] = d2_nd_total[n] + numb_density[n]*diamter**2
  
  knudsen_number = 1.0/( np.sqrt(2.0)*np.pi*d2_nd_total*length)
  atmosphere_dict[KEY_KN] = knudsen_number

  return atmosphere_dict


def get_atmosphere_property(altitude, atmosphere_dict):

  # Interpolate atmosphere data from altitude data of satellite
  # 補間器は initial_settings_atmosphere で構築済みのものを評価するだけ

  altitude_atm    = atmosphere_dict[KEY_Height]
  density_atm     = atmosphere_dict[KEY_Mass_density]
  temperature_atm = atmosphere_dict[KEY_Temperature_neutral]
  knudsen_atm     = atmosphere_dict[KEY_KN]

  if altitude < altitude_atm[0] :
    # 下端より下は端の値で止める（高度 0 以下ではソルバーが計算を打ち切る）
    density     = density_atm[0]
    temperature = temperature_atm[0]
    knudsen     = knudsen_atm[0]

  elif altitude > altitude_atm[-1] :
    scale_height = atmosphere_dict[KEY_SCALE_HEIGHT]
    if scale_height is None :
      density     = density_atm[-1]
      temperature = temperature_atm[-1]
      knudsen     = knudsen_atm[-1]
    else :
      # 等温・拡散平衡を仮定した気圧式で外挿する。
      # 温度は等温領域なので端の値のまま。Kn は数密度に反比例するので逆に増える。
      #
      # 指数には上限を置く。上端から離れすぎると exp がアンダーフローして密度 0 に
      # なり、Kn = 端の値/0 が Inf になって出力ファイルに Inf が並ぶ
      # （0 除算の RuntimeWarning つき）。物理としては真空・完全な自由分子流の極限で、
      # 空力係数もそこで飽和しているので、EXPONENT_MAXIMUM より遠くは同じ値で扱う。
      exponent    = min( (altitude - altitude_atm[-1])/scale_height, EXPONENT_MAXIMUM )
      factor      = np.exp( -exponent )
      density     = density_atm[-1]*factor
      temperature = temperature_atm[-1]
      knudsen     = knudsen_atm[-1]/factor

  else :
    interpolator_atm = atmosphere_dict[KEY_INTERP]
    density       = interpolator_atm[KEY_Mass_density](altitude)
    temperature   = interpolator_atm[KEY_Temperature_neutral](altitude)
    knudsen       = interpolator_atm[KEY_KN](altitude)

  return density, temperature, knudsen


