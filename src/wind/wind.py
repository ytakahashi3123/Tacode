#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/10
#
# 風。
#
# 従来このコードは「大気は地球と剛体的に共回転する」と仮定し、ECEF 速度を
# そのまま対気速度として使っていた（force_term.py の抗力、solver.get_aerodynamic_state）。
# 風はその仮定を緩めるだけで、定式化上の矛盾は生じない:
#
#   v_air = v_ecef - v_wind(ECEF)
#
# **対気速度を使うのは空力だけ**（3 自由度の抗力、6 自由度の力とモーメント、動圧、
# 迎角・横滑り角、姿勢振動の周期見積もり）。重力・コリオリ力・遠心力は ECEF 速度の
# ままで、**角速度にも触らない**（風は並進の量であって回転ではない。空力減衰の
# omega_relative は ECEF 基準のまま）。
#
# 既定は無効。config に wind セクションが無い（あるいは flag_wind: False）なら
# initial_settings_wind は None を返し、出力を含めて従来と完全に同じ挙動になる。
# 姿勢（attitude_dict）・エポック（epoch_dict）と同じ設計で、None かどうかだけで分岐する。
#
# 風の向きの規約:
#
#   [東, 北, 上] の**地心**ローカル系、m/s。
#
# これは初期速度（initial_settings.velocity）とまったく同じ規約で、同じ変換
# （coordinate_system.convert_polar_carteasian）を通す。気象データの u/v は測地の
# 水平面で定義されているので最大 0.19 度の傾きがあるが、既存の初期速度・出力と
# 規約を揃える方を採った（README の "Reference frame of the initial velocity" と同じ話）。
# 水平風 60 m/s に対して混入する鉛直成分は 0.2 m/s ほど。
#
# テーブルは最大 4 次元（時刻 x 経度 x 緯度 x 高度）。**節点が 1 つの軸は落とす**ので、
# 単一スナップショットの 3 次元場も、鉛直プロファイル 1 本も、同じ形式・同じ経路で扱える。
#
# wind.velocity_factor は風の場に一様に掛かる係数（既定 1.0）。モンテカルロで
# テーブルの風を振るための入口で、initial_settings.density_factor と同じ流儀。

import sys as sys
import numpy as np
import os as os
import scipy.interpolate   # set_interpolator のコメントに残した
from bisect import bisect_left
import coordinate_system.coordinate_system as coordinate_system
import epoch.epoch as epoch_module
from general.general import get_setting, get_database_directory

# Dict key
KEY_MODEL     = 'kind_wind_model'
KEY_VELOCITY  = 'velocity'
KEY_FACTOR    = 'velocity_factor'
KEY_TIME      = 'Time'
KEY_LONGITUDE = 'Longitude'
KEY_LATITUDE  = 'Latitude'
KEY_ALTITUDE  = 'Altitude'
KEY_WIND      = 'Wind'
KEY_INTERP    = 'Interpolator'
KEY_NODE      = 'Node'          # 補間器が持つ軸の節点（Python のリスト）
KEY_VALUE     = 'Value'         # 補間器が持つ風の値（縮退した軸を落としたもの）
KEY_AXIS      = 'Axis_active'
KEY_EPOCH     = 'Epoch'
KEY_OFFSET    = 'Time_offset'
KEY_FILE      = 'Filename'
KEY_EXTRAP    = 'kind_extrapolation'
KEY_WARNED    = 'Warned_out_of_range'

MODEL_CONSTANT = 'constant'
MODEL_FILEREAD = 'fileread'

EXTRAPOLATION_ZERO  = 'zero'
EXTRAPOLATION_CLAMP = 'clamp'

# テーブルのデータ行の列数
# --6 列: 経度, 緯度, 高度, 東, 北, 上（時刻軸なし）
# --7 列: 時刻, 経度, 緯度, 高度, 東, 北, 上（時刻はテーブルのエポックからの秒）
NUM_COLUMN_SNAPSHOT = 6
NUM_COLUMN_SERIES   = 7

# 軸の並び（時刻, 経度, 緯度, 高度）
INDEX_TIME      = 0
INDEX_LONGITUDE = 1
INDEX_LATITUDE  = 2
INDEX_ALTITUDE  = 3

# 経度軸が全球かどうかの判定（1 格子分を足して 360 度に届くか）
TOLERANCE_GLOBAL_LONGITUDE = 1.e-3

# エポックの食い違いを警告する閾値（秒）。6 時間格子に対して 1 時間
TOLERANCE_EPOCH = 3600.0

# 範囲内外の判定に使う相対の余裕。節点そのものを引くとき、deg -> rad -> deg の
# 往復で最下位ビットが動いて「範囲外」と判定されるのを防ぐ
# （境界は測度 0 の集合なので、素の比較では必ずこの事故が起きる）
TOLERANCE_RANGE = 1.e-9

# 単位換算（orbital をここから import すると循環参照になるので自前で持つ）
RAD2DEG = 180.0/np.pi
M2KM    = 1.e-3


def initial_settings_wind(config, epoch_dict=None):
  #
  # config の wind セクションを読む。
  # 戻り値は wind_dict、既定（セクション無し / flag_wind: False）では None。
  # epoch_dict はテーブルの時刻軸を絶対時刻に結びつけるのに使う。
  #
  section = get_setting(config, 'wind', None)
  if not bool( get_setting(section, 'flag_wind', False) ) :
    return None

  kind_wind_model = get_setting(section, KEY_MODEL, MODEL_CONSTANT)

  print('Setting the wind...')

  if kind_wind_model == MODEL_CONSTANT :
    velocity_wind = np.array( get_setting(section, KEY_VELOCITY, [0.0, 0.0, 0.0]), dtype=float )
    if velocity_wind.shape != (3,) :
      print('wind.velocity must have three components [East, North, Up] in m/s.')
      print('Program stopped.')
      sys.exit(1)
    print('--Model: constant')
    print('--Wind velocity [East, North, Up] (m/s):', ', '.join(['{:g}'.format(value) for value in velocity_wind]))
    wind_dict = {KEY_MODEL: kind_wind_model, KEY_VELOCITY: velocity_wind}

  elif kind_wind_model == MODEL_FILEREAD :
    print('--Model: fileread')
    wind_dict = read_wind_file(config)
    wind_dict = set_interpolator(wind_dict)
    wind_dict = set_time_offset(wind_dict, epoch_dict)

  else :
    print('kind_wind_model in config is incorrect:', kind_wind_model)
    print('--Available: constant, fileread')
    print('Program stopped.')
    sys.exit(1)

  wind_dict[KEY_FACTOR] = get_velocity_factor(section)
  if wind_dict[KEY_FACTOR] != 1.0 :
    print('--Wind velocity factor: {:g}'.format(wind_dict[KEY_FACTOR]))

  return wind_dict


def get_velocity_factor(section):
  #
  # 風のスケール係数（既定 1.0）。風の場そのものに一様に掛かる。
  #
  # モンテカルロがテーブルの風を振るための入口である: driver は制御ファイルの
  # **数値**を書き換える仕組みで、テーブルはファイル名で選ぶので、テーブルを
  # 振るにはこういう数値が要る。initial_settings.density_factor とまったく
  # 同じ流儀（1 要素のリスト）で書くのはそのため（driver がリストの要素ごとに
  # 置換するので、スカラーのままでは振れない）。
  #
  # 既定の 1.0 では掛け算が恒等（IEEE 754 で 1.0*x は厳密に x）なので、
  # 出力は 1 バイトも変わらない。
  #
  factor = get_setting(section, KEY_FACTOR, 1.0)

  if isinstance(factor, (list, tuple, np.ndarray)) :
    if len(factor) != 1 :
      print('wind.velocity_factor must be a single value (a one-element list).')
      print('Program stopped.')
      sys.exit(1)
    factor = factor[0]

  try:
    factor = float(factor)
  except (TypeError, ValueError):
    print('wind.velocity_factor must be a number:', factor)
    print('Program stopped.')
    sys.exit(1)

  return factor


def read_wind_file(config):
  #
  # 風のテーブルを読む。パスの解決は大気・空力と同じ規約
  # （default/auto ならスクリプト位置基準、manual なら config の相対パス）。
  #
  # 形式は 1 種類だけで、1 行 1 格子点。列数で時刻軸の有無が決まる:
  #
  #   6 列: 経度[deg.] 緯度[deg.] 高度[km] 東[m/s] 北[m/s] 上[m/s]
  #   7 列: 時刻[s] 経度[deg.] 緯度[deg.] 高度[km] 東[m/s] 北[m/s] 上[m/s]
  #
  # 時刻はテーブルのエポック（"# Epoch (UTC):" の行）からの秒。
  # **1 次元の鉛直プロファイルは、経度・緯度の節点が 1 つだけの場として書けばよい**
  # （多次元の縮退として同じ形式で扱う）。
  #
  section = config['wind']

  directory_path = get_database_directory(section, 'wind', 'wind', 'directory_wind', 'database/wind')

  filename_tmp = directory_path + '/' + get_setting(section, 'filename_wind', 'windmodel.txt')
  print('Reading wind model...:', filename_tmp)

  if not os.path.exists(filename_tmp) :
    print('Wind table not found:', filename_tmp)
    print('--Generate one with database/wind/generate_wind_table.py, or point')
    print('--wind.directory_wind / wind.filename_wind at an existing table.')
    print('Program stopped.')
    sys.exit(1)

  string_epoch = None
  row_list     = []
  num_column   = None
  with open(filename_tmp) as f:
    for line in f:
      stripped = line.strip()
      if len(stripped) == 0 :
        continue
      if stripped.startswith('#') :
        # 生成時刻だけは機械的に読む（時刻軸の基準と、config のエポックとの照合に使う）
        comment = stripped.lstrip('#').strip()
        if comment.lower().startswith('epoch') and ':' in comment :
          string_epoch = comment.split(':', 1)[1].strip()
        continue
      word_list = stripped.split()
      if len(word_list) not in (NUM_COLUMN_SNAPSHOT, NUM_COLUMN_SERIES) :
        continue
      try:
        value_list = [float(word) for word in word_list]
      except ValueError:
        continue
      if num_column is None :
        num_column = len(value_list)
      elif len(value_list) != num_column :
        print('The wind table mixes rows of {:d} and {:d} numbers:'.format(num_column, len(value_list)),
              filename_tmp)
        print('--Every data row must have the same number of columns (6 without a time axis, 7 with one).')
        print('Program stopped.')
        sys.exit(1)
      row_list.append(value_list)

  if len(row_list) == 0 :
    print('No data row was found in the wind table:', filename_tmp)
    print('--Each row must hold six numbers (longitude, latitude, altitude, east, north, up)')
    print('--or seven with the time from the epoch of the table in front.')
    print('Program stopped.')
    sys.exit(1)

  wind_dict = build_grid(np.array(row_list), filename_tmp, num_column == NUM_COLUMN_SERIES)

  wind_dict[KEY_MODEL]  = MODEL_FILEREAD
  wind_dict[KEY_FILE]   = filename_tmp
  wind_dict[KEY_EPOCH]  = string_epoch
  wind_dict[KEY_WARNED] = False

  kind_extrapolation = get_setting(section, KEY_EXTRAP, EXTRAPOLATION_ZERO)
  if kind_extrapolation not in (EXTRAPOLATION_ZERO, EXTRAPOLATION_CLAMP) :
    print('wind.kind_extrapolation in config is incorrect:', kind_extrapolation)
    print('--Available: zero, clamp')
    print('Program stopped.')
    sys.exit(1)
  wind_dict[KEY_EXTRAP] = kind_extrapolation

  print('--Grid: {:d} times x {:d} longitudes x {:d} latitudes x {:d} altitudes'.format(
        len(wind_dict[KEY_TIME]), len(wind_dict[KEY_LONGITUDE]),
        len(wind_dict[KEY_LATITUDE]), len(wind_dict[KEY_ALTITUDE])))
  print('--Longitude (deg.): {:g} to {:g}'.format(wind_dict[KEY_LONGITUDE][0], wind_dict[KEY_LONGITUDE][-1]))
  print('--Latitude  (deg.): {:g} to {:g}'.format(wind_dict[KEY_LATITUDE][0], wind_dict[KEY_LATITUDE][-1]))
  print('--Altitude  (km)  : {:g} to {:g}'.format(wind_dict[KEY_ALTITUDE][0], wind_dict[KEY_ALTITUDE][-1]))
  if len(wind_dict[KEY_TIME]) > 1 :
    print('--Time      (s)   : {:g} to {:g} from the epoch of the table'.format(
          wind_dict[KEY_TIME][0], wind_dict[KEY_TIME][-1]))
  print('--Above and below the altitude range:', kind_extrapolation)
  if string_epoch is not None :
    print('--Epoch of the table (UTC):', string_epoch)

  return wind_dict


def build_grid(rows, filename_tmp, flag_time):
  #
  # 「1 行 1 格子点」の行から (時刻 x 経度 x 緯度 x 高度 x 3) の配列を組み立てる。
  # 行の順序には依存しない。格子が埋まっていなければ停止する
  # （欠けた点を黙って 0 にすると、風が弱いのかデータが無いのか区別できなくなる）。
  #
  if flag_time :
    time      = rows[:,0]
    longitude = normalize_longitude(rows[:,1])
    latitude  = rows[:,2]
    altitude  = rows[:,3]
    velocity  = rows[:,4:7]
  else :
    time      = np.zeros(len(rows))
    longitude = normalize_longitude(rows[:,0])
    latitude  = rows[:,1]
    altitude  = rows[:,2]
    velocity  = rows[:,3:6]

  axis_time      = np.unique(time)
  axis_longitude = np.unique(longitude)
  axis_latitude  = np.unique(latitude)
  axis_altitude  = np.unique(altitude)

  num_expected = len(axis_time)*len(axis_longitude)*len(axis_latitude)*len(axis_altitude)
  if len(rows) != num_expected :
    print('The wind table does not fill its grid:', filename_tmp)
    print('--{:d} rows for {:d} x {:d} x {:d} x {:d} = {:d} grid points'.format(
          len(rows), len(axis_time), len(axis_longitude), len(axis_latitude),
          len(axis_altitude), num_expected))
    print('--(times x longitudes x latitudes x altitudes)')
    print('Program stopped.')
    sys.exit(1)

  wind = np.full((len(axis_time), len(axis_longitude), len(axis_latitude), len(axis_altitude), 3), np.nan)
  wind[np.searchsorted(axis_time,      time),
       np.searchsorted(axis_longitude, longitude),
       np.searchsorted(axis_latitude,  latitude),
       np.searchsorted(axis_altitude,  altitude), :] = velocity

  if np.isnan(wind).any() :
    print('The wind table has duplicated or missing grid points:', filename_tmp)
    print('--Longitudes are folded into [-180, 180), so -180 and 180 deg. are the same')
    print('--meridian and must not both appear.')
    print('Program stopped.')
    sys.exit(1)

  return {KEY_TIME: axis_time, KEY_LONGITUDE: axis_longitude,
          KEY_LATITUDE: axis_latitude, KEY_ALTITUDE: axis_altitude, KEY_WIND: wind}


def normalize_longitude(longitude):
  # 経度を [-180, 180) に畳む。気象データは 0-360 deg. で配られるが、
  # このコードの測地経度は -180 から 180 なので合わせる
  return ( np.asarray(longitude, dtype=float) + 180.0 ) % 360.0 - 180.0


def is_global_longitude(axis_longitude):
  # 経度軸が全球を覆っているか（覆っているなら継ぎ目を跨いだ内挿ができるようにする）
  if len(axis_longitude) < 2 :
    return False
  spacing = np.median( np.diff(axis_longitude) )
  return ( axis_longitude[-1] - axis_longitude[0] + spacing ) >= 360.0 - TOLERANCE_GLOBAL_LONGITUDE


def get_axis_all(wind_dict):
  # 軸を (時刻, 経度, 緯度, 高度) の並びで返す
  return (wind_dict[KEY_TIME], wind_dict[KEY_LONGITUDE], wind_dict[KEY_LATITUDE], wind_dict[KEY_ALTITUDE])


def set_interpolator(wind_dict):
  #
  # 補間器を初期化時に 1 度だけ構築する（大気・空力と同じ方針。
  # 毎ステップ再構築すると計算時間を支配する）。
  #
  # **節点が 1 つしかない軸は落とす。** 単一スナップショット（時刻軸が 1 点）も、
  # 1 次元の鉛直プロファイル（経度・緯度が 1 点）も、そうやって多次元の縮退として扱う。
  # RegularGridInterpolator が 1 点の軸を受けないという実装上の理由もある。
  #
  axis_longitude = wind_dict[KEY_LONGITUDE]
  wind           = wind_dict[KEY_WIND]

  # 全球なら経度の継ぎ目に折り返しの節点を足す（-180 の値を +180 にも置く）
  if is_global_longitude(axis_longitude) :
    axis_longitude = np.append(axis_longitude, axis_longitude[0] + 360.0)
    wind           = np.concatenate((wind, wind[:,0:1,:,:,:]), axis=INDEX_LONGITUDE)
    wind_dict[KEY_LONGITUDE] = axis_longitude
    wind_dict[KEY_WIND]      = wind

  axis_all    = get_axis_all(wind_dict)
  axis_active = [index for index, axis in enumerate(axis_all) if len(axis) > 1]

  if len(axis_active) == 0 :
    # 1 点だけの表。定数風と同じ
    wind_dict[KEY_INTERP] = None
    wind_dict[KEY_AXIS]   = axis_active
    return wind_dict

  values = wind
  for index in reversed(range(0,4)):
    if index not in axis_active :
      values = values.take(0, axis=index)

  # 補間器が持つのは軸の節点と値だけで、引くのは evaluate_multilinear。
  # scipy の RegularGridInterpolator は使わない: 1 点を引くのに 3 次元で 39 us かかり、
  # 風のテーブルを引くケースでは計算時間の半分を占めていた（手書きは 4 us）。
  # 空力表（satellite.set_interpolator）と同じ理由で、同じ手を当ててある。
  #
  # **値は RGI とビット単位で一致する**（evaluate_multilinear のコメントと
  # test_wind.TestMultilinearMatchesScipy）。
  #
  # 節点を ndarray ではなく Python のリストで持つのは速度のため。理由は同じ関数に書いた。
  #
  # 戻すときはここを次の 3 行に替え、get_wind_table の evaluate_multilinear の
  # 呼び出しを interpolator(np.array([query]))[0] に戻せばよい:
  #
  #   wind_dict[KEY_INTERP] = scipy.interpolate.RegularGridInterpolator(
  #                             tuple([axis_all[index] for index in axis_active]), values,
  #                             method='linear', bounds_error=False, fill_value=None)
  #
  wind_dict[KEY_INTERP] = { KEY_NODE : [axis_all[index].tolist() for index in axis_active],
                            KEY_VALUE: values }
  wind_dict[KEY_AXIS]   = axis_active

  return wind_dict


def set_time_offset(wind_dict, epoch_dict):
  #
  # テーブルの時刻軸を、計算の経過時間に結びつける。
  #
  #   テーブルの時刻 = (config のエポック - テーブルのエポック) + 経過時間
  #
  # 時刻軸を持つテーブルには config の epoch が必要（絶対時刻が無いと置けない）。
  # 単一スナップショットなら、エポックの食い違いを警告するだけで進む。
  #
  flag_series  = len(wind_dict[KEY_TIME]) > 1
  string_epoch = wind_dict.get(KEY_EPOCH)

  wind_dict[KEY_OFFSET] = 0.0

  if flag_series and string_epoch is None :
    print('The wind table has a time axis but no epoch.')
    print('--Add a line "# Epoch (UTC): <ISO 8601>" to the table, so that its time')
    print('--axis can be placed in absolute time.')
    print('Program stopped.')
    sys.exit(1)

  if flag_series and epoch_dict is None :
    print('The wind table has a time axis but the epoch section of config is off.')
    print('--Set epoch.flag_epoch: True and give epoch.datetime, so that the run can')
    print('--be placed on the time axis of the table.')
    print('Program stopped.')
    sys.exit(1)

  if string_epoch is None or epoch_dict is None :
    return wind_dict

  time_table = epoch_module.parse_datetime(string_epoch)
  difference = ( epoch_module.get_datetime(epoch_dict, 0.0) - time_table ).total_seconds()
  wind_dict[KEY_OFFSET] = difference

  if flag_series :
    axis_time = wind_dict[KEY_TIME]
    print('--The run starts at {:g} s on the time axis of the table'.format(difference))
    if difference < axis_time[0] or difference > axis_time[-1] :
      print('Warning: the epoch of the run is outside the time range of the wind table.')
      print('--Table: {:g} to {:g} s from {}'.format(axis_time[0], axis_time[-1], string_epoch))
      print('--The nearest snapshot is used (the time axis is clamped).')
    return wind_dict

  if abs(difference) > TOLERANCE_EPOCH :
    print('Warning: the wind table was generated for a different time than the epoch.')
    print('--Table:', string_epoch)
    print('--Epoch:', epoch_module.get_string(epoch_dict, 0.0))
    print('--Difference: {:.1f} hours. The wind is used as it is.'.format(difference/3600.0))

  return wind_dict


def get_wind_local(coordinate_geodetic, wind_dict, time_elapsed=0.0):
  #
  # 現在位置の風を地心ローカル系 [東, 北, 上] (m/s) で返す。
  # coordinate_geodetic: [経度, 緯度, 高度]（rad, rad, m）。
  # time_elapsed: 計算開始からの経過時間（s）。時刻軸を持つテーブルだけが使う。
  #
  if wind_dict[KEY_MODEL] == MODEL_CONSTANT :
    wind_local = wind_dict[KEY_VELOCITY]
  else :
    wind_local = get_wind_table(coordinate_geodetic, wind_dict, time_elapsed)

  # スケール係数（wind.velocity_factor）。既定の 1.0 では厳密に恒等なので、
  # 掛けても結果は 1 ビットも動かない
  return get_setting(wind_dict, KEY_FACTOR, 1.0)*wind_local


def is_outside(value, axis):
  #
  # 軸の範囲から出ているか。節点が 1 つだけの軸（縮退した方向）は常に範囲内とする。
  # 端の比較には TOLERANCE_RANGE の余裕を持たせる（上記）。
  #
  if len(axis) < 2 :
    return False

  margin = TOLERANCE_RANGE*max( abs(axis[0]), abs(axis[-1]), axis[-1] - axis[0] )

  return value < axis[0] - margin or value > axis[-1] + margin


def get_wind_table(coordinate_geodetic, wind_dict, time_elapsed=0.0):
  #
  # テーブルを内挿して風を返す（評価するだけ。補間器は初期化時に構築済み）。
  #
  # 範囲外の扱い:
  #   水平（経度・緯度）と時刻は**常にクランプ**する。領域や時間窓の端の値を使う
  #     ということで、全球のテーブルなら継ぎ目は折り返し節点で繋いであるので端に当たらない。
  #   高度は wind.kind_extrapolation:
  #     zero  (既定) — 範囲外では風なし。大気が共回転する従来の仮定に戻る。
  #                    NCEP は 31-48 km までしか無いので、その上へクランプすると
  #                    成層圏のジェットを熱圏まで引き延ばすことになる
  #     clamp        — 端の値で止める
  # いずれの場合も、最初に範囲を出たところで 1 度だけ警告する。
  #
  point = [ wind_dict[KEY_OFFSET] + time_elapsed,
            normalize_longitude( coordinate_geodetic[0]*RAD2DEG ),
            coordinate_geodetic[1]*RAD2DEG,
            coordinate_geodetic[2]*M2KM ]

  axis_all = get_axis_all(wind_dict)

  altitude_outside = is_outside(point[INDEX_ALTITUDE], axis_all[INDEX_ALTITUDE])
  other_outside    = is_outside(point[INDEX_TIME],      axis_all[INDEX_TIME])      \
                  or is_outside(point[INDEX_LONGITUDE], axis_all[INDEX_LONGITUDE]) \
                  or is_outside(point[INDEX_LATITUDE],  axis_all[INDEX_LATITUDE])

  if ( altitude_outside or other_outside ) and not wind_dict[KEY_WARNED] :
    wind_dict[KEY_WARNED] = True
    print('Caution: the trajectory leaves the range of the wind table.')
    print('--Table: longitude {:g} to {:g} deg., latitude {:g} to {:g} deg., altitude {:g} to {:g} km'.format(
          axis_all[INDEX_LONGITUDE][0], axis_all[INDEX_LONGITUDE][-1],
          axis_all[INDEX_LATITUDE][0], axis_all[INDEX_LATITUDE][-1],
          axis_all[INDEX_ALTITUDE][0], axis_all[INDEX_ALTITUDE][-1]))
    print('--Now  : longitude {:g} deg., latitude {:g} deg., altitude {:g} km'.format(
          point[INDEX_LONGITUDE], point[INDEX_LATITUDE], point[INDEX_ALTITUDE]))
    if len(axis_all[INDEX_TIME]) > 1 :
      print('--Time : {:g} s on a table axis of {:g} to {:g} s'.format(
            point[INDEX_TIME], axis_all[INDEX_TIME][0], axis_all[INDEX_TIME][-1]))
    print('--The horizontal directions and the time are clamped to the edge; outside the')
    print('--altitude range the wind is "{}" (wind.kind_extrapolation).'.format(wind_dict[KEY_EXTRAP]))

  if altitude_outside and wind_dict[KEY_EXTRAP] == EXTRAPOLATION_ZERO :
    return np.zeros(3)

  interpolator = wind_dict[KEY_INTERP]
  if interpolator is None :
    # 節点が 1 つだけの表
    return wind_dict[KEY_WIND][0,0,0,0,:]

  # 軸の端でクランプしてから引く（範囲外の扱いは上記）。素の float に落とすのは
  # 速度のため。np.float64 のままだと算術 1 つごとに numpy のスカラー経路を通る
  # （値は変わらない。float() は倍精度の値をそのまま取り出すだけ）
  query = []
  for index in wind_dict[KEY_AXIS]:
    axis = axis_all[index]
    query.append( float( min( max( point[index], axis[0] ), axis[-1] ) ) )

  return evaluate_multilinear(interpolator[KEY_NODE], interpolator[KEY_VALUE], query)


def evaluate_multilinear(axis_node, values, point):
  #
  # 最大 4 次元（時刻 x 経度 x 緯度 x 高度）の線形補間。軸は昇順、点は軸の内側にある
  # こと（呼び出し側でクランプ済み）。values の末尾の軸は風の 3 成分。
  #
  # scipy.interpolate.RegularGridInterpolator の代わりに手で書いてある。理由と
  # 戻し方は set_interpolator のコメントを見ること。空力表の
  # satellite.evaluate_bilinear の N 次元版だが、次元数も末尾の成分数も違うので
  # 共通化していない（共通化すると、この関数の目的である速さが出ない）。
  #
  # **値は RGI とビット単位で一致する。** RGI の _evaluate_linear は超立方体の頂点を
  # itertools.product の順（**最初の軸が最も外側、各軸は下側の節点が先**）に回り、
  # **重みを軸の順に左から掛け**、0 から順に足す。ここもそのとおりに書いてある。
  # 3 次元なら:
  #
  #   value = 0 + v[i,j,k]*((1-y0)*(1-y1)*(1-y2)) + v[i,j,k+1]*((1-y0)*(1-y1)*y2)
  #             + v[i,j+1,k]*((1-y0)*y1*(1-y2))   + ... + v[i+1,j+1,k+1]*(y0*y1*y2)
  #
  # 区間の決め方（searchsorted の既定 side='left' から 1 を引き、両端で丸める）も
  # RGI に合わせてある。まとめ方を変えると最下位桁が動き、バイト一致が崩れる。
  #
  # 節点を ndarray ではなく Python のリストで持ち、bisect で引くのは速度のため。
  # np.searchsorted は 1 点を引くだけでも ndarray を作って返す
  # （3 次元で 6.7 us -> 4.2 us）。区間の決め方は bisect_left でも同じ。
  #
  vertex = [((), 1.0)]
  for node, coordinate in zip(axis_node, point):
    index    = min( max(bisect_left(node, coordinate) - 1, 0), len(node) - 2 )
    distance = (coordinate - node[index])/(node[index+1] - node[index])
    vertex   = [ (index_vertex + (index_node,), weight*weight_node)
                 for index_vertex, weight in vertex
                 for index_node, weight_node in ((index, 1.0 - distance), (index + 1, distance)) ]

  wind_east = wind_north = wind_up = 0.0
  for index_vertex, weight in vertex:
    value      = values[index_vertex]
    wind_east  = wind_east  + value[0]*weight
    wind_north = wind_north + value[1]*weight
    wind_up    = wind_up    + value[2]*weight

  return np.array([wind_east, wind_north, wind_up])


def get_wind_velocity(config, coordinate, coordinate_geodetic, wind_dict, time_elapsed=0.0):
  #
  # 現在位置の風を ECEF 直交成分 (m/s) で返す。
  # 初期速度とまったく同じ経路（地心の極角 -> convert_polar_carteasian）を通す。
  #
  wind_local = get_wind_local(coordinate_geodetic, wind_dict, time_elapsed)

  coordinate_polar = coordinate_system.set_angle_polar(config, coordinate)
  latitude  = coordinate_polar[1]
  longitude = coordinate_polar[2]

  return np.array( coordinate_system.convert_polar_carteasian(config, wind_local, longitude, latitude) )


def get_relative_velocity(config, coordinate, coordinate_geodetic, velocity, wind_dict, time_elapsed=0.0):
  #
  # 対気速度（ECEF 成分, m/s）。wind_dict が None なら ECEF 速度をそのまま返す
  # （引き算を通さないので、風を切ったときの結果はビット単位で従来と同じになる）。
  #
  if wind_dict is None :
    return velocity

  return np.array(velocity) - get_wind_velocity(config, coordinate, coordinate_geodetic, wind_dict, time_elapsed)
