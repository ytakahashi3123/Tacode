#!/usr/bin/env python3
#
# 風のテーブルを作る。
#
# Tacode のソルバー（src/）はテキストのテーブルしか読まない。NetCDF や GRIB を
# 読むライブラリはここに閉じ込めてある（計算環境に重い依存を持ち込まないため。
# matplotlib を src/ で import しないのと同じ理由）。
#
# 出力形式は 1 種類だけで、1 行 1 格子点:
#
#   経度[deg.] 緯度[deg.] 高度[km] 東[m/s] 北[m/s] 上[m/s]
#
# **鉛直プロファイル（1 次元）は、経度・緯度の節点が 1 つだけの場として書く。**
# 多次元の場と同じ形式・同じ読み取り経路になる。
#
# 出力は最大 4 次元（時刻 x 経度 x 緯度 x 高度）。--duration を与えると時刻軸が付き、
# 行の先頭にテーブルのエポックからの秒が入る（7 列）。
#
# データ源:
#
#   ncep    NCEP/NCAR Reanalysis 1（過去の実データ）を NOAA PSL の OPeNDAP から
#           必要な範囲だけ取ってくる。**標準ライブラリと numpy だけで済む**
#           （.ascii の応答を読むので netCDF4 も cfgrib も要らない）。
#           気圧面 1000-10 hPa = 地表から約 31 km。それより上は持っていない
#   hwm14   HWM14（NRL の経験モデル、0-500 km）を呼ぶ。**HWM14 は Fortran と NRL の
#           係数バイナリなので、このスクリプトの依存にはしない**。使うときだけ
#           遅延 import し、無ければ導入方法を示して止まる
#   profile 高度と風速の対を並べたテキストから 1 次元プロファイルを作る。
#           オフラインで使えるので、形式の確認や感度解析に向く
#   merge   2 つのテーブルを高度方向に合成する。**NCEP（下層）と HWM14（上層）を
#           遷移層で繋ぐ**のがこれ。境目で風速が飛ばないようにする
#   calm    どこでも風速 0（テーブルを読む経路そのものの確認用）
#
# 使い方:
#   # 2024-01-01 00 UTC、南米上空の領域
#   python3 generate_wind_table.py --source ncep --datetime 2024-01-01T00:00:00Z \
#           --longitude -125 -115 --latitude -20 -10 -o wind_ncep_20240101.txt
#
#   # 同じ領域・同じ時刻を HWM14 で 20-500 km。--like で水平・時刻の軸を揃える
#   python3 generate_wind_table.py --source hwm14 --datetime 2024-01-01T00:00:00Z \
#           --like wind_ncep_20240101.txt --altitude 20 500 --altitude-step 20 \
#           --hwm14-data ~/hwm14 -o wind_hwm14_20240101.txt
#
#   # 20-30 km を遷移層にして合成する
#   python3 generate_wind_table.py --source merge \
#           -i wind_ncep_20240101.txt wind_hwm14_20240101.txt --transition 20 30 \
#           -o wind_merged_20240101.txt
#
#   # 6 時間おきの時刻軸を持つテーブル（時間内挿つき）
#   python3 generate_wind_table.py --source ncep --datetime 2024-01-01T00:00:00Z \
#           --duration 43200 --longitude -125 -115 --latitude -20 -10 -o wind_series.txt
#
#   # 高度 [km] と 東・北 [m/s] を並べたファイルから鉛直プロファイル
#   python3 generate_wind_table.py --source profile -i profile.txt -o wind_profile.txt
#
# 注意:
#   - 鉛直風は書かない（0 を入れる）。総観規模の鉛直流は cm/s の桁で、
#     軌道計算には効かない。omega [Pa/s] からの換算には密度が要る
#   - NCEP の風は気圧面上にある。高度は同じ時刻・同じ格子の**ジオポテンシャル高度
#     （hgt）から作る**。hgt は場所によって違うので、領域平均のジオポテンシャル高度を
#     共通の高度節点にとり、各鉛直列をその節点へ内挿してから書き出す
#   - ジオポテンシャル高度 -> 幾何高度の換算 z = R H/(R - H) を掛けてある
#     （30 km で 0.14 km の差）

import argparse
import datetime
import numpy as np
import os
import sys
import urllib.request

# テーブルの形式を 2 か所に持たないため、読み取りはソルバーの実装を使う
# （src_helper が src/attitude/attitude.py を import しているのと同じ理由）。
DIRECTORY_SRC = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'src')
sys.path.insert(0, os.path.normpath(DIRECTORY_SRC))
import wind.wind as wind_module

# NOAA PSL の OPeNDAP。年ごとに 1 ファイル、6 時間毎・2.5 度格子・17 気圧面
URL_NCEP_BASE = 'https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/pressure'

# NCEP の時間軸の起点（"hours since 1800-01-01 00:00:0.0"）
DATETIME_NCEP_ZERO = datetime.datetime(1800, 1, 1, tzinfo=datetime.timezone.utc)
SECOND_PER_HOUR    = 3600.0

# ジオポテンシャル高度から幾何高度への換算に使う地球半径, m
RADIUS_GEOPOTENTIAL = 6356766.0

TIMEOUT_REQUEST = 60


def parse_datetime(string_datetime):
  string_tmp = string_datetime.strip()
  if string_tmp.endswith('Z') or string_tmp.endswith('z') :
    string_tmp = string_tmp[:-1] + '+00:00'
  time_tmp = datetime.datetime.fromisoformat(string_tmp)
  if time_tmp.tzinfo is None :
    time_tmp = time_tmp.replace(tzinfo=datetime.timezone.utc)
  return time_tmp.astimezone(datetime.timezone.utc)


def fetch_opendap(url):
  # OPeNDAP の .ascii 応答を取ってくる
  print('--Fetching:', url)
  try:
    with urllib.request.urlopen(url, timeout=TIMEOUT_REQUEST) as response:
      return response.read().decode('utf-8', 'replace')
  except Exception as error:
    print('Failed to fetch the data:', type(error).__name__, error)
    print('--Check the network, or download the NetCDF file and convert it yourself.')
    sys.exit(1)


def parse_opendap_variable(text, name):
  #
  # .ascii 応答から 1 変数を取り出す。応答は
  #
  #   <name>[n]
  #   v1, v2, ...
  #
  # という節（座標変数）と
  #
  #   <name>.<name>[t][z][y][x]
  #   [0][0][0], v, v, v
  #
  # という節（グリッド変数、最後の次元が行内に並ぶ）でできている。
  #
  line_list = text.splitlines()
  value_list = []
  index = 0
  found = False
  while index < len(line_list):
    stripped = line_list[index].strip()
    if stripped.startswith(name + '[') or stripped.startswith(name + '.' + name + '[') :
      found = True
      index = index + 1
      while index < len(line_list):
        row = line_list[index].strip()
        if len(row) == 0 :
          break
        if '[' in row and ']' in row and ',' in row :
          row = row.split(',', 1)[1]
        elif '[' in row :
          break
        value_list.extend([float(word) for word in row.split(',') if len(word.strip()) > 0])
        index = index + 1
      break
    index = index + 1

  if not found or len(value_list) == 0 :
    print('The variable "{}" was not found in the OPeNDAP response.'.format(name))
    sys.exit(1)

  return np.array(value_list)


def index_range(axis, value_min, value_max):
  # 指定した範囲を含む最小の添字範囲（端が範囲外なら内挿できないので 1 格子広げる）
  order = 1 if axis[-1] >= axis[0] else -1
  if order > 0 :
    index_low  = max(0, int(np.searchsorted(axis, value_min, side='right')) - 1)
    index_high = min(len(axis) - 1, int(np.searchsorted(axis, value_max, side='left')))
  else :
    reversed_axis = axis[::-1]
    index_low  = max(0, int(np.searchsorted(reversed_axis, value_min, side='right')) - 1)
    index_high = min(len(axis) - 1, int(np.searchsorted(reversed_axis, value_max, side='left')))
    index_low, index_high = len(axis) - 1 - index_high, len(axis) - 1 - index_low
  return index_low, index_high


def geometric_altitude(height_geopotential):
  # ジオポテンシャル高度 [m] -> 幾何高度 [m]
  return RADIUS_GEOPOTENTIAL*height_geopotential/(RADIUS_GEOPOTENTIAL - height_geopotential)


def build_ncep(argument):
  #
  # NCEP/NCAR Reanalysis 1 から領域を切り出してテーブルにする。
  #
  time_target = parse_datetime(argument.datetime)
  year = time_target.year

  url_uwnd = '{}/uwnd.{:d}.nc'.format(URL_NCEP_BASE, year)
  url_vwnd = '{}/vwnd.{:d}.nc'.format(URL_NCEP_BASE, year)
  url_hgt  = '{}/hgt.{:d}.nc'.format(URL_NCEP_BASE, year)

  # 座標軸をまず取る
  text_axis = fetch_opendap(url_uwnd + '.ascii?level,lat,lon,time')
  axis_level = parse_opendap_variable(text_axis, 'level')
  axis_lat   = parse_opendap_variable(text_axis, 'lat')
  axis_lon   = parse_opendap_variable(text_axis, 'lon')
  axis_time  = parse_opendap_variable(text_axis, 'time')

  # 時刻の添字。--duration が 0 なら最近傍の 1 枚、正なら [t, t+duration] を覆う範囲
  hour_target = ( time_target - DATETIME_NCEP_ZERO ).total_seconds()/SECOND_PER_HOUR
  index_first = int(np.argmin(np.abs(axis_time - hour_target)))
  if axis_time[index_first] > hour_target and index_first > 0 :
    index_first = index_first - 1
  if argument.duration > 0.0 :
    hour_end   = hour_target + argument.duration/SECOND_PER_HOUR
    index_last = int(np.searchsorted(axis_time, hour_end, side='left'))
    index_last = min(len(axis_time) - 1, max(index_first, index_last))
  else :
    index_first = int(np.argmin(np.abs(axis_time - hour_target)))
    index_last  = index_first

  time_first = DATETIME_NCEP_ZERO + datetime.timedelta(hours=float(axis_time[index_first]))
  print('--Requested time:', time_target.isoformat())
  print('--Samples used  : {} to {} ({:d} snapshots)'.format(
        time_first.isoformat(),
        (DATETIME_NCEP_ZERO + datetime.timedelta(hours=float(axis_time[index_last]))).isoformat(),
        index_last - index_first + 1))
  if index_first == index_last and abs(float(axis_time[index_first]) - hour_target) > 3.0 :
    print('--Caution: the nearest sample is more than 3 hours away.')

  # 経度は 0-360 deg. で並んでいるので、要求範囲もその向きに直す
  longitude_min = argument.longitude[0] % 360.0
  longitude_max = argument.longitude[1] % 360.0
  if longitude_min > longitude_max :
    print('The longitude range crosses the 0/360 deg. seam of the NCEP grid.')
    print('--Give the range without crossing it, or use the whole globe (-180 180).')
    sys.exit(1)

  index_lon = index_range(axis_lon, longitude_min, longitude_max)
  index_lat = index_range(axis_lat, argument.latitude[0], argument.latitude[1])
  index_lev = (0, len(axis_level) - 1)

  selector = '[{0:d}:1:{1:d}][{2:d}:1:{3:d}][{4:d}:1:{5:d}][{6:d}:1:{7:d}]'.format(
             index_first, index_last, index_lev[0], index_lev[1],
             index_lat[0], index_lat[1], index_lon[0], index_lon[1])

  shape = (index_last - index_first + 1,
           index_lev[1] - index_lev[0] + 1,
           index_lat[1] - index_lat[0] + 1,
           index_lon[1] - index_lon[0] + 1)

  uwnd = parse_opendap_variable(fetch_opendap(url_uwnd + '.ascii?uwnd' + selector), 'uwnd').reshape(shape)
  vwnd = parse_opendap_variable(fetch_opendap(url_vwnd + '.ascii?vwnd' + selector), 'vwnd').reshape(shape)
  hgt  = parse_opendap_variable(fetch_opendap(url_hgt  + '.ascii?hgt'  + selector), 'hgt' ).reshape(shape)

  longitude = axis_lon[index_lon[0]:index_lon[1]+1]
  latitude  = axis_lat[index_lat[0]:index_lat[1]+1]
  time_second = ( axis_time[index_first:index_last+1] - axis_time[index_first] )*SECOND_PER_HOUR
  flag_time = len(time_second) > 1

  # 気圧面 -> 高度。ジオポテンシャル高度は場所（と時刻）ごとに違うので、
  # 領域・時間平均を共通の節点にとり、各鉛直列をそこへ内挿する
  altitude_node = geometric_altitude( hgt.mean(axis=(0,2,3)) )/1000.0
  altitude_node = np.sort(altitude_node)
  if argument.altitude_max is not None :
    altitude_node = altitude_node[altitude_node <= argument.altitude_max]
  print('--Pressure levels: {:g} to {:g} hPa'.format(axis_level[0], axis_level[-1]))
  print('--Altitude nodes (km): {:g} to {:g} ({:d} levels)'.format(
        altitude_node[0], altitude_node[-1], len(altitude_node)))

  row_list = []
  for index_t, value_time in enumerate(time_second):
    for index_i, value_lon in enumerate(longitude):
      for index_j, value_lat in enumerate(latitude):
        altitude_column = geometric_altitude( hgt[index_t,:,index_j,index_i] )/1000.0
        sort_column = np.argsort(altitude_column)
        east  = np.interp(altitude_node, altitude_column[sort_column], uwnd[index_t,sort_column,index_j,index_i])
        north = np.interp(altitude_node, altitude_column[sort_column], vwnd[index_t,sort_column,index_j,index_i])
        for index_k, value_alt in enumerate(altitude_node):
          row = [value_lon, value_lat, value_alt, east[index_k], north[index_k], 0.0]
          row_list.append([value_time] + row if flag_time else row)

  header = ['Source: NCEP/NCAR Reanalysis 1, pressure-level uwnd/vwnd/hgt',
            '--Obtained from ' + URL_NCEP_BASE,
            'Epoch (UTC): ' + time_first.strftime('%Y-%m-%dT%H:%M:%S') + '.000Z',
            'Region: longitude {:g} to {:g} deg. (0-360 convention), latitude {:g} to {:g} deg.'.format(
              longitude[0], longitude[-1], latitude[0], latitude[-1]),
            'Vertical: geometric altitude from the region- and time-mean geopotential height of each pressure level',
            'The vertical wind is not modelled and is written as zero']
  if flag_time :
    header.append('Time axis: {:d} snapshots, 0 to {:g} s from the epoch'.format(
                  len(time_second), time_second[-1]))

  return row_list, header


def import_hwm14():
  #
  # HWM14 を遅延 import する。
  #
  # HWM14 は Fortran のソースと NRL の係数バイナリ（hwm123114.bin / dwm07b104i.dat /
  # gd2qd.dat）でできており、pip で入るものではない。**使うときにだけ確かめ、
  # 無ければ導入方法を示して止める**（ソルバーはもちろんこれに依存しない）。
  #
  module = None
  try:
    import hwm14 as module
  except ImportError:
    try:
      from pyhwm2014 import hwm14 as module
    except ImportError:
      module = None

  if module is None :
    print('HWM14 is not available in this Python environment.')
    print('--HWM14 is Fortran plus the NRL coefficient files, so it is not a pip install.')
    print('--Obtain hwm14.f90 and the data files (hwm123114.bin, dwm07b104i.dat, gd2qd.dat)')
    print('--from NRL, then build the extension module in the directory holding them:')
    print('--    python3 -m numpy.f2py -c hwm14.f90 -m hwm14')
    print('--(numpy 1.26 and later drive f2py through meson, so meson and ninja must be')
    print('-- installed as well, along with a Fortran compiler such as gfortran.)')
    print('--Then run this script with that directory on PYTHONPATH and pass it as')
    print('--the --hwm14-data option, since the model opens its data files from the')
    print('--current directory.')
    print('Program stopped.')
    sys.exit(1)

  return module


def build_hwm14(argument):
  #
  # HWM14（NRL の経験モデル、0-500 km）を格子上で評価する。
  #
  # HWM14 は「子午線方向（北）」と「帯状方向（東）」の風を返す。鉛直風は持たない。
  # 引数は年通日 iyd = (年 - 1900 か 2000)*1000 + 通日、世界時の秒、高度 km、
  # 測地の緯度経度、地方時 stl（-1 で経度と時刻から内部計算）、f107a/f107（静穏成分は
  # 使わないので -1）、ap（2 要素目が擾乱成分の 3 時間 ap。-1 なら静穏のみ）。
  #
  hwm14 = import_hwm14()

  time_epoch = parse_datetime(argument.datetime)

  axis_time, longitude, latitude = axis_from_argument(argument)
  altitude = np.arange(argument.altitude[0], argument.altitude[1] + 0.5*argument.altitude_step,
                       argument.altitude_step)
  if len(altitude) < 1 :
    print('The altitude range is empty:', argument.altitude)
    sys.exit(1)

  print('--Grid: {:d} times x {:d} longitudes x {:d} latitudes x {:d} altitudes'.format(
        len(axis_time), len(longitude), len(latitude), len(altitude)))
  print('--Altitude (km): {:g} to {:g}'.format(altitude[0], altitude[-1]))
  print('--ap (3-hour index for the disturbance model): {:g}'.format(argument.ap))

  ap = np.array([-1.0, float(argument.ap)])
  flag_time = len(axis_time) > 1

  # HWM14 は係数ファイルをカレントディレクトリから開くので、その間だけ移動する
  directory_return = os.getcwd()
  if argument.hwm14_data is not None :
    if not os.path.isdir(argument.hwm14_data) :
      print('The HWM14 data directory does not exist:', argument.hwm14_data)
      sys.exit(1)
    os.chdir(argument.hwm14_data)

  row_list = []
  try:
    for value_time in axis_time:
      time_now = time_epoch + datetime.timedelta(seconds=float(value_time))
      year = time_now.year
      iyd = ( year - (2000 if year > 1999 else 1900) )*1000 + time_now.timetuple().tm_yday
      second = float(time_now.hour*3600 + time_now.minute*60 + time_now.second)
      for value_lon in longitude:
        for value_lat in latitude:
          for value_alt in altitude:
            result = hwm14.hwm14(iyd, second, float(value_alt), float(value_lat), float(value_lon),
                                 -1.0, -1.0, -1.0, ap)
            # result = [meridional (north), zonal (east)]
            row = [float(value_lon), float(value_lat), float(value_alt),
                   float(result[1]), float(result[0]), 0.0]
            row_list.append([float(value_time)] + row if flag_time else row)
  finally:
    os.chdir(directory_return)

  header = ['Source: HWM14 (Horizontal Wind Model 2014, NRL)',
            '--Drob et al., Earth and Space Science 2, 301-319 (2015)',
            'Epoch (UTC): ' + time_epoch.strftime('%Y-%m-%dT%H:%M:%S') + '.000Z',
            'ap for the disturbance component: {:g} (-1 means the quiet component only)'.format(argument.ap),
            'The vertical wind is not modelled and is written as zero']
  if flag_time :
    header.append('Time axis: {:d} steps, 0 to {:g} s from the epoch'.format(len(axis_time), axis_time[-1]))

  return row_list, header


def axis_from_argument(argument):
  #
  # 時刻・経度・緯度の軸を決める。--like を与えれば既にあるテーブルから引き継ぐ
  # （NCEP と HWM14 を合成するには軸が揃っている必要がある）。
  #
  if argument.like is not None :
    axis_time, longitude, latitude, altitude, wind, string_epoch = read_table(argument.like)
    print('--Taking the time and horizontal axes from', argument.like)
    return axis_time, longitude, latitude

  if argument.duration > 0.0 :
    axis_time = np.arange(0.0, argument.duration + 0.5*argument.time_step, argument.time_step)
  else :
    axis_time = np.array([0.0])

  longitude = np.arange(argument.longitude[0], argument.longitude[1] + 0.5*argument.longitude_step,
                        argument.longitude_step) \
              if argument.longitude[1] > argument.longitude[0] else np.array([argument.longitude[0]])
  latitude  = np.arange(argument.latitude[0], argument.latitude[1] + 0.5*argument.latitude_step,
                        argument.latitude_step) \
              if argument.latitude[1] > argument.latitude[0] else np.array([argument.latitude[0]])

  return axis_time, longitude, latitude


def read_table(filename):
  #
  # 書かれたテーブルを読む。**形式を 2 か所に持たないため、ソルバーの実装を使う**
  # （src/wind/wind.py）。合成のときと --like のときに要る。
  #
  if not os.path.exists(filename) :
    print('File not found:', filename)
    sys.exit(1)

  config = {'wind': {'flag_wind': True,
                     'kind_wind_model': 'fileread',
                     'directory_path_specify': 'manual',
                     'directory_wind': os.path.dirname(os.path.abspath(filename)),
                     'filename_wind': os.path.basename(filename)}}
  wind_dict = wind_module.read_wind_file(config)

  return (wind_dict[wind_module.KEY_TIME], wind_dict[wind_module.KEY_LONGITUDE],
          wind_dict[wind_module.KEY_LATITUDE], wind_dict[wind_module.KEY_ALTITUDE],
          wind_dict[wind_module.KEY_WIND], wind_dict[wind_module.KEY_EPOCH])


def build_merge(argument):
  #
  # 2 つのテーブルを高度方向に合成する。**NCEP（下層、地表から約 31 km）と
  # HWM14（上層、0-500 km）を遷移層で繋ぐ**のがこれ。
  #
  # 境目でそのまま切り替えると風速が飛ぶ（対流圏界面のジェットと中間圏の風は
  # 連続しない）ので、--transition A1 A2 の間で線形に重みを移す。
  #
  # 時刻・経度・緯度の軸は 2 つのテーブルで一致していなければならない
  # （HWM14 側を --like で揃えて作る）。高度節点は両者の和集合。
  #
  if argument.input is None or len(argument.input) != 2 :
    print('--source merge needs exactly two input tables (-i low.txt high.txt).')
    sys.exit(1)

  transition_low, transition_high = argument.transition
  if transition_high < transition_low :
    print('--transition must be given as the lower altitude first:', argument.transition)
    sys.exit(1)

  table_low  = read_table(argument.input[0])
  table_high = read_table(argument.input[1])

  for index, name in ((0, 'time'), (1, 'longitude'), (2, 'latitude')):
    if len(table_low[index]) != len(table_high[index]) or \
       not np.allclose(table_low[index], table_high[index], rtol=0.0, atol=1.e-6) :
      print('The two tables do not share the same {} axis.'.format(name))
      print('--Generate the upper one with --like pointing at the lower one, so that')
      print('--the axes match.')
      sys.exit(1)

  altitude_low  = table_low[3]
  altitude_high = table_high[3]
  for name, altitude in (('lower', altitude_low), ('upper', altitude_high)):
    if transition_low < altitude[0] - 1.e-6 or transition_high > altitude[-1] + 1.e-6 :
      print('The transition layer {:g}-{:g} km is not inside the {} table ({:g}-{:g} km).'.format(
            transition_low, transition_high, name, altitude[0], altitude[-1]))
      sys.exit(1)

  axis_time = table_low[0]
  longitude = table_low[1]
  latitude  = table_low[2]
  altitude_node = np.unique(np.concatenate((altitude_low, altitude_high)))
  flag_time = len(axis_time) > 1

  print('--Merging {} (low) and {} (high)'.format(
        os.path.basename(argument.input[0]), os.path.basename(argument.input[1])))
  print('--Transition layer: {:g} to {:g} km'.format(transition_low, transition_high))
  print('--Altitude nodes (km): {:g} to {:g} ({:d} levels)'.format(
        altitude_node[0], altitude_node[-1], len(altitude_node)))

  row_list = []
  for index_t, value_time in enumerate(axis_time):
    for index_i, value_lon in enumerate(longitude):
      for index_j, value_lat in enumerate(latitude):
        for value_alt in altitude_node:
          weight = blend_weight(value_alt, transition_low, transition_high)
          inside_low  = altitude_low[0]  - 1.e-6 <= value_alt <= altitude_low[-1]  + 1.e-6
          inside_high = altitude_high[0] - 1.e-6 <= value_alt <= altitude_high[-1] + 1.e-6
          if not inside_low :
            weight = 1.0
          if not inside_high :
            weight = 0.0
          value = np.zeros(3)
          if weight < 1.0 :
            value = value + (1.0 - weight)*column_value(table_low, index_t, index_i, index_j, value_alt)
          if weight > 0.0 :
            value = value + weight*column_value(table_high, index_t, index_i, index_j, value_alt)
          row = [float(value_lon), float(value_lat), float(value_alt),
                 float(value[0]), float(value[1]), float(value[2])]
          row_list.append([float(value_time)] + row if flag_time else row)

  header = ['Source: merged from {} (below) and {} (above)'.format(
              os.path.basename(argument.input[0]), os.path.basename(argument.input[1])),
            'Transition layer: {:g} to {:g} km, linear in the altitude'.format(transition_low, transition_high)]
  if table_low[5] is not None :
    header.append('Epoch (UTC): ' + table_low[5])
  if table_high[5] is not None and table_high[5] != table_low[5] :
    header.append('--Caution: the upper table carries a different epoch, ' + table_high[5])
  if flag_time :
    header.append('Time axis: {:d} steps, 0 to {:g} s from the epoch'.format(len(axis_time), axis_time[-1]))

  return row_list, header


def blend_weight(altitude, transition_low, transition_high):
  # 上層側の重み。遷移層の下で 0、上で 1
  if altitude <= transition_low :
    return 0.0
  if altitude >= transition_high :
    return 1.0
  if transition_high <= transition_low :
    return 1.0
  return ( altitude - transition_low )/( transition_high - transition_low )


def column_value(table, index_time, index_longitude, index_latitude, altitude):
  # 1 本の鉛直列を高度方向に線形内挿する（範囲外は端でクランプ）
  altitude_axis = table[3]
  wind = table[4][index_time, index_longitude, index_latitude, :, :]
  return np.array([np.interp(altitude, altitude_axis, wind[:,component]) for component in range(0,3)])


def build_profile(argument):
  #
  # 高度と風速の対を並べたテキストから 1 次元プロファイルを作る。
  #
  #   高度[km] 東[m/s] 北[m/s] [上[m/s]]
  #
  if argument.input is None or len(argument.input) != 1 :
    print('--source profile needs exactly one input file (-i).')
    sys.exit(1)
  filename_input = argument.input[0]
  if not os.path.exists(filename_input) :
    print('File not found:', filename_input)
    sys.exit(1)

  row_list = []
  with open(filename_input) as f:
    for line in f:
      stripped = line.strip()
      if len(stripped) == 0 or stripped.startswith('#') :
        continue
      word_list = [float(word) for word in stripped.split()]
      if len(word_list) == 3 :
        word_list.append(0.0)
      if len(word_list) != 4 :
        print('Each row of the profile must hold altitude, east, north and optionally up:', stripped)
        sys.exit(1)
      row_list.append([argument.centre[0], argument.centre[1],
                       word_list[0], word_list[1], word_list[2], word_list[3]])

  if len(row_list) == 0 :
    print('No data row was found in', filename_input)
    sys.exit(1)

  header = ['Source: vertical profile read from ' + os.path.basename(filename_input),
            'The longitude and the latitude have one node each, so the wind depends on the altitude alone']
  if argument.datetime is not None :
    header.append('Epoch (UTC): ' + parse_datetime(argument.datetime).strftime('%Y-%m-%dT%H:%M:%S') + '.000Z')

  return row_list, header


def build_calm(argument):
  # どこでも風速 0。テーブルを読む経路そのものの確認用
  row_list = []
  for altitude in (0.0, argument.altitude_max if argument.altitude_max is not None else 100.0):
    row_list.append([argument.centre[0], argument.centre[1], altitude, 0.0, 0.0, 0.0])

  return row_list, ['Source: calm air (zero wind everywhere), for checking the file path only']


def write_table(filename, row_list, header_list):
  #
  # 行が 6 個なら時刻軸なし、7 個なら先頭がテーブルのエポックからの秒。
  #
  flag_time = len(row_list[0]) == 7

  print('Writing the wind table...:', filename)
  with open(filename, 'w') as f:
    f.write('# Tacode wind table' + '\n')
    for line in header_list:
      f.write('# ' + line + '\n')
    f.write('# Generated: ' + datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%S')
            + 'Z by generate_wind_table.py' + '\n')
    f.write('#' + '\n')
    f.write('# The wind is given in the local geocentric horizon [East, North, Up]' + '\n')
    if flag_time :
      f.write('# Time[s] Longitude[deg.] Latitude[deg.] Altitude[km] East[m/s] North[m/s] Up[m/s]' + '\n')
      f.write('# --Time is measured from the epoch above' + '\n')
    else :
      f.write('# Longitude[deg.] Latitude[deg.] Altitude[km] East[m/s] North[m/s] Up[m/s]' + '\n')
    for row in row_list:
      if flag_time :
        f.write('  {:12.3f} {:12.5f} {:11.5f} {:11.5f} {:12.5f} {:12.5f} {:12.5f}'.format(*row) + '\n')
      else :
        f.write('  {:12.5f} {:11.5f} {:11.5f} {:12.5f} {:12.5f} {:12.5f}'.format(*row) + '\n')

  print('--{:d} grid points{}'.format(len(row_list), ' (with a time axis)' if flag_time else ''))

  return


def main():
  parser = argparse.ArgumentParser(description='Generate a wind table for Tacode.')
  parser.add_argument('--source', type=str, default='ncep',
                      choices=('ncep', 'hwm14', 'profile', 'merge', 'calm'))
  parser.add_argument('-o', '--output', type=str, default='windmodel.txt')
  parser.add_argument('-i', '--input', type=str, nargs='+', default=None,
                      help='input file for --source profile, or the two tables for --source merge')
  parser.add_argument('--datetime', type=str, default=None,
                      help='UTC in ISO 8601, e.g. 2024-01-01T00:00:00Z (required by --source ncep)')
  parser.add_argument('--longitude', type=float, nargs=2, default=[-180.0, 180.0],
                      metavar=('MIN', 'MAX'))
  parser.add_argument('--latitude', type=float, nargs=2, default=[-90.0, 90.0],
                      metavar=('MIN', 'MAX'))
  parser.add_argument('--altitude-max', type=float, default=None,
                      help='drop the levels above this altitude, km')
  parser.add_argument('--centre', type=float, nargs=2, default=[0.0, 0.0],
                      metavar=('LONGITUDE', 'LATITUDE'),
                      help='the single horizontal node of a one-dimensional profile, deg.')
  parser.add_argument('--duration', type=float, default=0.0,
                      help='length of the time axis, s (0 means a single snapshot)')
  parser.add_argument('--time-step', type=float, default=3600.0,
                      help='spacing of the time axis for --source hwm14, s')
  parser.add_argument('--altitude', type=float, nargs=2, default=[0.0, 500.0],
                      metavar=('MIN', 'MAX'), help='altitude range for --source hwm14, km')
  parser.add_argument('--altitude-step', type=float, default=10.0,
                      help='altitude spacing for --source hwm14, km')
  parser.add_argument('--longitude-step', type=float, default=2.5,
                      help='longitude spacing for --source hwm14, deg.')
  parser.add_argument('--latitude-step', type=float, default=2.5,
                      help='latitude spacing for --source hwm14, deg.')
  parser.add_argument('--ap', type=float, default=-1.0,
                      help='3-hour ap index for the HWM14 disturbance model (-1: quiet only)')
  parser.add_argument('--hwm14-data', type=str, default=None,
                      help='directory holding the HWM14 coefficient files')
  parser.add_argument('--like', type=str, default=None,
                      help='take the time and horizontal axes from this existing table')
  parser.add_argument('--transition', type=float, nargs=2, default=[20.0, 30.0],
                      metavar=('LOW', 'HIGH'), help='transition layer for --source merge, km')
  argument = parser.parse_args()

  # 出力先は、HWM14 のためにカレントディレクトリを移す前に絶対パスへ直しておく
  argument.output = os.path.abspath(argument.output)

  if argument.source == 'ncep' :
    if argument.datetime is None :
      print('--source ncep needs --datetime.')
      sys.exit(1)
    row_list, header_list = build_ncep(argument)
  elif argument.source == 'hwm14' :
    if argument.datetime is None :
      print('--source hwm14 needs --datetime.')
      sys.exit(1)
    row_list, header_list = build_hwm14(argument)
  elif argument.source == 'profile' :
    row_list, header_list = build_profile(argument)
  elif argument.source == 'merge' :
    row_list, header_list = build_merge(argument)
  else :
    row_list, header_list = build_calm(argument)

  write_table(argument.output, row_list, header_list)

  return


if __name__ == '__main__':
  main()
