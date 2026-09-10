#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/10
#
# モンテカルロ（src/tacode-montecarlo.py）が作ったケース群の軌跡を、
# 分散円（共分散楕円）と一緒にアニメーションにする。
# 計算はしない。各ケースの出力（Tecplot）を読むだけ。
#
#   python3 src_helper/montecarlo_animation/montecarlo_animation.py work_montecarlo_wind \
#       --reference ../work_reentry/output_result/tecplot.dat -o dispersion.mp4
#
# 設定はカレントディレクトリの config_helper.yml の montecarlo_animation セクションに
# 書ける（コマンドラインが優先。--save-config でいまの設定を書き出せる）。
#
# 図は 3 枚。--view で左側が変わる:
#
#   左（--view flat, 既定）— 上が地上軌跡（経度・緯度）。100 ケースは全行程の縮尺では
#          1 本に重なるので、ここは「いまどこを飛んでいるか」の見取り図
#   左（--view 3d）— **経度・緯度・高度の 3 次元**で軌跡そのもの（絶対座標）を描く。
#          分散円は床（高度 0）に、着地点の経度・緯度に置く。100 ケースはこの縮尺では
#          1 本に重なるので束として 1 色。箱の縦横比は見やすさで決めている（軸は等尺ではない）
#   左（--view 3d-relative）— 同じ箱を**基準ケースから見た東・北・高度**で描く。
#          降りながら束がほどけて楕円になるのが見えるのはこちら（散らばりの育ち方が読める）
#   左（--view globe）— **ECEF の絶対座標**（X, Y, Z [km]）で地球ごと描く。地球半径
#          6378 km に対して高度は 150 km なので、--exaggerate で高度だけ引き伸ばせる
#          （1.0 が実寸）。分散円は着地点の接平面に置くので、この縮尺では点にしかならない
#   左（--view follow）— **カメラが基準ケースを追いかける**。地球ごと入れると
#          散らばり（数十 km）が地球半径（6378 km）に埋もれて読めないので、基準ケースの
#          まわりの立方体だけを見る。ECEF の軸のまま、**各ケースの「基準ケースからのずれ」を
#          いまの基準位置に運んで**描く（1 フレームで機体は 90 km 進むので、絶対座標の
#          軌跡は 1 点を残して窓から出てしまう。ずれなら「どこで離れ始めたか」が残り、
#          先頭は本当の現在位置と一致する）。窓は散らばりに合わせて広がる
#          （--window で固定。縮む方向には動かさない）。地表は窓に入ったときだけ、
#          現在位置のまわりを毎フレーム描き直す（animate_trajectory の追従カメラと同じ作り）。
#          3d-relative との違いは、軸が東・北・高度ではなく ECEF のままで、
#          縮尺が 3 方向とも等しく、地表が入ること
#   左下 — 高度の履歴と、散らばり（1 sigma）の時間履歴
#   右   — **基準ケースから見た散らばり**（東・北, km）。各ケースの現在位置を
#          基準ケースの同時刻の位置から測り、共分散楕円（1/2/3 sigma）と CEP 50% を
#          毎フレーム引き直す。雲が 1 点から楕円へ育っていくのが本題
#
# 散らばりを「終端だけ」でなく時々刻々見るのは、風の効きが効いている高度帯を
# 目で読めるようにするため（この設定では 80-150 km と、32 km 以下の終端降下）。
#
# 色は各ケースの東向きの風（Tecplot の WindE 列の初期値）。扇の並び方が
# そのまま風の強弱になるので、ばらつきの向きが読める。

import argparse
import os
import sys

import numpy as np

# Tecplot の読み取りと、ローカル水平系・共分散楕円の計算は既存のものを使う
# （形式と規約を二重に持たないため）
HELPER_DIR = os.path.dirname(os.path.abspath(__file__))
READER_DIR = os.path.normpath(os.path.join(HELPER_DIR, '..', 'animate_trajectory'))
DISPERSION_DIR = os.path.normpath(os.path.join(HELPER_DIR, '..', 'montecarlo_dispersion'))
GENERAL_DIR = os.path.normpath(os.path.join(HELPER_DIR, '..', 'general'))
for directory in (READER_DIR, DISPERSION_DIR, GENERAL_DIR) :
  if directory not in sys.path :
    sys.path.insert(0, directory)

import tecplot_reader as tecplot_reader          # noqa: E402
import montecarlo_dispersion as dispersion       # noqa: E402
import helper_config as helper_config            # noqa: E402

# 設定ファイルの中でこのツールが読むセクション
NAME_SECTION = 'montecarlo_animation'

# 1 枚の HTML にフレームを埋め込んで出力する拡張子（animate_trajectory と同じ）
LIST_EXTENSION_HTML = ('.html', '.htm')

# 表示の色
COLOR_TRACK   = '#9e9e9e'
COLOR_MEAN    = '#d62728'
COLOR_ELLIPSE = '#1f77b4'
COLOR_CEP     = '#d62728'
COLOR_MARK    = ['#2ca02c', '#9467bd', '#8c564b', '#e377c2']
# 地球ごと描くときの軌跡の色。100 本が 1 本に重なって見える縮尺なので、
# 風で色分けしても読めない（淡い中間色は地表に埋もれる）。束として 1 色で描く
COLOR_BUNDLE  = '#d95f02'

# 散らばりの窓（km）に持たせる余白の割合
MARGIN_WINDOW = 0.18

# 追従カメラ（--view follow）の窓。散らばりの最大距離に掛ける余白と、
# 窓の半幅の下限 [km]（散らばりが 0 の突入直後に窓が潰れないようにする）
MARGIN_FOLLOW         = 1.25
WINDOW_FOLLOW_MINIMUM = 2.0

# 地球の赤道半径 [km]（描画用。animate_trajectory の同名の定数と同じ値。あちらは
# 読み込むだけで matplotlib を要求するので、ここでは値だけ持つ。
# 一致は test_montecarlo.py が検査する）
RADIUS_PLANET = 6378.137


def argument():
  parser = argparse.ArgumentParser(
    description='Animate the trajectories of a Monte-Carlo run together with the '
                'dispersion ellipses of the cases.')
  parser.add_argument('directory', nargs='?', default=None,
                      help='the Monte-Carlo working directory (montecarlo.work_dir), '
                           'which holds case0001, case0002, ... and case_template')
  parser.add_argument('-o', '--output', default='montecarlo.mp4',
                      help='output animation (.gif, .html, or .mp4 if ffmpeg is available), '
                           'or an image when --snapshot is given')
  parser.add_argument('--reference', default=None,
                      help='Tecplot output of a reference case (for instance the same '
                           'entry without wind). The dispersion is then measured from it '
                           'instead of from the mean of the cases')
  parser.add_argument('--mark', action='append', default=None, metavar='LABEL=PATH',
                      help='Tecplot output of another run to follow as a star; it is left '
                           'out of the statistics. May be given more than once')
  parser.add_argument('--tecplot', default=dispersion.PATH_TECPLOT_DEFAULT,
                      help='path of the Tecplot output inside a case directory '
                           '(default: %(default)s)')
  parser.add_argument('--view', default='flat',
                      choices=('flat', '3d', '3d-relative', 'globe', 'follow'),
                      help='flat (default): a ground track, the altitude and the dispersion '
                           'plane. 3d: the trajectories themselves in a box of longitude, '
                           'latitude and altitude. 3d-relative: the same box, but measured '
                           'from the reference case at the same time, which is what shows '
                           'the dispersion growing. globe: the trajectories in ECEF, drawn '
                           'with the Earth. follow: the same ECEF trajectories, but with the '
                           'camera tracking the reference case, so that the dispersion is not '
                           'lost against the size of the Earth')
  parser.add_argument('--window', type=float, default=0.0,
                      help='half width of the follow window, km. 0 (the default) widens it '
                           'with the dispersion itself; a value keeps it fixed')
  parser.add_argument('--ground', type=float, default=1.15,
                      help='keep the ground inside the follow window: the half width is at '
                           'least this factor times the altitude of the reference, so the '
                           'surface is in the frame from the first breath of the entry. '
                           '0 leaves the window to the dispersion alone')
  parser.add_argument('--scale-bar', action='store_true',
                      help='draw a scale bar with ticks in the follow view. The window is '
                           'written in the state box in any case')
  parser.add_argument('--exaggerate', type=float, default=1.0,
                      help='stretch the altitude by this factor in the globe view, so that '
                           'a 150 km descent is visible against a 6378 km radius '
                           '(1.0, the default, is the true scale)')
  parser.add_argument('--elevation', type=float, default=None,
                      help='camera elevation of the 3D view, deg (default: 20 for the box, '
                           'and the latitude of the trajectory for the globe)')
  parser.add_argument('--azimuth', type=float, default=None,
                      help='camera azimuth of the 3D view at the first frame, deg (default: '
                           '-62 for the box, and the longitude of the trajectory for the '
                           'globe, so that the track faces the camera). In the follow view '
                           'it is the offset from the longitude the camera is tracking '
                           '(default -55)')
  parser.add_argument('--spin', type=float, default=35.0,
                      help='how far the 3D camera turns over the animation, deg')
  parser.add_argument('--altitude-max', type=float, default=0.0,
                      help='upper limit of the altitude axis of the 3D view, km. '
                           '0 (the default) covers the whole descent')
  parser.add_argument('--frames', type=int, default=180,
                      help='number of frames (the run is resampled to this)')
  parser.add_argument('--fps', type=int, default=20, help='frames per second')
  parser.add_argument('--tail', type=float, default=0.0,
                      help='length of the trail behind each case, s. 0 keeps the whole path')
  parser.add_argument('--dpi', type=int, default=110, help='resolution of the output')
  parser.add_argument('--embed-limit', type=float, default=512.0,
                      help='size limit of the frames embedded in an .html output, MB')
  parser.add_argument('--snapshot', type=float, default=None,
                      help='write a single frame at this time (s) instead of an animation. '
                           'A negative value takes the last frame')
  helper_config.add_argument(parser)

  return helper_config.get_setting(parser, NAME_SECTION)


def read_case(filename):
  #
  # 1 ケース分の軌跡を読む。位置は ECEF [km]、時刻は [s]。
  #
  data = tecplot_reader.read_tecplot(filename)

  # 経度は +-180 度をまたぐと不連続になるので、読んだ時点で連続化しておく
  # （100 ケースの平均を取るときに、+179 と -179 が混じると軌跡が壊れる）。
  # 図の目盛りは normalize_longitude で [-180, 180) に戻して書く
  case = {'time': data['Time'],
          'position': tecplot_reader.position_array(data),
          'longitude': np.degrees(np.unwrap(np.radians(data['Long']))),
          'latitude': data['Lati'],
          'altitude': data['Alti']}

  # 色に使う風（無風のケース、あるいは風を切った計算では 0 とみなす）
  case['wind_east'] = float(data['WindE'][0]) if 'WindE' in data else 0.0

  return case


def resample(case, time_grid):
  #
  # 共通の時刻格子に載せ直す。np.interp は範囲外を端の値で止めるので、
  # 先に着地したケースはその点に留まる（着地後も点として残す）。
  #
  sampled = {'wind_east': case['wind_east']}
  for name in ('longitude', 'latitude', 'altitude') :
    sampled[name] = np.interp(time_grid, case['time'], case[name])
  sampled['position'] = np.stack([np.interp(time_grid, case['time'], case['position'][:, i])
                                  for i in range(0, 3)], axis=1)
  sampled['time_end'] = float(case['time'][-1])

  return sampled


def get_offset(position, position_reference):
  #
  # 基準ケースから見たずれ（東, 北 [km]）。ローカル水平系は基準ケースの
  # **その時刻の位置**で作る。地心の [東, 北, 上] で、ソルバーの入出力と同じ規約。
  #
  number_frame = position_reference.shape[0]
  east  = np.zeros((position.shape[0], number_frame))
  north = np.zeros((position.shape[0], number_frame))

  for n in range(0, number_frame) :
    unit_east, unit_north, _ = dispersion.get_local_horizon(position_reference[n])
    offset     = position[:, n, :] - position_reference[n]
    east[:, n]  = offset.dot(unit_east)
    north[:, n] = offset.dot(unit_north)

  return east, north


def get_ellipse_point(center_east, center_north, semi_major, semi_minor, angle, number=181):
  #
  # 共分散楕円を点の列にする。3 次元の図では床（高度 0）に線として置くので、
  # matplotlib の Ellipse パッチではなく座標が要る。
  #
  radian = np.radians(angle)
  parameter = np.linspace(0.0, 2.0*np.pi, number)
  x = semi_major*np.cos(parameter)
  y = semi_minor*np.sin(parameter)

  return (center_east + x*np.cos(radian) - y*np.sin(radian),
          center_north + x*np.sin(radian) + y*np.cos(radian))


def stretch_altitude(position, factor):
  #
  # 高度だけ factor 倍した位置（ECEF の描画用）。地球半径 6378 km に対して高度は
  # 150 km しかないので、そのままだと軌道が地表に貼り付いて見える。
  # factor = 1.0 なら元のまま（実寸）。
  #
  if factor == 1.0 :
    return position

  radius = np.linalg.norm(position, axis=-1, keepdims=True)
  radius = np.where(radius > 0.0, radius, 1.0)
  altitude = radius - RADIUS_PLANET

  return position*(RADIUS_PLANET + factor*altitude)/radius


def normalize_longitude(longitude):
  # 連続化した経度を目盛りの表示用に [-180, 180) へ戻す
  return (longitude + 180.0) % 360.0 - 180.0


def get_follow_window(position, position_reference, altitude_reference,
                      window_fixed=0.0, factor_ground=0.0):
  #
  # 追従カメラの窓（中心と半幅 [km]）をフレームごとに返す。
  #
  # 中心は基準ケースとケース群の重心の中点に置く。基準ケースにぴったり載せると、
  # ケースは風で片側（風下）に寄るので窓の半分が空くことになる。
  #
  # 半幅は 2 つの要求の大きいほうを取る:
  #
  #   1. 散らばり — 中心からいちばん離れているケースまでの距離に余白を掛けたもの。
  #      基準ケースは別に見なくてよい: 中心が中点なので基準ケースの距離は
  #      |重心 - 基準|/2 であり、最遠のケースの距離を超えられない（重心のずれは
  #      各ケースのずれの平均だから）。**狭める方向には動かさない**（累積の最大）。
  #      散らばりが一時的に縮むたびに寄っては引いてを繰り返すと読めなくなる
  #   2. 地表 — 基準ケースの高度の factor_ground 倍。散らばりだけで決めると、
  #      突入直後は窓が数 km なのに地表は 150 km 下にあり、**何も無い空間に点が
  #      1 つ**という絵になる。高度に追わせておけば、地表が最初から入り、
  #      降りるにつれて窓が閉じ、そのあと散らばりで開いていく
  #
  # --window を与えれば両方を無視して固定する。
  #
  centre = 0.5*( position_reference + np.mean(position, axis=0) )

  if window_fixed > 0.0 :
    return centre, np.full(position_reference.shape[0], window_fixed)

  distance = np.linalg.norm(position - centre[np.newaxis, :, :], axis=2)
  half     = np.maximum.accumulate( MARGIN_FOLLOW*np.max(distance, axis=0) )

  if factor_ground > 0.0 :
    half = np.maximum( half, factor_ground*np.asarray(altitude_reference, dtype=float) )

  return centre, np.maximum(half, WINDOW_FOLLOW_MINIMUM)


def get_scale_length(half_width):
  #
  # ものさしの長さ（km）。窓の半幅に収まる 1, 2, 5 x 10^n のうちいちばん大きいもの。
  # follow では軸目盛りを切っている（地球を入れると枠が邪魔になる）ので、
  # 代わりにこの 1 本で縮尺を示す。
  #
  if half_width <= 0.0 :
    return 0.0

  exponent = np.floor( np.log10(half_width) )
  for factor in (5.0, 2.0, 1.0) :
    length = factor*10.0**exponent
    if length <= half_width :
      return length

  return 10.0**exponent


def get_scale_bar(origin, length, unit_along, unit_up, ratio_tick=0.05):
  #
  # ものさしの折れ線を返す。両端に目盛りを立て、間を NaN で切って 1 本の線にする
  # （3 次元の線は NaN のところで途切れる。artist を増やさないため）。
  #
  tick = ratio_tick*length*unit_up
  end  = origin + length*unit_along
  gap  = np.full(3, np.nan)

  return np.stack([origin + tick, origin - tick, gap,
                   origin, end, gap,
                   end + tick, end - tick])


def clip_to_window(line, centre, half):
  #
  # 窓から外れる点を NaN にして落とす。
  #
  # matplotlib の 3 次元は線を箱で切ってくれないので、そのまま渡すと窓の外の
  # 点まで線が伸びて図が汚れる（animate_trajectory と同じ扱い）。
  #
  line    = np.array(line, dtype=float)
  if line.size == 0 :
    return line

  outside = np.any( np.abs(line - centre) > half, axis=1 )
  line[outside, :] = np.nan

  return line


def get_window(east, north):
  #
  # 散らばりの図の範囲。最後まで入るように全フレームから決め、原点も必ず含める。
  #
  limit = max(float(np.max(np.abs(east))), float(np.max(np.abs(north))), 1.0)

  center_east  = 0.5*(float(np.max(east)) + float(np.min(east)))
  center_north = 0.5*(float(np.max(north)) + float(np.min(north)))
  half = (1.0 + MARGIN_WINDOW)*max(float(np.max(east)) - center_east,
                                   float(np.max(north)) - center_north,
                                   abs(center_east), abs(center_north), 0.5*limit)

  return (center_east - half, center_east + half), (center_north - half, center_north + half)


def main():

  args = argument()

  if helper_config.save_file(args, NAME_SECTION) :
    return

  helper_config.require(args, 'directory', NAME_SECTION)

  if not os.path.isdir(args.directory) :
    print('Directory not found: ', args.directory)
    sys.exit(1)

  # ---------------------------------------------------------------- 読み込み
  case_path_list = dispersion.find_case_directory(args.directory)
  case_list      = []
  name_list      = []
  for case_path in case_path_list :
    filename = os.path.join(case_path, args.tecplot)
    if not os.path.exists(filename) :
      print('--Caution: the case has no result file, skipped: ', filename)
      continue
    case_list.append(read_case(filename))
    name_list.append(os.path.basename(case_path))

  if len(case_list) == 0 :
    print('No case with a result file was found in: ', args.directory)
    sys.exit(1)
  print('Cases read: ', len(case_list))

  # 時刻格子。いちばん長いケースまで覆う
  time_end   = max(float(case['time'][-1]) for case in case_list)
  time_grid  = np.linspace(0.0, time_end, max(args.frames, 2))
  case_list  = [resample(case, time_grid) for case in case_list]

  position  = np.stack([case['position'] for case in case_list], axis=0)
  altitude  = np.stack([case['altitude'] for case in case_list], axis=0)
  longitude = np.stack([case['longitude'] for case in case_list], axis=0)
  latitude  = np.stack([case['latitude'] for case in case_list], axis=0)
  wind_east = np.array([case['wind_east'] for case in case_list])

  # 基準（--reference なら参照ケース、そうでなければケースの平均）
  if args.reference is not None :
    reference = resample(read_case(args.reference), time_grid)
    position_reference = reference['position']
    altitude_reference = reference['altitude']
    text_reference     = dispersion.shorten_path(args.reference)
  else :
    position_reference = np.mean(position, axis=0)
    altitude_reference = np.mean([case['altitude'] for case in case_list], axis=0)
    text_reference     = 'the mean of the cases'

  east, north = get_offset(position, position_reference)

  # モンテカルロの外の計算（統計には入れない）
  mark_list = []
  for text in args.mark or [] :
    label, path = text.split('=', 1) if '=' in text else (os.path.basename(text), text)
    if not os.path.exists(path) :
      print('Mark file not found: ', path)
      sys.exit(1)
    mark = resample(read_case(path), time_grid)
    mark_east, mark_north = get_offset(mark['position'][np.newaxis, :, :], position_reference)
    mark_list.append((label.strip(), mark_east[0], mark_north[0], mark['altitude'],
                      mark['position'], mark['longitude'], mark['latitude']))

  # ---------------------------------------------------------------- 図の骨組み
  try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import animation
    from matplotlib.patches import Ellipse
  except ImportError:
    print('matplotlib is not available, so nothing can be drawn.')
    sys.exit(1)

  flag_follow   = (args.view == 'follow')
  flag_globe    = (args.view == 'globe')
  # ECEF の絶対座標で描く 2 つの view（地球ごとの globe と、追従カメラの follow）
  flag_ecef     = flag_globe or flag_follow
  flag_relative = (args.view == '3d-relative')
  flag_box      = args.view in ('3d', '3d-relative', 'globe', 'follow')

  # 追従窓の中心と半幅（follow のときだけ使う）
  centre_window = half_window = None
  if flag_follow :
    centre_window, half_window = get_follow_window(position, position_reference,
                                                   altitude_reference, args.window, args.ground)

  # カメラ
  # --globe: 既定で軌道の真上に置く（そうしないと地球の裏側を飛ぶことになる）
  # --follow: 仰角は斜め上から（高度差が読めるように）。方位は追いかけている点の
  #   経度からの相対で毎フレーム作るので、ここでは基準だけ持つ
  azimuth_offset = -55.0
  if flag_follow :
    elevation_view = 20.0
    azimuth_view   = 0.0
    if args.azimuth is not None :
      azimuth_offset = args.azimuth
  elif flag_globe :
    centre = position_reference[len(time_grid)//2]
    elevation_view = np.degrees(np.arcsin(centre[2]/np.linalg.norm(centre)))
    azimuth_view   = np.degrees(np.arctan2(centre[1], centre[0]))
  else :
    elevation_view = 20.0
    azimuth_view   = -62.0
  if args.elevation is not None :
    elevation_view = args.elevation
  if args.azimuth is not None and not flag_follow :
    azimuth_view = args.azimuth

  figure = plt.figure(figsize=(13.0, 6.8))
  if flag_box :
    # 3 次元の図を主役にして、平面の散布図と履歴を右に並べる
    grid    = figure.add_gridspec(2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.55, 1.0],
                                  wspace=0.16, hspace=0.34)
    ax_map  = None
    ax_box  = figure.add_subplot(grid[:, 0], projection='3d')
    ax_disp = figure.add_subplot(grid[0, 1])
    ax_alti = figure.add_subplot(grid[1, 1])
  else :
    grid    = figure.add_gridspec(2, 2, width_ratios=[1.0, 1.15], height_ratios=[2.0, 1.0],
                                  wspace=0.38, hspace=0.34)
    ax_box  = None
    ax_map  = figure.add_subplot(grid[0, 0])
    ax_alti = figure.add_subplot(grid[1, 0])
    ax_disp = figure.add_subplot(grid[:, 1])

  colour_map = plt.get_cmap('coolwarm')
  if float(np.ptp(wind_east)) > 0.0 :
    normalize = plt.Normalize(vmin=float(np.min(wind_east)), vmax=float(np.max(wind_east)))
  else :
    normalize = plt.Normalize(vmin=-1.0, vmax=1.0)
  colour = [colour_map(normalize(value)) for value in wind_east]

  # ---- 左上（flat）: 地上軌跡 -----------------------------------------------
  longitude_mean = np.mean([case['longitude'] for case in case_list], axis=0)
  latitude_mean  = np.mean([case['latitude'] for case in case_list], axis=0)
  longitude_reference = (reference['longitude'] if args.reference is not None
                         else longitude_mean)
  latitude_reference  = (reference['latitude'] if args.reference is not None
                         else latitude_mean)

  # 散らばりの図が見る範囲（--mark も入れて決める。散らばりの外に落ちることがある）
  east_window  = np.vstack([east] + [mark[1][np.newaxis, :] for mark in mark_list])
  north_window = np.vstack([north] + [mark[2][np.newaxis, :] for mark in mark_list])
  window_east, window_north = get_window(east_window, north_window)

  line_map = point_map = None
  if ax_map is not None :
    ax_map.plot(longitude_reference, latitude_reference, color=COLOR_TRACK, linewidth=1.4,
                linestyle='-', zorder=2, label='reference ground track')
    ax_map.plot(longitude_reference[0], latitude_reference[0], marker='o', markersize=6,
                color='black', linestyle='none', zorder=4, label='entry')
    line_map, = ax_map.plot([], [], color=COLOR_MEAN, linewidth=2.2, zorder=5)
    point_map, = ax_map.plot([], [], marker='o', markersize=8, color=COLOR_MEAN,
                             markeredgecolor='black', markeredgewidth=0.6,
                             linestyle='none', zorder=6, label='cases (mean)')

    # 散らばりの図が見ている場所を指す（40 km 四方は、この縮尺では点にしかならない）
    ax_map.annotate('impact area\n(right panel)',
                    xy=(longitude_reference[-1], latitude_reference[-1]),
                    xytext=(-70, 34), textcoords='offset points', fontsize=8,
                    ha='center', color=COLOR_ELLIPSE,
                    arrowprops=dict(arrowstyle='->', color=COLOR_ELLIPSE, linewidth=1.0))

    ax_map.xaxis.set_major_formatter(
      matplotlib.ticker.FuncFormatter(
        lambda value, position: '{:g}'.format(normalize_longitude(value))))
    ax_map.set_xlabel('Longitude [deg.]')
    ax_map.set_ylabel('Latitude [deg.]')
    ax_map.set_title('Ground track')
    ax_map.grid(True, linestyle=':', linewidth=0.5)
    ax_map.legend(loc='best', fontsize=8, framealpha=0.9)

  # ---- 左下: 高度 -----------------------------------------------------------
  ax_alti.plot(time_grid, altitude_reference, color='#d9d9d9', linewidth=1.6, zorder=2)
  line_alti, = ax_alti.plot([], [], color='#7f7f7f', linewidth=2.0, zorder=3)
  point_alti, = ax_alti.plot([], [], marker='o', markersize=7, color='#7f7f7f',
                             markeredgecolor='black', markeredgewidth=0.6,
                             linestyle='none', zorder=4)
  ax_alti.set_xlabel('Time [s]')
  ax_alti.set_ylabel('Altitude [km]', color='#7f7f7f')
  ax_alti.tick_params(axis='y', colors='#7f7f7f')
  ax_alti.set_xlim(0.0, time_end)
  ax_alti.set_ylim(0.0, 1.05*float(np.max(altitude_reference)))
  ax_alti.grid(True, linestyle=':', linewidth=0.5)

  # 散らばりの育ち方を同じ時間軸に重ねる。**どの高度帯で風が効いたか**が読める
  # （この設定では 80-150 km を通るあいだと、32 km 以下の終端降下）
  sigma_east  = np.std(east, axis=0, ddof=1) if len(case_list) > 1 else np.zeros(len(time_grid))
  sigma_north = np.std(north, axis=0, ddof=1) if len(case_list) > 1 else np.zeros(len(time_grid))
  ax_sigma = ax_alti.twinx()
  line_sigma_east,  = ax_sigma.plot([], [], color='#b2182b', linewidth=1.8, label='sigma East')
  line_sigma_north, = ax_sigma.plot([], [], color='#2166ac', linewidth=1.8, label='sigma North')
  ax_sigma.set_ylabel('Dispersion, 1 sigma [km]')
  ax_sigma.set_ylim(0.0, 1.15*max(float(np.max(sigma_east)), float(np.max(sigma_north)), 0.1))
  ax_sigma.legend(loc='upper center', fontsize=8, ncol=2, framealpha=0.9)

  # ---- 右: 散らばりと分散円 --------------------------------------------------
  ax_disp.axhline(0.0, color='black', linewidth=0.5, zorder=1)
  ax_disp.axvline(0.0, color='black', linewidth=0.5, zorder=1)
  ax_disp.plot(0.0, 0.0, marker='+', markersize=14, markeredgewidth=1.6, color='black',
               linestyle='none', zorder=6, label='reference')

  trail_list = []
  for i in range(0, len(case_list)) :
    trail, = ax_disp.plot([], [], color=colour[i], linewidth=0.7, alpha=0.55, zorder=3)
    trail_list.append(trail)
  # 点は最初のフレームの位置で作る（空で作ると色の数と合わずに落ちる）
  scatter = ax_disp.scatter(east[:, 0], north[:, 0], s=26, facecolor=colour,
                            edgecolor='black', linewidth=0.4, zorder=7)

  ellipse_list = []
  for k in (1, 2, 3) :
    ellipse = Ellipse((0.0, 0.0), 0.0, 0.0, angle=0.0, facecolor=COLOR_ELLIPSE,
                      alpha=0.16 if k == 1 else 0.08, edgecolor=COLOR_ELLIPSE,
                      linestyle='--', linewidth=1.0, zorder=2,
                      label='1, 2, 3 sigma ellipse' if k == 1 else None)
    ax_disp.add_patch(ellipse)
    ellipse_list.append(ellipse)

  line_cep, = ax_disp.plot([], [], color=COLOR_CEP, linewidth=1.4, zorder=5,
                           label='CEP 50%')
  point_mean, = ax_disp.plot([], [], marker='x', markersize=12, markeredgewidth=2.0,
                             color=COLOR_MEAN, linestyle='none', zorder=8,
                             label='mean of the cases')

  mark_point = []
  mark_trail = []
  for i, label in enumerate([mark[0] for mark in mark_list]) :
    trail, = ax_disp.plot([], [], color=COLOR_MARK[i % len(COLOR_MARK)], linewidth=1.2,
                          alpha=0.8, zorder=8)
    point, = ax_disp.plot([], [], marker='*', markersize=18, linestyle='none',
                          color=COLOR_MARK[i % len(COLOR_MARK)], markeredgecolor='black',
                          markeredgewidth=0.5, zorder=9, label=label)
    mark_trail.append(trail)
    mark_point.append(point)

  ax_disp.set_xlim(*window_east)
  ax_disp.set_ylim(*window_north)
  ax_disp.set_aspect('equal', adjustable='box')
  ax_disp.set_xlabel('East from the reference [km]')
  ax_disp.set_ylabel('North from the reference [km]')
  ax_disp.set_title('Dispersion of the cases')
  ax_disp.grid(True, linestyle=':', linewidth=0.5)
  if not flag_box :
    # 3 次元のときは箱のほうに凡例を置く（散布図が小さくなり、覆ってしまうため）
    ax_disp.legend(loc='upper left', fontsize=9, framealpha=0.9)

  if float(np.ptp(wind_east)) > 0.0 :
    mappable = plt.cm.ScalarMappable(norm=normalize, cmap=colour_map)
    colour_bar = figure.colorbar(mappable, ax=ax_disp, fraction=0.045, pad=0.02)
    colour_bar.set_label('Wind, east component [m/s]')

  # ---- 左（3d）: 東・北・高度の箱 -------------------------------------------
  # 軌跡そのものを 3 次元で描き、**分散円は床（高度 0）に落として**描く。
  # 水平は数十 km、鉛直は 150 km と桁が違うので、箱の縦横比は見やすさで決める
  # （set_box_aspect。軸は等尺ではない）。
  trail_box  = []
  point_box  = []
  floor_box  = []
  mark_box_trail = []
  mark_box_point = []
  shadow_box = None
  line_cep_box = None
  if ax_box is not None :
    for i in range(0, len(case_list)) :
      if flag_relative or flag_follow :
        # 窓を絞っているので 1 本ずつ見分けられる。色は東向きの風
        trail, = ax_box.plot([], [], [], color=colour[i], linewidth=0.8, alpha=0.6)
      else :
        # 絶対座標では 100 本が 1 本に重なるので、束として 1 色で描く
        # （淡い中間色は地表や背景に埋もれて見えない）
        trail, = ax_box.plot([], [], [], color=COLOR_BUNDLE, linewidth=1.6, alpha=0.5,
                             label='cases ({:d})'.format(len(case_list)) if i == 0 else None)
      point, = ax_box.plot([], [], [], marker='o', markersize=4, linestyle='none',
                           color=colour[i], markeredgecolor='black', markeredgewidth=0.3)
      trail_box.append(trail)
      point_box.append(point)

    # 床に落とした影。3 次元では奥行きが読みにくいので、高さの手がかりに置く
    shadow_box, = ax_box.plot([], [], [], marker='.', markersize=3, linestyle='none',
                              color='#c0c0c0')

    for k in (1, 2, 3) :
      floor, = ax_box.plot([], [], [], color=COLOR_ELLIPSE, linewidth=1.0, linestyle='--',
                           label='1, 2, 3 sigma ellipse' if k == 1 else None)
      floor_box.append(floor)
    line_cep_box, = ax_box.plot([], [], [], color=COLOR_CEP, linewidth=1.4, label='CEP 50%')

    for i, label in enumerate([mark[0] for mark in mark_list]) :
      trail, = ax_box.plot([], [], [], color=COLOR_MARK[i % len(COLOR_MARK)],
                           linewidth=1.4 if (flag_relative or flag_follow) else 1.0)
      point, = ax_box.plot([], [], [], marker='*', markersize=14, linestyle='none',
                           color=COLOR_MARK[i % len(COLOR_MARK)], markeredgecolor='black',
                           markeredgewidth=0.5, label=label)
      mark_box_trail.append(trail)
      mark_box_point.append(point)

    if flag_follow :
      # ---- ECEF（絶対座標）＋ 追従カメラ -------------------------------------
      # 地球ごと入れると散らばりが埋もれるので、基準ケースのまわりの立方体だけを見る。
      # 地表と等尺の箱は animate_trajectory のものを使う（形を二重に持たない）。
      # matplotlib を要求する import なので、この分岐に入ってから読む
      import animate_trajectory as animate_trajectory

      ax_box.computed_zorder = False
      for artist in trail_box + point_box + floor_box + mark_box_trail + mark_box_point :
        artist.set_zorder(4)
      shadow_box.set_visible(False)
      line_cep_box.set_zorder(4)

      # 基準ケースの現在位置。軌跡を基準からのずれで描くので、基準そのものは点になる
      point_reference, = ax_box.plot([], [], [], marker='o', markersize=8, linestyle='none',
                                     color=COLOR_TRACK, markeredgecolor='black',
                                     markeredgewidth=0.6, zorder=5, label='reference')
      # 地表は窓に入ったときだけ、現在位置のまわりに作り直す
      ground_follow = {'surface': None, 'graticule': None}

      # ものさし（--scale-bar）。軸目盛りを切っているので、頼まれたらこれで示す。
      # 目盛りと数字を付けないと、ただの線が何を意味するのか分からない
      line_scale = None
      text_scale = {'artist': None}
      if args.scale_bar :
        line_scale, = ax_box.plot([], [], [], color='black', linewidth=1.8, zorder=6)

      # 平行投影にする。透視投影だと窓の中の同じ距離が奥と手前で違って見えるので、
      # 散らばりを読むのに向かない。zoom は箱を画面いっぱいに寄せるため
      ax_box.set_proj_type('ortho')
      ax_box.set_box_aspect((1.0, 1.0, 1.0), zoom=2.2)
      ax_box.set_axis_off()
      ax_box.set_title('Departure from the reference, in ECEF axes\n'
                       '(the camera follows the reference; the window is in the box)',
                       fontsize=10)

    elif flag_globe :
      # ---- ECEF（絶対座標）。地球ごと描く -----------------------------------
      # 地表の点と等尺の箱は animate_trajectory のものを使う（形を二重に持たない）。
      # matplotlib を要求する import なので、この分岐に入ってから読む
      import animate_trajectory as animate_trajectory
      # 地表は不透明にする（半透明だと球の裏側が透けて奥行きが読めない）。
      # 奥行きの順序は computed_zorder = False と zorder で決める
      ax_box.computed_zorder = False
      surface = animate_trajectory.surface_points((-180.0, 180.0), (-90.0, 90.0), 96, 48)
      ax_box.plot_surface(surface[0], surface[1], surface[2], color='#9ecae1',
                          alpha=1.0, linewidth=0.0, edgecolor='none', antialiased=True,
                          shade=True, zorder=0)
      for longitude in np.arange(-180.0, 180.1, 30.0) :
        latitude = np.radians(np.linspace(-90.0, 90.0, 91))
        angle    = np.radians(longitude)
        ax_box.plot(RADIUS_PLANET*np.cos(latitude)*np.cos(angle),
                    RADIUS_PLANET*np.cos(latitude)*np.sin(angle),
                    RADIUS_PLANET*np.sin(latitude),
                    color='#6baed6', linewidth=0.4, zorder=1)
      for latitude in np.arange(-60.0, 60.1, 30.0) :
        longitude = np.radians(np.linspace(-180.0, 180.0, 181))
        angle     = np.radians(latitude)
        ax_box.plot(RADIUS_PLANET*np.cos(angle)*np.cos(longitude),
                    RADIUS_PLANET*np.cos(angle)*np.sin(longitude),
                    RADIUS_PLANET*np.sin(angle)*np.ones_like(longitude),
                    color='#6baed6', linewidth=0.4, zorder=1)

      # 基準ケースの軌道（あれば）を灰色で通しておく
      if args.reference is not None :
        line = stretch_altitude(position_reference, args.exaggerate)
        ax_box.plot(line[:, 0], line[:, 1], line[:, 2], color=COLOR_TRACK,
                    linewidth=1.0, zorder=3)

      for artist in trail_box + point_box + floor_box + mark_box_trail + mark_box_point :
        artist.set_zorder(4)
      shadow_box.set_visible(False)
      line_cep_box.set_zorder(4)

      # 突入点。どちら向きに飛んでいるかの手がかり
      entry = stretch_altitude(position_reference[0], args.exaggerate)
      ax_box.plot([entry[0]], [entry[1]], [entry[2]], marker='o', markersize=6,
                  color='black', linestyle='none', zorder=5, label='entry')

      animate_trajectory.set_equal_box(ax_box, RADIUS_PLANET*1.02)
      # 目盛りは切る。地球を入れると軸の枠が邪魔になるだけで、数字は他の図で読める
      ax_box.set_axis_off()
      title = 'Trajectories in ECEF (graticule every 30 deg.)'
      if args.exaggerate != 1.0 :
        title = title + '\n(the altitude is stretched {:g} times)'.format(args.exaggerate)
      ax_box.set_title(title, fontsize=10)

    else :
      altitude_top = args.altitude_max if args.altitude_max > 0.0 \
                     else 1.02*float(np.max(altitude))
      # matplotlib の 3 次元は線を箱で切ってくれないので、超える点は NaN にして落とす
      altitude_draw = np.where(altitude > altitude_top, np.nan, altitude)
      mark_altitude_draw = [np.where(mark[3] > altitude_top, np.nan, mark[3])
                            for mark in mark_list]
      ax_box.set_zlim(0.0, altitude_top)
      ax_box.set_zlabel('Altitude [km]')

      if flag_relative :
        # ---- 基準ケースから見た東・北・高度 ----------------------------------
        ax_box.set_xlim(*window_east)
        ax_box.set_ylim(*window_north)
        ax_box.set_box_aspect((1.0, 1.0, 0.78))
        ax_box.set_xlabel('East from the reference [km]')
        ax_box.set_ylabel('North from the reference [km]')
        ax_box.set_title('Trajectories relative to the reference\n'
                         '(the ellipses lie on the ground; the axes are not to scale)',
                         fontsize=10)

      else :
        # ---- 経度・緯度・高度（絶対座標） ------------------------------------
        # 軌道そのもの。100 ケースはこの縮尺では 1 本に重なるので、束として描く
        # （散らばりは右の図で読む。床の分散円も点にしかならない）
        margin_longitude = 0.04*max(float(np.ptp(longitude)), 1.0)
        margin_latitude  = 0.04*max(float(np.ptp(latitude)), 1.0)
        ax_box.set_xlim(float(np.min(longitude)) - margin_longitude,
                        float(np.max(longitude)) + margin_longitude)
        ax_box.set_ylim(float(np.min(latitude)) - margin_latitude,
                        float(np.max(latitude)) + margin_latitude)
        ax_box.set_box_aspect((1.0, 1.0, 0.6))
        ax_box.xaxis.set_major_formatter(
          matplotlib.ticker.FuncFormatter(
            lambda value, position: '{:g}'.format(normalize_longitude(value))))
        ax_box.set_xlabel('Longitude [deg.]')
        ax_box.set_ylabel('Latitude [deg.]')
        ax_box.set_title('Trajectories in longitude, latitude and altitude\n'
                         '(the ellipses lie on the ground; the axes are not to scale)',
                         fontsize=10)

        # 基準ケースの軌道と突入点
        if args.reference is not None :
          ax_box.plot(reference['longitude'], reference['latitude'],
                      np.where(reference['altitude'] > altitude_top, np.nan,
                               reference['altitude']),
                      color=COLOR_TRACK, linewidth=1.0, zorder=2)
        ax_box.plot([longitude[0, 0]], [latitude[0, 0]],
                    [min(altitude[0, 0], altitude_top)], marker='o', markersize=6,
                    color='black', linestyle='none', zorder=5, label='entry')

    ax_box.legend(loc='upper left', fontsize=8, framealpha=0.9)
    ax_box.view_init(elev=elevation_view, azim=azimuth_view)

  # 状態の表示は、広いほうの図の隅に置く
  if flag_box :
    # computed_zorder = False の図では、zorder を上げないと軌跡や分散円が文字に重なる
    text_state = ax_box.text2D(0.99, 0.99, '', transform=ax_box.transAxes, fontsize=10,
                               ha='right', va='top', family='monospace', zorder=20,
                               bbox=dict(boxstyle='round', facecolor='white',
                                         edgecolor='lightgray', alpha=0.92))
  else :
    text_state = ax_disp.text(0.98, 0.02, '', transform=ax_disp.transAxes, fontsize=10,
                              ha='right', va='bottom', family='monospace',
                              bbox=dict(boxstyle='round', facecolor='white',
                                        edgecolor='lightgray', alpha=0.92))

  figure.suptitle('Monte-Carlo dispersion: {:d} cases,  reference: {}'
                  .format(len(case_list), text_reference), fontsize=12)

  # ---------------------------------------------------------------- フレーム
  def update(index):

    east_now  = east[:, index]
    north_now = north[:, index]
    east_mean  = float(np.mean(east_now))
    north_mean = float(np.mean(north_now))

    index_first = 0
    if args.tail > 0.0 :
      index_first = max(0, index - int(args.tail/(time_grid[1] - time_grid[0])))

    for i, trail in enumerate(trail_list) :
      trail.set_data(east[i, index_first:index+1], north[i, index_first:index+1])
    scatter.set_offsets(np.column_stack([east_now, north_now]))
    point_mean.set_data([east_mean], [north_mean])

    # 分散円: 共分散楕円（1/2/3 sigma）と CEP 50%
    angle_circle = np.linspace(0.0, 2.0*np.pi, 181)
    semi_major = semi_minor = angle = 0.0
    if len(case_list) > 2 :
      semi_major, semi_minor, angle = dispersion.get_covariance_ellipse(east_now, north_now)
      for k, ellipse in enumerate(ellipse_list, start=1) :
        ellipse.set_center((east_mean, north_mean))
        ellipse.set_width(2.0*k*semi_major)
        ellipse.set_height(2.0*k*semi_minor)
        ellipse.set_angle(angle)

      radius = float(np.median(np.hypot(east_now - east_mean, north_now - north_mean)))
      line_cep.set_data(east_mean + radius*np.cos(angle_circle),
                        north_mean + radius*np.sin(angle_circle))
    else :
      radius = 0.0

    for i, point in enumerate(mark_point) :
      mark_trail[i].set_data(mark_list[i][1][index_first:index+1],
                             mark_list[i][2][index_first:index+1])
      point.set_data([mark_list[i][1][index]], [mark_list[i][2][index]])

    # 3 次元の図（--view 3d / globe / follow のときだけ）
    if ax_box is not None and flag_ecef :
      # 追従カメラでは高度の引き伸ばしは効かせない（窓が数十 km なので、
      # 引き伸ばすと窓の中の位置関係そのものが歪む）
      factor_exaggerate = 1.0 if flag_follow else args.exaggerate

      # 窓の中心（basis は基準ケース。中心はケース群との中点）
      centre_follow = centre_window[index] if flag_follow else None
      half_follow   = float(half_window[index]) if flag_follow else 0.0
      if flag_follow :
        unit_east_follow, unit_north_follow, unit_up_follow \
          = dispersion.get_local_horizon(position_reference[index])

      # ECEF の絶対座標。基準ケースからの差ではなく軌道そのもの
      # --follow だけは軌跡を**基準ケースからのずれ**にして、いまの基準位置に運ぶ。
      #   窓は数十 km、1 フレームで機体は 90 km も進むので、絶対座標の軌跡は
      #   1 点を残して窓から出てしまう。ずれで描けば「どこで離れ始めたか」が
      #   そのまま残る（先頭は本当の現在位置と一致する）
      for i in range(0, len(case_list)) :
        if flag_follow :
          offset = (position[i, index_first:index+1, :]
                    - position_reference[index_first:index+1, :])
          line   = position_reference[index] + offset
        else :
          line = stretch_altitude(position[i, index_first:index+1, :], factor_exaggerate)
        # 窓から外れる点は落とす（3 次元の線は箱で切ってくれない）
        line_draw = clip_to_window(line, centre_follow, half_follow) if flag_follow else line
        trail_box[i].set_data(line_draw[:, 0], line_draw[:, 1])
        trail_box[i].set_3d_properties(line_draw[:, 2])
        point_box[i].set_data([line[-1, 0]], [line[-1, 1]])
        point_box[i].set_3d_properties([line[-1, 2]])

      if flag_follow :
        # 基準ケースはずれ 0 なので点だけ。軌跡の束の起点がこれ
        point_reference.set_data([position_reference[index][0]], [position_reference[index][1]])
        point_reference.set_3d_properties([position_reference[index][2]])

      # 分散円は接平面に置く。globe では着地点の接平面（この縮尺では点にしかならない
      # が位置は正しい）、follow では基準ケースの現在位置（雲がそこに在るので）
      if len(case_list) > 2 :
        unit_east, unit_north, unit_up = dispersion.get_local_horizon(position_reference[index])
        origin = position_reference[index] if flag_follow else RADIUS_PLANET*unit_up
        for k, floor in enumerate(floor_box, start=1) :
          point_east, point_north = get_ellipse_point(east_mean, north_mean,
                                                      k*semi_major, k*semi_minor, angle)
          line = (origin + np.outer(point_east, unit_east) + np.outer(point_north, unit_north))
          floor.set_data(line[:, 0], line[:, 1])
          floor.set_3d_properties(line[:, 2])
        line = (origin
                + np.outer(east_mean + radius*np.cos(angle_circle), unit_east)
                + np.outer(north_mean + radius*np.sin(angle_circle), unit_north))
        line_cep_box.set_data(line[:, 0], line[:, 1])
        line_cep_box.set_3d_properties(line[:, 2])

      for i in range(0, len(mark_list)) :
        if flag_follow :
          line = (position_reference[index]
                  + mark_list[i][4][index_first:index+1, :]
                  - position_reference[index_first:index+1, :])
        else :
          line = stretch_altitude(mark_list[i][4][index_first:index+1, :], factor_exaggerate)
        line_draw = clip_to_window(line, centre_follow, half_follow) if flag_follow else line
        mark_box_trail[i].set_data(line_draw[:, 0], line_draw[:, 1])
        mark_box_trail[i].set_3d_properties(line_draw[:, 2])
        mark_box_point[i].set_data([line[-1, 0]], [line[-1, 1]])
        mark_box_point[i].set_3d_properties([line[-1, 2]])

      if flag_follow :
        # 窓を動かす。地表は窓に入ったときだけ、現在位置のまわりに描き直す
        ax_box.set_xlim(centre_follow[0]-half_follow, centre_follow[0]+half_follow)
        ax_box.set_ylim(centre_follow[1]-half_follow, centre_follow[1]+half_follow)
        ax_box.set_zlim(centre_follow[2]-half_follow, centre_follow[2]+half_follow)

        longitude_follow = np.degrees(np.arctan2(centre_follow[1], centre_follow[0]))
        latitude_follow  = np.degrees(np.arcsin(centre_follow[2]/np.linalg.norm(centre_follow)))

        for key in ('surface', 'graticule') :
          if ground_follow[key] is not None :
            ground_follow[key].remove()
            ground_follow[key] = None
        if line_scale is not None :
          # ものさしは窓の下手前に置く（雲と重ならない場所）。
          # 窓は立方体なので、中心からの**距離**が半幅を超えると窓の外に出る
          # （下へ 0.8、南へ 0.8 と取ると距離は 1.13 倍になり、実際に消えていた）。
          # 下へ 0.5・南へ 0.5、長さは半幅の 0.6 までに収めると距離は 0.77 倍で収まる
          length_scale = get_scale_length(0.6*half_follow)
          origin_scale = (centre_follow - 0.50*half_follow*unit_up_follow
                          - 0.50*half_follow*unit_north_follow
                          - 0.5*length_scale*unit_east_follow)
          line = get_scale_bar(origin_scale, length_scale,
                               unit_east_follow, unit_up_follow)
          line_scale.set_data(line[:, 0], line[:, 1])
          line_scale.set_3d_properties(line[:, 2])

          # 数字は毎フレーム作り直す（Text3D の位置更新は matplotlib の版に依るため）
          if text_scale['artist'] is not None :
            text_scale['artist'].remove()
          label = origin_scale + 0.5*length_scale*unit_east_follow \
                  - 0.11*length_scale*unit_up_follow
          text_scale['artist'] = ax_box.text(label[0], label[1], label[2],
                                             '{:g} km'.format(length_scale),
                                             color='black', fontsize=9, zorder=6,
                                             ha='center', va='top')

        if np.linalg.norm(centre_follow) - half_follow < RADIUS_PLANET :
          # 経緯線の間隔は窓の大きさに合わせる（数十 km の窓に 10 度刻みでは 1 本も入らない）
          interval = max(0.05, round(np.degrees(half_follow/RADIUS_PLANET), 2))
          ground_follow['surface'], ground_follow['graticule'] = \
            animate_trajectory.make_planet(ax_box, longitude_follow, latitude_follow,
                                           half_follow, interval)

    elif ax_box is not None :
      # 東・北（基準ケースから見た差）か、経度・緯度（絶対座標）か
      if flag_relative :
        coordinate_x, coordinate_y = east, north
        mark_x = [mark[1] for mark in mark_list]
        mark_y = [mark[2] for mark in mark_list]
      else :
        coordinate_x, coordinate_y = longitude, latitude
        mark_x = [mark[5] for mark in mark_list]
        mark_y = [mark[6] for mark in mark_list]

      for i in range(0, len(case_list)) :
        trail_box[i].set_data(coordinate_x[i, index_first:index+1],
                              coordinate_y[i, index_first:index+1])
        trail_box[i].set_3d_properties(altitude_draw[i, index_first:index+1])
        point_box[i].set_data([coordinate_x[i, index]], [coordinate_y[i, index]])
        point_box[i].set_3d_properties([altitude_draw[i, index]])

      shadow_box.set_data(coordinate_x[:, index], coordinate_y[:, index])
      shadow_box.set_3d_properties(np.zeros(len(case_list)))

      if len(case_list) > 2 :
        # 床の分散円。絶対座標では、着地点の経度・緯度に km を度へ直して置く
        if flag_relative :
          origin_x, origin_y = 0.0, 0.0
          scale_x = scale_y = 1.0
        else :
          origin_x = float(reference['longitude'][index]) if args.reference is not None \
                     else float(np.mean(longitude[:, index]))
          origin_y = float(reference['latitude'][index]) if args.reference is not None \
                     else float(np.mean(latitude[:, index]))
          scale_y  = np.degrees(1.0/RADIUS_PLANET)
          scale_x  = scale_y/max(np.cos(np.radians(origin_y)), 1.e-6)

        for k, floor in enumerate(floor_box, start=1) :
          point_east, point_north = get_ellipse_point(east_mean, north_mean,
                                                      k*semi_major, k*semi_minor, angle)
          floor.set_data(origin_x + scale_x*point_east, origin_y + scale_y*point_north)
          floor.set_3d_properties(np.zeros(len(point_east)))
        line_cep_box.set_data(
          origin_x + scale_x*(east_mean + radius*np.cos(angle_circle)),
          origin_y + scale_y*(north_mean + radius*np.sin(angle_circle)))
        line_cep_box.set_3d_properties(np.zeros(len(angle_circle)))

      for i in range(0, len(mark_list)) :
        mark_box_trail[i].set_data(mark_x[i][index_first:index+1],
                                   mark_y[i][index_first:index+1])
        mark_box_trail[i].set_3d_properties(mark_altitude_draw[i][index_first:index+1])
        mark_box_point[i].set_data([mark_x[i][index]], [mark_y[i][index]])
        mark_box_point[i].set_3d_properties([mark_altitude_draw[i][index]])

    # カメラをゆっくり回す（立体感の手がかり）。
    # follow では、追いかけている点の経度に張り付けたうえで回す
    if ax_box is not None :
      if len(time_grid) > 1 :
        progress = index/float(len(time_grid)-1)
        azimuth_now = azimuth_view
        if flag_follow :
          azimuth_now = np.degrees(np.arctan2(position_reference[index][1],
                                              position_reference[index][0])) + azimuth_offset
        ax_box.view_init(elev=elevation_view, azim=azimuth_now + args.spin*progress)

    # 地上軌跡と高度（ケースの平均で代表させる。この縮尺では 1 本に重なる）
    if ax_map is not None :
      line_map.set_data(longitude_mean[:index+1], latitude_mean[:index+1])
      point_map.set_data([longitude_mean[index]], [latitude_mean[index]])

    line_alti.set_data(time_grid[:index+1], altitude_reference[:index+1])
    point_alti.set_data([time_grid[index]], [altitude_reference[index]])
    line_sigma_east.set_data(time_grid[:index+1], sigma_east[:index+1])
    line_sigma_north.set_data(time_grid[:index+1], sigma_north[:index+1])

    message = ('t = {:6.0f} s   h = {:6.1f} km\n'
               'mean  E {:+7.2f}  N {:+7.2f} km\n'
               'sigma E {:6.2f}  N {:6.2f} km\n'
               'CEP50 {:6.2f} km'
               .format(time_grid[index], altitude_reference[index], east_mean, north_mean,
                       float(np.std(east_now, ddof=1)) if len(case_list) > 1 else 0.0,
                       float(np.std(north_now, ddof=1)) if len(case_list) > 1 else 0.0,
                       radius))
    if flag_follow :
      # 追従窓は毎フレーム広がるので、いまの縮尺を書いておく
      # 縮尺は窓の半幅で示す（軸目盛りを切っているので、これが唯一の数字になる）
      message = message + '\nwindow +-{:6.2f} km'.format(float(half_window[index]))
    text_state.set_text(message)

    return []

  # ---------------------------------------------------------------- 保存
  if args.snapshot is not None :
    index = len(time_grid)-1 if args.snapshot < 0.0 \
            else int(np.argmin(np.abs(time_grid - args.snapshot)))
    update(index)
    print('Writing snapshot at t = {:.1f} s...: '.format(time_grid[index]), args.output)
    figure.savefig(args.output, dpi=args.dpi)
    print('Done.: {:.2f} MB'.format(os.path.getsize(args.output)/1024.0**2))
    return

  movie = animation.FuncAnimation(figure, update, frames=len(time_grid),
                                  interval=1000.0/float(args.fps), blit=False)

  extension = os.path.splitext(args.output)[1].lower()
  if extension == '.mp4' :
    if not animation.FFMpegWriter.isAvailable() :
      print('ffmpeg was not found, so an mp4 cannot be written.')
      print('--Give a .gif or .html file name instead, or install ffmpeg.')
      sys.exit(1)
    writer = animation.FFMpegWriter(fps=args.fps, bitrate=3200)
  elif extension == '.gif' :
    writer = animation.PillowWriter(fps=args.fps)
  elif extension in LIST_EXTENSION_HTML :
    writer = animation.HTMLWriter(fps=args.fps, embed_frames=True,
                                  default_mode='loop', embed_limit=args.embed_limit)
  else :
    print('Unknown output format:', extension)
    print('--Use .gif or .html, or .mp4 if ffmpeg is installed.')
    sys.exit(1)

  print('Writing animation...: ', args.output, '({:d} frames)'.format(len(time_grid)))
  movie.save(args.output, writer=writer, dpi=args.dpi)

  # HTMLWriter は上限に達しても例外を出さず、残りのフレームを黙って捨てる
  if getattr(writer, '_hit_limit', False) :
    print('--The embedded frames reached the limit of {:g} MB,'.format(args.embed_limit),
          'so the animation is truncated.')
    print('  Raise --embed-limit, or lower --dpi and --frames.')

  print('Done.: {:.1f} MB'.format(os.path.getsize(args.output)/1024.0**2))

  return


if __name__ == '__main__':
  main()
