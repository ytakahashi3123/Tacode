#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/06
#
# Tacode の出力（Tecplot 形式）を読んで、軌道と姿勢のアニメーションを作る。
#
#   python3 src_helper/animate_trajectory/animate_trajectory.py <tecplot.dat> -o attitude.gif
#
# 計算はしない。src/ からは姿勢の座標変換だけを import して使う
# （クォータニオンの規約を二重に持たないため）。
#
# 画面は 4 つに分かれる。
#   左   : ECEF の全体図。地球・軌道・機体（拡大して描く）
#   右上 : ローカル水平系での機体の姿勢。速度ベクトルと機体軸を重ねるので迎角が読める
#   右中 : 高度の履歴
#   右下 : 迎角・横滑り角の履歴

import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

# src/ はリポジトリ直下。このスクリプトは src_helper/<tool>/ に置くので 2 つ上がる
HELPER_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.normpath(os.path.join(HELPER_DIR, '..', '..', 'src'))
if SRC_DIR not in sys.path :
  sys.path.insert(0, SRC_DIR)

import attitude.attitude as attitude   # noqa: E402

sys.path.insert(0, HELPER_DIR)
import tecplot_reader as tecplot_reader  # noqa: E402
import vehicle_shape as vehicle_shape    # noqa: E402

# 地球の赤道半径 [km]（描画用。config.yml の planet.radius と同じ値）
RADIUS_PLANET = 6378.137

# 表示の色
COLOR_PLANET   = '#9ecae1'
COLOR_GRATICULE = '#6baed6'
COLOR_PATH     = '#bdbdbd'
COLOR_TRAIL    = '#d95f02'
COLOR_VELOCITY = '#111111'
COLOR_AXIS     = ['#d62728', '#2ca02c', '#1f77b4']   # 機体 x, y, z

# 1 枚の HTML にフレームを埋め込んで出力する拡張子（matplotlib の HTMLWriter が受け付ける綴り）
LIST_EXTENSION_HTML = ('.html', '.htm')

# NED -> 表示用の ENU（上が上になるように）
MATRIX_NED_TO_ENU = np.array([[0.0, 1.0, 0.0],
                              [1.0, 0.0, 0.0],
                              [0.0, 0.0, -1.0]])


def argument():
  parser = argparse.ArgumentParser(
    description='Animate a Tacode trajectory, with the attitude if the file has it.')
  parser.add_argument('filename', type=str,
                      help='Tecplot file written by Tacode (output_result/tecplot.dat)')
  parser.add_argument('-o', '--output', type=str, default='trajectory.gif',
                      help='output animation (.gif, .html, or .mp4 if ffmpeg is available)')
  parser.add_argument('-s', '--shape', type=str, default=vehicle_shape.KIND_CAPSULE,
                      choices=vehicle_shape.LIST_KIND, help='vehicle shape to draw')
  parser.add_argument('--frames', type=int, default=150,
                      help='number of frames (the trajectory is subsampled to this)')
  parser.add_argument('--fps', type=int, default=20, help='frames per second')
  parser.add_argument('--scale', type=float, default=500.0,
                      help='size of the vehicle drawn on the ground view, km (exaggerated)')
  parser.add_argument('--window', type=float, default=1800.0,
                      help='half width of the ground view, km. 0 shows the whole trajectory')
  parser.add_argument('--fixed-view', action='store_true',
                      help='do not let the camera follow the vehicle')
  parser.add_argument('--elevation', type=float, default=22.0, help='camera elevation, deg')
  parser.add_argument('--dpi', type=int, default=100, help='resolution of the animation')
  parser.add_argument('--embed-limit', type=float, default=512.0,
                      help='size limit of the frames embedded in an .html output, MB')
  parser.add_argument('--snapshot', type=float, default=None,
                      help='write a single frame at this time (s) instead of an animation')
  return parser.parse_args()


def geodetic_range(position, margin_deg=8.0):
  #
  # 軌道が通る経度・緯度の範囲（deg）。地表をどこまで描くかを決めるのに使う。
  # 経度が半周以上に及ぶときは全球を描く。
  #
  radius = np.linalg.norm(position, axis=1)
  longitude = np.degrees(np.arctan2(position[:, 1], position[:, 0]))
  latitude = np.degrees(np.arcsin(position[:, 2]/radius))

  # 経度の不連続（+-180 度）をまたぐ場合に備えて連続化する
  longitude = np.degrees(np.unwrap(np.radians(longitude)))

  span_longitude = longitude.max() - longitude.min()
  if span_longitude > 180.0 :
    return (-180.0, 180.0), (-90.0, 90.0), True

  range_longitude = (longitude.min() - margin_deg, longitude.max() + margin_deg)
  range_latitude = (max(-90.0, latitude.min() - margin_deg),
                    min(90.0, latitude.max() + margin_deg))

  return range_longitude, range_latitude, False


def surface_points(range_longitude, range_latitude, num_longitude, num_latitude):
  longitude = np.radians(np.linspace(range_longitude[0], range_longitude[1], num_longitude))
  latitude = np.radians(np.linspace(range_latitude[0], range_latitude[1], num_latitude))
  grid_lon, grid_lat = np.meshgrid(longitude, latitude)

  return (RADIUS_PLANET*np.cos(grid_lat)*np.cos(grid_lon),
          RADIUS_PLANET*np.cos(grid_lat)*np.sin(grid_lon),
          RADIUS_PLANET*np.sin(grid_lat))


def make_planet(ax, longitude_center, latitude_center, half_width_km, interval_graticule=10.0):
  #
  # 現在位置のまわりの地表と緯度経度線を描き、作った artist を返す。
  #
  # 軌道全体の地表を一度に描くと、経度方向に何十度も広がって球の裏側まで回り込み、
  # 地平線のあたりが読めなくなる。追従窓のまわりだけを毎フレーム描き直す。
  #
  margin_latitude = np.degrees(half_width_km/RADIUS_PLANET) + 3.0
  margin_longitude = margin_latitude/max(0.15, np.cos(np.radians(latitude_center)))

  range_latitude = (max(-90.0, latitude_center - margin_latitude),
                    min(90.0, latitude_center + margin_latitude))
  range_longitude = (longitude_center - margin_longitude, longitude_center + margin_longitude)

  # 不透明にする。半透明だと面の継ぎ目が縞に見え、球の裏側も透けてしまう。
  # 奥行きの順序は ax.computed_zorder = False と zorder で決めている
  x, y, z = surface_points(range_longitude, range_latitude, 40, 32)
  surface = ax.plot_surface(x, y, z, color=COLOR_PLANET, alpha=1.0, linewidth=0.0,
                            edgecolor='none', antialiased=True, shade=False, zorder=0)

  num_point = 40
  segment = []
  for longitude in np.arange(np.ceil(range_longitude[0]/interval_graticule)*interval_graticule,
                             range_longitude[1] + 0.01, interval_graticule):
    latitude = np.radians(np.linspace(range_latitude[0], range_latitude[1], num_point))
    angle = np.radians(longitude)
    segment.append(np.stack([RADIUS_PLANET*np.cos(latitude)*np.cos(angle),
                             RADIUS_PLANET*np.cos(latitude)*np.sin(angle),
                             RADIUS_PLANET*np.sin(latitude)], axis=1))
  for latitude in np.arange(np.ceil(range_latitude[0]/interval_graticule)*interval_graticule,
                            range_latitude[1] + 0.01, interval_graticule):
    longitude = np.radians(np.linspace(range_longitude[0], range_longitude[1], num_point))
    angle = np.radians(latitude)
    segment.append(np.stack([RADIUS_PLANET*np.cos(angle)*np.cos(longitude),
                             RADIUS_PLANET*np.cos(angle)*np.sin(longitude),
                             RADIUS_PLANET*np.sin(angle)*np.ones_like(longitude)], axis=1))

  graticule = None
  if len(segment) > 0 :
    graticule = Line3DCollection(segment, colors=COLOR_GRATICULE, linewidths=0.5, alpha=0.7)
    ax.add_collection3d(graticule)

  return surface, graticule


def make_globe(ax, position):
  # 全球のインセット。どのあたりを飛んでいるかの手掛かり
  x, y, z = surface_points((-180.0, 180.0), (-90.0, 90.0), 48, 24)
  ax.plot_surface(x, y, z, color=COLOR_PLANET, alpha=0.35, linewidth=0.0,
                  edgecolor='none', antialiased=True, shade=False)
  ax.plot(position[:, 0], position[:, 1], position[:, 2], color=COLOR_TRAIL, linewidth=1.2)
  set_equal_box(ax, RADIUS_PLANET*1.05)
  ax.set_axis_off()


def set_equal_box(ax, limit):
  ax.set_xlim(-limit, limit)
  ax.set_ylim(-limit, limit)
  ax.set_zlim(-limit, limit)
  try:
    ax.set_box_aspect((1.0, 1.0, 1.0))
  except AttributeError:
    pass


def matrix_body_to_ecef(quaternion):
  # 機体軸 -> ECEF
  return attitude.quaternion_to_matrix(attitude.quaternion_normalize(quaternion)).T


def matrix_ecef_to_enu(position):
  # 位置（ECEF）から地心ローカル系 [東, 北, 上] への変換
  radius = np.linalg.norm(position)
  longitude = np.arctan2(position[1], position[0])
  latitude = np.arcsin(position[2]/radius)
  return np.dot(MATRIX_NED_TO_ENU, attitude.matrix_ecef_to_ned(longitude, latitude))


def main():
  args = argument()

  data = tecplot_reader.read_tecplot(args.filename)
  flag_attitude = tecplot_reader.has_attitude(data)

  print('Reading:', args.filename)
  print('--Steps: {:d}, time: {:.1f} - {:.1f} s'.format(
        len(data['Time']), data['Time'][0], data['Time'][-1]))
  if flag_attitude :
    print('--The file has the attitude columns; the vehicle will be oriented.')
  else :
    print('--No attitude columns found; only the trajectory will be drawn.')
    print('  (run a case with attitude.flag_attitude: True to get them)')

  position = tecplot_reader.position_array(data)
  index_frame = np.unique(np.linspace(0, len(data['Time'])-1, args.frames).astype(int))

  vertex_body, color_body = vehicle_shape.get_vertices_and_colors(args.shape)

  # ------------------------------------------------------------------ 画面
  figure = plt.figure(figsize=(14.0, 7.6))
  grid = figure.add_gridspec(4, 3, width_ratios=[1.0, 1.0, 0.95],
                             left=0.02, right=0.97, top=0.94, bottom=0.07,
                             wspace=0.25, hspace=0.55)

  ax_orbit = figure.add_subplot(grid[:, 0:2], projection='3d')
  ax_body = figure.add_subplot(grid[0:2, 2], projection='3d')
  ax_altitude = figure.add_subplot(grid[2, 2])
  ax_angle = figure.add_subplot(grid[3, 2])

  # --- 左: 軌道の俯瞰図（地表は現在位置のまわりだけを毎フレーム描く）
  range_longitude, range_latitude, flag_global = geodetic_range(position, margin_deg=2.0)
  interval_graticule = 30.0 if flag_global else 10.0

  # 描画順を自前で決める（既定では奥行きで自動に並べ替えられ、
  # 地表が軌道や機体を覆ってしまう）
  ax_orbit.computed_zorder = False

  ax_orbit.plot(position[:, 0], position[:, 1], position[:, 2],
                color=COLOR_PATH, linewidth=1.2, zorder=1)
  line_trail, = ax_orbit.plot([], [], [], color=COLOR_TRAIL, linewidth=2.6, zorder=2)
  # 遠景では稜線を描かない（小さく描くと線で潰れて灰色の塊になる）
  mesh_orbit = Poly3DCollection([], facecolors=color_body, edgecolors='none',
                                linewidths=0.0, zorder=5)
  ax_orbit.add_collection3d(mesh_orbit)

  # 表示範囲。既定では機体のまわりだけを切り出して追いかける
  # （軌道全体を入れると機体が点になってしまう）
  if args.window > 0.0 :
    half_view = args.window
  else :
    half_view = 0.55*np.max(position.max(axis=0) - position.min(axis=0)) + RADIUS_PLANET*0.05
  center_fixed = 0.5*(position.max(axis=0) + position.min(axis=0))
  try:
    ax_orbit.set_box_aspect((1.0, 1.0, 1.0))
  except AttributeError:
    pass
  ax_orbit.set_axis_off()
  ax_orbit.set_title('Trajectory over the ground (ECEF)\nvehicle drawn {:.0f} km across (exaggerated)'.format(
                     args.scale), fontsize=10)

  ground = {'surface': None, 'graticule': None}

  # 全球のインセット（いまどのあたりか）
  ax_globe = figure.add_axes([0.015, 0.60, 0.17, 0.30], projection='3d')
  make_globe(ax_globe, position)
  marker_globe, = ax_globe.plot([], [], [], 'o', color='#d62728', markersize=5)

  # --- 右上: ローカル水平系での姿勢
  mesh_body = Poly3DCollection([], facecolors=color_body, edgecolors='#333333',
                               linewidths=0.3)
  ax_body.add_collection3d(mesh_body)
  horizon = np.array([[-1.25, -1.25, 0.0], [1.25, -1.25, 0.0],
                      [1.25, 1.25, 0.0], [-1.25, 1.25, 0.0]])
  ax_body.add_collection3d(Poly3DCollection([horizon], facecolors='#dfe6ee',
                                            edgecolors='#9fb0c0', alpha=0.35))
  line_velocity, = ax_body.plot([], [], [], color=COLOR_VELOCITY, linewidth=1.8,
                                linestyle=(0, (6, 3)))
  head_velocity, = ax_body.plot([], [], [], 'o', color=COLOR_VELOCITY, markersize=4)
  # 機体軸は機体より少し長く伸ばす（機首の向きが胴体に隠れないように）
  line_axis = [ax_body.plot([], [], [], color=color, linewidth=2.0)[0] for color in COLOR_AXIS]
  set_equal_box(ax_body, 1.25)
  ax_body.set_axis_off()
  ax_body.view_init(elev=18.0, azim=65.0)
  ax_body.set_title('Attitude in the local horizon frame', fontsize=10)
  for label, vector in (('E', [1.15, 0.0, 0.0]), ('N', [0.0, 1.15, 0.0]), ('Up', [0.0, 0.0, 1.15])):
    ax_body.text(vector[0], vector[1], vector[2], label, fontsize=8, color='#555555')
  ax_body.text2D(0.01, 0.02,
                 'body axes: x red, y green, z blue      velocity: dashed',
                 transform=ax_body.transAxes, fontsize=7, color='#555555')

  # --- 右中: 高度
  ax_altitude.plot(data['Time'], data['Alti'], color='#2166ac', linewidth=1.2)
  cursor_altitude = ax_altitude.axvline(data['Time'][0], color='#d95f02', linewidth=1.0)
  marker_altitude, = ax_altitude.plot([], [], 'o', color='#d95f02', markersize=4)
  ax_altitude.set_ylabel('Altitude [km]', fontsize=9)
  ax_altitude.tick_params(labelsize=8)
  ax_altitude.grid(alpha=0.3)

  # --- 右下: 空力角
  if flag_attitude :
    ax_angle.plot(data['Time'], data['AoA'], color='#1b7837', linewidth=1.0, label='AoA')
    ax_angle.plot(data['Time'], data['Sideslip'], color='#762a83', linewidth=1.0, label='Sideslip')
    ax_angle.legend(fontsize=7, loc='upper right', ncol=2, framealpha=0.7)
    ax_angle.set_ylabel('Angle [deg]', fontsize=9)
  else :
    ax_angle.plot(data['Time'], data['VelplAbs'], color='#1b7837', linewidth=1.0)
    ax_angle.set_ylabel('Velocity [m/s]', fontsize=9)
  cursor_angle = ax_angle.axvline(data['Time'][0], color='#d95f02', linewidth=1.0)
  ax_angle.set_xlabel('Time [s]', fontsize=9)
  ax_angle.tick_params(labelsize=8)
  ax_angle.grid(alpha=0.3)

  text_info = figure.text(0.020, 0.030, '', fontsize=9.5, family='monospace',
                          verticalalignment='bottom',
                          bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                                    edgecolor='#cccccc', alpha=0.85))

  # ------------------------------------------------------------------ 更新
  def update(frame_index):
    n = index_frame[frame_index]
    position_now = position[n]

    line_trail.set_data(position[:n+1, 0], position[:n+1, 1])
    line_trail.set_3d_properties(position[:n+1, 2])

    marker_globe.set_data([position_now[0]], [position_now[1]])
    marker_globe.set_3d_properties([position_now[2]])

    if flag_attitude :
      quaternion = np.array([data['q0'][n], data['q1'][n], data['q2'][n], data['q3'][n]])
      matrix_eb = matrix_body_to_ecef(quaternion)

      # 左: ECEF の中に拡大して置く
      face_orbit = [position_now + args.scale*np.dot(vertex, matrix_eb.T) for vertex in vertex_body]
      mesh_orbit.set_verts(face_orbit)

      # 右上: ローカル水平系（表示は東・北・上）
      matrix_le = matrix_ecef_to_enu(position_now)
      matrix_lb = np.dot(matrix_le, matrix_eb)
      mesh_body.set_verts([np.dot(vertex, matrix_lb.T) for vertex in vertex_body])

      velocity = np.array([data['Upl'][n], data['Vpl'][n], data['Wpl'][n]])
      velocity = velocity/np.linalg.norm(velocity)

      # 速度方向の斜め前から見る（進行方向が変わってもいつも同じ向きに見える）。
      # 65 度ずらすと横からに近くなり、ピッチ（迎角）が読みやすい
      ax_body.view_init(elev=18.0,
                        azim=np.degrees(np.arctan2(velocity[1], velocity[0])) + 65.0)
      line_velocity.set_data([-1.1*velocity[0], 1.1*velocity[0]],
                             [-1.1*velocity[1], 1.1*velocity[1]])
      line_velocity.set_3d_properties([-1.1*velocity[2], 1.1*velocity[2]])
      head_velocity.set_data([1.1*velocity[0]], [1.1*velocity[1]])
      head_velocity.set_3d_properties([1.1*velocity[2]])

      for index in range(0, 3):
        unit = np.zeros(3)
        unit[index] = 1.0
        direction = np.dot(matrix_lb, unit)
        line_axis[index].set_data([0.0, 1.15*direction[0]], [0.0, 1.15*direction[1]])
        line_axis[index].set_3d_properties([0.0, 1.15*direction[2]])
    else :
      face_orbit = [position_now + args.scale*vertex for vertex in vertex_body]
      mesh_orbit.set_verts(face_orbit)

    cursor_altitude.set_xdata([data['Time'][n], data['Time'][n]])
    cursor_angle.set_xdata([data['Time'][n], data['Time'][n]])
    marker_altitude.set_data([data['Time'][n]], [data['Alti'][n]])

    longitude = np.degrees(np.arctan2(position_now[1], position_now[0]))
    latitude = np.degrees(np.arcsin(position_now[2]/np.linalg.norm(position_now)))
    ax_globe.view_init(elev=20.0, azim=longitude - 20.0)

    # 地表を現在位置のまわりに作り直す
    for key in ('surface', 'graticule'):
      if ground[key] is not None :
        ground[key].remove()
    ground['surface'], ground['graticule'] = make_planet(
      ax_orbit, longitude, latitude, half_view, interval_graticule)

    center_view = position_now if args.window > 0.0 else center_fixed
    ax_orbit.set_xlim(center_view[0]-half_view, center_view[0]+half_view)
    ax_orbit.set_ylim(center_view[1]-half_view, center_view[1]+half_view)
    ax_orbit.set_zlim(center_view[2]-half_view, center_view[2]+half_view)

    if not args.fixed_view :
      ax_orbit.view_init(elev=args.elevation, azim=longitude - 55.0)

    message = 'time {:8.1f} s   alt {:8.2f} km   vel {:8.1f} m/s'.format(
              data['Time'][n], data['Alti'][n], data['VelplAbs'][n])
    if flag_attitude :
      message += '\nyaw {:7.2f}   pitch {:7.2f}   roll {:7.2f}  [deg]'.format(
                 data['Yaw'][n], data['Pitch'][n], data['Roll'][n])
      message += '\nAoA {:7.2f}   sides {:7.2f}   total {:7.2f}  [deg]'.format(
                 data['AoA'][n], data['Sideslip'][n], data['AoAtotal'][n])
      message += '\np   {:7.2f}   q     {:7.2f}   r     {:7.2f}  [deg/s]'.format(
                 data['P'][n], data['Q'][n], data['R'][n])
    text_info.set_text(message)

    return []

  # ------------------------------------------------------------------ 保存
  # 静止画: 論文や報告書の図に使う 1 コマ
  if args.snapshot is not None :
    index_snapshot = int(np.argmin(np.abs(data['Time'] - args.snapshot)))
    index_frame = np.array([index_snapshot])
    update(0)
    print('Writing snapshot at t = {:.1f} s...:'.format(data['Time'][index_snapshot]), args.output)
    figure.savefig(args.output, dpi=args.dpi)
    print('Done.')
    return

  movie = animation.FuncAnimation(figure, update, frames=len(index_frame),
                                  interval=1000.0/float(args.fps), blit=False)

  extension = os.path.splitext(args.output)[1].lower()
  if extension == '.mp4' :
    if not animation.FFMpegWriter.isAvailable() :
      print('ffmpeg was not found, so an mp4 cannot be written.')
      print('--Give a .gif or .html file name instead, or install ffmpeg.')
      sys.exit(1)
    writer = animation.FFMpegWriter(fps=args.fps, bitrate=2400)
  elif extension == '.gif' :
    writer = animation.PillowWriter(fps=args.fps)
  elif extension in LIST_EXTENSION_HTML :
    # ブラウザで再生する 1 枚の HTML。ffmpeg が要らず、GIF より画質が落ちない。
    # embed_frames でフレームを base64 の PNG として埋め込む（別ディレクトリを作らない）。
    # 既定の埋め込み上限は 20 MB で、超えたフレームは黙って捨てられるので明示的に上げる
    writer = animation.HTMLWriter(fps=args.fps, embed_frames=True,
                                  default_mode='loop', embed_limit=args.embed_limit)
  else :
    print('Unknown output format:', extension)
    print('--Use .gif or .html, or .mp4 if ffmpeg is installed.')
    sys.exit(1)

  print('Writing animation...:', args.output, '({:d} frames)'.format(len(index_frame)))
  movie.save(args.output, writer=writer, dpi=args.dpi)

  # HTMLWriter は上限に達しても例外を出さず、残りのフレームを黙って捨てる
  if getattr(writer, '_hit_limit', False) :
    print('--The embedded frames reached the limit of {:g} MB,'.format(args.embed_limit),
          'so the animation is truncated.')
    print('  Raise --embed-limit, or lower --dpi and --frames.')

  print('Done.: {:.1f} MB'.format(os.path.getsize(args.output)/1024.0**2))


if __name__ == '__main__':
  main()
