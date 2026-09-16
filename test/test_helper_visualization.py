#!/usr/bin/env python3
"""
src_helper/animate_trajectory/ の可視化ツールのテスト。

ソルバーには触らないツールなので、検査するのは
「出力ファイルを正しく読めるか」「機体形状が壊れていないか」
「CLI が実際に絵を書けるか」の 3 点。

matplotlib は Tacode 本体の依存ではないので、絵を書くテストは
入っていなければ skip する。
"""

import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from context import ROOT_DIR

HELPER_DIR = os.path.join(ROOT_DIR, 'src_helper', 'animate_trajectory')
if HELPER_DIR not in sys.path:
    sys.path.insert(0, HELPER_DIR)

GENERAL_DIR = os.path.join(ROOT_DIR, 'src_helper', 'general')
if GENERAL_DIR not in sys.path:
    sys.path.insert(0, GENERAL_DIR)

import tecplot_reader  # noqa: E402
import vehicle_shape   # noqa: E402
import coastline       # noqa: E402

TECPLOT_3DOF = os.path.join(ROOT_DIR, 'tutorial', 'work', 'output_result', 'tecplot.dat')
TECPLOT_6DOF = os.path.join(ROOT_DIR, 'tutorial', 'work_reentry_6dof', 'output_result', 'tecplot.dat')

try:
    import matplotlib  # noqa: F401
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def animation_writer_is_available(name):
    # ffmpeg は環境によって無いので、あるときだけ mp4 を検査する
    if not HAS_MATPLOTLIB:
        return False
    from matplotlib import animation
    return animation.writers.is_available(name)


class TestTecplotReader(unittest.TestCase):

    def test_reads_the_three_degree_of_freedom_output(self):
        data = tecplot_reader.read_tecplot(TECPLOT_3DOF)

        for name in ('Time', 'X', 'Y', 'Z', 'Long', 'Lati', 'Alti', 'Dens', 'Temp', 'Kn'):
            self.assertIn(name, data)
        self.assertFalse(tecplot_reader.has_attitude(data))
        self.assertEqual(data['Time'][0], 0.0)
        self.assertTrue(np.all(data['Alti'] > 0.0))

    def test_reads_the_six_degree_of_freedom_output(self):
        data = tecplot_reader.read_tecplot(TECPLOT_6DOF)

        self.assertTrue(tecplot_reader.has_attitude(data))
        for name in tecplot_reader.COLUMN_ATTITUDE:
            self.assertIn(name, data)

        quaternion = tecplot_reader.quaternion_array(data)
        self.assertEqual(quaternion.shape[1], 4)
        norm = np.linalg.norm(quaternion, axis=1)
        np.testing.assert_allclose(norm, np.ones_like(norm), atol=1.e-9)

    def test_units_are_stripped_from_the_column_names(self):
        self.assertEqual(tecplot_reader.strip_unit('Alti[km]'), 'Alti')
        self.assertEqual(tecplot_reader.strip_unit(' Long[deg.] '), 'Long')
        self.assertEqual(tecplot_reader.strip_unit('q0'), 'q0')

    def test_position_and_velocity_have_three_components(self):
        data = tecplot_reader.read_tecplot(TECPLOT_6DOF)
        self.assertEqual(tecplot_reader.position_array(data).shape[1], 3)
        self.assertEqual(tecplot_reader.velocity_geocentric_array(data).shape[1], 3)

    def test_a_file_of_several_zones_stops_the_program(self):
        """
        モンテカルロのまとめ出力（1 ケース = 1 ゾーン）を 1 本の軌跡として読むと、
        別ケースの行が繋がって終端点が「最後のケースの最後の点」になる。
        黙って繋げずに止めること。
        """
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'gathered.dat')
            with open(path, 'w') as f:
                f.write('Variables = Time[s],Alti[km]\n')
                f.write('zone t="case0001" i= 2 f=point\n')
                f.write('0.0 100.0\n1.0 90.0\n')
                f.write('zone t="case0002" i= 3 f=point\n')
                f.write('0.0 100.0\n1.0 80.0\n2.0 70.0\n')
            with self.assertRaises(SystemExit):
                tecplot_reader.read_tecplot(path)

            zone_list = tecplot_reader.read_tecplot_zone(path)
            self.assertEqual(len(zone_list), 2)
            self.assertEqual(len(zone_list[0]['Time']), 2)
            self.assertEqual(len(zone_list[1]['Time']), 3)
            self.assertEqual(zone_list[1]['Alti'][-1], 70.0)

    def test_a_single_zone_file_is_read_as_one_zone(self):
        zone_list = tecplot_reader.read_tecplot_zone(TECPLOT_3DOF)
        self.assertEqual(len(zone_list), 1)
        np.testing.assert_array_equal(zone_list[0]['Time'],
                                      tecplot_reader.read_tecplot(TECPLOT_3DOF)['Time'])

    def test_a_missing_file_stops_the_program(self):
        with self.assertRaises(SystemExit):
            tecplot_reader.read_tecplot(os.path.join(ROOT_DIR, 'no_such_file.dat'))

    def test_a_file_without_a_variables_line_stops_the_program(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'broken.dat')
            with open(path, 'w') as f:
                f.write('1.0 2.0 3.0\n')
            with self.assertRaises(SystemExit):
                tecplot_reader.read_tecplot(path)


class TestVehicleShape(unittest.TestCase):

    def test_every_shape_is_built(self):
        for kind in vehicle_shape.LIST_KIND:
            vertex, color = vehicle_shape.get_vertices_and_colors(kind)

            self.assertGreater(len(vertex), 10)
            self.assertEqual(len(vertex), len(color))

            for face in vertex:
                self.assertGreaterEqual(face.shape[0], 3)
                self.assertEqual(face.shape[1], 3)
                self.assertTrue(np.all(np.isfinite(face)))

    def test_shapes_are_normalised_to_about_one(self):
        # 描画側が拡大率を掛けるので、素の大きさは 1 前後に揃っていること
        for kind in vehicle_shape.LIST_KIND:
            vertex, color = vehicle_shape.get_vertices_and_colors(kind)
            point = np.vstack(vertex)
            size = np.max(point.max(axis=0) - point.min(axis=0))
            self.assertGreater(size, 0.5)
            self.assertLess(size, 2.5)

    def mirror_distance(self, point, axis):
        # 鏡像を作り、元の点群からどれだけ離れるかを測る。
        # 対称な形状なら 0 に近くなる
        mirrored = point.copy()
        mirrored[:, axis] = -mirrored[:, axis]
        distance = np.sqrt(((mirrored[:, None, :] - point[None, :, :])**2).sum(axis=2))
        return float(distance.min(axis=1).max())

    def test_shapes_distinguish_front_from_back(self):
        # 前後が見分けられない形状だと、機首の向きが読めない
        for kind in vehicle_shape.LIST_KIND:
            vertex, color = vehicle_shape.get_vertices_and_colors(kind)
            point = np.vstack(vertex)
            with self.subTest(shape=kind):
                self.assertGreater(self.mirror_distance(point, 0), 0.05)

    def test_roll_is_visible_on_every_shape(self):
        #
        # 機体軸まわりに 180 度回した姿を、形か色のどちらかで見分けられること。
        # 回転体（カプセル）は形だけでは区別できないので、色の基準線を入れてある。
        #
        for kind in vehicle_shape.LIST_KIND:
            vertex, color = vehicle_shape.get_vertices_and_colors(kind)
            centroid = np.array([face.mean(axis=0) for face in vertex])

            rotated = centroid.copy()
            rotated[:, 1] = -rotated[:, 1]
            rotated[:, 2] = -rotated[:, 2]

            distance = np.sqrt(((rotated[:, None, :] - centroid[None, :, :])**2).sum(axis=2))
            nearest = distance.argmin(axis=1)

            flag_distinguishable = False
            for index in range(0, len(centroid)):
                if distance[index, nearest[index]] > 0.05 :
                    flag_distinguishable = True     # 形で分かる
                    break
                if color[index] != color[nearest[index]] :
                    flag_distinguishable = True     # 色で分かる
                    break
            with self.subTest(shape=kind):
                self.assertTrue(flag_distinguishable,
                                '{:s} is unchanged by a 180 deg roll'.format(kind))

    def test_an_unknown_shape_is_rejected(self):
        with self.assertRaises(SystemExit):
            vehicle_shape.get_shape('banana')


@unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
class TestCoastline(unittest.TestCase):
    """
    3 次元の図に入れる海岸線（src_helper/general/coastline.py と同梱のテキスト）。

    cartopy / shapely / pyproj を持ち込まないためにテキスト 1 枚に落としてあるので、
    ここで見るのは「読めること」「窓で正しく切り出せること」
    「Line3DCollection に渡せる形になること」「地球の裏側を落とせること」。
    """

    @classmethod
    def setUpClass(cls):
        cls.segment = coastline.read_coastline()

    def test_the_shipped_file_is_read(self):
        self.assertGreater(len(self.segment), 100)
        self.assertGreater(sum(len(line) for line in self.segment), 5000)
        for line in self.segment:
            self.assertEqual(line.shape[1], 2)
            self.assertGreater(len(line), 1)

    def test_the_coordinates_are_degrees_of_longitude_and_latitude(self):
        longitude = np.concatenate([line[:, 0] for line in self.segment])
        latitude = np.concatenate([line[:, 1] for line in self.segment])
        self.assertGreaterEqual(float(np.min(longitude)), -180.0)
        self.assertLessEqual(float(np.max(longitude)), 180.0)
        self.assertGreaterEqual(float(np.min(latitude)), -90.0)
        self.assertLessEqual(float(np.max(latitude)), 90.0)

    def test_a_window_holds_the_coast_it_should(self):
        # 日本のあたりと、太平洋の真ん中（着地点の近く）
        near_japan = coastline.select_range(self.segment, (135.0, 145.0), (30.0, 40.0))
        self.assertGreater(len(near_japan), 0)

        empty = coastline.select_range(self.segment, (-160.0, -159.0), (24.0, 25.0))
        self.assertEqual(len(empty), 0)

    def test_a_polyline_is_cut_where_it_leaves_the_window(self):
        # 切らずに端の点を結ぶと、窓を横切る直線が 1 本引かれてしまう
        line = np.array([[0.0, 0.0], [1.0, 0.0], [50.0, 0.0], [51.0, 0.0],
                         [2.0, 0.0], [3.0, 0.0]])
        selected = coastline.select_range([line], (-5.0, 5.0), (-5.0, 5.0))
        # 断片は 2 つ。それぞれ窓の中の点と、外へ出る 1 点だけ
        self.assertEqual([len(piece) for piece in selected], [3, 3])
        np.testing.assert_allclose(selected[0][:, 0], [0.0, 1.0, 50.0])
        np.testing.assert_allclose(selected[1][:, 0], [51.0, 2.0, 3.0])

    def test_the_window_may_cross_the_180th_meridian(self):
        # 窓は中心 ± 余白で作られるので 180 度を越えることがある
        line = np.array([[179.0, 0.0], [-179.0, 0.0]])
        selected = coastline.select_range([line], (175.0, 185.0), (-5.0, 5.0))
        self.assertEqual(len(selected), 1)
        np.testing.assert_allclose(selected[0][:, 0], [179.0, 181.0])

    def test_the_folding_keeps_the_longitude_inside_the_window(self):
        np.testing.assert_allclose(coastline.fold_longitude(np.array([-179.0]), (175.0, 185.0)),
                                   [181.0])
        np.testing.assert_allclose(coastline.fold_longitude(np.array([179.0]), (-185.0, -175.0)),
                                   [-181.0])

    def test_the_segments_are_pairs_of_points_on_the_sphere(self):
        # Line3DCollection は点数の揃ったセグメントしか受けない（不揃いだと ValueError）
        segment = coastline.get_segment_ecef(self.segment[0:3])
        self.assertEqual(segment.ndim, 3)
        self.assertEqual(segment.shape[1:], (2, 3))
        radius = np.linalg.norm(segment.reshape(-1, 3), axis=1)
        np.testing.assert_allclose(radius, coastline.RADIUS_PLANET, rtol=1.e-12)

        # 折れ線の点の数より 1 つ少ないセグメントに割れること
        self.assertEqual(len(segment), sum(len(line) - 1 for line in self.segment[0:3]))

    def test_nothing_in_the_window_gives_no_segment(self):
        self.assertEqual(len(coastline.get_segment_ecef([])), 0)

    def test_the_far_side_of_the_globe_is_dropped(self):
        # 不透明な地表でも線は透けるので、落とさないと大陸が鏡像で重なる
        radius = coastline.RADIUS_PLANET
        segment = np.array([[[radius, 0.0, 0.0], [radius, 1.0, 0.0]],
                            [[-radius, 0.0, 0.0], [-radius, 1.0, 0.0]]])
        direction = coastline.get_view_direction(0.0, 0.0)
        visible = coastline.select_visible(segment, direction)

        self.assertEqual(len(visible), 1)
        self.assertGreater(float(visible[0, 0, 0]), 0.0)

    def test_the_view_direction_follows_the_camera(self):
        np.testing.assert_allclose(coastline.get_view_direction(0.0, 0.0), [1.0, 0.0, 0.0],
                                   atol=1.e-12)
        np.testing.assert_allclose(coastline.get_view_direction(0.0, 90.0), [0.0, 1.0, 0.0],
                                   atol=1.e-12)
        np.testing.assert_allclose(coastline.get_view_direction(90.0, 0.0), [0.0, 0.0, 1.0],
                                   atol=1.e-12)

    def test_the_planet_radius_agrees_with_the_drawing_tools(self):
        # 二重に持っている値なので、ずれたら困る
        import animate_trajectory
        self.assertEqual(coastline.RADIUS_PLANET, animate_trajectory.RADIUS_PLANET)


@unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
class TestAnimationScript(unittest.TestCase):
    """CLI が実際に絵を書けること。1 コマだけ書かせて確かめる。"""

    SCRIPT = os.path.join(ROOT_DIR, 'src_helper', 'animate_trajectory', 'animate_trajectory.py')

    def run_script(self, arguments):
        return subprocess.run([sys.executable, self.SCRIPT] + arguments,
                              capture_output=True, text=True)

    def test_snapshot_of_a_six_degree_of_freedom_case(self):
        with tempfile.TemporaryDirectory() as directory:
            output = os.path.join(directory, 'snapshot.png')
            completed = self.run_script([TECPLOT_6DOF, '-o', output,
                                         '--snapshot', '300', '--dpi', '50'])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertTrue(os.path.exists(output))
            self.assertGreater(os.path.getsize(output), 1000)
            self.assertIn('attitude columns', completed.stdout)

    def test_a_three_degree_of_freedom_file_is_reported_and_still_drawn(self):
        with tempfile.TemporaryDirectory() as directory:
            output = os.path.join(directory, 'snapshot.png')
            completed = self.run_script([TECPLOT_3DOF, '-o', output,
                                         '--snapshot', '100', '--dpi', '50', '--window', '0'])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertTrue(os.path.exists(output))
            self.assertIn('No attitude columns', completed.stdout)

    def test_the_coastline_can_be_left_off(self):
        with tempfile.TemporaryDirectory() as directory:
            output = os.path.join(directory, 'snapshot.png')
            completed = self.run_script([TECPLOT_6DOF, '-o', output, '--no-coastline',
                                         '--snapshot', '300', '--dpi', '50'])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertGreater(os.path.getsize(output), 1000)

    def test_an_html_animation_is_written_as_a_single_file(self):
        #
        # .html は ffmpeg の要らない動画出力。フレームを base64 の PNG として
        # 埋め込むので、別ディレクトリを作らず 1 ファイルで完結すること。
        #
        with tempfile.TemporaryDirectory() as directory:
            output = os.path.join(directory, 'anim.html')
            completed = self.run_script([TECPLOT_6DOF, '-o', output,
                                         '--frames', '3', '--dpi', '40'])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertEqual(os.listdir(directory), ['anim.html'])

            with open(output) as f:
                page = f.read()
            self.assertEqual(page.count('data:image/png;base64'), 3)

    def test_html_frames_beyond_the_embed_limit_are_reported(self):
        # 上限を超えると HTMLWriter は黙ってフレームを捨てるので、必ず知らせること
        with tempfile.TemporaryDirectory() as directory:
            output = os.path.join(directory, 'anim.html')
            completed = self.run_script([TECPLOT_6DOF, '-o', output, '--frames', '3',
                                         '--dpi', '40', '--embed-limit', '0.01'])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertIn('truncated', completed.stdout)

    @unittest.skipUnless(HAS_MATPLOTLIB and animation_writer_is_available('ffmpeg'),
                         'ffmpeg is not installed')
    def test_an_mp4_animation_is_written(self):
        with tempfile.TemporaryDirectory() as directory:
            output = os.path.join(directory, 'anim.mp4')
            completed = self.run_script([TECPLOT_6DOF, '-o', output,
                                         '--frames', '3', '--dpi', '40'])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertEqual(os.listdir(directory), ['anim.mp4'])
            self.assertGreater(os.path.getsize(output), 1000)

    def test_an_unknown_output_format_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            completed = self.run_script([TECPLOT_6DOF, '-o', os.path.join(directory, 'a.xyz'),
                                         '--frames', '2'])
            self.assertNotEqual(completed.returncode, 0)


if __name__ == '__main__':
    unittest.main()
