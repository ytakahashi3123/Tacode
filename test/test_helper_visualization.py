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

import tecplot_reader  # noqa: E402
import vehicle_shape   # noqa: E402

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
