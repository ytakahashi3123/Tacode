#!/usr/bin/env python3
"""
6 自由度（姿勢込み）チュートリアルの参照出力との回帰テスト。

test_regression_3dof.py が「姿勢オフ = 従来どおり」を守るのに対し、こちらは
「姿勢オンの結果が動いていない」ことを守る。tutorial/work_reentry_6dof の
出力 3 ファイル（tecplot.dat / restart.dat / geodetic.kml）はリポジトリに
入っているので、それを毎回作り直して突き合わせる。

姿勢まわりの単体テストは解析解・保存量・対称性を見るので、そこを通っていても
「結合の係数がわずかに変わった」「出力の列がずれた」といった退行は捕まらない。
数値そのものを 6 自由度側で見るのはこのテストだけである。

判定は 3 自由度側と同じ context.reference_mismatch（列の最大絶対値を分母に
した相対誤差 1e-9）。クォータニオンも角速度もゼロを横切るので、要素の値を
分母にすると使えない。

このケースは 1000 s / 20000 ステップを実際に走らせるので、1 回で 13 秒ほど
かかる。テスト 3 つで setUpClass の 1 回に集約してある。
"""

import copy
import os
import re
import shutil
import tempfile
import unittest

import numpy as np

from context import ROOT_DIR, load_config, quiet, reference_mismatch

import atmosphere.atmosphere as atmosphere
import output_gpsdata.output_gpsdata as output_gpsdata
import satellite.satellite as satellite
import solver.solver as solver
from orbital.orbital import orbital

RELATIVE_TOLERANCE = 1.e-9

CASE = 'work_reentry_6dof'
REFERENCE_DIR = os.path.join(ROOT_DIR, 'tutorial', CASE)

# Tecplot の列数（時刻・座標・大気量に姿勢 13 列が続く）
NUMBER_COLUMN_TECPLOT = 27
# restart の列数（座標 3・速度 3・クォータニオン 4・角速度 3）
NUMBER_COLUMN_RESTART = 13


def read_tecplot(path):
    """Tecplot ファイルを (Variables 行, zone 行, 数値の配列) にばらす。"""
    variables = None
    zone = None
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('Variables'):
                variables = line.strip()
            elif line.startswith('zone'):
                zone = line.strip()
            elif not line.startswith('#') and line.split():
                rows.append([float(word) for word in line.split()])
    return variables, zone, np.array(rows)


def read_restart(path):
    """restart ファイルを (コメント行, 数値の配列) にばらす。"""
    comment = []
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                comment.append(line.strip())
            elif line.split():
                rows.append([float(word) for word in line.split()])
    return comment, np.array(rows)


def read_kml_coordinates(path):
    """KML の <coordinates> を (経度, 緯度, 高度) の配列にする。"""
    with open(path) as f:
        text = f.read()
    match = re.search(r'<coordinates>(.*?)</coordinates>', text, re.DOTALL)
    if match is None:
        raise AssertionError('<coordinates> が見つからない: ' + path)
    return np.array([[float(word) for word in point.split(',')]
                     for point in match.group(1).split()])


def run_case(config, directory):
    """6 自由度で 1 ケース走らせ、出力 3 ファイルを directory に書く。"""
    config = copy.deepcopy(config)
    config['post_process']['directory_output'] = os.path.join(directory, 'output_result')
    config['restart_process']['directory_output'] = os.path.join(directory, 'output_restart')

    with quiet():
        atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
        aerodynamic_dict = satellite.initial_settings_satellite(config)

        orb = orbital()
        orb.make_directory_output(config)

        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            orb.initial_settings(config)
        attitude_dict = orb.initial_settings_attitude(config, coordinate_dict, velocity_dict)

        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            solver.solve_equation_motion(config, iteration, time_elapsed,
                                         coordinate_dict, velocity_dict, trajectory_dict,
                                         atmosphere_dict, aerodynamic_dict, attitude_dict)

        output_gpsdata.output_routine(config, iteration, coordinate_dict, velocity_dict)
        orb.output_tecplot(config, iteration, time_elapsed,
                           coordinate_dict, velocity_dict, trajectory_dict, attitude_dict)
        orb.output_restart(config, iteration, time_elapsed,
                           coordinate_dict['cartesian'], velocity_dict['cartesian'], attitude_dict)

    return attitude_dict, iteration


class TestReferenceOutputs(unittest.TestCase):
    """work_reentry_6dof の参照出力を再現できること。"""

    @classmethod
    def setUpClass(cls):
        cls.workdir = tempfile.mkdtemp(prefix='tacode_regression_6dof_')
        config = load_config(os.path.join(REFERENCE_DIR, 'config.yml'))
        cls.attitude_dict, cls.iteration = run_case(config, cls.workdir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.workdir, ignore_errors=True)

    def assert_matches_reference(self, rows, rows_ref, label):
        message = reference_mismatch(rows, rows_ref, label, RELATIVE_TOLERANCE)
        if message is not None:
            self.fail(message)

    def test_the_attitude_was_actually_solved(self):
        # 姿勢を頼んだのだから辞書が作られ、座標と同じ長さで揃っていること
        self.assertIsNotNone(self.attitude_dict)
        for key in (orbital.KEY_ATTITUDE_QUATERNION, orbital.KEY_ATTITUDE_OMEGA):
            self.assertEqual(len(self.attitude_dict[key]), self.iteration+1)

    def test_tecplot_matches_the_reference(self):
        variables, zone, rows = read_tecplot(
            os.path.join(self.workdir, 'output_result', 'tecplot.dat'))
        variables_ref, zone_ref, rows_ref = read_tecplot(
            os.path.join(REFERENCE_DIR, 'output_result', 'tecplot.dat'))

        # 列の構成と点数（zone ヘッダの i=）が変わっていないこと
        self.assertEqual(variables, variables_ref)
        self.assertEqual(zone, zone_ref)
        self.assertEqual(rows.shape, rows_ref.shape)
        self.assertEqual(rows.shape[1], NUMBER_COLUMN_TECPLOT)

        # 時刻は演算列が同じなので厳密一致でなければならない
        np.testing.assert_allclose(rows[:, 0], rows_ref[:, 0], rtol=0.0, atol=0.0)
        self.assert_matches_reference(rows, rows_ref, 'tecplot.dat')

    def test_restart_matches_the_reference(self):
        comment, rows = read_restart(
            os.path.join(self.workdir, 'output_restart', 'restart.dat'))
        comment_ref, rows_ref = read_restart(
            os.path.join(REFERENCE_DIR, 'output_restart', 'restart.dat'))

        # 姿勢の列を書いている旨の見出しと、出力間隔の行まで含めて一致すること
        self.assertEqual(comment, comment_ref)
        self.assertEqual(rows.shape, rows_ref.shape)
        self.assertEqual(rows.shape[1], NUMBER_COLUMN_RESTART)

        self.assert_matches_reference(rows, rows_ref, 'restart.dat')

    def test_kml_matches_the_reference(self):
        point = read_kml_coordinates(
            os.path.join(self.workdir, 'output_result', 'geodetic.kml'))
        point_ref = read_kml_coordinates(
            os.path.join(REFERENCE_DIR, 'output_result', 'geodetic.kml'))

        self.assertEqual(point.shape, point_ref.shape)

        # 経度・緯度は普段どおりの相対誤差で見る
        self.assert_matches_reference(point[:, 0:2], point_ref[:, 0:2], 'geodetic.kml')

        # 高度は output_kml が int() で切り捨てるので、最下位ビットの差が
        # そのまま 1 m の飛びになりうる。ここだけ絶対 1 m で見る
        self.assertTrue(np.all(np.abs(point[:, 2] - point_ref[:, 2]) <= 1.0),
                        'geodetic.kml の高度が参照から 1 m 以上ずれている')


if __name__ == '__main__':
    unittest.main()
