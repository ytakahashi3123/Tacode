#!/usr/bin/env python3
"""
姿勢計算をオフにしたときに、従来の 3 自由度質点計算と同じ結果になることの回帰テスト。

v2.3.0 で 6 自由度を足したときの受け入れ条件は「3 自由度の出力が 1 バイトも変わらない」
だった。それを一度確かめて終わりにせず、リポジトリに入っている参照出力
（tutorial/work と tutorial/work_reentry）と毎回突き合わせる。

数値そのものを見る唯一のテストなので、ソルバー・空力テーブルの読み込み・出力
ルーチンのいずれかが 3 自由度側の結果を動かしたらここで落ちる。

浮動小数点の最下位ビットは numpy/scipy の版で変わりうるので、比較は
バイト単位ではなく相対誤差 1e-9 で行う（構造とロジックの退行はこれで十分捕まる）。
分母には要素の値ではなく**その列の代表スケール（最大絶対値）**を使う。速度成分は
ゼロを横切るので、値そのものを分母にすると、ほぼ 0 の点で丸め誤差が相対誤差として
無限に効いてしまう。バイト単位の一致を見たいときは CLAUDE.md の手順で直接 diff すること。
"""

import copy
import os
import shutil
import tempfile
import unittest

import numpy as np

from context import ROOT_DIR, load_config, quiet, reference_mismatch

import atmosphere.atmosphere as atmosphere
import satellite.satellite as satellite
import solver.solver as solver
from orbital.orbital import orbital

RELATIVE_TOLERANCE = 1.e-9


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


def run_case(config, directory):
    """3 自由度で 1 ケース走らせ、Tecplot と restart を directory に書く。"""
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

        # 姿勢を渡さない = 従来の質点 3 自由度
        attitude_dict = orb.initial_settings_attitude(config, coordinate_dict, velocity_dict)

        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            solver.solve_equation_motion(config, iteration, time_elapsed,
                                         coordinate_dict, velocity_dict, trajectory_dict,
                                         atmosphere_dict, aerodynamic_dict, attitude_dict)

        orb.output_tecplot(config, iteration, time_elapsed,
                           coordinate_dict, velocity_dict, trajectory_dict, attitude_dict)
        orb.output_restart(config, iteration, time_elapsed,
                           coordinate_dict['cartesian'], velocity_dict['cartesian'], attitude_dict)

    return attitude_dict, iteration


class TestReferenceOutputs(unittest.TestCase):
    """リポジトリに入っている参照出力を再現できること。"""

    CASES = ('work', 'work_reentry')

    def assert_matches_reference(self, rows, rows_ref, label):
        """参照出力と一致することを見る（判定は context.reference_mismatch と共有）。"""
        message = reference_mismatch(rows, rows_ref, label, RELATIVE_TOLERANCE)
        if message is not None:
            self.fail(message)

    def check_case(self, case):
        reference_directory = os.path.join(ROOT_DIR, 'tutorial', case)
        config = load_config(os.path.join(reference_directory, 'config.yml'))

        workdir = tempfile.mkdtemp(prefix='tacode_regression_')
        try:
            attitude_dict, iteration = run_case(config, workdir)

            # 姿勢を頼んでいないので、辞書は作られていないこと
            self.assertIsNone(attitude_dict)

            variables, zone, rows = read_tecplot(
                os.path.join(workdir, 'output_result', 'tecplot.dat'))
            variables_ref, zone_ref, rows_ref = read_tecplot(
                os.path.join(reference_directory, 'output_result', 'tecplot.dat'))

            # 列の構成が増えていないこと（姿勢の列が紛れ込んでいないこと）
            self.assertEqual(variables, variables_ref)
            self.assertEqual(zone, zone_ref)
            self.assertEqual(rows.shape, rows_ref.shape)

            np.testing.assert_allclose(rows[:, 0], rows_ref[:, 0], rtol=0.0, atol=0.0)
            self.assert_matches_reference(rows, rows_ref, 'tecplot.dat')

            # restart も行数と列数を見る（姿勢を書き足していないこと）
            with open(os.path.join(workdir, 'output_restart', 'restart.dat')) as f:
                lines = [line for line in f.read().splitlines() if not line.startswith('#')]
            with open(os.path.join(reference_directory, 'output_restart', 'restart.dat')) as f:
                lines_ref = [line for line in f.read().splitlines() if not line.startswith('#')]

            self.assertEqual(len(lines), len(lines_ref))
            self.assertTrue(all(len(line.split()) == 6 for line in lines))
            self.assert_matches_reference(
                np.array([[float(word) for word in line.split()] for line in lines]),
                np.array([[float(word) for word in line.split()] for line in lines_ref]),
                'restart.dat')
        finally:
            shutil.rmtree(workdir, ignore_errors=True)

    def test_orbital_case(self):
        self.check_case('work')

    def test_reentry_case(self):
        self.check_case('work_reentry')


class TestAttitudeSectionSwitchedOff(unittest.TestCase):
    """
    attitude セクションを持っていても flag_attitude が False なら、
    セクションが無い config と厳密に同じ軌道になること。
    """

    TIME_MAX = 100.0

    def run_trajectory(self, config):
        config = copy.deepcopy(config)
        config['computational_setup']['time_elapsed_maximum'] = self.TIME_MAX
        config['post_process']['tecplot']['flag_output'] = False

        with quiet():
            atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
            aerodynamic_dict = satellite.initial_settings_satellite(config)

            orb = orbital()
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                orb.initial_settings(config)
            attitude_dict = orb.initial_settings_attitude(config, coordinate_dict, velocity_dict)

            solver.solve_equation_motion(config, iteration, time_elapsed,
                                         coordinate_dict, velocity_dict, trajectory_dict,
                                         atmosphere_dict, aerodynamic_dict, attitude_dict)

        return attitude_dict, np.array(coordinate_dict['cartesian']), np.array(velocity_dict['cartesian'])

    def test_identical_to_a_configuration_without_the_section(self):
        config_plain = load_config()

        config_section = copy.deepcopy(config_plain)
        config_section['attitude'] = {
            'flag_attitude': False,
            'inertia_tensor': {'Ixx': 0.75, 'Iyy': 0.50, 'Izz': 0.50},
            'damping_coefficient': {'Clp': 0.0, 'Cmq': -0.5, 'Cnr': -0.5},
            'static_stability_derivative': -0.5,
        }
        config_section['initial_settings']['attitude'] = [90.0, 20.0, 0.0]
        config_section['initial_settings']['angular_velocity'] = [1.0, 2.0, 3.0]

        attitude_plain, coordinate_plain, velocity_plain = self.run_trajectory(config_plain)
        attitude_section, coordinate_section, velocity_section = self.run_trajectory(config_section)

        self.assertIsNone(attitude_plain)
        self.assertIsNone(attitude_section)

        # 丸めではなく完全一致でなければならない（同じ演算列を通るはずなので）
        np.testing.assert_array_equal(coordinate_section, coordinate_plain)
        np.testing.assert_array_equal(velocity_section, velocity_plain)


if __name__ == '__main__':
    unittest.main()
