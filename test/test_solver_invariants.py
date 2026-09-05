#!/usr/bin/env python3
"""
ソルバーと出力の不変条件のテスト。

いずれも過去に実在したバグの回帰テストである。
- B-2  時間ループが 1 ステップ余分に回る
- B-3  軌道量が 1 ステップずれて出力される / 配列長が揃わない
- D-11 Tecplot の zone ヘッダの点数が実際と合わない
- D-13 出力で最終ステップが欠落する
- A-2  constant モデルが実行できない
"""

import os
import shutil
import tempfile
import unittest

import numpy as np

from context import ROOT_DIR, SRC_DIR, load_config, quiet

import atmosphere.atmosphere as atmosphere
import satellite.satellite as satellite
import solver.solver as solver
from orbital.orbital import orbital

TIME_MAX = 30.0
TIMESTEP = 1.0


def run_solver(config):
    """初期化からソルバーまでを回し、辞書ごと返す。"""
    with quiet():
        atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
        aerodynamic_dict = satellite.initial_settings_satellite(config)

        orb = orbital()
        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            orb.initial_settings(config)

        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            solver.solve_equation_motion(
                config, iteration, time_elapsed,
                coordinate_dict, velocity_dict, trajectory_dict,
                atmosphere_dict, aerodynamic_dict)

    return {
        'iteration': iteration,
        'time_elapsed': time_elapsed,
        'coordinate': coordinate_dict,
        'velocity': velocity_dict,
        'trajectory': trajectory_dict,
        'atmosphere': atmosphere_dict,
    }


def short_config():
    config = load_config()
    config['computational_setup']['time_elapsed_maximum'] = TIME_MAX
    config['time_integration']['timestep_constant'] = TIMESTEP
    return config


class TestTimeLoop(unittest.TestCase):
    """B-2: ループが指定時間ちょうどで止まること。"""

    @classmethod
    def setUpClass(cls):
        cls.result = run_solver(short_config())

    def test_iteration_count(self):
        self.assertEqual(self.result['iteration'], int(TIME_MAX / TIMESTEP),
                         '時間ループのステップ数が指定と合わない')

    def test_elapsed_time_is_returned_and_correct(self):
        """A-3: solve_equation_motion が経過時間を返すこと。"""
        self.assertAlmostEqual(self.result['time_elapsed'], TIME_MAX, places=9)

    def test_elapsed_time_has_no_accumulated_error(self):
        """経過時間を累積加算せずに求めていること。"""
        config = short_config()
        config['time_integration']['timestep_constant'] = 0.1
        config['computational_setup']['time_elapsed_maximum'] = 100.0
        result = run_solver(config)
        # 0.1 を 1000 回足すと 1e-13 程度ずれる。オフセット計算なら厳密に一致する
        self.assertEqual(result['time_elapsed'], 100.0)


class TestArrayLengths(unittest.TestCase):
    """B-3 / D-13: 返り値の配列長がすべて iteration+1 で揃うこと。"""

    @classmethod
    def setUpClass(cls):
        cls.result = run_solver(short_config())

    def test_all_arrays_have_the_same_length(self):
        expected = self.result['iteration'] + 1
        lengths = {
            'coordinate/cartesian': len(self.result['coordinate']['cartesian']),
            'coordinate/geodetic': len(self.result['coordinate']['geodetic']),
            'velocity/cartesian': len(self.result['velocity']['cartesian']),
            'velocity/polar': len(self.result['velocity']['polar']),
            'trajectory/density': len(self.result['trajectory']['density']),
            'trajectory/temperature': len(self.result['trajectory']['temperature']),
            'trajectory/knudsen': len(self.result['trajectory']['knudsen']),
        }
        for name, length in lengths.items():
            self.assertEqual(length, expected, '%s の長さが iteration+1 と違う' % name)


class TestTrajectoryAlignment(unittest.TestCase):
    """
    B-3: n 行目の大気量が n 行目の位置に対応していること。

    出力の高度から補間器で密度を再計算し、格納された値と突き合わせる。
    RK4 の 4 段目の値を出力していた頃は 0.4% ほどずれていた。
    """

    @classmethod
    def setUpClass(cls):
        cls.result = run_solver(short_config())

    def test_density_matches_the_altitude_of_the_same_index(self):
        atm = self.result['atmosphere']
        geodetic = self.result['coordinate']['geodetic']
        density = self.result['trajectory']['density']

        worst = 0.0
        for n in range(len(density)):
            altitude_km = geodetic[n][2] * orbital.m2km
            expected, _, _ = atmosphere.get_atmosphere_property(altitude_km, atm)
            worst = max(worst, abs(float(density[n]) - float(expected)) / float(expected))

        self.assertLess(worst, 1.0e-12,
                        '大気量が位置と 1 ステップずれている（最大相対差 %.3e）' % worst)


class TestOutputFiles(unittest.TestCase):
    """D-11 / D-13: 出力ファイルの点数とヘッダが実態と合うこと。"""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix='tacode_test_')
        cls.cwd = os.getcwd()
        os.chdir(cls.tmpdir)

        config = short_config()
        cls.result = run_solver(config)

        with quiet():
            orb = orbital()
            orb.make_directory_output(config)
            orb.output_tecplot(config, cls.result['iteration'], cls.result['time_elapsed'],
                               cls.result['coordinate'], cls.result['velocity'],
                               cls.result['trajectory'])
        cls.tecplot = os.path.join(config['post_process']['directory_output'],
                                   config['post_process']['tecplot']['filename_output'])
        cls.frequency = config['post_process']['tecplot']['frequency_output']

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.cwd)
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def _read_rows(self):
        rows = []
        with open(self.tecplot) as f:
            for line in f:
                words = line.split()
                if len(words) >= 14 and not line.startswith(('#', 'Variables', 'zone')):
                    rows.append([float(w) for w in words])
        return rows

    def test_zone_header_matches_the_number_of_rows(self):
        with open(self.tecplot) as f:
            header = [line for line in f if line.startswith('zone')][0]

        declared = header.split('i=')[1].split()[0]
        self.assertEqual(int(declared), len(self._read_rows()),
                         'zone ヘッダの点数と実際の行数が違う')
        self.assertNotIn('.', declared, 'zone ヘッダの点数が整数でない')

    def test_final_step_is_present(self):
        rows = self._read_rows()
        self.assertAlmostEqual(rows[0][0], 0.0, places=9)
        self.assertAlmostEqual(rows[-1][0], TIME_MAX, places=9,
                               msg='出力に最終ステップが含まれていない')
        self.assertEqual(len(rows), self.result['iteration'] // self.frequency + 1)


class TestModelCombinations(unittest.TestCase):
    """A-2: 大気・空力モデルの 4 通りの組み合わせがすべて動くこと。"""

    def test_all_four_combinations_run(self):
        for kind_atm in ('fileread', 'constant'):
            for kind_aero in ('fileread', 'constant'):
                with self.subTest(atmosphere=kind_atm, aerodynamic=kind_aero):
                    config = short_config()
                    config['atmosphere']['kind_atmosphere_model'] = kind_atm
                    config['satellite']['kind_aerodynamic_model'] = kind_aero

                    result = run_solver(config)
                    self.assertEqual(result['iteration'], int(TIME_MAX / TIMESTEP))

    def test_constant_model_uses_the_config_values(self):
        config = short_config()
        config['atmosphere']['kind_atmosphere_model'] = 'constant'
        config['atmosphere']['density'] = 1.0e-11
        config['atmosphere']['temperature'] = 250.0
        config['atmosphere']['knudsen'] = 12.5

        result = run_solver(config)

        self.assertEqual(set(result['trajectory']['density']), {1.0e-11})
        self.assertEqual(set(result['trajectory']['temperature']), {250.0})
        self.assertEqual(set(result['trajectory']['knudsen']), {12.5})


if __name__ == '__main__':
    unittest.main()
