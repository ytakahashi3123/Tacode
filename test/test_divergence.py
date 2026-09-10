#!/usr/bin/env python3
"""
発散・非物理状態の検知のテスト（レビュー A-9）。

ソルバーの打ち切り条件は「高度が負」だけだった。そのため時間刻みが粗すぎる計算は
高度 6.6e6 km・速度 1.4e6 m/s といった状態をそのまま Tecplot・restart・KML に書き、
終了コード 0 を返していた。数値が並んでいるので一見それらしく、後処理ツールもそのまま読む。

ここで確かめるのは 3 点:

  1. 状態が非物理になったら **止まる**（終了コード 1。出力を書き切らない）
  2. 健全な計算には一切触らない。既定の上限（地心距離の 10 倍・脱出速度の 10 倍）は
     チュートリアルのどのケースでも引っかからない
  3. 上限と on/off が config から効くこと
"""

import os
import unittest

import numpy as np

from context import ROOT_DIR, load_config, quiet

import solver.solver as solver

TIME_MAX = 30.0
TIMESTEP = 1.0

# 発散させる時間刻み（レビューの再現条件）
TIMESTEP_COARSE = 400.0


def config_short(timestep=TIMESTEP, time_maximum=TIME_MAX):
    config = load_config()
    config['computational_setup']['time_elapsed_maximum'] = time_maximum
    config['time_integration']['timestep_constant'] = timestep
    return config


def run_solver(config):
    """初期化からソルバーまでを回す（test_solver_invariants と同じ手順）。"""
    import atmosphere.atmosphere as atmosphere
    import satellite.satellite as satellite
    from orbital.orbital import orbital

    atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
    aerodynamic_dict = satellite.initial_settings_satellite(config)

    orb = orbital()
    iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
        orb.initial_settings(config)

    return solver.solve_equation_motion(config, iteration, time_elapsed,
                                        coordinate_dict, velocity_dict, trajectory_dict,
                                        atmosphere_dict, aerodynamic_dict)


class TestTheLimitsComeFromTheInitialState(unittest.TestCase):
    """上限は初期状態から作る。絶対値で決めると軌道の高さで意味が変わる。"""

    def setUp(self):
        self.config = load_config()
        self.coordinate = np.array([7.0e6, 0.0, 0.0])

    def test_the_default_is_ten_times_the_initial_radius(self):
        flag, radius_maximum, velocity_maximum = solver.set_divergence_limit(self.config,
                                                                            self.coordinate)
        self.assertTrue(flag)
        self.assertAlmostEqual(radius_maximum, 10.0*np.linalg.norm(self.coordinate), places=6)

    def test_the_default_is_ten_times_the_escape_velocity(self):
        planet = self.config['planet']
        expected = np.sqrt(2.0*planet['gravitational_constant']*planet['mass']
                           / np.linalg.norm(self.coordinate))
        flag, radius_maximum, velocity_maximum = solver.set_divergence_limit(self.config,
                                                                            self.coordinate)
        self.assertAlmostEqual(velocity_maximum, 10.0*expected, places=6)

    def test_the_factors_are_read_from_the_configuration(self):
        self.config['computational_setup']['factor_radius_maximum'] = 2.0
        self.config['computational_setup']['factor_velocity_maximum'] = 3.0
        flag, radius_maximum, velocity_maximum = solver.set_divergence_limit(self.config,
                                                                            self.coordinate)
        self.assertAlmostEqual(radius_maximum, 2.0*np.linalg.norm(self.coordinate), places=6)

        planet = self.config['planet']
        escape = np.sqrt(2.0*planet['gravitational_constant']*planet['mass']
                         / np.linalg.norm(self.coordinate))
        self.assertAlmostEqual(velocity_maximum, 3.0*escape, places=6)

    def test_the_check_can_be_switched_off(self):
        self.config['computational_setup']['flag_check_divergence'] = False
        flag, radius_maximum, velocity_maximum = solver.set_divergence_limit(self.config,
                                                                            self.coordinate)
        self.assertFalse(flag)


class TestTheStateIsChecked(unittest.TestCase):
    """check_divergence そのもの。止まるときは終了コード 1 で止まる。"""

    RADIUS_MAXIMUM = 7.0e7
    VELOCITY_MAXIMUM = 1.0e5

    def check(self, coordinate, velocity, quaternion=None, angular_velocity=None):
        with quiet():
            return solver.check_divergence(np.array(coordinate), np.array(velocity),
                                           self.RADIUS_MAXIMUM, self.VELOCITY_MAXIMUM, 1.0,
                                           quaternion, angular_velocity)

    def expect_a_stop(self, coordinate, velocity, quaternion=None, angular_velocity=None):
        with self.assertRaises(SystemExit) as raised:
            self.check(coordinate, velocity, quaternion, angular_velocity)
        self.assertEqual(raised.exception.code, 1)

    def test_a_sound_state_passes(self):
        self.assertIsNone(self.check([7.0e6, 0.0, 0.0], [0.0, 7.5e3, 0.0]))

    def test_a_state_just_inside_the_limits_passes(self):
        self.assertIsNone(self.check([self.RADIUS_MAXIMUM*(1.0 - 1.0e-12), 0.0, 0.0],
                                     [self.VELOCITY_MAXIMUM*(1.0 - 1.0e-12), 0.0, 0.0]))

    def test_a_distance_beyond_the_limit_stops_the_run(self):
        self.expect_a_stop([self.RADIUS_MAXIMUM*1.1, 0.0, 0.0], [0.0, 7.5e3, 0.0])

    def test_a_velocity_beyond_the_limit_stops_the_run(self):
        self.expect_a_stop([7.0e6, 0.0, 0.0], [0.0, self.VELOCITY_MAXIMUM*1.1, 0.0])

    def test_a_position_which_is_not_a_number_stops_the_run(self):
        self.expect_a_stop([np.nan, 0.0, 0.0], [0.0, 7.5e3, 0.0])

    def test_a_velocity_which_is_infinite_stops_the_run(self):
        self.expect_a_stop([7.0e6, 0.0, 0.0], [np.inf, 0.0, 0.0])

    def test_an_attitude_which_is_not_a_number_stops_the_run(self):
        # 姿勢を解いているときは、クォータニオンと角速度も見る
        self.expect_a_stop([7.0e6, 0.0, 0.0], [0.0, 7.5e3, 0.0],
                           quaternion=[1.0, 0.0, np.nan, 0.0],
                           angular_velocity=[0.0, 0.0, 0.0])

    def test_an_angular_velocity_which_is_infinite_stops_the_run(self):
        self.expect_a_stop([7.0e6, 0.0, 0.0], [0.0, 7.5e3, 0.0],
                           quaternion=[1.0, 0.0, 0.0, 0.0],
                           angular_velocity=[np.inf, 0.0, 0.0])

    def test_a_sound_attitude_passes(self):
        self.assertIsNone(self.check([7.0e6, 0.0, 0.0], [0.0, 7.5e3, 0.0],
                                     quaternion=[1.0, 0.0, 0.0, 0.0],
                                     angular_velocity=[0.0, 0.01, 0.0]))


class TestTheSolverStops(unittest.TestCase):
    """時間ループの中で効くこと。"""

    def test_a_coarse_timestep_stops_the_run(self):
        # レビューの再現条件。それまでは高度 6.6e6 km を書いて終了コード 0 だった
        with self.assertRaises(SystemExit) as raised:
            with quiet():
                run_solver(config_short(timestep=TIMESTEP_COARSE, time_maximum=5000.0))
        self.assertEqual(raised.exception.code, 1)

    def test_the_same_run_goes_through_with_the_check_off(self):
        # 逃げ道が要る（上限のほうが間違っている場合のため）
        config = config_short(timestep=TIMESTEP_COARSE, time_maximum=5000.0)
        config['computational_setup']['flag_check_divergence'] = False
        with quiet():
            iteration = run_solver(config)[0]
        self.assertGreater(iteration, 0)

    def test_a_sound_run_is_untouched(self):
        with quiet():
            iteration, time_elapsed = run_solver(config_short())[0:2]
        self.assertEqual(iteration, int(TIME_MAX/TIMESTEP))
        self.assertAlmostEqual(time_elapsed, TIME_MAX, places=9)

    def test_the_message_names_the_reason_and_the_way_out(self):
        with quiet() as printed:
            with self.assertRaises(SystemExit):
                run_solver(config_short(timestep=TIMESTEP_COARSE, time_maximum=5000.0))
        text = printed.getvalue()
        self.assertIn('The solution has diverged', text)
        self.assertIn('timestep_constant', text)
        self.assertIn('flag_check_divergence', text)


if __name__ == '__main__':
    unittest.main()
