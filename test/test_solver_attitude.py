#!/usr/bin/env python3
"""
6 自由度計算をソルバーごと走らせるテスト。

姿勢と並進が双方向に結合しているので、単体テストだけでは
「結合の向き」や「配列長」の取り違えが見つからない。ここでは
解析解（振動周期）と保存量（角運動量）、対称性（面内運動）で確かめる。
"""

import copy
import os
import unittest

import numpy as np

from context import ROOT_DIR, quiet

import atmosphere.atmosphere as atmosphere
import attitude.attitude as attitude
import force_term.force_term as force_term
import moment_term.moment_term as moment_term
import satellite.satellite as satellite
import solver.solver as solver
from orbital.orbital import orbital

import yaml

CONFIG_6DOF = os.path.join(ROOT_DIR, 'tutorial', 'work_reentry_6dof', 'config.yml')


def base_config(time_max, timestep):
    with open(CONFIG_6DOF) as f:
        config = yaml.safe_load(f)
    config['computational_setup']['time_elapsed_maximum'] = time_max
    config['time_integration']['timestep_constant'] = timestep
    config['post_process']['kml']['flag_output'] = False
    config['post_process']['tecplot']['flag_output'] = False
    return config


def constant_model_config(time_max, timestep, density=1.e-8, drag=1.0, stability=-0.5):
    """空力を定数モデルにして、振動周期が手で書けるようにした config。"""
    config = base_config(time_max, timestep)
    config['atmosphere']['kind_atmosphere_model'] = 'constant'
    config['atmosphere']['density'] = density
    config['atmosphere']['temperature'] = 1000.0
    config['atmosphere']['knudsen'] = 1.0
    config['satellite']['kind_aerodynamic_model'] = 'constant'
    config['satellite']['drag_coefficient'] = drag
    config['attitude']['static_stability_derivative'] = stability
    config['attitude']['flag_moment_damping'] = False
    config['attitude']['flag_moment_gravity_gradient'] = False
    return config


def run_solver(config):
    """初期化からソルバーまでを回し、辞書ごと返す。"""
    with quiet():
        atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
        aerodynamic_dict = satellite.initial_settings_satellite(config)

        orb = orbital()
        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            orb.initial_settings(config)
        attitude_dict = orb.initial_settings_attitude(config, coordinate_dict, velocity_dict)

        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            solver.solve_equation_motion(config, iteration, time_elapsed,
                                         coordinate_dict, velocity_dict, trajectory_dict,
                                         atmosphere_dict, aerodynamic_dict, attitude_dict)

    return iteration, coordinate_dict, velocity_dict, attitude_dict


def aerodynamic_angle_history(coordinate_dict, velocity_dict, attitude_dict, iteration):
    """各ステップの迎角・横滑り角（deg）を出力ルーチンと同じ手順で作り直す。"""
    alpha_list = []
    beta_list = []
    for n in range(0, iteration+1):
        matrix_be = attitude.quaternion_to_matrix(attitude_dict[orbital.KEY_ATTITUDE_QUATERNION][n])
        velocity_body = np.dot(matrix_be, velocity_dict['cartesian'][n])
        alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(velocity_body)
        alpha_list.append(alpha*orbital.rad2deg)
        beta_list.append(beta*orbital.rad2deg)
    return np.array(alpha_list), np.array(beta_list)


class TestAttitudeStateArrays(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.iteration, cls.coordinate, cls.velocity, cls.attitude = \
            run_solver(constant_model_config(10.0, 0.1))

    def test_attitude_arrays_have_the_same_length_as_the_coordinates(self):
        # 座標・速度と同じく iteration+1 個でなければ出力で添字がずれる
        self.assertEqual(len(self.attitude[orbital.KEY_ATTITUDE_QUATERNION]), self.iteration+1)
        self.assertEqual(len(self.attitude[orbital.KEY_ATTITUDE_OMEGA]), self.iteration+1)
        self.assertEqual(len(self.coordinate['cartesian']), self.iteration+1)

    def test_quaternion_stays_of_unit_length(self):
        for quaternion in self.attitude[orbital.KEY_ATTITUDE_QUATERNION]:
            self.assertAlmostEqual(float(np.linalg.norm(quaternion)), 1.0, places=12)

    def test_initial_attitude_matches_the_configuration(self):
        config = constant_model_config(10.0, 0.1)
        quaternion = self.attitude[orbital.KEY_ATTITUDE_QUATERNION][0]
        coordinate_polar = self.coordinate['polar'][0]
        yaw, pitch, roll = attitude.get_euler_angle(quaternion, coordinate_polar[2], coordinate_polar[1])

        expected = np.array(config['initial_settings']['attitude'])*orbital.deg2rad
        np.testing.assert_allclose([yaw, pitch, roll], expected, atol=1.e-10)


class TestTorqueFreeMotion(unittest.TestCase):
    """トルクを全て切った状態。角運動量が保存し、姿勢は一定速度で回る。"""

    @classmethod
    def setUpClass(cls):
        config = constant_model_config(20.0, 0.05, density=0.0, drag=0.0)
        config['attitude']['flag_moment_aerodynamic'] = False
        config['initial_settings']['angular_velocity'] = [5.0, 0.0, 0.0]
        cls.config = config
        cls.iteration, cls.coordinate, cls.velocity, cls.attitude = run_solver(config)

    def test_angular_momentum_magnitude_is_conserved(self):
        inertia = np.diag([self.config['attitude']['inertia_tensor']['Ixx'],
                           self.config['attitude']['inertia_tensor']['Iyy'],
                           self.config['attitude']['inertia_tensor']['Izz']])
        momentum = [np.linalg.norm(np.dot(inertia, omega))
                    for omega in self.attitude[orbital.KEY_ATTITUDE_OMEGA]]
        self.assertLess(abs(momentum[-1]/momentum[0] - 1.0), 1.e-12)

    def test_spin_about_a_principal_axis_stays_on_that_axis(self):
        # 慣性系基準の角速度には地球自転の分だけ横成分が残る。それは軸対称体の
        # 歳差として y-z 面内を回るので、軸成分と横成分の大きさは変わらない
        omega_last = self.attitude[orbital.KEY_ATTITUDE_OMEGA][-1]
        omega_first = self.attitude[orbital.KEY_ATTITUDE_OMEGA][0]

        self.assertAlmostEqual(omega_last[0], omega_first[0], places=12)
        self.assertAlmostEqual(float(np.linalg.norm(omega_last[1:3])),
                               float(np.linalg.norm(omega_first[1:3])), places=12)

    def test_the_body_turns_by_the_expected_angle(self):
        # ECEF に対する回転角 = |omega_rel| * t（地球自転の分だけ僅かにずれる）
        quaternion_first = self.attitude[orbital.KEY_ATTITUDE_QUATERNION][0]
        quaternion_last = self.attitude[orbital.KEY_ATTITUDE_QUATERNION][-1]

        matrix_relative = np.dot(attitude.quaternion_to_matrix(quaternion_last),
                                 attitude.quaternion_to_matrix(quaternion_first).T)
        cosine = 0.5*(np.trace(matrix_relative) - 1.0)
        angle = np.arccos(min(1.0, max(-1.0, cosine)))*orbital.rad2deg

        # 5 deg/s * 20 s = 100 deg（回転角なので 0-180 度に収まる範囲を選んである）
        # 地球自転（0.0042 deg/s）の分だけずれるので、その大きさを許容する
        self.assertAlmostEqual(angle, 100.0, delta=0.1)


class TestPitchOscillation(unittest.TestCase):
    """
    静安定微係数を与えた定数モデルでは、微小迎角のピッチ振動の周期が
      T = 2 pi sqrt( Iyy / (q S L |Cm_alpha|) )
    になる。復元モーメントの符号・大きさ・結合の向きがまとめて確かめられる。
    """

    DENSITY = 1.e-8
    STABILITY = -0.5
    TIME_MAX = 60.0
    TIMESTEP = 0.02

    @classmethod
    def setUpClass(cls):
        config = constant_model_config(cls.TIME_MAX, cls.TIMESTEP,
                                       density=cls.DENSITY, stability=cls.STABILITY)
        # 微小振動にするため初期迎角は 2 度だけ与える
        config['initial_settings']['attitude'] = [90.0, 2.0, 0.0]
        cls.config = config
        cls.iteration, cls.coordinate, cls.velocity, cls.attitude = run_solver(config)
        cls.alpha, cls.beta = aerodynamic_angle_history(cls.coordinate, cls.velocity,
                                                        cls.attitude, cls.iteration)

    def period_analytic(self):
        velocity = np.linalg.norm(self.velocity['cartesian'][0])
        dynamic_pressure = 0.5*self.DENSITY*velocity**2
        stiffness = dynamic_pressure \
                  * self.config['satellite']['characteristic_area'] \
                  * self.config['satellite']['characteristic_length'] \
                  * abs(self.STABILITY)
        return 2.0*np.pi*np.sqrt(self.config['attitude']['inertia_tensor']['Iyy']/stiffness)

    def test_period_matches_the_analytic_value(self):
        # 迎角の符号が変わる時刻を線形補間で拾い、その間隔から周期を求める
        crossing = []
        for n in range(1, len(self.alpha)):
            if self.alpha[n-1]*self.alpha[n] < 0.0:
                fact = self.alpha[n-1]/(self.alpha[n-1] - self.alpha[n])
                crossing.append((n - 1 + fact)*self.TIMESTEP)

        self.assertGreaterEqual(len(crossing), 3)
        period = 2.0*(crossing[-1] - crossing[0])/float(len(crossing) - 1)

        self.assertAlmostEqual(period/self.period_analytic(), 1.0, delta=0.01)

    def test_oscillation_is_bounded_by_the_initial_amplitude(self):
        # 減衰を切ってあるので振幅は増えも減りもしない
        self.assertLess(np.max(np.abs(self.alpha)), 2.0*1.02)
        self.assertGreater(np.max(np.abs(self.alpha)), 2.0*0.98)

    def test_motion_stays_in_the_plane(self):
        # 面内に置いた初期姿勢なら横滑りは立たない
        self.assertLess(np.max(np.abs(self.beta)), 0.05)


class TestDampingAndTable(unittest.TestCase):

    def test_damping_reduces_the_amplitude(self):
        config = constant_model_config(60.0, 0.02)
        config['initial_settings']['attitude'] = [90.0, 5.0, 0.0]
        config['attitude']['flag_moment_damping'] = True
        config['attitude']['damping_coefficient'] = {'Clp': 0.0, 'Cmq': -5.0, 'Cnr': -5.0}
        iteration, coordinate, velocity, attitude_dict = run_solver(config)
        alpha, beta = aerodynamic_angle_history(coordinate, velocity, attitude_dict, iteration)

        half = len(alpha)//2
        self.assertLess(np.max(np.abs(alpha[half:])), np.max(np.abs(alpha[:half])))

    def test_table_model_restores_the_attitude(self):
        # 迎角依存の係数表（球円錐）でも迎角 0 のまわりに戻る
        config = base_config(60.0, 0.05)
        config['initial_settings']['attitude'] = [90.0, 10.0, 0.0]
        iteration, coordinate, velocity, attitude_dict = run_solver(config)
        alpha, beta = aerodynamic_angle_history(coordinate, velocity, attitude_dict, iteration)

        self.assertLess(np.min(alpha), 0.0)          # 0 度を通り越して戻ってくる
        self.assertLess(np.max(np.abs(alpha)), 10.5)  # 発散しない
        self.assertLess(np.max(np.abs(beta)), 0.05)


class TestCouplingWithTheTranslation(unittest.TestCase):
    """迎角 0 では 6 自由度の空力が 3 自由度の抗力と一致しなければならない。"""

    def test_aerodynamic_force_matches_the_three_degree_of_freedom_drag(self):
        config = constant_model_config(1.0, 0.1, density=1.e-6, drag=1.2)

        coordinate = np.array([4.0e6, 3.0e6, 3.5e6])
        velocity = np.array([1000.0, -6000.0, 4000.0])

        # 機体 x 軸を速度方向に合わせたクォータニオンを作る
        axis_x = velocity/np.linalg.norm(velocity)
        axis_y = np.cross([0.0, 0.0, 1.0], axis_x)
        axis_y = axis_y/np.linalg.norm(axis_y)
        axis_z = np.cross(axis_x, axis_y)
        quaternion = attitude.matrix_to_quaternion(np.array([axis_x, axis_y, axis_z]))

        with quiet():
            aerodynamic_dict = satellite.initial_settings_satellite(config)
            moment, property_dict = moment_term.moment_initialsettings(config)

        force_aerodynamic, moment_total, omega_relative = solver.get_aerodynamic_state(
            config, property_dict, aerodynamic_dict, 'constant',
            coordinate, velocity, quaternion, np.zeros(3),
            config['satellite']['mass'], config['satellite']['characteristic_area'],
            config['satellite']['characteristic_length'],
            config['atmosphere']['density'], 1.0, 1.0,
            config['satellite']['drag_coefficient'], config['attitude']['static_stability_derivative'],
            config['planet']['rotation_rate'], moment)

        force = force_term.force_initialsettings(config)
        force = force_term.force_routine(config, coordinate, velocity,
                                         config['satellite']['mass'],
                                         config['satellite']['characteristic_area'],
                                         config['satellite']['drag_coefficient'], 1.0,
                                         config['atmosphere']['density'], force)

        np.testing.assert_allclose(force_aerodynamic, force[4, :], rtol=1.e-12)
        # 迎角 0 なので復元モーメントも立たない
        np.testing.assert_allclose(moment[1, :], np.zeros(3), atol=1.e-12)


if __name__ == '__main__':
    unittest.main()
