#!/usr/bin/env python3
"""
モーメント項と Euler の運動方程式のテスト。

力側（test_force_term.py）と同じ方針で、解析解と保存量に突き合わせる。
"""

import copy
import unittest

import numpy as np

from context import load_config, quiet

import attitude.attitude as attitude
import moment_term.moment_term as moment_term


def attitude_config(inertia=None, center_of_gravity=None, damping=None, flags=None):
    """姿勢セクションを付けた config を返す。"""
    config = copy.deepcopy(load_config())
    section = {'inertia_tensor': inertia if inertia is not None
                                 else {'Ixx': 0.75, 'Iyy': 0.50, 'Izz': 0.50}}
    if center_of_gravity is not None:
        section['center_of_gravity'] = center_of_gravity
    if damping is not None:
        section['damping_coefficient'] = damping
    if flags is not None:
        section.update(flags)
    config['attitude'] = section
    return config


def integrate_torque_free(property_dict, omega, time_max, timestep):
    """トルクが無い状態で Euler の運動方程式を RK4 で積分する。"""
    moment = np.zeros(3)
    num_step = int(round(time_max/timestep))
    history = [np.array(omega)]
    for step in range(0, num_step):
        k1 = moment_term.solve_angular_acceleration(property_dict, omega, moment)
        k2 = moment_term.solve_angular_acceleration(property_dict, omega + 0.5*timestep*k1, moment)
        k3 = moment_term.solve_angular_acceleration(property_dict, omega + 0.5*timestep*k2, moment)
        k4 = moment_term.solve_angular_acceleration(property_dict, omega + timestep*k3, moment)
        omega = omega + timestep/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
        history.append(np.array(omega))
    return np.array(history)


class TestInertiaSettings(unittest.TestCase):

    def test_tensor_is_assembled_with_the_products_of_inertia(self):
        inertia = {'Ixx': 1.0, 'Iyy': 2.0, 'Izz': 3.0, 'Ixy': 0.1, 'Iyz': 0.2, 'Izx': 0.3}
        with quiet():
            moment, property_dict = moment_term.moment_initialsettings(attitude_config(inertia))

        expected = np.array([[1.0, 0.1, 0.3],
                             [0.1, 2.0, 0.2],
                             [0.3, 0.2, 3.0]])
        np.testing.assert_allclose(property_dict[moment_term.KEY_INERTIA], expected)
        np.testing.assert_allclose(
            np.dot(property_dict[moment_term.KEY_INERTIA], property_dict[moment_term.KEY_INERTIA_INV]),
            np.eye(3), atol=1.e-12)
        self.assertEqual(moment.shape, (4, 3))

    def test_singular_tensor_is_rejected(self):
        with self.assertRaises(SystemExit):
            with quiet():
                moment_term.moment_initialsettings(attitude_config({'Ixx': 0.0, 'Iyy': 0.0, 'Izz': 0.0}))

    def test_missing_tensor_is_rejected(self):
        config = copy.deepcopy(load_config())
        config['attitude'] = {'flag_attitude': True}
        with self.assertRaises(SystemExit):
            with quiet():
                moment_term.moment_initialsettings(config)


class TestEulerEquation(unittest.TestCase):

    def setUp(self):
        with quiet():
            self.moment, self.property_asymmetric = moment_term.moment_initialsettings(
                attitude_config({'Ixx': 1.0, 'Iyy': 2.0, 'Izz': 3.0}))
            _, self.property_axisymmetric = moment_term.moment_initialsettings(
                attitude_config({'Ixx': 0.75, 'Iyy': 0.50, 'Izz': 0.50}))

    def test_angular_momentum_and_energy_are_conserved(self):
        inertia = self.property_asymmetric[moment_term.KEY_INERTIA]
        history = integrate_torque_free(self.property_asymmetric,
                                        np.array([0.3, 0.2, -0.5]), 5.0, 1.e-3)

        momentum = [np.linalg.norm(np.dot(inertia, omega)) for omega in history]
        energy = [0.5*np.dot(omega, np.dot(inertia, omega)) for omega in history]

        self.assertLess(abs(momentum[-1]/momentum[0] - 1.0), 1.e-10)
        self.assertLess(abs(energy[-1]/energy[0] - 1.0), 1.e-10)

    def test_rotation_about_a_principal_axis_is_steady(self):
        for axis in range(0, 3):
            omega = np.zeros(3)
            omega[axis] = 0.7
            acceleration = moment_term.solve_angular_acceleration(
                self.property_asymmetric, omega, np.zeros(3))
            np.testing.assert_allclose(acceleration, np.zeros(3), atol=1.e-15)

    def test_axisymmetric_precession_matches_the_analytic_rate(self):
        #
        # 軸対称体（Iyy = Izz）のトルクフリー運動では、機体軸から見た横方向の
        # 角速度ベクトルが一定の速さ
        #   lambda = p (Ixx - Iyy)/Iyy
        # で回る。
        #
        inertia = self.property_axisymmetric[moment_term.KEY_INERTIA]
        spin = 2.0
        rate_analytic = spin*(inertia[0, 0] - inertia[1, 1])/inertia[1, 1]

        timestep = 2.e-4
        time_max = 2.0
        history = integrate_torque_free(self.property_axisymmetric,
                                        np.array([spin, 0.1, 0.0]), time_max, timestep)

        # 横成分の位相の進み方から回転速度を求める
        angle = np.unwrap(np.arctan2(history[:, 2], history[:, 1]))
        rate_computed = (angle[-1] - angle[0])/time_max

        self.assertAlmostEqual(rate_computed, rate_analytic, places=6)

        # スピン成分と横成分の大きさは変わらない
        self.assertAlmostEqual(history[-1, 0], spin, places=10)
        self.assertAlmostEqual(np.linalg.norm(history[-1, 1:3]), 0.1, places=10)


class TestGravityGradientMoment(unittest.TestCase):

    def setUp(self):
        self.config = load_config()
        self.coordinate = np.array([6.7e6, 0.0, 0.0])

    def moment_of(self, inertia, matrix_be):
        with quiet():
            moment, property_dict = moment_term.moment_initialsettings(
                attitude_config(inertia, flags={'flag_moment_aerodynamic': False,
                                                'flag_moment_damping': False}))
        moment = moment_term.moment_routine(self.config, property_dict, self.coordinate, matrix_be,
                                            np.zeros(3), np.zeros(3), np.zeros(3),
                                            0.0, 1.0, 1.0, 0.0, moment)
        return moment[3, :].copy(), property_dict

    def test_vanishes_for_an_isotropic_body(self):
        moment, property_dict = self.moment_of({'Ixx': 1.0, 'Iyy': 1.0, 'Izz': 1.0},
                                               attitude.euler_to_matrix(0.3, 0.4, 0.5))
        np.testing.assert_allclose(moment, np.zeros(3), atol=1.e-12)

    def test_vanishes_when_a_principal_axis_points_at_the_planet(self):
        # 機体軸 x が地心方向を向いていれば u x I u = 0
        moment, property_dict = self.moment_of({'Ixx': 1.0, 'Iyy': 2.0, 'Izz': 3.0}, np.eye(3))
        np.testing.assert_allclose(moment, np.zeros(3), atol=1.e-12)

    def test_matches_the_analytic_expression(self):
        inertia = {'Ixx': 1.0, 'Iyy': 2.0, 'Izz': 3.0}
        matrix_be = attitude.euler_to_matrix(0.6, -0.35, 1.1)
        moment, property_dict = self.moment_of(inertia, matrix_be)

        gravitational_parameter = self.config['planet']['gravitational_constant'] \
                                * self.config['planet']['mass']
        radius = np.linalg.norm(self.coordinate)
        unit_body = np.dot(matrix_be, self.coordinate/radius)
        expected = 3.0*gravitational_parameter/radius**3 \
                 * np.cross(unit_body, np.dot(property_dict[moment_term.KEY_INERTIA], unit_body))

        np.testing.assert_allclose(moment, expected, rtol=1.e-12)
        # 大きさの桁の確認（低軌道の重力傾斜トルクは 1e-6 N m 程度）
        self.assertLess(np.linalg.norm(moment), 1.e-4)

    def test_can_be_switched_off(self):
        with quiet():
            moment, property_dict = moment_term.moment_initialsettings(
                attitude_config({'Ixx': 1.0, 'Iyy': 2.0, 'Izz': 3.0},
                                flags={'flag_moment_gravity_gradient': False}))
        moment = moment_term.moment_routine(self.config, property_dict, self.coordinate,
                                            attitude.euler_to_matrix(0.3, 0.4, 0.5),
                                            np.zeros(3), np.zeros(3), np.zeros(3),
                                            0.0, 1.0, 1.0, 0.0, moment)
        np.testing.assert_allclose(moment[3, :], np.zeros(3), atol=0.0)


class TestAerodynamicMoment(unittest.TestCase):

    def setUp(self):
        self.config = load_config()
        self.coordinate = np.array([6.7e6, 0.0, 0.0])

    def build(self, center_of_gravity=None, damping=None):
        with quiet():
            return moment_term.moment_initialsettings(
                attitude_config({'Ixx': 0.75, 'Iyy': 0.50, 'Izz': 0.50},
                                center_of_gravity=center_of_gravity, damping=damping,
                                flags={'flag_moment_gravity_gradient': False}))

    def test_static_moment_follows_the_coefficient(self):
        moment, property_dict = self.build()
        dynamic_pressure, area, length = 100.0, 0.8, 0.5
        coefficient_moment = np.array([0.0, -0.03, 0.0])

        moment = moment_term.moment_routine(self.config, property_dict, self.coordinate, np.eye(3),
                                            np.zeros(3), coefficient_moment, np.zeros(3),
                                            dynamic_pressure, area, length, 1000.0, moment)
        np.testing.assert_allclose(moment[1, :],
                                   dynamic_pressure*area*length*coefficient_moment, rtol=1.e-12)

    def test_centre_of_gravity_offset_adds_the_expected_moment(self):
        # 係数の基準点より重心が後ろ（-x 側）にあると、軸力が重心まわりに
        # モーメントを作らないことと、法線力がモーメントを作ることを確かめる
        offset = np.array([-0.2, 0.0, 0.0])
        moment, property_dict = self.build(center_of_gravity=list(offset))

        force_axial = np.array([-50.0, 0.0, 0.0])
        moment = moment_term.moment_routine(self.config, property_dict, self.coordinate, np.eye(3),
                                            np.zeros(3), np.zeros(3), force_axial,
                                            0.0, 1.0, 1.0, 1000.0, moment)
        np.testing.assert_allclose(moment[1, :], np.zeros(3), atol=1.e-12)

        force_normal = np.array([0.0, 0.0, -50.0])
        moment = moment_term.moment_routine(self.config, property_dict, self.coordinate, np.eye(3),
                                            np.zeros(3), np.zeros(3), force_normal,
                                            0.0, 1.0, 1.0, 1000.0, moment)
        np.testing.assert_allclose(moment[1, :], np.cross(-offset, force_normal), rtol=1.e-12)

    def test_damping_removes_rotational_energy(self):
        damping = {'Clp': -0.2, 'Cmq': -0.5, 'Cnr': -0.5}
        moment, property_dict = self.build(damping=damping)

        omega_relative = np.array([0.05, -0.3, 0.2])
        moment = moment_term.moment_routine(self.config, property_dict, self.coordinate, np.eye(3),
                                            omega_relative, np.zeros(3), np.zeros(3),
                                            100.0, 0.8, 0.5, 1000.0, moment)

        # 減衰モーメントは角速度と逆を向く（回転エネルギーを減らす）
        self.assertLess(float(np.dot(moment[2, :], omega_relative)), 0.0)
        for index in range(0, 3):
            self.assertLess(moment[2, index]*omega_relative[index], 0.0)

    def test_damping_vanishes_without_velocity(self):
        moment, property_dict = self.build(damping={'Clp': -0.2, 'Cmq': -0.5, 'Cnr': -0.5})
        moment = moment_term.moment_routine(self.config, property_dict, self.coordinate, np.eye(3),
                                            np.array([0.1, 0.1, 0.1]), np.zeros(3), np.zeros(3),
                                            0.0, 0.8, 0.5, 0.0, moment)
        np.testing.assert_allclose(moment[2, :], np.zeros(3), atol=0.0)

    def test_total_is_the_sum_of_the_terms(self):
        moment, property_dict = self.build(center_of_gravity=[-0.1, 0.0, 0.05],
                                           damping={'Clp': -0.2, 'Cmq': -0.5, 'Cnr': -0.5})
        moment = moment_term.moment_routine(self.config, property_dict, self.coordinate,
                                            attitude.euler_to_matrix(0.2, 0.3, -0.4),
                                            np.array([0.05, -0.3, 0.2]),
                                            np.array([0.0, -0.02, 0.01]),
                                            np.array([-40.0, 2.0, -5.0]),
                                            100.0, 0.8, 0.5, 1000.0, moment)
        np.testing.assert_allclose(moment[0, :], moment[1, :] + moment[2, :] + moment[3, :], rtol=1.e-14)


if __name__ == '__main__':
    unittest.main()
