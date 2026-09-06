#!/usr/bin/env python3
"""
姿勢（6 自由度）のキネマティクスのテスト。

座標系の取り違えと符号の反転は 6 自由度計算で最も入りやすいバグなので、
変換は往復と既存の 3 自由度側の規約との突き合わせで確かめる。
"""

import unittest

import numpy as np

from context import load_config, quiet

import attitude.attitude as attitude
import coordinate_system.coordinate_system as coordinate_system


def rotation_matrix_passive(axis, angle):
    """軸 axis まわりに角 angle だけ座標系を回す受動回転行列（解析解）。"""
    unit = np.array(axis, dtype=float)
    unit = unit/np.linalg.norm(unit)
    cross = np.array([[0.0, -unit[2], unit[1]],
                      [unit[2], 0.0, -unit[0]],
                      [-unit[1], unit[0], 0.0]])
    return np.eye(3) - np.sin(angle)*cross + (1.0 - np.cos(angle))*np.dot(cross, cross)


def integrate_quaternion(quaternion, omega, time_max, timestep):
    """一定角速度でクォータニオンを RK4 で積分する。"""
    num_step = int(round(time_max/timestep))
    for step in range(0, num_step):
        k1 = attitude.quaternion_derivative(quaternion, omega)
        k2 = attitude.quaternion_derivative(quaternion + 0.5*timestep*k1, omega)
        k3 = attitude.quaternion_derivative(quaternion + 0.5*timestep*k2, omega)
        k4 = attitude.quaternion_derivative(quaternion + timestep*k3, omega)
        quaternion = attitude.quaternion_normalize(
            quaternion + timestep/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4))
    return quaternion


class TestRotationConversion(unittest.TestCase):

    EULER_CASES = [(0.1, 0.2, 0.3), (1.0, -1.2, 2.0), (-2.0, 0.0, 0.5), (3.0, 1.4, -3.0)]

    def test_euler_matrix_round_trip(self):
        for yaw, pitch, roll in self.EULER_CASES:
            matrix = attitude.euler_to_matrix(yaw, pitch, roll)
            result = attitude.matrix_to_euler(matrix)
            np.testing.assert_allclose(result, (yaw, pitch, roll), atol=1.e-12)

    def test_euler_matrix_is_orthonormal(self):
        for yaw, pitch, roll in self.EULER_CASES:
            matrix = attitude.euler_to_matrix(yaw, pitch, roll)
            np.testing.assert_allclose(np.dot(matrix, matrix.T), np.eye(3), atol=1.e-12)
            self.assertAlmostEqual(np.linalg.det(matrix), 1.0, places=12)

    def test_quaternion_matrix_round_trip(self):
        for yaw, pitch, roll in self.EULER_CASES:
            matrix = attitude.euler_to_matrix(yaw, pitch, roll)
            quaternion = attitude.matrix_to_quaternion(matrix)
            np.testing.assert_allclose(attitude.quaternion_to_matrix(quaternion), matrix, atol=1.e-12)
            self.assertAlmostEqual(np.linalg.norm(quaternion), 1.0, places=12)

    def test_quaternion_matches_analytic_rotation(self):
        # 各軸まわりの単一回転をクォータニオンと解析解で比べる
        for axis in ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 2.0, -0.5]):
            for angle in (0.3, 1.7, -2.4):
                matrix = rotation_matrix_passive(axis, angle)
                quaternion = attitude.matrix_to_quaternion(matrix)
                np.testing.assert_allclose(attitude.quaternion_to_matrix(quaternion), matrix, atol=1.e-12)

    def test_gimbal_lock_is_handled(self):
        # ピッチ 90 度ではロールとヨーが縮退する。落ちずに姿勢そのものは再現できること
        matrix = attitude.euler_to_matrix(0.4, 0.5*np.pi, 0.0)
        yaw, pitch, roll = attitude.matrix_to_euler(matrix)
        self.assertAlmostEqual(pitch, 0.5*np.pi, places=9)
        np.testing.assert_allclose(attitude.euler_to_matrix(yaw, pitch, roll), matrix, atol=1.e-7)


class TestLocalHorizonFrame(unittest.TestCase):
    """
    ローカル水平系が既存の 3 自由度側の規約（初期速度の [東, 北, 上]）と
    同じものを指していることを確かめる。
    """

    def setUp(self):
        self.config = load_config()

    def test_ned_matches_the_existing_east_north_up(self):
        vector = np.array([1.0, 2.0, 3.0])
        for longitude, latitude in [(0.0, 0.0), (0.7, -0.4), (-2.5, 1.2), (np.pi, 0.0)]:
            east_north_up = np.array(
                coordinate_system.convert_carteasian_polar(self.config, vector, longitude, latitude))
            north_east_down = np.dot(attitude.matrix_ecef_to_ned(longitude, latitude), vector)
            np.testing.assert_allclose(
                north_east_down,
                [east_north_up[1], east_north_up[0], -east_north_up[2]], atol=1.e-12)

    def test_euler_angle_round_trip_through_the_local_frame(self):
        for longitude, latitude in [(0.0, 0.0), (0.7, -0.4), (-2.5, 1.2)]:
            for yaw, pitch, roll in [(0.3, 0.2, -0.1), (1.6, -0.9, 2.2)]:
                quaternion = attitude.get_quaternion_from_euler(yaw, pitch, roll, longitude, latitude)
                result = attitude.get_euler_angle(quaternion, longitude, latitude)
                np.testing.assert_allclose(result, (yaw, pitch, roll), atol=1.e-10)

    def test_zero_attitude_points_north_and_level(self):
        # ヨー・ピッチ・ロールがすべて 0 なら機体の x 軸は北、z 軸は地心向き
        longitude, latitude = 0.6, 0.3
        quaternion = attitude.get_quaternion_from_euler(0.0, 0.0, 0.0, longitude, latitude)
        matrix_be = attitude.quaternion_to_matrix(quaternion)

        vector_north = np.array([-np.sin(latitude)*np.cos(longitude),
                                 -np.sin(latitude)*np.sin(longitude),
                                  np.cos(latitude)])
        vector_down = np.array([-np.cos(latitude)*np.cos(longitude),
                                -np.cos(latitude)*np.sin(longitude),
                                -np.sin(latitude)])

        np.testing.assert_allclose(np.dot(matrix_be, vector_north), [1.0, 0.0, 0.0], atol=1.e-12)
        np.testing.assert_allclose(np.dot(matrix_be, vector_down), [0.0, 0.0, 1.0], atol=1.e-12)


class TestQuaternionKinematics(unittest.TestCase):

    def test_constant_rate_reproduces_the_analytic_rotation(self):
        omega = np.array([0.3, -0.2, 0.5])
        time_max = 1.0
        quaternion = integrate_quaternion(np.array([1.0, 0.0, 0.0, 0.0]), omega, time_max, 1.e-4)

        matrix_analytic = rotation_matrix_passive(omega, np.linalg.norm(omega)*time_max)
        np.testing.assert_allclose(attitude.quaternion_to_matrix(quaternion), matrix_analytic, atol=1.e-9)

    def test_norm_is_preserved(self):
        quaternion = integrate_quaternion(np.array([1.0, 0.0, 0.0, 0.0]),
                                          np.array([0.1, 0.4, -0.7]), 5.0, 1.e-3)
        self.assertAlmostEqual(np.linalg.norm(quaternion), 1.0, places=12)

    def test_derivative_is_orthogonal_to_the_quaternion(self):
        # 単位長を保つ条件 q . dq/dt = 0
        quaternion = attitude.quaternion_normalize(np.array([0.3, -0.5, 0.7, 0.2]))
        derivative = attitude.quaternion_derivative(quaternion, np.array([0.2, -0.3, 0.9]))
        self.assertAlmostEqual(float(np.dot(quaternion, derivative)), 0.0, places=14)

    def test_earth_rate_is_removed_for_the_relative_rate(self):
        rotation_rate = 7.292115e-5
        quaternion = attitude.quaternion_normalize(np.array([0.3, -0.5, 0.7, 0.2]))
        matrix_be = attitude.quaternion_to_matrix(quaternion)

        omega_inertial = attitude.get_earth_rate_body(rotation_rate, matrix_be)
        omega_relative = attitude.get_omega_relative(omega_inertial, rotation_rate, matrix_be)
        np.testing.assert_allclose(omega_relative, np.zeros(3), atol=1.e-18)

        # 地球と一緒に回っているだけなら ECEF から見た姿勢は動かない
        np.testing.assert_allclose(attitude.quaternion_derivative(quaternion, omega_relative),
                                   np.zeros(4), atol=1.e-18)


class TestAerodynamicAngle(unittest.TestCase):

    def test_angles_of_a_known_velocity(self):
        # 迎角 20 度、横滑り 0（機体軸 x-z 面内）
        angle = 20.0*np.pi/180.0
        velocity = 100.0*np.array([np.cos(angle), 0.0, np.sin(angle)])
        alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(velocity)
        self.assertAlmostEqual(alpha, angle, places=12)
        self.assertAlmostEqual(beta, 0.0, places=12)
        self.assertAlmostEqual(alpha_total, angle, places=12)
        self.assertAlmostEqual(phi_aero, 0.0, places=12)

    def test_sideslip_only(self):
        angle = 15.0*np.pi/180.0
        velocity = 100.0*np.array([np.cos(angle), np.sin(angle), 0.0])
        alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(velocity)
        self.assertAlmostEqual(alpha, 0.0, places=12)
        self.assertAlmostEqual(beta, angle, places=12)
        self.assertAlmostEqual(alpha_total, angle, places=12)
        self.assertAlmostEqual(phi_aero, 0.5*np.pi, places=12)

    def test_total_angle_is_never_negative(self):
        for velocity in ([1.0, 0.0, -1.0], [-1.0, 0.5, 0.3], [0.0, 0.0, -2.0]):
            alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(np.array(velocity))
            self.assertGreaterEqual(alpha_total, 0.0)
            self.assertLessEqual(alpha_total, np.pi)

    def test_zero_velocity_does_not_divide_by_zero(self):
        result = attitude.get_aerodynamic_angle(np.zeros(3))
        np.testing.assert_allclose(result, np.zeros(4))

    def test_roll_matrix_maps_the_table_plane_onto_the_cross_flow(self):
        # 係数表は横流れが +z 側にある面で作られている。
        # 実際の横流れ方向へ回す行列であることを確かめる
        for velocity in ([1.0, 0.4, 0.3], [1.0, -0.6, 0.2], [0.5, 0.0, -0.8]):
            velocity = np.array(velocity)
            alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(velocity)
            matrix_roll = attitude.matrix_aerodynamic_roll(phi_aero)

            cross_flow = np.array([0.0, velocity[1], velocity[2]])
            cross_flow = cross_flow/np.linalg.norm(cross_flow)
            np.testing.assert_allclose(np.dot(matrix_roll, [0.0, 0.0, 1.0]), cross_flow, atol=1.e-12)

            # 固有回転であること（モーメントの擬ベクトルにも同じ変換を使うため）
            np.testing.assert_allclose(np.dot(matrix_roll, matrix_roll.T), np.eye(3), atol=1.e-12)
            self.assertAlmostEqual(np.linalg.det(matrix_roll), 1.0, places=12)

    def test_axial_force_is_recovered_at_zero_incidence(self):
        # 迎角 0 では表の [CFx, 0, 0] がそのまま機体軸の軸力になる
        matrix_roll = attitude.matrix_aerodynamic_roll(0.0)
        np.testing.assert_allclose(np.dot(matrix_roll, [1.2, 0.0, 0.0]), [1.2, 0.0, 0.0], atol=1.e-14)


if __name__ == '__main__':
    unittest.main()
