#!/usr/bin/env python3
"""
力の項のテスト。

中心となるのは「重力ポテンシャル U を数値微分した -grad U が、
force_term が返す重力ベクトルと一致するか」という検査である。
U はテスト側で README の定義から独立に実装しているので、
J 項の係数や符号を取り違えるとこのテストが落ちる。
"""

import unittest

import numpy as np

from context import load_config, two_body_config, gravitational_parameter

import force_term.force_term as force_term

# force_routine が返す配列の行の意味
IDX_TOTAL = 0
IDX_GRAVITY = 1
IDX_CORIOLIS = 2
IDX_CENTRIFUGAL = 3
IDX_AERO = 4


def potential(config, coord):
    """
    単位質量あたりの重力ポテンシャル U(x, y, z)。

    README の定義をそのまま書き下したもの。beta は極座標の緯度、
    alpha は極座標の経度。P_n はルジャンドル多項式。
    """
    gm = gravitational_parameter(config)
    radius_equat = config['planet']['radius']
    factor = config['planet']['potential_factor']

    x, y, z = coord
    r = np.sqrt(x * x + y * y + z * z)
    sin_b = z / r
    cos_b = np.sqrt(x * x + y * y) / r
    alpha = np.arctan2(y, x)

    alpha22 = factor['Lambda22'] * np.pi / 180.0
    ratio = radius_equat / r

    bracket = 1.0
    bracket -= ratio ** 2 * factor['J2'] * (3.0 * sin_b ** 2 - 1.0) / 2.0
    bracket -= ratio ** 2 * factor['J22'] * 3.0 * cos_b ** 2 * np.cos(2.0 * (alpha + alpha22))
    bracket -= ratio ** 3 * factor['J3'] * (5.0 * sin_b ** 3 - 3.0 * sin_b) / 2.0
    bracket -= ratio ** 4 * factor['J4'] * (35.0 * sin_b ** 4 - 30.0 * sin_b ** 2 + 3.0) / 8.0

    return -gm / r * bracket


def gravity_from_code(config, coord):
    """force_routine から重力成分だけを取り出す（空力は密度 0 で無効化）。"""
    force = force_term.force_initialsettings(config)
    force = force_term.force_routine(
        config, np.array(coord, dtype=float), np.zeros(3),
        mass_satellite=1.0, area_satellite=1.0,
        cdmean_aerodynamic=0.0, density_factor=0.0, density=0.0, force=force)
    return np.array(force[IDX_GRAVITY, :])


class TestGravityMatchesPotential(unittest.TestCase):
    """-grad U と force_term の重力が一致すること。"""

    def setUp(self):
        self.config = load_config()

    def _numerical_gradient(self, coord, step=200.0):
        """
        5 点公式による勾配。打ち切り誤差 O(h^4)。

        2 点の中心差分だと h を小さくしたときの桁落ちが J22/J30 起因の
        微小成分（主成分の 10^-7 倍）を埋めてしまうため、高次の公式を使う。
        """
        grad = np.zeros(3)
        for i in range(3):
            def shifted(delta):
                point = np.array(coord, dtype=float)
                point[i] += delta
                return potential(self.config, point)

            grad[i] = (-shifted(2.0 * step) + 8.0 * shifted(step)
                       - 8.0 * shifted(-step) + shifted(-2.0 * step)) / (12.0 * step)
        return grad

    def test_matches_at_several_points(self):
        points = [
            [6.8e6, 0.0, 0.0],
            [0.0, 6.9e6, 0.0],
            [3.0e6, -4.0e6, 4.5e6],
            [-5.0e6, 2.0e6, -3.5e6],
            [-6.8e6, 0.0, 1.0e6],      # 経度 180 度
            [1.0e5, 1.0e5, 7.0e6],     # 極の近く
        ]
        for coord in points:
            expected = -self._numerical_gradient(coord)
            obtained = gravity_from_code(self.config, coord)
            np.testing.assert_allclose(
                obtained, expected, rtol=1.0e-5, atol=1.0e-9,
                err_msg='重力が -grad U と一致しない: coord=%s' % (coord,))

    def test_point_mass_limit(self):
        """J 項をすべて 0 にすると -GM/r^2 * r_hat になること。"""
        config = two_body_config()
        gm = gravitational_parameter(config)

        for coord in ([7.0e6, 0.0, 0.0], [0.0, -6.9e6, 1.0e6], [2.0e6, 3.0e6, -5.0e6]):
            coord = np.array(coord, dtype=float)
            r = np.linalg.norm(coord)
            expected = -gm / r ** 3 * coord

            force = force_term.force_initialsettings(config)
            force = force_term.force_routine(
                config, coord, np.zeros(3), 1.0, 1.0, 0.0, 0.0, 0.0, force)

            # 厳密に 0 になる成分があるので atol を併用する
            np.testing.assert_allclose(force[IDX_GRAVITY, :], expected, rtol=1.0e-12, atol=1.0e-12)

    def test_j2_makes_equatorial_pull_stronger_than_polar(self):
        """
        J2 > 0（扁平）なら、同じ半径では赤道方向の引力が極方向より強い。
        符号を取り違えるとこの関係が逆転する。
        """
        radius = 7.0e6
        equator = gravity_from_code(self.config, [radius, 0.0, 0.0])
        pole = gravity_from_code(self.config, [0.0, 0.0, radius])

        self.assertGreater(self.config['planet']['potential_factor']['J2'], 0.0)
        self.assertGreater(np.linalg.norm(equator), np.linalg.norm(pole))


class TestInertialAndAeroTerms(unittest.TestCase):

    def setUp(self):
        self.config = load_config()

    def test_coriolis_sign(self):
        """コリオリ項が -2 omega x v であること。"""
        omega = np.array([0.0, 0.0, self.config['planet']['rotation_rate']])
        velocity = np.array([1.0e3, -2.0e3, 5.0e2])
        coord = np.array([7.0e6, 0.0, 0.0])

        force = force_term.force_initialsettings(self.config)
        force = force_term.force_routine(
            self.config, coord, velocity, 1.0, 1.0, 0.0, 0.0, 0.0, force)

        np.testing.assert_allclose(force[IDX_CORIOLIS, :], -2.0 * np.cross(omega, velocity),
                                   rtol=1.0e-12, atol=1.0e-15)

    def test_centrifugal_sign(self):
        """遠心力項が -omega x (omega x x) であること。"""
        omega = np.array([0.0, 0.0, self.config['planet']['rotation_rate']])
        coord = np.array([4.0e6, -5.0e6, 2.0e6])

        force = force_term.force_initialsettings(self.config)
        force = force_term.force_routine(
            self.config, coord, np.zeros(3), 1.0, 1.0, 0.0, 0.0, 0.0, force)

        expected = -np.cross(omega, np.cross(omega, coord))
        np.testing.assert_allclose(force[IDX_CENTRIFUGAL, :], expected, rtol=1.0e-12, atol=1.0e-15)

    def test_drag_opposes_velocity(self):
        """抗力が速度と逆向きで、大きさが 1/2 rho Cd S |v|^2 / m であること。"""
        velocity = np.array([7.0e3, 1.0e3, -5.0e2])
        coord = np.array([7.0e6, 0.0, 0.0])
        density, cd, area, mass, factor = 1.0e-11, 2.2, 0.5, 4.0, 1.0

        force = force_term.force_initialsettings(self.config)
        force = force_term.force_routine(
            self.config, coord, velocity, mass, area, cd, factor, density, force)

        drag = np.array(force[IDX_AERO, :])
        speed = np.linalg.norm(velocity)

        # 向き: 速度の逆
        np.testing.assert_allclose(drag / np.linalg.norm(drag), -velocity / speed, rtol=1.0e-12)
        # 大きさ
        self.assertAlmostEqual(np.linalg.norm(drag),
                               0.5 * density * cd * area * speed ** 2 / mass,
                               delta=1.0e-18)

    def test_total_is_the_sum_of_the_parts(self):
        force = force_term.force_initialsettings(self.config)
        force = force_term.force_routine(
            self.config, np.array([5.0e6, 3.0e6, 2.0e6]), np.array([1.0e3, 2.0e3, 3.0e3]),
            4.0, 0.5, 2.2, 1.0, 1.0e-11, force)

        expected = force[IDX_GRAVITY] + force[IDX_CORIOLIS] + force[IDX_CENTRIFUGAL] + force[IDX_AERO]
        np.testing.assert_allclose(force[IDX_TOTAL, :], expected, rtol=1.0e-14)


if __name__ == '__main__':
    unittest.main()
