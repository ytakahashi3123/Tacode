#!/usr/bin/env python3
"""
力の項のテスト。

中心となるのは「重力ポテンシャル U を数値微分した -grad U が、
force_term が返す重力ベクトルと一致するか」という検査である。
U はテスト側で README の定義から独立に実装しているので、
J 項の係数や符号を取り違えるとこのテストが落ちる。
"""

import copy
import os
import unittest

import numpy as np

from context import ROOT_DIR, load_config, quiet, two_body_config, gravitational_parameter

import atmosphere.atmosphere as atmosphere
import force_term.force_term as force_term
import satellite.satellite as satellite
import solver.solver as solver
from orbital.orbital import orbital

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

    def test_lift_is_off_by_default(self):
        """lift_coefficient が無い config では、空力は従来の抗力そのもの。"""
        velocity = np.array([7.0e3, 1.0e3, -5.0e2])
        coord = np.array([7.0e6, 0.0, 0.0])
        self.assertNotIn('lift_coefficient', self.config['satellite'])

        force = force_term.force_initialsettings(self.config)
        force = force_term.force_routine(
            self.config, coord, velocity, 4.0, 0.5, 2.2, 1.0, 1.0e-11, force)
        aero_without = np.array(force[IDX_AERO, :])

        self.config['satellite']['lift_coefficient'] = 0.0
        force = force_term.force_initialsettings(self.config)
        force = force_term.force_routine(
            self.config, coord, velocity, 4.0, 0.5, 2.2, 1.0, 1.0e-11, force)

        # ビット単位で同じ（CL = 0 は足し算すら通らない）
        np.testing.assert_array_equal(force[IDX_AERO, :], aero_without)

    def _aero_with_lift(self, coefficient_lift, angle_bank, velocity, coord,
                        velocity_air=None, density=1.0e-5, cd=1.2, area=12.0, mass=5000.0):
        config = copy.deepcopy(self.config)
        config['satellite']['lift_coefficient'] = coefficient_lift
        config['satellite']['bank_angle'] = angle_bank
        force = force_term.force_initialsettings(config)
        force = force_term.force_routine(
            config, coord, velocity, mass, area, cd, 1.0, density, force,
            velocity_air=velocity_air)
        aero = np.array(force[IDX_AERO, :])

        # 抗力だけの分を引くと揚力が残る
        config['satellite']['lift_coefficient'] = 0.0
        force = force_term.force_initialsettings(config)
        force = force_term.force_routine(
            config, coord, velocity, mass, area, cd, 1.0, density, force,
            velocity_air=velocity_air)
        return aero - np.array(force[IDX_AERO, :]), (density, area, mass)

    def test_lift_is_perpendicular_to_the_velocity_and_has_the_right_magnitude(self):
        velocity = np.array([7.0e3, 1.0e3, -5.0e2])
        coord = np.array([6.5e6, 1.0e6, 2.0e6])
        coefficient_lift = 0.4

        lift, (density, area, mass) = self._aero_with_lift(coefficient_lift, 0.0, velocity, coord)
        speed = np.linalg.norm(velocity)

        self.assertAlmostEqual(np.dot(lift, velocity)/(np.linalg.norm(lift)*speed), 0.0, places=12)
        np.testing.assert_allclose(np.linalg.norm(lift),
                                   0.5*density*coefficient_lift*area*speed**2/mass,
                                   rtol=1.0e-12)

    def test_the_bank_angle_turns_the_lift_about_the_velocity(self):
        velocity = np.array([7.0e3, 1.0e3, -5.0e2])
        coord = np.array([6.5e6, 1.0e6, 2.0e6])
        direction_up = coord/np.linalg.norm(coord)
        direction_velocity = velocity/np.linalg.norm(velocity)
        direction_right = np.cross(direction_velocity, direction_up)
        direction_right = direction_right/np.linalg.norm(direction_right)

        lift_up, _ = self._aero_with_lift(0.4, 0.0, velocity, coord)
        lift_down, _ = self._aero_with_lift(0.4, 180.0, velocity, coord)
        lift_right, _ = self._aero_with_lift(0.4, 90.0, velocity, coord)

        # バンク 0 は上向き側、180 度はその真逆
        self.assertGreater(np.dot(lift_up, direction_up), 0.0)
        np.testing.assert_allclose(lift_down, -lift_up, rtol=1.0e-12)
        # バンク 90 度は鉛直成分を持たず、進行方向の右を向く
        self.assertAlmostEqual(np.dot(lift_right, direction_up)/np.linalg.norm(lift_right),
                               0.0, places=12)
        self.assertGreater(np.dot(lift_right, direction_right), 0.0)
        # 大きさはバンク角によらない
        np.testing.assert_allclose(np.linalg.norm(lift_right), np.linalg.norm(lift_up), rtol=1.0e-12)

    def test_the_lift_follows_the_air_relative_velocity(self):
        # 風があるときは対気速度に直交する（抗力と同じ扱い）
        velocity = np.array([7.0e3, 1.0e3, -5.0e2])
        velocity_air = velocity - np.array([0.0, 300.0, 0.0])
        coord = np.array([6.5e6, 1.0e6, 2.0e6])

        lift, _ = self._aero_with_lift(0.4, 0.0, velocity, coord, velocity_air=velocity_air)

        self.assertAlmostEqual(np.dot(lift, velocity_air)/(np.linalg.norm(lift)*np.linalg.norm(velocity_air)),
                               0.0, places=12)
        self.assertGreater(abs(np.dot(lift, velocity)), 0.0)

    def test_the_six_degree_of_freedom_path_is_untouched(self):
        # 姿勢を解いているときは呼び出し側が空力加速度を作るので、CL は使われない
        config = copy.deepcopy(self.config)
        config['satellite']['lift_coefficient'] = 0.4
        given = np.array([1.0, -2.0, 3.0])
        force = force_term.force_initialsettings(config)
        force = force_term.force_routine(
            config, np.array([6.5e6, 1.0e6, 2.0e6]), np.array([7.0e3, 1.0e3, -5.0e2]),
            5000.0, 12.0, 1.2, 1.0, 1.0e-5, force, force_aerodynamic=given)
        np.testing.assert_array_equal(force[IDX_AERO, :], given)

    def test_the_bank_angle_table_is_interpolated_and_clamped(self):
        """satellite.bank_angle_table は時刻について線形内挿し、両端の外は端の値。"""
        config = copy.deepcopy(self.config)
        config['satellite']['bank_angle_table'] = [[0.0, 0.0], [10.0, 90.0], [20.0, 90.0]]

        self.assertAlmostEqual(force_term.get_bank_angle(config, 0.0), 0.0)
        self.assertAlmostEqual(force_term.get_bank_angle(config, 5.0), 45.0)
        self.assertAlmostEqual(force_term.get_bank_angle(config, 10.0), 90.0)
        self.assertAlmostEqual(force_term.get_bank_angle(config, 15.0), 90.0)
        # 範囲外は端の値
        self.assertAlmostEqual(force_term.get_bank_angle(config, -100.0), 0.0)
        self.assertAlmostEqual(force_term.get_bank_angle(config, 1.0e4), 90.0)

    def test_the_table_wins_over_the_constant(self):
        config = copy.deepcopy(self.config)
        config['satellite']['bank_angle'] = 30.0
        self.assertAlmostEqual(force_term.get_bank_angle(config, 0.0), 30.0)
        config['satellite']['bank_angle_table'] = [[0.0, 60.0], [10.0, 60.0]]
        self.assertAlmostEqual(force_term.get_bank_angle(config, 5.0), 60.0)

    def test_a_malformed_table_stops_the_run(self):
        for table in ([[0.0, 0.0]],                      # 1 点だけ
                      [[0.0, 0.0], [0.0, 90.0]],         # 時刻が増えない
                      [[10.0, 0.0], [0.0, 90.0]],        # 時刻が逆
                      [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]):   # 列数が違う
            config = copy.deepcopy(self.config)
            config['satellite']['bank_angle_table'] = table
            with self.assertRaises(SystemExit):
                with quiet():
                    force_term.get_bank_angle(config, 1.0)

    def test_the_lift_uses_the_bank_handed_in(self):
        # solver は段の時刻で引いた値を渡す。渡された値が使われること
        velocity = np.array([7.0e3, 1.0e3, -5.0e2])
        coord = np.array([6.5e6, 1.0e6, 2.0e6])
        config = copy.deepcopy(self.config)
        config['satellite']['lift_coefficient'] = 0.4
        config['satellite']['bank_angle'] = 0.0

        def aero(angle_bank, coefficient_lift=0.4):
            config_tmp = copy.deepcopy(config)
            config_tmp['satellite']['lift_coefficient'] = coefficient_lift
            force = force_term.force_initialsettings(config_tmp)
            force = force_term.force_routine(
                config_tmp, coord, velocity, 5000.0, 12.0, 1.2, 1.0, 1.0e-5, force,
                angle_bank=angle_bank)
            return np.array(force[IDX_AERO, :])

        drag = aero(0.0, coefficient_lift=0.0)
        # 引数 180 deg は 0 deg とは逆向きの揚力になる
        np.testing.assert_allclose(aero(180.0) - drag, -(aero(0.0) - drag), rtol=1.0e-12)
        # 引数を与えなければ config の値（0 deg）が使われる
        force = force_term.force_initialsettings(config)
        force = force_term.force_routine(
            config, coord, velocity, 5000.0, 12.0, 1.2, 1.0, 1.0e-5, force)
        np.testing.assert_allclose(np.array(force[IDX_AERO, :]), aero(0.0), rtol=1.0e-14)

    def test_total_is_the_sum_of_the_parts(self):
        force = force_term.force_initialsettings(self.config)
        force = force_term.force_routine(
            self.config, np.array([5.0e6, 3.0e6, 2.0e6]), np.array([1.0e3, 2.0e3, 3.0e3]),
            4.0, 0.5, 2.2, 1.0, 1.0e-11, force)

        expected = force[IDX_GRAVITY] + force[IDX_CORIOLIS] + force[IDX_CENTRIFUGAL] + force[IDX_AERO]
        np.testing.assert_allclose(force[IDX_TOTAL, :], expected, rtol=1.0e-14)


class TestTheLiftThroughTheSolver(unittest.TestCase):
    """
    solver との結線。

    バンク角は**各 Runge-Kutta 段の時刻で**引かれる（風と同じ扱い）。定数で与えたときと
    同じ値の表を与えたときが一致すること、表で向きを変えると軌道が変わることを見る。
    """

    TIME_MAX = 120.0

    @classmethod
    def setUpClass(cls):
        config = load_config(os.path.join(ROOT_DIR, 'tutorial', 'work_reentry', 'config.yml'))
        config['computational_setup']['time_elapsed_maximum'] = cls.TIME_MAX
        config['satellite']['kind_aerodynamic_model'] = 'constant'
        config['satellite']['drag_coefficient'] = 1.0
        cls.config = config

    def _run(self, **setting):
        config = copy.deepcopy(self.config)
        config['satellite'].update(setting)
        with quiet():
            atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
            aerodynamic_dict = satellite.initial_settings_satellite(config)
            orb = orbital()
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                orb.initial_settings(config)
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                solver.solve_equation_motion(
                    config, iteration, time_elapsed, coordinate_dict, velocity_dict,
                    trajectory_dict, atmosphere_dict, aerodynamic_dict)
        return np.array(coordinate_dict['cartesian'][-1])

    def test_a_constant_table_reproduces_the_scalar_bank(self):
        position_scalar = self._run(lift_coefficient=0.4, bank_angle=45.0)
        position_table = self._run(lift_coefficient=0.4,
                                   bank_angle_table=[[0.0, 45.0], [1000.0, 45.0]])
        np.testing.assert_array_equal(position_table, position_scalar)

    def test_no_lift_is_the_earlier_result(self):
        position_without = self._run()
        position_zero = self._run(lift_coefficient=0.0, bank_angle=45.0)
        np.testing.assert_array_equal(position_zero, position_without)

    def test_the_bank_changes_the_trajectory(self):
        position_up = self._run(lift_coefficient=0.4, bank_angle=0.0)
        position_down = self._run(lift_coefficient=0.4, bank_angle=180.0)
        position_none = self._run()
        # 揚力を上に向けたほうが地心距離が大きい
        self.assertGreater(np.linalg.norm(position_up), np.linalg.norm(position_none))
        self.assertLess(np.linalg.norm(position_down), np.linalg.norm(position_none))

    def test_a_table_which_switches_is_between_the_two_constants(self):
        # 途中でリフトアップからリフトダウンへ倒す表。結果は両者の間に入る
        position_up = self._run(lift_coefficient=0.4, bank_angle=0.0)
        position_down = self._run(lift_coefficient=0.4, bank_angle=180.0)
        position_switch = self._run(lift_coefficient=0.4,
                                    bank_angle_table=[[0.0, 0.0], [59.0, 0.0],
                                                      [61.0, 180.0], [1000.0, 180.0]])
        radius = [np.linalg.norm(item) for item in (position_down, position_switch, position_up)]
        self.assertLess(radius[0], radius[1])
        self.assertLess(radius[1], radius[2])


if __name__ == '__main__':
    unittest.main()
