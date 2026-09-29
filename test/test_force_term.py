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

# acceleration_routine が返す配列の行の意味
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
    bracket -= ratio ** 2 * factor['J22'] * 3.0 * cos_b ** 2 * np.cos(2.0 * (alpha - alpha22))
    bracket -= ratio ** 3 * factor['J3'] * (5.0 * sin_b ** 3 - 3.0 * sin_b) / 2.0
    bracket -= ratio ** 4 * factor['J4'] * (35.0 * sin_b ** 4 - 30.0 * sin_b ** 2 + 3.0) / 8.0

    return -gm / r * bracket


def gravity_from_code(config, coord):
    """acceleration_routine から重力成分だけを取り出す（空力は密度 0 で無効化）。"""
    acceleration = force_term.acceleration_initialsettings(config)
    acceleration = force_term.acceleration_routine(
        config, np.array(coord, dtype=float), np.zeros(3),
        mass_satellite=1.0, area_satellite=1.0,
        cdmean_aerodynamic=0.0, density_factor=0.0, density=0.0, acceleration=acceleration)
    return np.array(acceleration[IDX_GRAVITY, :])


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

            acceleration = force_term.acceleration_initialsettings(config)
            acceleration = force_term.acceleration_routine(
                config, coord, np.zeros(3), 1.0, 1.0, 0.0, 0.0, 0.0, acceleration)

            # 厳密に 0 になる成分があるので atol を併用する
            np.testing.assert_allclose(acceleration[IDX_GRAVITY, :], expected, rtol=1.0e-12, atol=1.0e-12)

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


class TestTheEarthJ22PointsTheRightWay(unittest.TestCase):
    """
    地球の J22 の長軸が西経 14.5 度にあること（EGM96 は西経 14.93 度）。

    README が挙げる Wagner（NASA TN D-3317 式 99、TN D-3557 付録 B）は cos 2(lon - Lambda22)
    で、Lambda22 は赤道の長軸の経度（J22 < 0）。
    2026-09-29 まで、コードと README は J22 の項を cos 2(lon + Lambda22) と書いていた。
    出典の J22 = -1.81222e-6, Lambda22 = -14.545 度をこれに入れると長軸は東経 14.5 度に
    来る（S22 の符号が逆）。TestGravityMatchesPotential は README の式とコードの一致しか
    見ないので、両方が同じ向きに間違っていると通ってしまう。そこで、コードの重力から
    C22 / S22 を逆算して、独立な重力場モデル（EGM96）と比べる。

    赤道上の重力の東向き成分は、帯状項（J2 / J3 / J4）からは出ず、
        g_east = 6 GM/r^2 (a/r)^2 (S22 cos 2lon - C22 sin 2lon)
    だけになる（C22, S22 は正規化していない係数）。
    """

    # EGM96 の完全正規化係数
    C22_NORMALISED = 0.243914352398e-5
    S22_NORMALISED = -0.140016683654e-5

    def test_the_long_axis_and_the_amplitude_match_egm96(self):
        config = load_config()
        gm = gravitational_parameter(config)
        radius_equat = config['planet']['radius']
        radius = radius_equat + 200.0e3

        longitude = np.linspace(-np.pi, np.pi, 72, endpoint=False)
        east = []
        for lon in longitude:
            coord = radius * np.array([np.cos(lon), np.sin(lon), 0.0])
            east.append(np.dot(gravity_from_code(config, coord), [-np.sin(lon), np.cos(lon), 0.0]))
        scale = 6.0 * gm / radius ** 2 * (radius_equat / radius) ** 2
        design = np.stack([np.cos(2.0 * longitude), -np.sin(2.0 * longitude)], axis=1) * scale
        s22, c22 = np.linalg.lstsq(design, np.array(east), rcond=None)[0]

        factor = np.sqrt(5.0 / 12.0)
        c22_model, s22_model = self.C22_NORMALISED * factor, self.S22_NORMALISED * factor
        axis = 0.5 * np.degrees(np.arctan2(s22, c22))
        axis_model = 0.5 * np.degrees(np.arctan2(s22_model, c22_model))

        # 出典の値は EGM96 から長軸で 0.39 度、振幅で 0.18 % ずれている（古い重力場）。
        # 符号を取り違えると長軸は東経 14.5 度に来て、29 度ずれる
        self.assertAlmostEqual(axis, -14.545, delta=1.0e-6)
        self.assertAlmostEqual(axis, axis_model, delta=0.5)
        self.assertAlmostEqual(np.hypot(c22, s22) / np.hypot(c22_model, s22_model), 1.0, delta=0.005)


class TestMarsGravityMatchesTheFieldModel(unittest.TestCase):
    """
    火星のチュートリアルの planet 節が、元の重力場モデルを再現すること。

    config の J2 / J22 / J3 / J4 / Lambda22 は JGMRO_120F（PDS の MRO 電波科学、
    基準半径 3396.0 km）の**完全正規化**係数を README の形に換算したもの。
    換算（正規化の係数、J_n = -C_n0 の符号、C22/S22 から J22/Lambda22 への角度）の
    どれかを取り違えると、コードが解く重力は元のモデルと食い違う。

    ここでは換算を経ずに、正規化係数のまま標準形のポテンシャル

        V = GM/r [ 1 + sum (R/r)^n Pbar_nm(sin lat) (Cbar_nm cos m lon + Sbar_nm sin m lon) ]

    を書き、その勾配とコードの重力を比べる。点質量の分を引いてから比べるので、
    主項の 1e-5 倍しかない C22/S22 の取り違えも見える。
    """

    CONFIG = os.path.join(ROOT_DIR, 'tutorial', 'work_reentry_mars', 'config.yml')

    # JGMRO_120F_SHA.TAB の先頭（完全正規化）と、ラベルにある火星単体の GM
    GM = 42828.3748574e9
    RADIUS = 3396.0e3
    C20 = -0.8750219819894000e-03
    C22 = -0.8463283575906001e-04
    S22 = 0.4893975901192000e-04
    C30 = -0.1189685487260000e-04
    C40 = 0.5129215056400000e-05

    def setUp(self):
        self.config = load_config(self.CONFIG)

    def potential_model(self, coord):
        x, y, z = coord
        r = np.sqrt(x * x + y * y + z * z)
        s = z / r
        c = np.sqrt(x * x + y * y) / r
        lon = np.arctan2(y, x)
        ratio = self.RADIUS / r
        # 完全正規化されたルジャンドル陪関数
        p20 = np.sqrt(5.0) * (3.0 * s ** 2 - 1.0) / 2.0
        p30 = np.sqrt(7.0) * (5.0 * s ** 3 - 3.0 * s) / 2.0
        p40 = 3.0 * (35.0 * s ** 4 - 30.0 * s ** 2 + 3.0) / 8.0
        p22 = np.sqrt(5.0 / 12.0) * 3.0 * c ** 2
        bracket = (1.0
                   + ratio ** 2 * p20 * self.C20
                   + ratio ** 2 * p22 * (self.C22 * np.cos(2.0 * lon) + self.S22 * np.sin(2.0 * lon))
                   + ratio ** 3 * p30 * self.C30
                   + ratio ** 4 * p40 * self.C40)
        return self.GM / r * bracket

    def gradient_model(self, coord, step=100.0):
        grad = np.zeros(3)
        for i in range(3):
            def shifted(delta):
                point = np.array(coord, dtype=float)
                point[i] += delta
                return self.potential_model(point)
            grad[i] = (-shifted(2.0 * step) + 8.0 * shifted(step)
                       - 8.0 * shifted(-step) + shifted(-2.0 * step)) / (12.0 * step)
        return grad

    def test_the_mass_times_the_constant_is_the_gm_of_mars(self):
        self.assertAlmostEqual(gravitational_parameter(self.config) / self.GM, 1.0, delta=1.0e-7)
        self.assertEqual(self.config['planet']['radius'], self.RADIUS)

    def test_the_gravity_is_the_gradient_of_the_field_model(self):
        rng = np.random.default_rng(20260929)
        for _ in range(12):
            lon = rng.uniform(-np.pi, np.pi)
            lat = rng.uniform(-1.4, 1.4)
            r = self.RADIUS + rng.uniform(0.0, 200.0e3)
            coord = r * np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])

            # 点質量の分はそれぞれの GM で引く（config の質量は 8 桁に丸めてあり、
            # 共通の GM で引くと丸めの 1.5e-9 が主項の 1e-3 倍の差を 6e-7 だけ汚す）
            expected = self.gradient_model(coord) + self.GM / r ** 3 * coord
            obtained = gravity_from_code(self.config, coord) \
                + gravitational_parameter(self.config) / r ** 3 * coord
            # 2026-09-29 の実測: 出荷値で相対 1e-8 の桁。Lambda22 を 1 度ずらすと 6e-3、
            # J2 を 0.1 % 変えると 1e-3、J22 / J3 / J4 の符号を反転すると 0.02 以上
            np.testing.assert_allclose(
                obtained, expected, rtol=0.0, atol=1.0e-5 * np.linalg.norm(expected),
                err_msg='火星の重力が JGMRO_120F と一致しない: lon=%.3f lat=%.3f' % (lon, lat))


class TestInertialAndAeroTerms(unittest.TestCase):

    def setUp(self):
        self.config = load_config()

    def test_coriolis_sign(self):
        """コリオリ項が -2 omega x v であること。"""
        omega = np.array([0.0, 0.0, self.config['planet']['rotation_rate']])
        velocity = np.array([1.0e3, -2.0e3, 5.0e2])
        coord = np.array([7.0e6, 0.0, 0.0])

        acceleration = force_term.acceleration_initialsettings(self.config)
        acceleration = force_term.acceleration_routine(
            self.config, coord, velocity, 1.0, 1.0, 0.0, 0.0, 0.0, acceleration)

        np.testing.assert_allclose(acceleration[IDX_CORIOLIS, :], -2.0 * np.cross(omega, velocity),
                                   rtol=1.0e-12, atol=1.0e-15)

    def test_centrifugal_sign(self):
        """遠心力項が -omega x (omega x x) であること。"""
        omega = np.array([0.0, 0.0, self.config['planet']['rotation_rate']])
        coord = np.array([4.0e6, -5.0e6, 2.0e6])

        acceleration = force_term.acceleration_initialsettings(self.config)
        acceleration = force_term.acceleration_routine(
            self.config, coord, np.zeros(3), 1.0, 1.0, 0.0, 0.0, 0.0, acceleration)

        expected = -np.cross(omega, np.cross(omega, coord))
        np.testing.assert_allclose(acceleration[IDX_CENTRIFUGAL, :], expected, rtol=1.0e-12, atol=1.0e-15)

    def test_drag_opposes_velocity(self):
        """抗力が速度と逆向きで、大きさが 1/2 rho Cd S |v|^2 / m であること。"""
        velocity = np.array([7.0e3, 1.0e3, -5.0e2])
        coord = np.array([7.0e6, 0.0, 0.0])
        density, cd, area, mass, factor = 1.0e-11, 2.2, 0.5, 4.0, 1.0

        acceleration = force_term.acceleration_initialsettings(self.config)
        acceleration = force_term.acceleration_routine(
            self.config, coord, velocity, mass, area, cd, factor, density, acceleration)

        drag = np.array(acceleration[IDX_AERO, :])
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

        acceleration = force_term.acceleration_initialsettings(self.config)
        acceleration = force_term.acceleration_routine(
            self.config, coord, velocity, 4.0, 0.5, 2.2, 1.0, 1.0e-11, acceleration)
        aero_without = np.array(acceleration[IDX_AERO, :])

        self.config['satellite']['lift_coefficient'] = 0.0
        acceleration = force_term.acceleration_initialsettings(self.config)
        acceleration = force_term.acceleration_routine(
            self.config, coord, velocity, 4.0, 0.5, 2.2, 1.0, 1.0e-11, acceleration)

        # ビット単位で同じ（CL = 0 は足し算すら通らない）
        np.testing.assert_array_equal(acceleration[IDX_AERO, :], aero_without)

    def _aero_with_lift(self, coefficient_lift, angle_bank, velocity, coord,
                        velocity_air=None, density=1.0e-5, cd=1.2, area=12.0, mass=5000.0):
        config = copy.deepcopy(self.config)
        config['satellite']['lift_coefficient'] = coefficient_lift
        config['satellite']['bank_angle'] = angle_bank
        acceleration = force_term.acceleration_initialsettings(config)
        acceleration = force_term.acceleration_routine(
            config, coord, velocity, mass, area, cd, 1.0, density, acceleration,
            velocity_air=velocity_air)
        aero = np.array(acceleration[IDX_AERO, :])

        # 抗力だけの分を引くと揚力が残る
        config['satellite']['lift_coefficient'] = 0.0
        acceleration = force_term.acceleration_initialsettings(config)
        acceleration = force_term.acceleration_routine(
            config, coord, velocity, mass, area, cd, 1.0, density, acceleration,
            velocity_air=velocity_air)
        return aero - np.array(acceleration[IDX_AERO, :]), (density, area, mass)

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
        acceleration = force_term.acceleration_initialsettings(config)
        acceleration = force_term.acceleration_routine(
            config, np.array([6.5e6, 1.0e6, 2.0e6]), np.array([7.0e3, 1.0e3, -5.0e2]),
            5000.0, 12.0, 1.2, 1.0, 1.0e-5, acceleration, acceleration_aerodynamic=given)
        np.testing.assert_array_equal(acceleration[IDX_AERO, :], given)

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
            acceleration = force_term.acceleration_initialsettings(config_tmp)
            acceleration = force_term.acceleration_routine(
                config_tmp, coord, velocity, 5000.0, 12.0, 1.2, 1.0, 1.0e-5, acceleration,
                angle_bank=angle_bank)
            return np.array(acceleration[IDX_AERO, :])

        drag = aero(0.0, coefficient_lift=0.0)
        # 引数 180 deg は 0 deg とは逆向きの揚力になる
        np.testing.assert_allclose(aero(180.0) - drag, -(aero(0.0) - drag), rtol=1.0e-12)
        # 引数を与えなければ config の値（0 deg）が使われる
        acceleration = force_term.acceleration_initialsettings(config)
        acceleration = force_term.acceleration_routine(
            config, coord, velocity, 5000.0, 12.0, 1.2, 1.0, 1.0e-5, acceleration)
        np.testing.assert_allclose(np.array(acceleration[IDX_AERO, :]), aero(0.0), rtol=1.0e-14)

    def test_total_is_the_sum_of_the_parts(self):
        acceleration = force_term.acceleration_initialsettings(self.config)
        acceleration = force_term.acceleration_routine(
            self.config, np.array([5.0e6, 3.0e6, 2.0e6]), np.array([1.0e3, 2.0e3, 3.0e3]),
            4.0, 0.5, 2.2, 1.0, 1.0e-11, acceleration)

        expected = acceleration[IDX_GRAVITY] + acceleration[IDX_CORIOLIS] + acceleration[IDX_CENTRIFUGAL] + acceleration[IDX_AERO]
        np.testing.assert_allclose(acceleration[IDX_TOTAL, :], expected, rtol=1.0e-14)


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
