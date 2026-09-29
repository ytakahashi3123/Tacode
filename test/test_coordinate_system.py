#!/usr/bin/env python3
"""座標変換のテスト。"""

import unittest

import numpy as np

from context import load_config

import coordinate_system.coordinate_system as cs

RAD2DEG = 180.0 / np.pi


class TestRoundTrip(unittest.TestCase):
    """直交 <-> 測地 <-> 直交 の往復で元に戻ること。"""

    def setUp(self):
        self.config = load_config()

    def test_cartesian_geodetic_cartesian(self):
        rng = np.random.default_rng(20260905)
        worst = 0.0
        for _ in range(500):
            direction = rng.normal(size=3)
            direction /= np.linalg.norm(direction)
            point = direction * rng.uniform(6.5e6, 7.5e6)

            geodetic = cs.convert_cartesian_geodetic(self.config, list(point))
            restored = cs.convert_geodetic_cartesian(self.config, geodetic)

            worst = max(worst, np.linalg.norm(np.array(restored) - point))

        # 変換は解析的だが平方根と三次方程式を経るので 1 um 程度の誤差は許容する
        self.assertLess(worst, 1.0e-6, 'cartesian -> geodetic -> cartesian の復元誤差が大きい')

    def test_geodetic_cartesian_geodetic(self):
        for lon_deg in (-179.0, -90.0, 0.0, 90.0, 179.0):
            for lat_deg in (-80.0, -45.0, 0.0, 45.0, 80.0):
                for alt in (0.0, 2.0e5, 5.0e5):
                    geodetic = [lon_deg / RAD2DEG, lat_deg / RAD2DEG, alt]
                    cartesian = cs.convert_geodetic_cartesian(self.config, geodetic)
                    restored = cs.convert_cartesian_geodetic(self.config, cartesian)

                    self.assertAlmostEqual(restored[0] * RAD2DEG, lon_deg, places=8)
                    self.assertAlmostEqual(restored[1] * RAD2DEG, lat_deg, places=8)
                    self.assertAlmostEqual(restored[2], alt, delta=1.0e-6)

    def test_polar_vector_transforms_are_inverse(self):
        """convert_carteasian_polar と convert_polar_carteasian が互いの逆変換であること。"""
        rng = np.random.default_rng(7)
        for _ in range(200):
            vec = rng.normal(size=3) * 1.0e3
            longitude = rng.uniform(-np.pi, np.pi)
            latitude = rng.uniform(-0.49 * np.pi, 0.49 * np.pi)

            polar = cs.convert_carteasian_polar(self.config, list(vec), longitude, latitude)
            restored = cs.convert_polar_carteasian(self.config, polar, longitude, latitude)

            np.testing.assert_allclose(restored, vec, rtol=1.0e-12, atol=1.0e-9)


class TestSingularities(unittest.TestCase):
    """かつて破綻していた特異点（レビュー B-4）の回帰テスト。"""

    def setUp(self):
        self.config = load_config()

    def test_longitude_180_degrees(self):
        """
        y = 0, x < 0 は経度 180 度。

        sign(y)*arccos(x/r) では sign(0) = 0 となり 0 度と誤算出されていた。
        """
        geodetic = cs.convert_cartesian_geodetic(self.config, [-6.6e6, 0.0, 1.0e6])
        self.assertAlmostEqual(abs(geodetic[0] * RAD2DEG), 180.0, places=9)

        polar = cs.set_angle_polar(self.config, [-6.6e6, 0.0, 1.0e6])
        self.assertAlmostEqual(abs(polar[2] * RAD2DEG), 180.0, places=9)

    def test_longitude_is_continuous_across_180(self):
        """経度 180 度をまたいでも値が飛ばないこと（絶対値で連続）。"""
        for y in (1.0e-3, 0.0, -1.0e-3):
            geodetic = cs.convert_cartesian_geodetic(self.config, [-6.6e6, y, 1.0e6])
            self.assertAlmostEqual(abs(geodetic[0] * RAD2DEG), 180.0, places=6)

    def test_pole_is_not_nan(self):
        """x = y = 0（極）で nan にならず、緯度 +-90 度と正しい高度を返すこと。"""
        radius_equat = self.config['planet']['radius']
        radius_polar = radius_equat * (1.0 - self.config['planet']['ellipticity'])

        for sign in (1.0, -1.0):
            z = sign * 6.4e6
            geodetic = cs.convert_cartesian_geodetic(self.config, [0.0, 0.0, z])

            self.assertTrue(np.isfinite(geodetic).all(), '極で nan が出ている')
            self.assertAlmostEqual(geodetic[1] * RAD2DEG, sign * 90.0, places=9)
            self.assertAlmostEqual(geodetic[2], abs(z) - radius_polar, delta=1.0e-6)

    def test_pole_matches_the_limit_from_nearby(self):
        """極の特別扱いが一般式の極限と連続していること。"""
        z = 6.4e6
        near = cs.convert_cartesian_geodetic(self.config, [1.0e-3, 0.0, z])
        at_pole = cs.convert_cartesian_geodetic(self.config, [0.0, 0.0, z])

        self.assertAlmostEqual(near[1], at_pole[1], places=9)
        self.assertAlmostEqual(near[2], at_pole[2], delta=1.0e-3)

    def test_equator_is_not_singular(self):
        """z を厳密に 0 にしても緯度 0 度・高度が正しいこと。"""
        radius_equat = self.config['planet']['radius']
        geodetic = cs.convert_cartesian_geodetic(self.config, [radius_equat + 2.0e5, 0.0, 0.0])

        self.assertTrue(np.isfinite(geodetic).all())
        self.assertAlmostEqual(geodetic[1] * RAD2DEG, 0.0, places=9)
        self.assertAlmostEqual(geodetic[2], 2.0e5, delta=1.0e-6)


class TestBelowTheSurfaceOnTheEquator(unittest.TestCase):
    """
    赤道面上 (z=0) で地表より下の点が、緯度 0・負の高度になること。

    B = sign(z) b と書いていた間は z=0 で B=0 になり、地表より下の点が
    緯度 90 度・高度 +0 と出た。ソルバーは高度 <= 0 で止まるので、赤道から
    真東に打った再突入は地表で止まらず、時間の上限まで地中を進んだ
    （2026-09-29、火星のケースの検証で見つけた。地球でも同じ）。
    """

    def test_the_point_is_below_the_surface(self):
        for config in (load_config(), {'planet': {'radius': 3396.0e3, 'ellipticity': 0.0}}):
            radius_equat = config['planet']['radius']
            for depth in (1.0, 100.0, 5.0e3):
                for point in ([radius_equat - depth, 0.0, 0.0], [0.0, -(radius_equat - depth), 0.0]):
                    geodetic = cs.convert_cartesian_geodetic(config, point)
                    self.assertAlmostEqual(geodetic[1], 0.0, places=12)
                    self.assertAlmostEqual(geodetic[2], -depth, delta=1.0e-6)


class TestASphericalPlanet(unittest.TestCase):
    """
    扁平率 0（球）でも変換が壊れないこと。

    火星のチュートリアルは ellipticity: 0.0 で走らせる（大気テーブルの高度が
    半径 3396 km の球からの高さで与えられているため）。一般式は扁平率で割る形を
    していないが、三次方程式の解き方が球の極限で退化しないかは確かめておく。
    球なら測地高度は r - a、測地緯度は地心緯度そのものになる。
    """

    def setUp(self):
        self.config = {'planet': {'radius': 3396.0e3, 'ellipticity': 0.0}}

    def test_the_altitude_is_the_distance_above_the_sphere(self):
        rng = np.random.default_rng(20260929)
        for _ in range(2000):
            direction = rng.normal(size=3)
            direction /= np.linalg.norm(direction)
            point = direction * rng.uniform(3.39e6, 3.70e6)
            radius = np.linalg.norm(point)

            geodetic = cs.convert_cartesian_geodetic(self.config, list(point))
            self.assertAlmostEqual(geodetic[2], radius - 3396.0e3, delta=1.0e-6)
            self.assertAlmostEqual(geodetic[1], np.arcsin(point[2] / radius), delta=1.0e-12)

            restored = cs.convert_geodetic_cartesian(self.config, geodetic)
            self.assertLess(np.linalg.norm(np.array(restored) - point), 1.0e-6)


class TestPolarAngles(unittest.TestCase):

    def setUp(self):
        self.config = load_config()

    def test_known_directions(self):
        radius = 7.0e6
        cases = [
            ([radius, 0.0, 0.0], 0.0, 0.0),
            ([0.0, radius, 0.0], 0.0, 90.0),
            ([0.0, -radius, 0.0], 0.0, -90.0),
            ([0.0, 0.0, radius], 90.0, 0.0),
        ]
        for coord, beta_deg, alpha_deg in cases:
            polar = cs.set_angle_polar(self.config, coord)
            self.assertAlmostEqual(polar[0], radius, delta=1.0e-6)
            self.assertAlmostEqual(polar[1] * RAD2DEG, beta_deg, places=9)
            self.assertAlmostEqual(polar[2] * RAD2DEG, alpha_deg, places=9)


if __name__ == '__main__':
    unittest.main()
