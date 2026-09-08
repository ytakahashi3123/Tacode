#!/usr/bin/env python3
"""大気・空力テーブルの読み込みと内挿のテスト。"""

import os
import tempfile
import unittest

import numpy as np

from context import load_config, quiet

import atmosphere.atmosphere as atmosphere
import satellite.satellite as satellite


TABLE_FILENAME = 'aerodynamic_test.txt'
COLUMN_SPARE = '\t'.join(['0.0']*6)


def write_aerodynamic_table(block, directory):
    """
    合成の空力係数表を書く。

    block は (迎角, 行の並び) の並びで、行は
    [Kn, CFx, CFy, CFz, CMx, CMy, CMz, Altitude]。迎角に None を渡すと
    "AOA" 行を書かない（従来形式のファイル）。標準偏差の 6 列は 0 で埋める。
    """
    path = os.path.join(directory, TABLE_FILENAME)
    with open(path, 'w') as f:
        f.write('test table\n')
        f.write('variables = Kn, CFx, CFy, CFz, CMx, CMy, CMz, SDV..., Altitude\n')
        for angle, rows in block:
            if angle is not None:
                f.write('AOA {:g}\n'.format(angle))
            for row in rows:
                f.write('\t'.join(['{:.18e}'.format(value) for value in row[0:7]])
                        + '\t' + COLUMN_SPARE + '\t{:.18e}\n'.format(row[7]))
    return path


def aerodynamic_config(directory):
    """合成テーブルを指す config を返す。"""
    config = load_config()
    config['satellite']['directory_path_specify'] = 'manual'
    config['satellite']['directory_aerodynamic'] = directory
    config['satellite']['filename_aerodynamic'] = TABLE_FILENAME
    return config


def read_aerodynamic_table(block):
    """合成テーブルを一時ディレクトリに書いて読み込む。"""
    with tempfile.TemporaryDirectory() as directory:
        write_aerodynamic_table(block, directory)
        with quiet():
            return satellite.initial_settings_satellite(aerodynamic_config(directory))


class TestAtmosphereTable(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.config = load_config()
        with quiet():
            cls.atm = atmosphere.initial_settings_atmosphere(cls.config)

    def test_table_is_loaded(self):
        height = self.atm[atmosphere.KEY_Height]
        self.assertEqual(len(height), self.atm[atmosphere.KEY_DATA])
        self.assertTrue(np.all(np.diff(height) > 0.0), '高度が単調増加でない')
        self.assertGreater(self.atm[atmosphere.KEY_Mass_density].min(), 0.0)
        self.assertGreater(self.atm[atmosphere.KEY_Temperature_neutral].min(), 0.0)
        self.assertGreater(self.atm[atmosphere.KEY_KN].min(), 0.0)

    def test_interpolator_is_prebuilt(self):
        """
        補間器が初期化時に構築されて dict に入っていること（レビュー C-1 の回帰）。

        評価のたびに構築し直す実装に戻ると、この鍵が無くなるか使われなくなる。
        """
        self.assertIn(atmosphere.KEY_INTERP, self.atm)
        for key in (atmosphere.KEY_Mass_density,
                    atmosphere.KEY_Temperature_neutral,
                    atmosphere.KEY_KN):
            self.assertIn(key, self.atm[atmosphere.KEY_INTERP])

    def _property_at(self, altitude_km, atm=None):
        return atmosphere.get_atmosphere_property(altitude_km, atm if atm is not None else self.atm)

    def test_interpolation_reproduces_table_nodes(self):
        """テーブルの格子点では元の値をそのまま返すこと。"""
        height = self.atm[atmosphere.KEY_Height]
        for index in (1, 50, 200, 399):
            density, temperature, knudsen = self._property_at(height[index])
            self.assertAlmostEqual(float(density) / self.atm[atmosphere.KEY_Mass_density][index],
                                   1.0, places=9)
            self.assertAlmostEqual(float(temperature) / self.atm[atmosphere.KEY_Temperature_neutral][index],
                                   1.0, places=9)
            self.assertAlmostEqual(float(knudsen) / self.atm[atmosphere.KEY_KN][index],
                                   1.0, places=9)

    def test_below_the_table_is_clamped(self):
        """下端より下は端の値で止まること。"""
        height = self.atm[atmosphere.KEY_Height]
        below = self._property_at(height[0] - 10.0)
        self.assertEqual(float(below[0]), self.atm[atmosphere.KEY_Mass_density][0])
        self.assertEqual(float(below[1]), self.atm[atmosphere.KEY_Temperature_neutral][0])
        self.assertEqual(float(below[2]), self.atm[atmosphere.KEY_KN][0])

    def test_density_decreases_with_altitude(self):
        """高度が上がると密度は下がること（物理的な健全性）。"""
        altitudes = [100.0, 150.0, 200.0, 300.0, 399.0]
        densities = [float(self._property_at(a)[0]) for a in altitudes]
        for lower, higher in zip(densities[:-1], densities[1:]):
            self.assertLess(higher, lower)
        self.assertTrue(all(d > 0.0 for d in densities), '内挿で密度が負になっている')

    def test_knudsen_scales_inversely_with_characteristic_length(self):
        """Kn は代表長さに反比例すること。"""
        config = load_config()
        config['satellite']['characteristic_length'] *= 2.0
        with quiet():
            doubled = atmosphere.initial_settings_atmosphere(config)

        np.testing.assert_allclose(doubled[atmosphere.KEY_KN],
                                   self.atm[atmosphere.KEY_KN] / 2.0, rtol=1.0e-12)


class TestExtrapolationAboveTheTable(unittest.TestCase):
    """
    テーブル上端より上の扱い（レビュー B-9）。

    熱圏上部は等温・拡散平衡なので密度は指数関数で減る。
    端の値で止める（クランプ）と 522 km で 11 倍以上の過大評価になっていた。
    """

    @classmethod
    def setUpClass(cls):
        with quiet():
            cls.atm = atmosphere.initial_settings_atmosphere(load_config())

        config_clamp = load_config()
        config_clamp['atmosphere']['kind_extrapolation'] = atmosphere.KIND_EXTRAPOLATION_CLAMP
        with quiet():
            cls.atm_clamp = atmosphere.initial_settings_atmosphere(config_clamp)

    def _at(self, altitude_km, atm):
        return atmosphere.get_atmosphere_property(altitude_km, atm)

    def test_scale_height_matches_the_table(self):
        """
        フィットしたスケールハイトが、上端の組成から求まる kT/(m g) と合うこと。

        400 km では原子状酸素が質量の 96% を占めるので H = kT/(m_O g) になる。
        """
        scale_height = self.atm[atmosphere.KEY_SCALE_HEIGHT]
        self.assertIsNotNone(scale_height)

        config = load_config()
        boltzmann = 1.380649e-23
        atomic_mass_unit = 1.66053907e-27
        gravity = (config['planet']['gravitational_constant'] * config['planet']['mass']
                   / (config['planet']['radius'] + self.atm[atmosphere.KEY_Height][-1] * 1.0e3) ** 2)
        expected = (boltzmann * self.atm[atmosphere.KEY_Temperature_neutral][-1]
                    / (16.0 * atomic_mass_unit * gravity) / 1.0e3)

        self.assertAlmostEqual(scale_height / expected, 1.0, delta=0.1,
                               msg='H = %.2f km, kT/(m_O g) = %.2f km' % (scale_height, expected))

    def test_density_decays_exponentially(self):
        """外挿区間で密度が exp(-dh/H) に従うこと。"""
        height = self.atm[atmosphere.KEY_Height]
        scale_height = self.atm[atmosphere.KEY_SCALE_HEIGHT]
        density_top = self.atm[atmosphere.KEY_Mass_density][-1]

        for delta in (20.0, 50.0, 100.0, 122.0, 200.0):
            density = float(self._at(height[-1] + delta, self.atm)[0])
            expected = density_top * np.exp(-delta / scale_height)
            self.assertAlmostEqual(density / expected, 1.0, places=9)

    def test_knudsen_grows_and_temperature_stays(self):
        """Kn は数密度に反比例して増え、温度は等温なので据え置きであること。"""
        height = self.atm[atmosphere.KEY_Height]
        delta = 100.0
        density, temperature, knudsen = self._at(height[-1] + delta, self.atm)

        ratio = float(density) / self.atm[atmosphere.KEY_Mass_density][-1]
        self.assertAlmostEqual(float(knudsen) / self.atm[atmosphere.KEY_KN][-1], 1.0 / ratio, places=9)
        self.assertEqual(float(temperature), self.atm[atmosphere.KEY_Temperature_neutral][-1])

    def test_it_is_continuous_at_the_top_of_the_table(self):
        """テーブル上端で内挿と外挿が滑らかにつながること。"""
        height = self.atm[atmosphere.KEY_Height]
        inside = float(self._at(height[-1] - 1.0e-6, self.atm)[0])
        outside = float(self._at(height[-1] + 1.0e-6, self.atm)[0])
        self.assertAlmostEqual(inside / outside, 1.0, places=6)

    def test_clamp_mode_reproduces_the_old_behaviour(self):
        """kind_extrapolation: clamp なら従来どおり端の値で止まること。"""
        height = self.atm_clamp[atmosphere.KEY_Height]
        self.assertIsNone(self.atm_clamp[atmosphere.KEY_SCALE_HEIGHT])

        above = self._at(height[-1] + 200.0, self.atm_clamp)
        self.assertEqual(float(above[0]), self.atm_clamp[atmosphere.KEY_Mass_density][-1])
        self.assertEqual(float(above[1]), self.atm_clamp[atmosphere.KEY_Temperature_neutral][-1])
        self.assertEqual(float(above[2]), self.atm_clamp[atmosphere.KEY_KN][-1])

    def test_extrapolation_is_far_below_the_clamped_value(self):
        """外挿はクランプより十分小さいこと（522 km で 10 倍以上の差）。"""
        height = self.atm[atmosphere.KEY_Height]
        altitude = 522.0
        if altitude <= height[-1]:
            self.skipTest('テーブルが %.0f km まであるので外挿が起きない' % height[-1])

        extrapolated = float(self._at(altitude, self.atm)[0])
        clamped = float(self._at(altitude, self.atm_clamp)[0])
        self.assertGreater(clamped / extrapolated, 10.0)

    def test_unknown_mode_is_rejected(self):
        config = load_config()
        config['atmosphere']['kind_extrapolation'] = 'bogus'
        with self.assertRaises(SystemExit):
            with quiet():
                atmosphere.initial_settings_atmosphere(config)

class TestAerodynamicTable(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.config = load_config()
        with quiet():
            cls.aero = satellite.initial_settings_satellite(cls.config)

    def test_interpolator_is_prebuilt(self):
        self.assertIn(satellite.KEY_INTERP, self.aero)
        self.assertIn(satellite.KEY_CD_MEAN, self.aero[satellite.KEY_INTERP])

    def _cd_at(self, knudsen):
        return satellite.get_aerodynamic_coefficient(
            knudsen,
            self.aero[satellite.KEY_KN],
            self.aero[satellite.KEY_CD_MEAN],
            self.aero[satellite.KEY_INTERP])

    def test_reproduces_table_nodes(self):
        kn_table = self.aero[satellite.KEY_KN]
        cd_table = self.aero[satellite.KEY_CD_MEAN]
        for index in range(1, len(kn_table) - 1):
            self.assertAlmostEqual(float(self._cd_at(kn_table[index])), cd_table[index], places=9)

    def test_values_are_clamped_outside_the_table(self):
        kn_table = self.aero[satellite.KEY_KN]
        cd_table = self.aero[satellite.KEY_CD_MEAN]

        self.assertEqual(float(self._cd_at(kn_table[0] * 0.5)), cd_table[0])
        self.assertEqual(float(self._cd_at(kn_table[-1] * 2.0)), cd_table[-1])

    def test_cd_is_positive_over_the_table_range(self):
        kn_table = self.aero[satellite.KEY_KN]
        for knudsen in np.linspace(kn_table[0], kn_table[-1], 50):
            self.assertGreater(float(self._cd_at(knudsen)), 0.0)


class TestAerodynamicTableAngleOfAttack(unittest.TestCase):
    """
    迎角依存の空力係数表（6 自由度計算用）。

    従来の「AOA 0 だけ」のファイルもそのまま読めること、
    迎角ブロックが複数あるときに (AOA, Kn) の 2 次元内挿になることを確かめる。
    """

    def test_table_without_an_aoa_line_is_read_as_zero(self):
        rows = [[1.0, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0],
                [10.0, 1.3, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0]]
        aerodynamic_dict = read_aerodynamic_table([(None, rows)])

        np.testing.assert_allclose(aerodynamic_dict[satellite.KEY_AOA], [0.0])
        np.testing.assert_allclose(aerodynamic_dict[satellite.KEY_CD_MEAN], [1.2, 1.3])

    def test_single_block_is_independent_of_the_angle_of_attack(self):
        rows = [[1.0, 1.2, 0.0, -0.1, 0.0, -0.02, 0.0, 100.0],
                [10.0, 1.3, 0.0, -0.1, 0.0, -0.02, 0.0, 200.0]]
        aerodynamic_dict = read_aerodynamic_table([(0.0, rows)])

        for angle in (0.0, 30.0, 90.0):
            force, moment = satellite.get_aerodynamic_coefficient_attitude(1.0, angle, aerodynamic_dict)
            np.testing.assert_allclose(force, [1.2, 0.0, -0.1], atol=1.e-12)
            np.testing.assert_allclose(moment, [0.0, -0.02, 0.0], atol=1.e-12)

    def test_two_dimensional_interpolation_reproduces_the_nodes(self):
        block = [(0.0,  [[1.0, 1.2, 0.0, 0.0, 0.0,  0.00, 0.0, 100.0],
                         [10.0, 1.3, 0.0, 0.0, 0.0,  0.00, 0.0, 200.0]]),
                 (20.0, [[1.0, 1.1, 0.0, 0.4, 0.0, -0.10, 0.0, 100.0],
                         [10.0, 1.15, 0.0, 0.5, 0.0, -0.12, 0.0, 200.0]])]
        aerodynamic_dict = read_aerodynamic_table(block)

        np.testing.assert_allclose(aerodynamic_dict[satellite.KEY_AOA], [0.0, 20.0])

        force, moment = satellite.get_aerodynamic_coefficient_attitude(1.0, 0.0, aerodynamic_dict)
        np.testing.assert_allclose(force, [1.2, 0.0, 0.0], atol=1.e-12)
        force, moment = satellite.get_aerodynamic_coefficient_attitude(10.0, 20.0, aerodynamic_dict)
        np.testing.assert_allclose(force, [1.15, 0.0, 0.5], atol=1.e-12)
        np.testing.assert_allclose(moment, [0.0, -0.12, 0.0], atol=1.e-12)

    def test_interpolates_linearly_between_the_blocks(self):
        block = [(0.0,  [[1.0, 1.2, 0.0, 0.0, 0.0,  0.00, 0.0, 100.0],
                         [10.0, 1.2, 0.0, 0.0, 0.0,  0.00, 0.0, 200.0]]),
                 (20.0, [[1.0, 1.0, 0.0, 0.4, 0.0, -0.10, 0.0, 100.0],
                         [10.0, 1.0, 0.0, 0.4, 0.0, -0.10, 0.0, 200.0]])]
        aerodynamic_dict = read_aerodynamic_table(block)

        force, moment = satellite.get_aerodynamic_coefficient_attitude(1.0, 10.0, aerodynamic_dict)
        np.testing.assert_allclose(force, [1.1, 0.0, 0.2], atol=1.e-12)
        np.testing.assert_allclose(moment, [0.0, -0.05, 0.0], atol=1.e-12)

    def test_values_are_clamped_outside_the_table(self):
        block = [(0.0,  [[1.0, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0],
                         [10.0, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0]]),
                 (20.0, [[1.0, 1.0, 0.0, 0.4, 0.0, -0.1, 0.0, 100.0],
                         [10.0, 1.0, 0.0, 0.4, 0.0, -0.1, 0.0, 200.0]])]
        aerodynamic_dict = read_aerodynamic_table(block)

        force_low, moment_low = satellite.get_aerodynamic_coefficient_attitude(0.01, -30.0, aerodynamic_dict)
        np.testing.assert_allclose(force_low, [1.2, 0.0, 0.0], atol=1.e-12)
        force_high, moment_high = satellite.get_aerodynamic_coefficient_attitude(1.e4, 90.0, aerodynamic_dict)
        np.testing.assert_allclose(force_high, [1.0, 0.0, 0.4], atol=1.e-12)

    def test_inconsistent_knudsen_numbers_are_rejected(self):
        block = [(0.0,  [[1.0, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0],
                         [10.0, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0]]),
                 (20.0, [[2.0, 1.0, 0.0, 0.4, 0.0, -0.1, 0.0, 100.0],
                         [10.0, 1.0, 0.0, 0.4, 0.0, -0.1, 0.0, 200.0]])]
        with self.assertRaises(SystemExit):
            read_aerodynamic_table(block)


class TestSphereConeTable(unittest.TestCase):
    """同梱の球円錐テーブル（6 自由度チュートリアル用）の性質。"""

    @classmethod
    def setUpClass(cls):
        config = load_config()
        # リポジトリ直下の database/ を見る（チュートリアルの作業ディレクトリには
        # そのケースが使うテーブルしか置いていないため）
        config['satellite']['directory_path_specify'] = 'default'
        config['satellite']['filename_aerodynamic'] = 'aerodynamic_spherecone_aoa.txt'
        with quiet():
            cls.aero = satellite.initial_settings_satellite(config)

    def test_covers_the_whole_range_of_the_angle_of_attack(self):
        aoa = self.aero[satellite.KEY_AOA]
        self.assertAlmostEqual(aoa[0], 0.0)
        self.assertAlmostEqual(aoa[-1], 180.0)

    def test_is_statically_stable_about_zero_incidence(self):
        # 復元モーメント: 迎角が正なら機首下げ（Cm < 0）
        for knudsen in (1.e-3, 1.0, 1.e4):
            for angle in (5.0, 10.0, 20.0, 40.0):
                force, moment = satellite.get_aerodynamic_coefficient_attitude(knudsen, angle, self.aero)
                self.assertLess(moment[1], 0.0)

    def test_is_symmetric_at_zero_incidence(self):
        for knudsen in (1.e-3, 1.0, 1.e4):
            force, moment = satellite.get_aerodynamic_coefficient_attitude(knudsen, 0.0, self.aero)
            self.assertGreater(force[0], 0.0)          # 軸力（= 迎角 0 での CD）は正
            self.assertAlmostEqual(force[1], 0.0, places=9)
            self.assertAlmostEqual(force[2], 0.0, places=9)
            np.testing.assert_allclose(moment, np.zeros(3), atol=1.e-9)

    def test_drag_is_larger_in_free_molecular_flow(self):
        force_continuum, moment_continuum = satellite.get_aerodynamic_coefficient_attitude(1.e-4, 0.0, self.aero)
        force_molecular, moment_molecular = satellite.get_aerodynamic_coefficient_attitude(1.e5, 0.0, self.aero)
        self.assertGreater(force_molecular[0], force_continuum[0])


class TestAxisymmetryOfTheAerodynamicTable(unittest.TestCase):
    """
    係数表の軸対称性の検査。

    6 自由度では表を全迎角だけで引き、attitude.matrix_aerodynamic_roll で実際の
    横流れ面へ回す。これが厳密なのは軸対称の機体だけで、そのとき表は全迎角・全 Kn で
    CFy = CMx = CMz = 0 になる。非軸対称の表を入れると、その 3 成分が面内の量として
    誤った向きへ回されるが、表に迎角以外の姿勢変数が無いので正しい向きは復元できない。
    黙って誤った答えを出さないよう、読み込み時に警告することにした（v2.3.1）。
    """

    def synthetic_table(self, cfy=0.0, cmx=0.0, cmz=0.0, cfx=1.2, cfz=0.4, cmy=-0.1):
        # 迎角 0 は軸対称（横流れが無い）、迎角 20 deg 側に非対称成分を入れる
        return [(0.0,  [[1.0,  cfx, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0],
                        [10.0, cfx, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0]]),
                (20.0, [[1.0,  cfx, cfy, cfz, cmx, cmy, cmz, 100.0],
                        [10.0, cfx, cfy, cfz, cmx, cmy, cmz, 200.0]])]

    def check(self, **keyword):
        aerodynamic_dict = read_aerodynamic_table(self.synthetic_table(**keyword))
        return satellite.check_axisymmetry(aerodynamic_dict)

    def names(self, violation):
        return [item['name'] for item in violation]

    def test_an_axisymmetric_table_passes(self):
        self.assertEqual(self.check(), [])

    def test_the_shipped_spherecone_table_is_axisymmetric(self):
        """6 自由度チュートリアルの表は解析モデルなので厳密に軸対称（残る 1e-17 は丸め誤差）。"""
        config = load_config()
        config['satellite']['directory_path_specify'] = 'default'
        config['satellite']['filename_aerodynamic'] = 'aerodynamic_spherecone_aoa.txt'
        with quiet():
            aerodynamic_dict = satellite.initial_settings_satellite(config)

        self.assertEqual(satellite.check_axisymmetry(aerodynamic_dict), [])
        self.assertGreater(np.abs(aerodynamic_dict[satellite.KEY_CM][:,:,1]).max(), 0.0)

    def test_the_shipped_egg_table_is_flagged_for_its_moments(self):
        """
        3 自由度用の EGG の表（DSMC+CFD）は計測のばらつきで CMx と CMz が残っている。

        3 自由度では CFx しか使わないので無害だが、6 自由度に持ち込むと
        ロール・ヨーが立つ。CFy は CFx の 4e-4 しかないので閾値には掛からない。
        """
        with quiet():
            aerodynamic_dict = satellite.initial_settings_satellite(load_config())

        self.assertEqual(self.names(satellite.check_axisymmetry(aerodynamic_dict)), ['CMx', 'CMz'])

    def test_a_side_force_is_detected(self):
        violation = self.check(cfy=0.05)

        self.assertEqual(self.names(violation), ['CFy'])
        self.assertAlmostEqual(violation[0]['magnitude'], 0.05)
        self.assertAlmostEqual(violation[0]['ratio'], 0.05/1.2)
        self.assertAlmostEqual(violation[0]['angle_of_attack'], 20.0)

    def test_a_rolling_and_a_yawing_moment_are_detected(self):
        self.assertEqual(self.names(self.check(cmx=0.01)), ['CMx'])
        self.assertEqual(self.names(self.check(cmz=0.01)), ['CMz'])
        self.assertEqual(self.names(self.check(cfy=0.05, cmx=0.01, cmz=0.01)), ['CFy', 'CMx', 'CMz'])

    def test_the_threshold_is_relative_to_the_in_plane_coefficients(self):
        tolerance = satellite.TOLERANCE_AXISYMMETRY
        self.assertEqual(self.check(cfy=0.5*tolerance*1.2), [])
        self.assertEqual(self.names(self.check(cfy=2.0*tolerance*1.2)), ['CFy'])
        # モーメントは CMy に対して見るので、CMy を大きくすれば同じ CMx が埋もれる
        self.assertEqual(self.names(self.check(cmx=1.e-3, cmy=-0.1)), ['CMx'])
        self.assertEqual(self.check(cmx=1.e-3, cmy=-10.0), [])

    def test_round_off_dust_is_ignored(self):
        """解析モデルの表に残る 1e-17 級の値で警告を出さない（面内が 0 のときも）。"""
        self.assertEqual(self.check(cfy=1.e-17, cmx=1.e-17, cmz=1.e-17), [])
        self.assertEqual(self.check(cfy=1.e-17, cmx=1.e-17, cmz=1.e-17, cmy=0.0), [])

    def test_a_component_is_detected_when_the_in_plane_moment_vanishes(self):
        """CMy が 0 の表でも、ロール・ヨーが立っていれば軸対称ではない。"""
        violation = self.check(cmx=1.e-6, cmy=0.0)

        self.assertEqual(self.names(violation), ['CMx'])
        self.assertEqual(violation[0]['ratio'], float('inf'))

    def test_the_warning_names_the_file_and_the_components(self):
        aerodynamic_dict = read_aerodynamic_table(self.synthetic_table(cfy=0.05, cmz=0.01))
        with quiet() as output:
            flag_warned = satellite.warn_if_not_axisymmetric(aerodynamic_dict)
        message = output.getvalue()

        self.assertTrue(flag_warned)
        self.assertIn('Warning', message)
        self.assertIn(TABLE_FILENAME, message)
        # 逸脱した成分だけが報告されること（説明文にも成分名が出るので行で見る）
        reported = [line.split()[0].lstrip('-') for line in message.splitlines() if 'reaches' in line]
        self.assertEqual(reported, ['CFy', 'CMz'])

    def test_the_warning_appears_only_when_the_attitude_is_solved(self):
        """3 自由度では CFy も CM も使わないので、警告は出さない。"""
        with tempfile.TemporaryDirectory() as directory:
            write_aerodynamic_table(self.synthetic_table(cfy=0.05), directory)

            config = aerodynamic_config(directory)
            with quiet() as output:
                satellite.initial_settings_satellite(config)
            self.assertNotIn('Warning', output.getvalue())

            config = aerodynamic_config(directory)
            config['attitude'] = {'flag_attitude': True}
            with quiet() as output:
                aerodynamic_dict = satellite.initial_settings_satellite(config)
            self.assertIn('Warning', output.getvalue())

            # 空力モデルが constant なら表そのものを使わないので黙る
            config = aerodynamic_config(directory)
            config['attitude'] = {'flag_attitude': True}
            config['satellite']['kind_aerodynamic_model'] = 'constant'
            with quiet() as output:
                satellite.initial_settings_satellite(config)
            self.assertNotIn('Warning', output.getvalue())

        # 警告であって停止ではない: 表はそのまま使える
        self.assertIn(satellite.KEY_INTERP, aerodynamic_dict)
        force, moment = satellite.get_aerodynamic_coefficient_attitude(1.0, 20.0, aerodynamic_dict)
        self.assertAlmostEqual(force[1], 0.05)


class TestFileFormats(unittest.TestCase):
    """
    大気モデルファイルの形式判定（CCMC / NRLMSISE-00 Fortran 版）。

    Fortran 版のファイルはリポジトリに含めていないので、
    NRLMSISE-00_readctl.FOR が書く形式をここで組み立てて読ませる。
    """

    HEADER_FORTRAN = [
        'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
        ' IDAY           0',
        ' UT   0.00000000    ',
        ' F107A   76.5999985    ',
        ' F107   76.5999985    ',
        ' APH   8.39999962       7.00000000    ',
        ' 1, ALTITUDE (KM)',
        ' 2, D(2) - O NUMBER DENSITY(CM-3)',
        ' 3, D(3) - N2 NUMBER DENSITY(CM-3)',
        ' 4, D(4) - O2 NUMBER DENSITY(CM-3)',
        ' 5, D(6) - TOTAL MASS DENSITY(GM/CM3)',
        ' 6, T(2) - TEMPERATURE AT ALT',
        ' 7, D(8) - N NUMBER DENSITY(CM-3)',
        '           1           2           3           4           5           6           7',
    ]

    # 高度[km], O, N2, O2, 質量密度[g/cm3], 温度[K], N
    ROWS = [
        (0.0,   0.0,      2.04e+19, 5.47e+18, 1.256e-03, 279.0, 0.0),
        (100.0, 4.30e+11, 9.49e+12, 2.20e+12, 5.60e-10,  190.0, 1.00e+07),
        (200.0, 3.20e+15, 2.60e+14, 1.00e+13, 1.87e-13,  742.0, 5.00e+10),
        (300.0, 3.60e+14, 8.90e+12, 2.10e+11, 9.20e-15,  780.0, 1.20e+10),
        (400.0, 5.50e+13, 6.10e+11, 8.90e+09, 8.95e-16,  783.0, 3.00e+09),
    ]

    def _write_fortran_file(self, directory, extra_lines=()):
        lines = list(self.HEADER_FORTRAN)
        for row in self.ROWS:
            lines.append('   %.6f       %.8E   %.8E   %.8E   %.8E   %.6f       %.8E    ' % row)
        lines.extend(extra_lines)

        path = os.path.join(directory, 'atmospheremodel.txt')
        with open(path, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        return path

    def _config_for(self, directory):
        config = load_config()
        config['atmosphere']['directory_path_specify'] = 'manual'
        config['atmosphere']['directory_atmosphere'] = directory
        config['atmosphere']['filename_atmosphere'] = 'atmospheremodel.txt'
        return config

    def test_fortran_format_is_detected_and_parsed(self):
        with tempfile.TemporaryDirectory() as directory:
            self._write_fortran_file(directory)
            with quiet():
                atm = atmosphere.read_atmosphere_file(self._config_for(directory))

        self.assertEqual(atm[atmosphere.KEY_DATA], len(self.ROWS))
        self.assertEqual(atm[atmosphere.KEY_ATM], 7)

        # 高度は km のまま、数密度は cm-3 -> m-3、質量密度は g/cm3 -> kg/m3
        np.testing.assert_allclose(atm[atmosphere.KEY_Height], [r[0] for r in self.ROWS])
        np.testing.assert_allclose(atm[atmosphere.KEY_O], [r[1] * 1.0e6 for r in self.ROWS])
        np.testing.assert_allclose(atm[atmosphere.KEY_Mass_density], [r[4] * 1.0e3 for r in self.ROWS])
        np.testing.assert_allclose(atm[atmosphere.KEY_Temperature_neutral], [r[5] for r in self.ROWS])
        np.testing.assert_allclose(atm[atmosphere.KEY_N], [r[6] * 1.0e6 for r in self.ROWS])

    def test_both_formats_produce_the_same_keys(self):
        with quiet():
            ccmc = atmosphere.read_atmosphere_file(load_config())
        with tempfile.TemporaryDirectory() as directory:
            self._write_fortran_file(directory)
            with quiet():
                fortran = atmosphere.read_atmosphere_file(self._config_for(directory))

        self.assertEqual(set(ccmc.keys()), set(fortran.keys()))

    def test_the_whole_pipeline_works_on_a_fortran_file(self):
        """Kn と補間器とスケールハイトまで通ること。"""
        with tempfile.TemporaryDirectory() as directory:
            self._write_fortran_file(directory)
            with quiet():
                atm = atmosphere.initial_settings_atmosphere(self._config_for(directory))

        self.assertGreater(atm[atmosphere.KEY_KN].min(), 0.0)
        self.assertIn(atmosphere.KEY_INTERP, atm)
        self.assertGreater(atm[atmosphere.KEY_SCALE_HEIGHT], 0.0)

    def test_trailing_junk_lines_are_ignored(self):
        """末尾に空行やゴミ行があっても落ちないこと（レビュー D-15）。"""
        with tempfile.TemporaryDirectory() as directory:
            self._write_fortran_file(directory, extra_lines=('', '   ', 'end of file', ''))
            with quiet():
                atm = atmosphere.read_atmosphere_file(self._config_for(directory))

        self.assertEqual(atm[atmosphere.KEY_DATA], len(self.ROWS))

    def test_unrecognized_format_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'atmospheremodel.txt')
            with open(path, 'w') as f:
                f.write('this file is not an atmosphere table\n1.0 2.0 3.0\n')
            with self.assertRaises(SystemExit):
                with quiet():
                    atmosphere.read_atmosphere_file(self._config_for(directory))

    def test_unknown_parameter_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'atmospheremodel.txt')
            lines = list(self.HEADER_FORTRAN)
            lines[14] = ' 2, D(2) - HE NUMBER DENSITY(CM-3)'
            with open(path, 'w') as f:
                f.write('\n'.join(lines) + '\n')
            with self.assertRaises(SystemExit):
                with quiet():
                    atmosphere.read_atmosphere_file(self._config_for(directory))



if __name__ == '__main__':
    unittest.main()
