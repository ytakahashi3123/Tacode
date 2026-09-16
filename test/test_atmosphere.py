#!/usr/bin/env python3
"""大気・空力テーブルの読み込みと内挿のテスト。"""

import glob
import os
import re
import tempfile
import unittest

import numpy as np
import scipy.interpolate

from context import ROOT_DIR, SRC_DIR, load_config, quiet

import atmosphere.atmosphere as atmosphere
import satellite.satellite as satellite


TABLE_FILENAME = 'aerodynamic_test.txt'
COLUMN_SPARE = '\t'.join(['0.0']*6)


def write_aerodynamic_table(block, directory, length_reference=None):
    """
    合成の空力係数表を書く。

    block は (迎角, 行の並び) の並びで、行は
    [Kn, CFx, CFy, CFz, CMx, CMy, CMz, Altitude]。迎角に None を渡すと
    "AOA" 行を書かない（従来形式のファイル）。標準偏差の 6 列は 0 で埋める。

    length_reference を渡すと "# Reference length: ..." の行を足す（文字列で
    渡せばそのまま書くので、壊れた行も作れる）。
    """
    path = os.path.join(directory, TABLE_FILENAME)
    with open(path, 'w') as f:
        f.write('test table\n')
        if length_reference is not None:
            if isinstance(length_reference, str):
                f.write('# {}{}\n'.format(satellite.MARKER_LENGTH_REFERENCE, length_reference))
            else:
                f.write('# {}: {:g} m\n'.format(satellite.MARKER_LENGTH_REFERENCE, length_reference))
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


def read_aerodynamic_table_with_length(block, length_reference, length_config):
    """代表長さを書いた合成テーブルを読み込み、(辞書, 出力) を返す。"""
    with tempfile.TemporaryDirectory() as directory:
        write_aerodynamic_table(block, directory, length_reference)
        config = aerodynamic_config(directory)
        config['satellite']['characteristic_length'] = length_config
        with quiet() as buffer:
            aerodynamic_dict = satellite.initial_settings_satellite(config)
        return aerodynamic_dict, buffer.getvalue()


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
        補間器は 3 つの量をまとめた**ベクトル値スプライン 1 本**で、
        成分の並びは KEYS_INTERPOLATED（密度・温度・Kn）。
        """
        self.assertIn(atmosphere.KEY_INTERP, self.atm)
        self.assertEqual(atmosphere.KEYS_INTERPOLATED,
                         [atmosphere.KEY_Mass_density,
                          atmosphere.KEY_Temperature_neutral,
                          atmosphere.KEY_KN])

        interpolator = self.atm[atmosphere.KEY_INTERP]
        altitude = self.atm[atmosphere.KEY_Height]
        value = interpolator(0.5*(altitude[0] + altitude[-1]))
        self.assertEqual(value.shape, (len(atmosphere.KEYS_INTERPOLATED),))

        # 1 本にまとめても、量ごとに作ったスプラインとビット一致すること
        # （まとめたのは速度のため。RK4 の各段で 3 回評価するのをやめられる）
        import scipy.interpolate
        for index, key in enumerate(atmosphere.KEYS_INTERPOLATED):
            spline = scipy.interpolate.make_interp_spline(altitude, self.atm[key], k=3)
            for fraction in (0.13, 0.5, 0.77):
                query = altitude[0] + fraction*(altitude[-1] - altitude[0])
                self.assertEqual(float(interpolator(query)[index]), float(spline(query)))

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

    チュートリアルの既定テーブルは 700 km まであるので、周回軌道（遠地点 522 km）は
    外挿区間に入らない。それでも外挿そのものは 700 km より上で効くので、ここでは
    **400 km で終わる `atmospheremodel.txt` を明示的に読んで**外挿を試験する。
    """

    @classmethod
    def setUpClass(cls):
        with quiet():
            cls.atm = atmosphere.initial_settings_atmosphere(cls.config_short())

        config_clamp = cls.config_short()
        config_clamp['atmosphere']['kind_extrapolation'] = atmosphere.KIND_EXTRAPOLATION_CLAMP
        with quiet():
            cls.atm_clamp = atmosphere.initial_settings_atmosphere(config_clamp)

    @staticmethod
    def config_short():
        """400 km で終わるテーブルを読む config。"""
        config = load_config(os.path.join(ROOT_DIR, 'tutorial', 'work_reentry', 'config.yml'))
        return config

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

    def test_the_exponent_is_capped_so_that_the_knudsen_number_stays_finite(self):
        """
        上端から遠く離れても Kn が Inf にならないこと。

        指数に上限を置かないと exp がアンダーフローして密度が 0 になり、
        Kn = 端の値/0 が Inf になる。発散した計算の出力に Inf が並び、
        0 除算の RuntimeWarning も出ていた（レビュー A-10）。
        """
        height = self.atm[atmosphere.KEY_Height]
        scale_height = self.atm[atmosphere.KEY_SCALE_HEIGHT]

        for delta in (scale_height*atmosphere.EXPONENT_MAXIMUM,
                      1.0e5, 1.0e7, 6.6e9):
            density, temperature, knudsen = self._at(height[-1] + delta, self.atm)
            self.assertTrue(np.isfinite(float(density)))
            self.assertTrue(np.isfinite(float(knudsen)))
            self.assertGreater(float(knudsen), 0.0)

    def test_the_cap_does_not_touch_the_extrapolation_in_range(self):
        """上限より内側では従来どおり exp(-dh/H) のままであること。"""
        height = self.atm[atmosphere.KEY_Height]
        scale_height = self.atm[atmosphere.KEY_SCALE_HEIGHT]
        density_top = self.atm[atmosphere.KEY_Mass_density][-1]

        delta = scale_height*(atmosphere.EXPONENT_MAXIMUM - 1.0)
        density = float(self._at(height[-1] + delta, self.atm)[0])
        self.assertAlmostEqual(density/(density_top*np.exp(-delta/scale_height)), 1.0, places=9)

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
        # 3 自由度の CD は np.interp で引くので補間器を持たない。
        # 辞書に入るのは迎角依存の表を読んだときの (AOA, Kn) の補間器だけ
        self.assertIn(satellite.KEY_INTERP, self.aero)
        self.assertNotIn(satellite.KEY_CD_MEAN, self.aero[satellite.KEY_INTERP])

    def _cd_at(self, knudsen):
        return satellite.get_aerodynamic_coefficient(
            knudsen,
            self.aero[satellite.KEY_KN],
            self.aero[satellite.KEY_CD_MEAN])

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

    def test_the_shipped_apollo_table_is_axisymmetric_and_statically_stable(self):
        """
        検証ケースのアポロ表（修正ニュートン流）。回転体なので厳密に軸対称で、
        軸上の基準点まわりでは迎角を増やすと頭上げ（Cm > 0）になる。
        機首下げのモーメントは重心を軸から外して初めて立つ（config の
        attitude.center_of_gravity）ので、この表そのものは静不安定に見えてよい。
        """
        config = load_config()
        config['satellite']['directory_path_specify'] = 'default'
        config['satellite']['filename_aerodynamic'] = 'aerodynamic_apollo_aoa.txt'
        with quiet():
            aerodynamic_dict = satellite.initial_settings_satellite(config)

        self.assertEqual(satellite.check_axisymmetry(aerodynamic_dict), [])

        # 迎角 0 では横流れが無いので Cm = 0、迎角が付くと軸上基準点まわりの Cm が立つ
        force_zero, moment_zero = satellite.get_aerodynamic_coefficient_attitude(1.e-4, 0.0, aerodynamic_dict)
        force_trim, moment_trim = satellite.get_aerodynamic_coefficient_attitude(1.e-4, 24.4, aerodynamic_dict)
        self.assertAlmostEqual(float(moment_zero[1]), 0.0, places=12)
        self.assertLess(float(moment_trim[1]), 0.0)

        # 連続流の軸力は公表されているトリム時の CD に近い（1.2891, NASA TN D-6725）
        angle = 24.4*np.pi/180.0
        drag = float(force_trim[0])*np.cos(angle) + float(force_trim[2])*np.sin(angle)
        self.assertAlmostEqual(drag, 1.2891, delta=0.05)

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


class TestConstantAerodynamicModel(unittest.TestCase):
    """
    kind_aerodynamic_model: constant は係数表を使わない。

    solver は config の drag_coefficient を読むので、表を読みに行くと、
    使いもしないファイルが無いだけで計算が止まる。
    """

    def test_the_table_is_not_read(self):
        config = load_config()
        config['satellite']['kind_aerodynamic_model'] = 'constant'
        config['satellite']['filename_aerodynamic'] = 'no_such_table.txt'
        with quiet():
            aerodynamic_dict = satellite.initial_settings_satellite(config)
        self.assertEqual(aerodynamic_dict, {})

    def test_the_table_is_still_read_for_fileread(self):
        config = load_config()
        self.assertEqual(config['satellite']['kind_aerodynamic_model'], 'fileread')
        with quiet():
            aerodynamic_dict = satellite.initial_settings_satellite(config)
        self.assertIn(satellite.KEY_INTERP, aerodynamic_dict)


class TestGeneratedAtmosphereTable(unittest.TestCase):
    """
    database/atmosphere/generate_atmosphere_table.py が書いた表を読めること。

    生成器は CCMC(VITMO) 形式で書く。読み取り側とずれると検証ケースが黙って
    別の大気で走るので、同梱の生成物（Apollo 4 の検証ケースの表）で往復を見る。
    """

    TABLE_DIRECTORY = os.path.join(ROOT_DIR, 'validation', 'apollo4', 'database', 'atmosphere')
    TABLE_FILENAME = 'atmospheremodel_apollo4.txt'

    @classmethod
    def setUpClass(cls):
        config = load_config()
        config['atmosphere']['directory_path_specify'] = 'manual'
        config['atmosphere']['directory_atmosphere'] = cls.TABLE_DIRECTORY
        config['atmosphere']['filename_atmosphere'] = cls.TABLE_FILENAME
        with quiet():
            cls.atmosphere_dict = atmosphere.initial_settings_atmosphere(config)

    def test_the_format_is_detected_as_ccmc(self):
        path = os.path.join(self.TABLE_DIRECTORY, self.TABLE_FILENAME)
        with open(path) as f:
            lines = [line.strip() for line in f.readlines()]
        self.assertEqual(atmosphere.detect_atmosphere_file_kind(lines, path),
                         atmosphere.KIND_FILE_CCMC)

    def test_it_covers_the_entry_altitudes(self):
        altitude = self.atmosphere_dict[atmosphere.KEY_Height]
        self.assertLessEqual(altitude[0], 0.0)
        self.assertGreaterEqual(altitude[-1], 123.5)

    def test_the_density_decreases_with_altitude(self):
        density = self.atmosphere_dict[atmosphere.KEY_Mass_density]
        self.assertTrue(np.all(np.diff(density) < 0.0))

    def test_every_property_is_finite_and_positive(self):
        for key in (atmosphere.KEY_Mass_density, atmosphere.KEY_Temperature_neutral,
                    atmosphere.KEY_KN):
            value = self.atmosphere_dict[key]
            self.assertTrue(np.all(np.isfinite(value)), key)
            self.assertTrue(np.all(value > 0.0), key)



class TestTheTableCoversTheTutorial(unittest.TestCase):
    """
    チュートリアルが自分の大気テーブルの範囲の中を飛ぶこと（レビュー B-9）。

    周回軌道のケースは遠地点 522 km で、400 km で終わるテーブルを読んでいたため
    飛行時間の 45 % を外挿で飛んでいた。外挿は止まりはしないが、実測の 700 km
    テーブルと比べると遠地点の密度を 31 % 小さく見積もる。テーブルを 700 km まで
    伸ばして直したので、**戻されたらここで落ちる**。
    """

    CASES = ('work', 'work_reentry', 'work_reentry_6dof')

    def test_every_case_stays_inside_its_table(self):
        for name in self.CASES:
            directory = os.path.join(ROOT_DIR, 'tutorial', name)
            path_result = os.path.join(directory, 'output_result', 'tecplot.dat')
            if not os.path.exists(path_result):
                continue

            with quiet():
                atmosphere_dict = atmosphere.initial_settings_atmosphere(
                    load_config(os.path.join(directory, 'config.yml')))
            top = atmosphere_dict[atmosphere.KEY_Height][-1]

            altitude = []
            with open(path_result) as stream:
                for line in stream:
                    item = line.split()
                    if len(item) > 7 and item[0][0].isdigit():
                        altitude.append(float(item[6]))
            self.assertGreater(len(altitude), 0, name)
            self.assertLessEqual(
                max(altitude), top,
                '%s reaches %.1f km but its table stops at %.1f km'
                % (name, max(altitude), top))


class TestBilinearMatchesScipy(unittest.TestCase):
    """
    手書きの双一次補間が RegularGridInterpolator とビット一致すること。

    RGI をやめたのは速度のため（1 点 22 us -> 5 us、6 自由度計算の 1 割）。
    速度のためだけの置き換えなので、検査することは「値が動かないこと」だけになる。
    足し込む順序を RGI の _evaluate_linear に合わせてあるので一致する
    （順序を変えると 1.5e-16 ずれ、チュートリアルと validation の出力が
    バイト一致しなくなる ── それが高速化の合格条件だった）。

    節点そのものも見る（格子の内側では区間の選び方によらず同じ値になるはずで、
    そこがずれるのは重みの作り方が違うとき）。
    """

    def setUp(self):
        self.rng = np.random.default_rng(20260916)

    def build(self, num_first, num_second, num_component):
        grid_first = np.sort(self.rng.uniform(-180.0, 180.0, num_first))
        grid_second = np.sort(np.exp(self.rng.uniform(-9.0, 9.0, num_second)))
        values = self.rng.normal(size=(num_first, num_second, num_component))
        interpolator = scipy.interpolate.RegularGridInterpolator(
            (grid_first, grid_second), values,
            method='linear', bounds_error=False, fill_value=None)
        return grid_first, grid_second, values, interpolator

    def test_it_matches_on_a_synthetic_grid(self):
        grid_first, grid_second, values, interpolator = self.build(37, 10, 6)

        query = [(float(a), float(k))
                 for a in self.rng.uniform(grid_first[0], grid_first[-1], 400)
                 for k in self.rng.uniform(grid_second[0], grid_second[-1], 5)]
        query += [(float(a), float(k)) for a in grid_first for k in grid_second]

        for first, second in query:
            reference = interpolator(np.array([[first, second]]))[0]
            np.testing.assert_array_equal(
                satellite.evaluate_bilinear(grid_first, grid_second, values, first, second),
                reference)

    def test_it_matches_on_the_shipped_tables(self):
        for case in ('work_reentry_6dof',):
            config = load_config(os.path.join(ROOT_DIR, 'tutorial', case, 'config.yml'))
            with quiet():
                aerodynamic_dict = satellite.initial_settings_satellite(config)

            aoa_table = aerodynamic_dict[satellite.KEY_AOA]
            knudsen_table = aerodynamic_dict[satellite.KEY_KN]
            values = aerodynamic_dict[satellite.KEY_INTERP][satellite.KEY_COEF]
            self.assertGreater(len(aoa_table), 1, case)

            interpolator = scipy.interpolate.RegularGridInterpolator(
                (aoa_table, knudsen_table), values,
                method='linear', bounds_error=False, fill_value=None)

            query = [(float(a), float(k))
                     for a in self.rng.uniform(aoa_table[0], aoa_table[-1], 200)
                     for k in self.rng.uniform(knudsen_table[0], knudsen_table[-1], 5)]
            query += [(float(a), float(k)) for a in aoa_table for k in knudsen_table]

            for first, second in query:
                np.testing.assert_array_equal(
                    satellite.evaluate_bilinear(aoa_table, knudsen_table, values, first, second),
                    interpolator(np.array([[first, second]]))[0])


# 呼び出しとしての interp1d（コメントや説明文の中は見ない）
PATTERN_INTERP1D = re.compile(r'\binterp1d\s*\(')


ROWS_AERODYNAMIC = [[1.0e-2, 1.10, 0.0, 0.0, 0.0, 0.0, 0.0, 80.0],
                    [1.0e+0, 1.20, 0.0, 0.0, 0.0, 0.0, 0.0, 110.0],
                    [1.0e+2, 1.30, 0.0, 0.0, 0.0, 0.0, 0.0, 160.0]]


class TestTheReferenceLengthOfTheAerodynamicTable(unittest.TestCase):
    """
    係数表が「どの代表長さで Kn 軸を作ったか」を宣言でき、食い違いが知らされること
    （CODE_REVIEW B-7）。

    Kn = lambda/L なので、表の Kn 軸は表を作ったときの L に紐づいている。別の L で
    引くと同じ高度に対して表の別の場所を読むことになり、CD が黙ってずれる。
    表が代表長さを書いていれば読み込み時に突き合わせられる。

    **停止はしない。**手元の表を別の機体に当てるのは近似と割り切れば実際に行うことで、
    同梱の再突入チュートリアルがまさにそれをしている（`aerodynamic.txt` は 0.8 m で
    作られているが、ケースの characteristic_length は 0.5 m）。軸対称でない表を
    使うときと同じく、知らせるだけにしてある。
    """

    def test_a_matching_length_says_nothing(self):
        aerodynamic_dict, output = read_aerodynamic_table_with_length(
            [(0.0, ROWS_AERODYNAMIC)], 0.8, 0.8)
        self.assertEqual(aerodynamic_dict[satellite.KEY_LENGTH_REF], 0.8)
        self.assertNotIn('Caution', output)

    def test_a_different_length_is_reported_with_both_values(self):
        aerodynamic_dict, output = read_aerodynamic_table_with_length(
            [(0.0, ROWS_AERODYNAMIC)], 0.8, 0.5)
        self.assertEqual(aerodynamic_dict[satellite.KEY_LENGTH_REF], 0.8)
        self.assertIn('Caution', output)
        self.assertIn('0.8', output)
        self.assertIn('0.5', output)
        # 表そのものは読めているので計算は続く
        np.testing.assert_allclose(aerodynamic_dict[satellite.KEY_CD_MEAN],
                                   [row[1] for row in ROWS_AERODYNAMIC])

    def test_a_table_without_the_line_is_not_checked(self):
        # 出所の分からない古いファイルが読めなくなっては困る
        aerodynamic_dict, output = read_aerodynamic_table_with_length(
            [(0.0, ROWS_AERODYNAMIC)], None, 0.5)
        self.assertIsNone(aerodynamic_dict[satellite.KEY_LENGTH_REF])
        self.assertNotIn('Caution', output)

    def test_a_line_without_a_number_stops(self):
        # 目印はあるのに値が読めない行は、書いたつもりの検査が黙って消えるので止める
        for broken in (': m', '', ': abc m'):
            with self.subTest(line=broken):
                with self.assertRaises(SystemExit) as raised:
                    read_aerodynamic_table_with_length(
                        [(0.0, ROWS_AERODYNAMIC)], broken, 0.5)
                self.assertEqual(raised.exception.code, 1)

    def test_every_shipped_table_declares_a_reference_length(self):
        for path in sorted(glob.glob(os.path.join(ROOT_DIR, 'database', 'aerodynamic', '*.txt'))):
            with self.subTest(table=os.path.basename(path)):
                with quiet():
                    _, _, length_reference = satellite.parse_aerodynamic_file(path)
                self.assertIsNotNone(length_reference)
                self.assertGreater(length_reference, 0.0)

    def test_the_length_of_the_egg_table_follows_from_its_altitude_column(self):
        """
        `aerodynamic.txt` の 0.8 m が、その表自身の Altitude 列から出てくること。

        この表は外から来た DSMC + CFD のデータで、代表長さがどこにも書かれていな
        かった。各行は Kn とその Kn に対応する高度を持つので、大気テーブルから
        「代表長さ 1 m での Kn」を引けば L = Kn(L=1)/Kn_table として復元できる。
        2026-09-16 にそうやって求めたところ、**8 行すべてで 0.8000 m** だった。
        表に書いた宣言はその値で、ここはその根拠をテストとして残しておくもの。

        解析モデルで作った表（球円錐・アポロ）の Altitude 列は Kn に対する線形内挿で
        求めた参考値でしかないので、この復元は効かない（生成器が代表長さを書くので、
        そちらは作りで保証されている）。`aerodynamic_fire2.txt` は Altitude 列が 0。
        """
        config = load_config()
        config['satellite']['characteristic_length'] = 1.0
        config['atmosphere']['directory_path_specify'] = 'default'
        config['atmosphere']['filename_atmosphere'] = 'atmospheremodel.txt'
        with quiet():
            atm = atmosphere.initial_settings_atmosphere(config)

        path = os.path.join(ROOT_DIR, 'database', 'aerodynamic', 'aerodynamic.txt')
        with quiet():
            _, block, length_reference = satellite.parse_aerodynamic_file(path)

        recovered = np.interp(block[0][:, 13], atm[atmosphere.KEY_Height],
                              atm[atmosphere.KEY_KN])/block[0][:, 0]
        np.testing.assert_allclose(recovered, length_reference, rtol=1.e-4)

    def test_the_cases_that_differ_from_their_table_are_the_known_ones(self):
        """
        同梱の config のうち、表の代表長さと食い違っているものを固定する。

        食い違っているのは**再突入のチュートリアル 4 つ**だけで、いずれも EGG の表
        （0.8 m で作られている）を 0.5 m の機体に当てている。CD の差は 100 km より上で
        最大 0.73 %、それより下では表の下端にクランプされて同じ値になるので、
        チュートリアルの軌道はこの食い違いで動かない（2026-09-16 実測）。

        新しいケースを足したときに黙って仲間が増えないよう、ここで並びを固定する。
        減らす（ケースを直す）ときは出力が動くので、参照出力の更新とセットになる。
        """
        expected = {'tutorial/template_wind', 'tutorial/template_wind_table',
                    'tutorial/work_montecarlo_wind',
                    'tutorial/work_reentry', 'tutorial/work_reentry_wind_table'}

        different = set()
        for path in sorted(glob.glob(os.path.join(ROOT_DIR, 'tutorial', '*', 'config.yml'))
                           + glob.glob(os.path.join(ROOT_DIR, 'validation', '*', 'config*.yml'))):
            config = load_config(path)
            section = config['satellite']
            if section.get('kind_aerodynamic_model') != 'fileread':
                continue
            directory = satellite.get_database_directory(section, 'satellite', 'aerodynamic',
                                                         'directory_aerodynamic')
            with quiet():
                _, _, length_reference = satellite.parse_aerodynamic_file(
                    os.path.join(directory, section['filename_aerodynamic']))
            if length_reference is None:
                continue
            if abs(length_reference - section['characteristic_length']) \
               > satellite.TOLERANCE_LENGTH_REFERENCE*length_reference:
                different.add(os.path.relpath(os.path.dirname(path), ROOT_DIR).replace(os.sep, '/'))

        self.assertEqual(different, expected)


class TestTheKnudsenNumberNeedsMolecularDensities(unittest.TestCase):
    """
    分子種の数密度を 1 つも持たないテーブルで止まること（CODE_REVIEW B-10）。

    Kn は sum( n d^2 ) から作る。その和を取るループは「その種が無い」を読み飛ばす
    ので、**4 種とも無いテーブルでは和が 0 のままになり、Kn が全高度で Inf になる**。
    Inf は係数表の上端にクランプされるので、警告も出ないまま自由分子流の係数で
    計算が進んでしまう。
    """

    HEADER = ['a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
              ' IDAY           0',
              ' UT   0.00000000    ',
              ' F107A   76.5999985    ',
              ' F107   76.5999985    ',
              ' APH   8.39999962       7.00000000    ',
              ' 1, ALTITUDE (KM)',
              ' 2, D(6) - TOTAL MASS DENSITY(GM/CM3)',
              ' 3, T(2) - TEMPERATURE AT ALT',
              '           1           2           3']

    ROWS = [(0.0, 1.256e-03, 279.0),
            (200.0, 1.870e-13, 742.0),
            (400.0, 8.950e-16, 783.0)]

    def test_a_table_without_any_species_stops(self):
        with tempfile.TemporaryDirectory() as directory:
            lines = list(self.HEADER)
            for row in self.ROWS:
                lines.append('   %.6f       %.8E   %.6f       ' % row)
            with open(os.path.join(directory, 'atmospheremodel.txt'), 'w') as f:
                f.write('\n'.join(lines) + '\n')

            config = load_config()
            config['atmosphere']['directory_path_specify'] = 'manual'
            config['atmosphere']['directory_atmosphere'] = directory
            config['atmosphere']['filename_atmosphere'] = 'atmospheremodel.txt'

            with quiet() as buffer:
                with self.assertRaises(SystemExit) as raised:
                    atmosphere.initial_settings_atmosphere(config)

        self.assertEqual(raised.exception.code, 1)
        self.assertIn('no molecular number density', buffer.getvalue())

    def test_the_shipped_tables_carry_all_four_species(self):
        config = load_config()
        for name in ('atmospheremodel.txt', 'atmospheremodel_2015_700km.txt',
                     'atmospheremodel_700km.txt'):
            with self.subTest(table=name):
                config['atmosphere']['directory_path_specify'] = 'default'
                config['atmosphere']['filename_atmosphere'] = name
                with quiet() as buffer:
                    atm = atmosphere.initial_settings_atmosphere(config)
                for key in atmosphere.LIST_MOLECULAR_KIND:
                    self.assertIn(key, atm, name)
                self.assertTrue(np.all(np.isfinite(atm[atmosphere.KEY_KN])), name)
                self.assertIn('Molecular species used', buffer.getvalue())


class TestTheSourceDoesNotUseLegacySciPy(unittest.TestCase):
    """
    scipy.interpolate.interp1d は SciPy 1.10 以降 legacy で、いつ消えてもおかしくない。

    大気（3 次スプライン）は make_interp_spline(k=3) に、空力の CD（1 次元線形）は
    np.interp に移した。どちらも interp1d の値とビット単位で一致することを確かめて
    入れ替えてある（CubicSpline は同じ not-a-knot でも評価が PPoly になり 4e-16 ずれる
    ので採らなかった）。呼び出しが戻ると将来の SciPy で落ちるので、ここで止める。
    """

    def test_no_module_calls_interp1d(self):
        offending = []
        for path in sorted(glob.glob(os.path.join(SRC_DIR, '**', '*.py'), recursive=True)):
            with open(path) as stream:
                for number, line in enumerate(stream, start=1):
                    if line.strip().startswith('#'):
                        continue
                    if PATTERN_INTERP1D.search(line):
                        offending.append('%s:%d: %s' % (os.path.relpath(path, SRC_DIR),
                                                        number, line.strip()))
        self.assertEqual(
            offending, [],
            'src/ が legacy の interp1d を呼んでいる:\n' + '\n'.join(offending))


if __name__ == '__main__':
    unittest.main()
