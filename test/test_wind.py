#!/usr/bin/env python3
"""
風のテスト。

風は「大気は地球と剛体的に共回転する」という従来の仮定を緩めるもので、
空力だけが対気速度 v_air = v_ecef - v_wind を使う。ここで確かめるのは 3 点:

  1. 既定（flag_wind: False）では計算に一切触れない。同梱の config もすべて無効
  2. 向きの規約。[東, 北, 上] の地心ローカル系で、初期速度とまったく同じ経路を通る
  3. **恒等式**: 風 w のもとで速度 v の空力は、無風で速度 v - w の空力に等しい。
     風が空力にしか効いていないこと（重力・コリオリ力・遠心力・角速度は不変）も同時に見る

符号を裏返す、風を足してしまう、対地速度と対気速度を取り違える、といった間違い方は
どれも 3 の恒等式で落ちる。
"""

import copy
import glob
import os
import unittest

import numpy as np
import scipy.interpolate

from context import ROOT_DIR, load_config, quiet

import atmosphere.atmosphere as atmosphere
import attitude.attitude as attitude
import coordinate_system.coordinate_system as coordinate_system
import force_term.force_term as force_term
import satellite.satellite as satellite
import solver.solver as solver
import wind.wind as wind
from orbital.orbital import orbital


CONFIG_6DOF = os.path.join(ROOT_DIR, 'tutorial', 'work_reentry_6dof', 'config.yml')
CONFIG_REENTRY = os.path.join(ROOT_DIR, 'tutorial', 'work_reentry', 'config.yml')

# 風を意図して有効にしている config。風のチュートリアル（モンテカルロとその
# テンプレート、実データのテーブルを読むケース）だけで、いずれも参照出力を
# 持たない（回帰テストの対象外）
CONFIG_WIND_ON = [os.path.join(ROOT_DIR, 'tutorial', 'work_montecarlo_wind', 'config.yml'),
                  os.path.join(ROOT_DIR, 'tutorial', 'work_reentry_wind_table', 'config.yml'),
                  os.path.join(ROOT_DIR, 'tutorial', 'template_wind', 'config.yml'),
                  os.path.join(ROOT_DIR, 'tutorial', 'template_wind_table', 'config.yml')]


def read_tecplot(path):
    """Tecplot の point 形式を列名と数値配列にして返す。"""
    variables = None
    rows = []
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.lower().startswith('variables'):
                variables = [name.strip() for name in stripped.split('=', 1)[1].split(',')]
                continue
            if not stripped or stripped.startswith('#') or stripped.lower().startswith('zone'):
                continue
            rows.append([float(value) for value in stripped.split()])
    return variables, np.array(rows)


def write_wind_table(path, rows, epoch=None, comment=None):
    """テーブルを 1 本書く。rows は (経度, 緯度, 高度, 東, 北, 上) の並び。"""
    with open(path, 'w') as f:
        f.write('# Tacode wind table\n')
        if comment is not None:
            f.write('# ' + comment + '\n')
        if epoch is not None:
            f.write('# Epoch (UTC): ' + epoch + '\n')
        if len(rows) > 0 and len(rows[0]) == 7:
            f.write('# Time[s] Longitude[deg.] Latitude[deg.] Altitude[km] East[m/s] North[m/s] Up[m/s]\n')
        else:
            f.write('# Longitude[deg.] Latitude[deg.] Altitude[km] East[m/s] North[m/s] Up[m/s]\n')
        for row in rows:
            f.write('  ' + ' '.join(['{:g}'.format(value) for value in row]) + '\n')


def grid_rows(longitude, latitude, altitude, function):
    """格子を全点埋めた行の並びを作る。function(lon, lat, alt) が [東, 北, 上] を返す。"""
    rows = []
    for value_lon in longitude:
        for value_lat in latitude:
            for value_alt in altitude:
                east, north, up = function(value_lon, value_lat, value_alt)
                rows.append([value_lon, value_lat, value_alt, east, north, up])
    return rows


def grid_rows_time(time, longitude, latitude, altitude, function):
    """時刻軸つきの行（7 列）。function(t, lon, lat, alt) が [東, 北, 上] を返す。"""
    rows = []
    for value_time in time:
        for value_lon in longitude:
            for value_lat in latitude:
                for value_alt in altitude:
                    east, north, up = function(value_time, value_lon, value_lat, value_alt)
                    rows.append([value_time, value_lon, value_lat, value_alt, east, north, up])
    return rows


def epoch_of(string_datetime):
    import epoch.epoch as epoch_module

    with quiet():
        return epoch_module.initial_settings_epoch(
            {'epoch': {'flag_epoch': True, 'datetime': string_datetime}})


def table_config(directory, filename, extrapolation='zero'):
    return {'wind': {'flag_wind': True,
                     'kind_wind_model': 'fileread',
                     'directory_path_specify': 'manual',
                     'directory_wind': directory,
                     'filename_wind': filename,
                     'kind_extrapolation': extrapolation}}


def make_table_wind(rows, extrapolation='zero', epoch=None, epoch_dict=None, comment=None):
    """テーブルを一時ディレクトリに書いて読み込ませ、wind_dict と出力を返す。"""
    import tempfile

    directory = tempfile.mkdtemp()
    path = os.path.join(directory, 'windmodel.txt')
    write_wind_table(path, rows, epoch=epoch, comment=comment)
    config = table_config(directory, 'windmodel.txt', extrapolation)
    with quiet() as buffer:
        wind_dict = wind.initial_settings_wind(config, epoch_dict)
    return wind_dict, buffer.getvalue(), directory


def local(wind_dict, longitude, latitude, altitude, time_elapsed=0.0):
    """度・km で引いて [東, 北, 上] を返す。"""
    return np.array(wind.get_wind_local([np.deg2rad(longitude), np.deg2rad(latitude),
                                         altitude*1000.0], wind_dict, time_elapsed))


def wind_config(config, velocity, kind='constant'):
    """config に風を入れて返す（元は壊さない）。"""
    config = copy.deepcopy(config)
    config['wind'] = {'flag_wind': True, 'kind_wind_model': kind, 'velocity': list(velocity)}
    return config


def make_wind(velocity, kind='constant'):
    with quiet():
        return wind.initial_settings_wind(wind_config({}, velocity, kind))


def run_solver(config, time_max, timestep, wind_velocity=None):
    """初期化からソルバーまでを 3 自由度で回す。"""
    config = copy.deepcopy(config)
    config['computational_setup']['time_elapsed_maximum'] = time_max
    config['time_integration']['timestep_constant'] = timestep
    if wind_velocity is not None:
        config = wind_config(config, wind_velocity)

    with quiet():
        atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
        aerodynamic_dict = satellite.initial_settings_satellite(config)
        wind_dict = wind.initial_settings_wind(config)

        orb = orbital()
        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            orb.initial_settings(config)

        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            solver.solve_equation_motion(config, iteration, time_elapsed,
                                         coordinate_dict, velocity_dict, trajectory_dict,
                                         atmosphere_dict, aerodynamic_dict, None, wind_dict)

    return iteration, coordinate_dict, velocity_dict


class TestTheWindIsOffByDefault(unittest.TestCase):

    def test_a_config_without_the_section_gives_none(self):
        self.assertIsNone(wind.initial_settings_wind({}))

    def test_an_explicit_false_gives_none(self):
        config = {'wind': {'flag_wind': False, 'kind_wind_model': 'constant', 'velocity': [10.0, 0.0, 0.0]}}
        self.assertIsNone(wind.initial_settings_wind(config))

    def test_the_shipped_configs_leave_it_off(self):
        # 参照出力を壊さないための条件。ここが True になった時点で回帰テストが落ちる。
        # 風のモンテカルロのチュートリアル（CONFIG_WIND_ON）だけは意図して有効で、
        # そちらは参照出力を持たない
        paths = sorted(glob.glob(os.path.join(ROOT_DIR, 'tutorial', '*', 'config.yml')))
        paths.append(os.path.join(ROOT_DIR, 'tutorial', 'template', 'config.yml'))
        paths.append(os.path.join(ROOT_DIR, 'src', 'config.yml'))
        paths = [path for path in paths if path not in CONFIG_WIND_ON]
        self.assertTrue(len(paths) >= 6)
        for path in paths:
            config = load_config(path)
            self.assertIn('wind', config, path)
            self.assertFalse(config['wind']['flag_wind'], path)
            self.assertIsNone(wind.initial_settings_wind(config), path)

    def test_only_the_wind_tutorials_turn_it_on(self):
        # 逆向きの検査。除外リストが実在の config を指していること、そこでは
        # 風が実際に組み立てられることを見る（config を消して除外だけが残ると気付けない）
        for path in CONFIG_WIND_ON:
            self.assertTrue(os.path.exists(path), path)
            config = load_config(path)
            self.assertTrue(config['wind']['flag_wind'], path)
            with quiet():
                wind_dict = wind.initial_settings_wind(config)
            self.assertIsNotNone(wind_dict, path)

            if config['wind']['kind_wind_model'] == wind.MODEL_CONSTANT:
                # ばらつかせる基準値になるので、水平成分がゼロでないこと
                self.assertGreater(abs(config['wind']['velocity'][0])
                                   + abs(config['wind']['velocity'][1]), 0.0, path)
            else:
                # テーブルを読むケース。同梱のテーブルが実在し、格子が組み上がること
                self.assertEqual(config['wind']['kind_wind_model'], wind.MODEL_FILEREAD, path)
                self.assertIn(wind.KEY_INTERP, wind_dict, path)
                self.assertGreater(len(wind_dict[wind.KEY_ALTITUDE]), 1, path)

    def test_the_relative_velocity_is_the_velocity_itself_when_off(self):
        # 引き算を通さないので、風を切ったときの結果はビット単位で従来と同じになる
        velocity = [1234.5, -678.9, 0.0]
        result = wind.get_relative_velocity({}, [1.0, 2.0, 3.0], [0.0, 0.0, 0.0], velocity, None)
        self.assertIs(result, velocity)


class TestTheConfigurationIsChecked(unittest.TestCase):

    def test_an_unknown_model_stops_the_run(self):
        with self.assertRaises(SystemExit):
            make_wind([0.0, 0.0, 0.0], kind='hwm14')

    def test_a_velocity_of_the_wrong_length_stops_the_run(self):
        with self.assertRaises(SystemExit):
            make_wind([10.0, 0.0])

    def test_a_missing_velocity_is_taken_as_calm(self):
        with quiet():
            wind_dict = wind.initial_settings_wind({'wind': {'flag_wind': True}})
        np.testing.assert_array_equal(wind_dict[wind.KEY_VELOCITY], np.zeros(3))


class TestTheFrameOfTheWind(unittest.TestCase):
    """
    [東, 北, 上] の地心ローカル系。初期速度と同じ規約で、同じ変換を通す。
    """

    def setUp(self):
        self.config = load_config()

    def test_east_at_the_prime_meridian_on_the_equator_is_plus_y(self):
        wind_dict = make_wind([10.0, 0.0, 0.0])
        result = wind.get_wind_velocity(self.config, np.array([6.4e6, 0.0, 0.0]), [0.0, 0.0, 0.0], wind_dict)
        np.testing.assert_allclose(result, [0.0, 10.0, 0.0], atol=1.e-9)

    def test_north_at_the_prime_meridian_on_the_equator_is_plus_z(self):
        wind_dict = make_wind([0.0, 10.0, 0.0])
        result = wind.get_wind_velocity(self.config, np.array([6.4e6, 0.0, 0.0]), [0.0, 0.0, 0.0], wind_dict)
        np.testing.assert_allclose(result, [0.0, 0.0, 10.0], atol=1.e-9)

    def test_up_at_the_prime_meridian_on_the_equator_is_plus_x(self):
        wind_dict = make_wind([0.0, 0.0, 10.0])
        result = wind.get_wind_velocity(self.config, np.array([6.4e6, 0.0, 0.0]), [0.0, 0.0, 0.0], wind_dict)
        np.testing.assert_allclose(result, [10.0, 0.0, 0.0], atol=1.e-9)

    def test_east_at_ninety_degrees_of_longitude_is_minus_x(self):
        wind_dict = make_wind([10.0, 0.0, 0.0])
        result = wind.get_wind_velocity(self.config, np.array([0.0, 6.4e6, 0.0]), [0.0, 0.0, 0.0], wind_dict)
        np.testing.assert_allclose(result, [-10.0, 0.0, 0.0], atol=1.e-9)

    def test_the_magnitude_survives_the_rotation(self):
        wind_dict = make_wind([13.0, -7.0, 2.0])
        coordinate = np.array([2.1e6, -3.4e6, 4.9e6])
        result = wind.get_wind_velocity(self.config, coordinate, [0.0, 0.0, 0.0], wind_dict)
        self.assertAlmostEqual(np.linalg.norm(result), np.linalg.norm([13.0, -7.0, 2.0]), places=9)

    def test_it_agrees_with_the_frame_the_initial_velocity_uses(self):
        # 風と初期速度が同じ規約であることを、初期速度の変換経路そのものと突き合わせる
        wind_dict = make_wind([11.0, -5.0, 3.0])
        coordinate = np.array([2.1e6, -3.4e6, 4.9e6])
        polar = coordinate_system.set_angle_polar(self.config, coordinate)
        reference = coordinate_system.convert_polar_carteasian(self.config, [11.0, -5.0, 3.0], polar[2], polar[1])
        np.testing.assert_allclose(wind.get_wind_velocity(self.config, coordinate, [0.0, 0.0, 0.0], wind_dict),
                                   reference, rtol=0.0, atol=0.0)

    def test_the_constant_model_does_not_depend_on_the_position(self):
        wind_dict = make_wind([11.0, -5.0, 3.0])
        for geodetic in ([0.0, 0.0, 0.0], [1.0, -0.5, 1.e5], [3.0, 1.2, 4.e5]):
            np.testing.assert_array_equal(wind.get_wind_local(geodetic, wind_dict), [11.0, -5.0, 3.0])


class TestTheVelocityFactor(unittest.TestCase):
    """
    wind.velocity_factor は風の場に一様に掛かる係数（既定 1.0）。

    これはテーブルの風をモンテカルロで振るための入口である。driver は制御ファイルの
    **数値**を書き換える仕組みで、テーブルはファイル名で選ぶので、テーブルそのものは
    振れない。initial_settings.density_factor と同じ流儀（1 要素のリスト）で書く。

    既定の 1.0 では厳密に恒等（IEEE 754）なので、係数を書かない config と
    ビット単位で同じ結果になる。ここはその 2 点を押さえる。
    """

    def setUp(self):
        self.config = load_config()
        self.coordinate = np.array([2.1e6, -3.4e6, 4.9e6])
        self.directories = []

    def tearDown(self):
        import shutil

        for directory in self.directories:
            shutil.rmtree(directory, ignore_errors=True)

    def build_wind(self, velocity, factor=None):
        config = wind_config({}, velocity)
        if factor is not None:
            config['wind'][wind.KEY_FACTOR] = factor
        with quiet():
            return wind.initial_settings_wind(config)

    def build_table(self, rows, factor=None):
        import tempfile

        directory = tempfile.mkdtemp()
        self.directories.append(directory)
        write_wind_table(os.path.join(directory, 'windmodel.txt'), rows)
        config = table_config(directory, 'windmodel.txt')
        if factor is not None:
            config['wind'][wind.KEY_FACTOR] = factor
        with quiet():
            return wind.initial_settings_wind(config)

    def test_the_default_is_one(self):
        wind_dict = self.build_wind([20.0, 10.0, 0.0])
        self.assertEqual(wind_dict[wind.KEY_FACTOR], 1.0)
        # 掛け算を通しても値は 1 ビットも動かない
        np.testing.assert_array_equal(wind.get_wind_local([0.1, 0.2, 1.e4], wind_dict),
                                      [20.0, 10.0, 0.0])

    def test_it_scales_the_constant_wind(self):
        wind_dict = self.build_wind([20.0, 10.0, -4.0], [0.5])
        np.testing.assert_allclose(wind.get_wind_local([0.1, 0.2, 1.e4], wind_dict),
                                   [10.0, 5.0, -2.0], rtol=0.0, atol=0.0)

    def test_a_bare_number_is_accepted_as_well(self):
        # config には 1 要素のリストで書くが、スカラーで与えても同じに読む
        wind_dict = self.build_wind([20.0, 10.0, 0.0], 2.0)
        np.testing.assert_allclose(wind.get_wind_local([0.0, 0.0, 0.0], wind_dict),
                                   [40.0, 20.0, 0.0], rtol=0.0, atol=0.0)

    def test_it_scales_the_table(self):
        rows = grid_rows([-10.0, 10.0], [0.0, 20.0], [0.0, 20.0],
                         lambda lon, lat, alt: [lon + alt, lat - alt, 1.0])
        wind_plain = self.build_table(rows)
        wind_scaled = self.build_table(rows, [3.0])
        for point in ([-10.0, 0.0, 0.0], [0.0, 10.0, 10.0], [5.0, 20.0, 20.0]):
            np.testing.assert_allclose(local(wind_scaled, *point),
                                       3.0*np.array(local(wind_plain, *point)), rtol=1.e-12)

    def test_it_scales_the_wind_and_nothing_else(self):
        # 係数 f の風 w は、係数なしの風 f*w と同じ。対気速度まで通して見る
        velocity = np.array([1200.0, 3400.0, -560.0])
        wind_scaled = self.build_wind([20.0, -10.0, 1.0], [0.25])
        wind_plain = self.build_wind([5.0, -2.5, 0.25])
        np.testing.assert_allclose(
            wind.get_relative_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0],
                                       velocity, wind_scaled),
            wind.get_relative_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0],
                                       velocity, wind_plain), rtol=0.0, atol=0.0)

    def test_zero_brings_back_the_co_rotating_atmosphere(self):
        # 係数 0 は「風なし」。対気速度が ECEF 速度に戻る
        velocity = np.array([1200.0, 3400.0, -560.0])
        wind_dict = self.build_wind([60.0, -30.0, 5.0], [0.0])
        np.testing.assert_allclose(
            wind.get_relative_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0],
                                       velocity, wind_dict), velocity, rtol=0.0, atol=0.0)

    def test_a_factor_of_the_wrong_length_stops_the_run(self):
        with self.assertRaises(SystemExit):
            self.build_wind([20.0, 10.0, 0.0], [1.0, 2.0])

    def test_a_factor_which_is_not_a_number_stops_the_run(self):
        with self.assertRaises(SystemExit):
            self.build_wind([20.0, 10.0, 0.0], ['strong'])

    def test_the_shipped_configs_leave_it_at_one(self):
        # 1.0 以外がまぎれ込むと、風を有効にしたケースの結果が黙って変わる
        paths = sorted(glob.glob(os.path.join(ROOT_DIR, 'tutorial', '*', 'config.yml')))
        paths.append(os.path.join(ROOT_DIR, 'tutorial', 'template', 'config.yml'))
        paths.append(os.path.join(ROOT_DIR, 'tutorial', 'template_wind', 'config.yml'))
        paths.append(os.path.join(ROOT_DIR, 'src', 'config.yml'))
        self.assertTrue(len(paths) >= 7)
        for path in paths:
            config = load_config(path)
            self.assertIn('wind', config, path)
            self.assertEqual(wind.get_velocity_factor(config['wind']), 1.0, path)


class TestTheRelativeVelocity(unittest.TestCase):

    def setUp(self):
        self.config = load_config()
        self.coordinate = np.array([2.1e6, -3.4e6, 4.9e6])
        self.velocity = np.array([1200.0, 3400.0, -560.0])

    def test_the_wind_is_subtracted(self):
        wind_dict = make_wind([20.0, -10.0, 1.0])
        velocity_wind = wind.get_wind_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0], wind_dict)
        result = wind.get_relative_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0], self.velocity, wind_dict)
        np.testing.assert_allclose(result, self.velocity - velocity_wind, rtol=0.0, atol=0.0)

    def test_a_calm_wind_leaves_the_velocity_unchanged(self):
        wind_dict = make_wind([0.0, 0.0, 0.0])
        result = wind.get_relative_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0], self.velocity, wind_dict)
        np.testing.assert_allclose(result, self.velocity, rtol=0.0, atol=0.0)

    def test_a_wind_equal_to_the_velocity_leaves_the_body_at_rest_in_the_air(self):
        # 空力がちょうど消える状況。符号を取り違えると 2 倍になって現れる
        wind_dict = make_wind([100.0, 0.0, 0.0])
        velocity = wind.get_wind_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0], wind_dict)
        result = wind.get_relative_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0], velocity, wind_dict)
        np.testing.assert_allclose(result, np.zeros(3), atol=1.e-9)


class TestTheWindEntersOnlyTheAerodynamics(unittest.TestCase):
    """
    恒等式: 風 w のもとで速度 v の空力は、無風で速度 v - w の空力に等しい。
    重力・コリオリ力・遠心力はそのいずれでも変わらない。
    """

    def setUp(self):
        self.config = load_config()
        self.coordinate = np.array([2.1e6, -3.4e6, 4.9e6])
        self.velocity = np.array([1200.0, 3400.0, -560.0])
        self.wind_dict = make_wind([40.0, -25.0, 3.0])
        self.velocity_air = wind.get_relative_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0],
                                                       self.velocity, self.wind_dict)

    def acceleration(self, velocity, velocity_air=None):
        acceleration = force_term.acceleration_initialsettings(self.config)
        return np.array(force_term.acceleration_routine(self.config, self.coordinate, velocity,
                                                        10.0, 0.785, 1.3, 1.0, 1.e-5,
                                                        acceleration, None, velocity_air))

    def test_the_drag_with_wind_equals_the_drag_at_the_relative_velocity(self):
        with_wind = self.acceleration(self.velocity, self.velocity_air)
        without = self.acceleration(self.velocity_air)
        np.testing.assert_allclose(with_wind[4], without[4], rtol=0.0, atol=0.0)

    def test_gravity_and_the_rotation_terms_do_not_see_the_wind(self):
        with_wind = self.acceleration(self.velocity, self.velocity_air)
        without = self.acceleration(self.velocity)
        for row, name in ((1, 'gravity'), (2, 'Coriolis'), (3, 'centrifugal')):
            np.testing.assert_allclose(with_wind[row], without[row], rtol=0.0, atol=0.0,
                                       err_msg=name)

    def test_the_coriolis_term_keeps_using_the_ground_velocity(self):
        # コリオリ力に対気速度を使ってしまう取り違えを捕まえる
        with_wind = self.acceleration(self.velocity, self.velocity_air)
        wrong = self.acceleration(self.velocity_air)
        self.assertFalse(np.allclose(with_wind[2], wrong[2]))

    def test_the_drag_opposes_the_relative_velocity(self):
        with_wind = self.acceleration(self.velocity, self.velocity_air)
        direction = with_wind[4]/np.linalg.norm(with_wind[4])
        np.testing.assert_allclose(direction, -self.velocity_air/np.linalg.norm(self.velocity_air), atol=1.e-12)

    def test_the_total_is_still_the_sum_of_the_terms(self):
        with_wind = self.acceleration(self.velocity, self.velocity_air)
        np.testing.assert_allclose(with_wind[0], with_wind[1:5].sum(axis=0), rtol=1.e-14, atol=0.0)


class TestTheWindInTheSixDegreeOfFreedomState(unittest.TestCase):
    """
    同じ恒等式を、姿勢に依存する空力（solver.get_aerodynamic_state）で確かめる。
    """

    def setUp(self):
        self.config = load_config(CONFIG_6DOF)
        with quiet():
            self.aerodynamic_dict = satellite.initial_settings_satellite(self.config)
        import moment_term.moment_term as moment_term
        with quiet():
            self.moment, self.attitude_property = moment_term.moment_initialsettings(self.config)
        self.coordinate = np.array([2.1e6, -3.4e6, 4.9e6])
        self.velocity = np.array([1200.0, 3400.0, -560.0])
        self.quaternion = attitude.quaternion_normalize(np.array([0.83, 0.21, -0.44, 0.27]))
        self.omega = np.array([0.01, -0.02, 0.03])
        self.wind_dict = make_wind([60.0, -30.0, 5.0])
        self.velocity_air = wind.get_relative_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0],
                                                       self.velocity, self.wind_dict)

    def state(self, velocity, velocity_air=None):
        return solver.get_aerodynamic_state(self.config, self.attitude_property, self.aerodynamic_dict,
                                            'fileread', self.coordinate, velocity,
                                            self.quaternion, self.omega,
                                            10.0, 0.785, 1.0,
                                            1.e-5, 1.0, 0.1, 1.3, -0.5,
                                            self.config['planet']['rotation_rate'],
                                            np.zeros(4*3).reshape(4, 3), velocity_air)

    def test_the_force_and_the_moment_match_the_relative_velocity(self):
        force_wind, moment_wind, omega_wind = self.state(self.velocity, self.velocity_air)
        force_ref, moment_ref, omega_ref = self.state(self.velocity_air)
        np.testing.assert_allclose(force_wind, force_ref, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(moment_wind, moment_ref, rtol=0.0, atol=0.0)

    def test_the_wind_actually_changes_the_force(self):
        force_wind, moment_wind, _ = self.state(self.velocity, self.velocity_air)
        force_calm, moment_calm, _ = self.state(self.velocity)
        self.assertFalse(np.allclose(force_wind, force_calm))
        self.assertFalse(np.allclose(moment_wind, moment_calm))

    def test_the_angular_velocity_is_not_touched_by_the_wind(self):
        # 風は並進の量。omega_relative は地球自転だけを引いたもののままでなければならない
        _, _, omega_wind = self.state(self.velocity, self.velocity_air)
        _, _, omega_calm = self.state(self.velocity)
        np.testing.assert_allclose(omega_wind, omega_calm, rtol=0.0, atol=0.0)


class TestTheSolverWithWind(unittest.TestCase):

    def test_a_calm_wind_reproduces_the_run_without_wind(self):
        # 引き算が入っても値が動かないことの確認（風を有効にしたまま風速 0）
        config = load_config(CONFIG_REENTRY)
        iteration_a, coordinate_a, velocity_a = run_solver(config, 300.0, 0.1)
        iteration_b, coordinate_b, velocity_b = run_solver(config, 300.0, 0.1, wind_velocity=[0.0, 0.0, 0.0])
        self.assertEqual(iteration_a, iteration_b)
        np.testing.assert_array_equal(np.array(coordinate_a['cartesian']), np.array(coordinate_b['cartesian']))
        np.testing.assert_array_equal(np.array(velocity_a['cartesian']), np.array(velocity_b['cartesian']))

    def test_a_tailwind_and_a_headwind_move_the_trajectory_the_opposite_ways(self):
        # 初速は東向き 7450 m/s。東風は下流へ、西風は上流へ動かす
        config = load_config(CONFIG_REENTRY)
        _, coordinate_calm, _ = run_solver(config, 300.0, 0.1, wind_velocity=[0.0, 0.0, 0.0])
        _, coordinate_east, _ = run_solver(config, 300.0, 0.1, wind_velocity=[60.0, 0.0, 0.0])
        _, coordinate_west, _ = run_solver(config, 300.0, 0.1, wind_velocity=[-60.0, 0.0, 0.0])
        longitude_calm = coordinate_calm['geodetic'][-1][0]
        longitude_east = coordinate_east['geodetic'][-1][0]
        longitude_west = coordinate_west['geodetic'][-1][0]
        self.assertGreater(longitude_east, longitude_calm)
        self.assertLess(longitude_west, longitude_calm)
        # 追い風と向かい風の効きはほぼ対称になる。厳密に対称でないのは抗力が
        # 速度の 2 乗に比例するため（追い風は対気速度を下げ、向かい風は上げる）
        ratio = (longitude_east - longitude_calm)/(longitude_calm - longitude_west)
        self.assertAlmostEqual(ratio, 1.0, delta=0.05)

    def test_a_crosswind_pushes_the_trajectory_sideways(self):
        config = load_config(CONFIG_REENTRY)
        _, coordinate_calm, _ = run_solver(config, 300.0, 0.1, wind_velocity=[0.0, 0.0, 0.0])
        _, coordinate_north, _ = run_solver(config, 300.0, 0.1, wind_velocity=[0.0, 60.0, 0.0])
        self.assertGreater(coordinate_north['geodetic'][-1][1], coordinate_calm['geodetic'][-1][1])

    def test_the_body_settles_at_the_same_speed_relative_to_the_air(self):
        #
        # 終端速度は対気の釣り合いで決まるので、水平な一様風のもとでも
        # 対気速度の大きさは無風のときの対地速度に一致する。
        # 対地速度そのものは風の分だけ増える。
        #
        config = load_config(CONFIG_REENTRY)
        _, coordinate_calm, velocity_calm = run_solver(config, 2700.0, 1.0, wind_velocity=[0.0, 0.0, 0.0])
        _, coordinate_wind, velocity_wind = run_solver(config, 2700.0, 1.0, wind_velocity=[20.0, 0.0, 0.0])

        with quiet():
            wind_dict = make_wind([20.0, 0.0, 0.0])
        velocity_air = wind.get_relative_velocity(config, coordinate_wind['cartesian'][-1],
                                                  coordinate_wind['geodetic'][-1],
                                                  velocity_wind['cartesian'][-1], wind_dict)

        speed_calm = np.linalg.norm(velocity_calm['cartesian'][-1])
        speed_air = np.linalg.norm(velocity_air)
        speed_ground = np.linalg.norm(velocity_wind['cartesian'][-1])

        # 終端では対気速度が無風の対地速度と一致する（残差は高度差と密度差の分）
        self.assertLess(abs(speed_air - speed_calm)/speed_calm, 0.02)
        # 対地速度は風の分だけ大きい。ここが一致してしまうなら風が効いていない
        self.assertAlmostEqual(speed_ground, np.sqrt(speed_air**2 + 20.0**2), delta=0.5)


class TestTheAttitudeOutputUsesTheRelativeVelocity(unittest.TestCase):
    """
    迎角・横滑り角は保存せず出力時に作り直すので、ソルバーが使ったのと同じ
    対気速度から作らなければ出力とソルバーが食い違う。
    """

    def setUp(self):
        self.config = load_config(CONFIG_6DOF)
        with quiet():
            self.orb = orbital()
        self.coordinate = np.array([2.1e6, -3.4e6, 4.9e6])
        self.velocity = np.array([120.0, 340.0, -56.0])
        self.quaternion = attitude.quaternion_normalize(np.array([0.83, 0.21, -0.44, 0.27]))
        self.omega = np.array([0.01, -0.02, 0.03])
        self.rotation_rate = self.config['planet']['rotation_rate']
        wind_dict = make_wind([60.0, -30.0, 5.0])
        self.velocity_air = wind.get_relative_velocity(self.config, self.coordinate, [0.0, 0.0, 0.0],
                                                       self.velocity, wind_dict)

    def output(self, velocity, velocity_air=None):
        return [float(value) for value in
                self.orb.get_attitude_output(self.config, self.coordinate, velocity,
                                             self.quaternion, self.omega, self.rotation_rate,
                                             velocity_air).split()]

    def test_the_aerodynamic_angles_match_the_relative_velocity(self):
        with_wind = self.output(self.velocity, self.velocity_air)
        reference = self.output(self.velocity_air)
        # 末尾 3 つが迎角・横滑り角・全迎角
        np.testing.assert_allclose(with_wind[-3:], reference[-3:], rtol=0.0, atol=0.0)

    def test_the_wind_actually_changes_the_aerodynamic_angles(self):
        with_wind = self.output(self.velocity, self.velocity_air)
        calm = self.output(self.velocity)
        self.assertFalse(np.allclose(with_wind[-3:], calm[-3:]))

    def test_the_quaternion_and_the_euler_angles_do_not_see_the_wind(self):
        # 姿勢そのものは風で変わらない。空力角だけが対気速度に依る
        with_wind = self.output(self.velocity, self.velocity_air)
        calm = self.output(self.velocity)
        np.testing.assert_allclose(with_wind[:10], calm[:10], rtol=0.0, atol=0.0)


class TestTheWindInTheOutput(unittest.TestCase):
    """
    出力の風の列。風は保存せず出力時に引き直しているので、ここが独立の検査になる。
    Upl/Vpl/Wpl と WindE/WindN/WindU はどちらも地心ローカル系 [東, 北, 上] なので、
    VelairAbs をその差から再計算して突き合わせられる。
    """

    @classmethod
    def setUpClass(cls):
        import tempfile

        config = copy.deepcopy(load_config(CONFIG_REENTRY))
        config['computational_setup']['time_elapsed_maximum'] = 60.0
        config['time_integration']['timestep_constant'] = 0.5
        config = wind_config(config, [45.0, -30.0, 2.0])

        cls.directory = tempfile.mkdtemp()
        config['post_process']['directory_output'] = cls.directory
        config['post_process']['kml']['flag_output'] = False
        config['post_process']['tecplot']['flag_output'] = True
        config['post_process']['tecplot']['frequency_output'] = 10

        with quiet():
            atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
            aerodynamic_dict = satellite.initial_settings_satellite(config)
            wind_dict = wind.initial_settings_wind(config)

            orb = orbital()
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                orb.initial_settings(config)
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                solver.solve_equation_motion(config, iteration, time_elapsed,
                                             coordinate_dict, velocity_dict, trajectory_dict,
                                             atmosphere_dict, aerodynamic_dict, None, wind_dict)
            orb.output_tecplot(config, iteration, time_elapsed, coordinate_dict, velocity_dict,
                               trajectory_dict, None, None, wind_dict)

        cls.variables, cls.rows = read_tecplot(os.path.join(cls.directory, 'tecplot.dat'))

    @classmethod
    def tearDownClass(cls):
        import shutil

        shutil.rmtree(cls.directory, ignore_errors=True)

    def column(self, name):
        return self.rows[:, self.variables.index(name)]

    def test_the_wind_columns_are_present(self):
        for name in ('WindE[m/s]', 'WindN[m/s]', 'WindU[m/s]', 'VelairAbs[m/s]'):
            self.assertIn(name, self.variables)

    def test_the_constant_wind_is_written_on_every_row(self):
        np.testing.assert_allclose(self.column('WindE[m/s]'), 45.0, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(self.column('WindN[m/s]'), -30.0, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(self.column('WindU[m/s]'), 2.0, rtol=0.0, atol=0.0)

    def test_the_air_relative_speed_is_consistent_with_the_other_columns(self):
        # 対地速度（地心ローカル系）から風を引いて独立に再計算する
        velocity_local = np.column_stack((self.column('Upl[m/s]'),
                                          self.column('Vpl[m/s]'),
                                          self.column('Wpl[m/s]')))
        wind_local = np.column_stack((self.column('WindE[m/s]'),
                                      self.column('WindN[m/s]'),
                                      self.column('WindU[m/s]')))
        reference = np.linalg.norm(velocity_local - wind_local, axis=1)
        np.testing.assert_allclose(self.column('VelairAbs[m/s]'), reference, rtol=1.e-12)

    def test_the_ground_speed_column_is_still_the_ground_speed(self):
        # 対気速度で上書きしてしまう取り違えを捕まえる
        velocity_local = np.column_stack((self.column('Upl[m/s]'),
                                          self.column('Vpl[m/s]'),
                                          self.column('Wpl[m/s]')))
        np.testing.assert_allclose(self.column('VelplAbs[m/s]'),
                                   np.linalg.norm(velocity_local, axis=1), rtol=1.e-12)
        self.assertFalse(np.allclose(self.column('VelplAbs[m/s]'), self.column('VelairAbs[m/s]')))

    def test_a_run_without_wind_has_no_wind_columns(self):
        import tempfile

        config = copy.deepcopy(load_config(CONFIG_REENTRY))
        config['computational_setup']['time_elapsed_maximum'] = 5.0
        config['time_integration']['timestep_constant'] = 0.5
        directory = tempfile.mkdtemp()
        config['post_process']['directory_output'] = directory
        config['post_process']['kml']['flag_output'] = False

        try:
            with quiet():
                atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
                aerodynamic_dict = satellite.initial_settings_satellite(config)
                orb = orbital()
                iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                    orb.initial_settings(config)
                iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                    solver.solve_equation_motion(config, iteration, time_elapsed,
                                                 coordinate_dict, velocity_dict, trajectory_dict,
                                                 atmosphere_dict, aerodynamic_dict)
                orb.output_tecplot(config, iteration, time_elapsed, coordinate_dict, velocity_dict,
                                   trajectory_dict)
            variables, rows = read_tecplot(os.path.join(directory, 'tecplot.dat'))
        finally:
            import shutil

            shutil.rmtree(directory, ignore_errors=True)

        self.assertEqual(len(variables), 14)
        for name in ('WindE[m/s]', 'VelairAbs[m/s]'):
            self.assertNotIn(name, variables)


class TestReadingAWindTable(unittest.TestCase):
    """
    形式は「1 行 1 格子点」の 1 種類だけ。**1 次元の鉛直プロファイルは、
    経度・緯度の節点が 1 つだけの場**として同じ形式・同じ経路で読む。
    """

    def tearDown(self):
        import shutil

        for directory in getattr(self, 'directories', []):
            shutil.rmtree(directory, ignore_errors=True)

    def build(self, rows, **keyword):
        wind_dict, output, directory = make_table_wind(rows, **keyword)
        self.directories = getattr(self, 'directories', []) + [directory]
        return wind_dict, output

    def test_a_three_dimensional_table_reproduces_its_nodes(self):
        longitude = [-10.0, 0.0, 10.0]
        latitude = [20.0, 30.0]
        altitude = [0.0, 10.0, 20.0]

        def field(value_lon, value_lat, value_alt):
            return [value_lon + value_alt, value_lat - value_alt, 0.0]

        wind_dict, _ = self.build(grid_rows(longitude, latitude, altitude, field))
        for value_lon in longitude:
            for value_lat in latitude:
                for value_alt in altitude:
                    np.testing.assert_allclose(local(wind_dict, value_lon, value_lat, value_alt),
                                               field(value_lon, value_lat, value_alt), atol=1.e-9)

    def test_the_interpolation_is_linear_between_the_nodes(self):
        rows = grid_rows([-10.0, 10.0], [0.0, 20.0], [0.0, 20.0],
                         lambda lo, la, al: [lo, la, al])
        wind_dict, _ = self.build(rows)
        np.testing.assert_allclose(local(wind_dict, 0.0, 10.0, 10.0), [0.0, 10.0, 10.0], atol=1.e-9)
        np.testing.assert_allclose(local(wind_dict, 5.0, 5.0, 15.0), [5.0, 5.0, 15.0], atol=1.e-9)

    def test_the_row_order_does_not_matter(self):
        rows = grid_rows([-10.0, 10.0], [0.0, 20.0], [0.0, 20.0],
                         lambda lo, la, al: [lo + la, al, 0.0])
        wind_direct, _ = self.build(rows)
        wind_shuffled, _ = self.build(list(reversed(rows)))
        np.testing.assert_array_equal(wind_direct[wind.KEY_WIND], wind_shuffled[wind.KEY_WIND])

    def test_a_one_dimensional_profile_depends_on_the_altitude_alone(self):
        # 経度・緯度の節点が 1 つ = 縮退した場。どこで引いても同じ鉛直分布になる
        altitude = [0.0, 10.0, 20.0, 30.0]
        rows = grid_rows([0.0], [40.0], altitude, lambda lo, la, al: [2.0*al, -al, 0.0])
        wind_dict, _ = self.build(rows)
        self.assertEqual(wind_dict[wind.KEY_AXIS], [3])
        for longitude, latitude in ((0.0, 40.0), (150.0, -60.0), (-179.0, 12.0)):
            np.testing.assert_allclose(local(wind_dict, longitude, latitude, 15.0),
                                       [30.0, -15.0, 0.0], atol=1.e-9)

    def test_a_table_with_one_altitude_is_constant_in_the_vertical(self):
        rows = grid_rows([-10.0, 10.0], [0.0, 20.0], [5.0], lambda lo, la, al: [lo, la, 0.0])
        wind_dict, _ = self.build(rows)
        self.assertEqual(wind_dict[wind.KEY_AXIS], [1, 2])
        np.testing.assert_allclose(local(wind_dict, 0.0, 10.0, 5.0), [0.0, 10.0, 0.0], atol=1.e-9)
        # 高度軸が縮退しているので範囲外の判定も起きない
        np.testing.assert_allclose(local(wind_dict, 0.0, 10.0, 400.0), [0.0, 10.0, 0.0], atol=1.e-9)

    def test_a_single_point_table_behaves_like_a_constant_wind(self):
        wind_dict, _ = self.build([[0.0, 0.0, 0.0, 7.0, -3.0, 1.0]])
        self.assertIsNone(wind_dict[wind.KEY_INTERP])
        np.testing.assert_allclose(local(wind_dict, 100.0, -40.0, 250.0), [7.0, -3.0, 1.0], atol=0.0)

    def test_the_longitude_is_folded_into_minus_180_to_180(self):
        # 気象データは 0-360 deg. で配られる。畳んだあとの軸で引ける必要がある
        rows = grid_rows([235.0, 240.0], [-15.0, -10.0], [0.0, 10.0],
                         lambda lo, la, al: [lo, la, 0.0])
        wind_dict, _ = self.build(rows)
        np.testing.assert_allclose(wind_dict[wind.KEY_LONGITUDE], [-125.0, -120.0], atol=1.e-9)
        np.testing.assert_allclose(local(wind_dict, -125.0, -15.0, 0.0)[0], 235.0, atol=1.e-9)
        np.testing.assert_allclose(local(wind_dict, 235.0, -15.0, 0.0)[0], 235.0, atol=1.e-9)

    def test_a_global_table_interpolates_across_the_seam(self):
        # 全球なら経度 180 度の継ぎ目で不連続になってはいけない
        longitude = list(np.arange(0.0, 360.0, 45.0))
        rows = grid_rows(longitude, [0.0], [0.0, 40.0],
                         lambda lo, la, al: [np.cos(np.deg2rad(lo)), np.sin(np.deg2rad(lo)), 0.0])
        wind_dict, _ = self.build(rows)
        self.assertTrue(wind.is_global_longitude(wind_dict[wind.KEY_LONGITUDE]))
        below = local(wind_dict, 179.999, 0.0, 10.0)
        above = local(wind_dict, -179.999, 0.0, 10.0)
        np.testing.assert_allclose(below, above, atol=1.e-4)
        # 継ぎ目でも節点の値そのものが出る
        np.testing.assert_allclose(local(wind_dict, 180.0, 0.0, 10.0), [-1.0, 0.0, 0.0], atol=1.e-9)

    def test_a_regional_table_is_clamped_horizontally(self):
        rows = grid_rows([-10.0, 10.0], [0.0, 20.0], [0.0, 20.0],
                         lambda lo, la, al: [lo, la, 0.0])
        wind_dict, _ = self.build(rows)
        np.testing.assert_allclose(local(wind_dict, 60.0, 10.0, 10.0), [10.0, 10.0, 0.0], atol=1.e-9)
        np.testing.assert_allclose(local(wind_dict, -60.0, -40.0, 10.0), [-10.0, 0.0, 0.0], atol=1.e-9)

    def test_the_altitude_range_is_reported_and_zeroed_by_default(self):
        rows = grid_rows([0.0], [0.0], [0.0, 30.0], lambda lo, la, al: [10.0, 5.0, 0.0])
        wind_dict, output = self.build(rows)
        self.assertIn('zero', output)
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 15.0), [10.0, 5.0, 0.0], atol=1.e-9)
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 150.0), np.zeros(3), atol=0.0)
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, -5.0), np.zeros(3), atol=0.0)

    def test_the_altitude_range_can_be_clamped_instead(self):
        rows = grid_rows([0.0], [0.0], [0.0, 30.0], lambda lo, la, al: [al, 0.0, 0.0])
        wind_dict, _ = self.build(rows, extrapolation='clamp')
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 150.0), [30.0, 0.0, 0.0], atol=1.e-9)
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, -5.0), [0.0, 0.0, 0.0], atol=1.e-9)

    def test_the_nodes_at_the_very_edge_are_inside(self):
        # deg -> rad -> deg の往復で最下位ビットが動くので、端の節点で
        # 範囲外と判定されると zero 外挿で風が消える
        rows = grid_rows([-125.0, -120.0], [-15.0, -10.0], [0.1344, 30.9994],
                         lambda lo, la, al: [11.0, -4.0, 0.0])
        wind_dict, _ = self.build(rows)
        for longitude in (-125.0, -120.0):
            for latitude in (-15.0, -10.0):
                for altitude in (0.1344, 30.9994):
                    np.testing.assert_allclose(local(wind_dict, longitude, latitude, altitude),
                                               [11.0, -4.0, 0.0], atol=1.e-9)
        self.assertFalse(wind_dict[wind.KEY_WARNED])

    def test_leaving_the_range_warns_only_once(self):
        rows = grid_rows([0.0], [0.0], [0.0, 30.0], lambda lo, la, al: [10.0, 0.0, 0.0])
        wind_dict, _ = self.build(rows)
        with quiet() as first:
            local(wind_dict, 0.0, 0.0, 150.0)
        with quiet() as second:
            local(wind_dict, 0.0, 0.0, 160.0)
        self.assertIn('Caution', first.getvalue())
        self.assertNotIn('Caution', second.getvalue())

    def test_the_interpolator_is_built_once(self):
        # 評価のたびに作り直すと計算時間を支配する（大気・空力と同じ方針）
        rows = grid_rows([-10.0, 10.0], [0.0, 20.0], [0.0, 20.0],
                         lambda lo, la, al: [lo, la, al])
        wind_dict, _ = self.build(rows)
        interpolator = wind_dict[wind.KEY_INTERP]
        local(wind_dict, 0.0, 10.0, 10.0)
        self.assertIs(wind_dict[wind.KEY_INTERP], interpolator)


class TestTheWindTableIsChecked(unittest.TestCase):

    def tearDown(self):
        import shutil

        for directory in getattr(self, 'directories', []):
            shutil.rmtree(directory, ignore_errors=True)

    def build(self, rows, **keyword):
        wind_dict, output, directory = make_table_wind(rows, **keyword)
        self.directories = getattr(self, 'directories', []) + [directory]
        return wind_dict, output

    def test_an_incomplete_grid_stops_the_run(self):
        # 欠けた点を黙って 0 にすると、風が弱いのかデータが無いのか区別できない
        rows = grid_rows([-10.0, 10.0], [0.0, 20.0], [0.0, 20.0],
                         lambda lo, la, al: [1.0, 2.0, 0.0])
        with self.assertRaises(SystemExit):
            self.build(rows[:-1])

    def test_a_duplicated_grid_point_stops_the_run(self):
        rows = grid_rows([-10.0, 10.0], [0.0, 20.0], [0.0, 20.0],
                         lambda lo, la, al: [1.0, 2.0, 0.0])
        rows[1] = list(rows[0])
        with self.assertRaises(SystemExit):
            self.build(rows)

    def test_writing_both_minus_180_and_180_stops_the_run(self):
        # 同じ子午線なので、畳んだあとに重複する
        rows = grid_rows([-180.0, 0.0, 180.0], [0.0], [0.0, 20.0],
                         lambda lo, la, al: [1.0, 0.0, 0.0])
        with self.assertRaises(SystemExit):
            self.build(rows)

    def test_an_empty_table_stops_the_run(self):
        with self.assertRaises(SystemExit):
            self.build([])

    def test_a_missing_file_stops_the_run(self):
        import tempfile

        directory = tempfile.mkdtemp()
        self.directories = getattr(self, 'directories', []) + [directory]
        with self.assertRaises(SystemExit):
            with quiet():
                wind.initial_settings_wind(table_config(directory, 'no_such_table.txt'))

    def test_an_unknown_extrapolation_stops_the_run(self):
        rows = grid_rows([0.0], [0.0], [0.0, 30.0], lambda lo, la, al: [1.0, 0.0, 0.0])
        with self.assertRaises(SystemExit):
            self.build(rows, extrapolation='exponential')

    def test_rows_with_the_wrong_number_of_columns_are_ignored(self):
        # 末尾の空行やゴミ行は無害（大気テーブルと同じ扱い）
        import tempfile

        directory = tempfile.mkdtemp()
        self.directories = getattr(self, 'directories', []) + [directory]
        path = os.path.join(directory, 'windmodel.txt')
        write_wind_table(path, grid_rows([0.0], [0.0], [0.0, 30.0],
                                         lambda lo, la, al: [7.0, 0.0, 0.0]))
        with open(path, 'a') as f:
            f.write('\n')
            f.write('  end of file\n')
            f.write('  1.0 2.0\n')
        with quiet():
            wind_dict = wind.initial_settings_wind(table_config(directory, 'windmodel.txt'))
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 10.0), [7.0, 0.0, 0.0], atol=1.e-9)


class TestTheEpochOfTheWindTable(unittest.TestCase):
    """
    テーブルの生成時刻と config のエポックの食い違いを警告する。
    大気テーブルで「生成条件が違うと 400 km で密度が 2 倍違う」という前例がある。
    """

    def tearDown(self):
        import shutil

        for directory in getattr(self, 'directories', []):
            shutil.rmtree(directory, ignore_errors=True)

    def build(self, epoch_table, epoch_config):
        import epoch.epoch as epoch_module

        rows = grid_rows([0.0], [0.0], [0.0, 30.0], lambda lo, la, al: [1.0, 0.0, 0.0])
        epoch_dict = None
        if epoch_config is not None:
            with quiet():
                epoch_dict = epoch_module.initial_settings_epoch(
                    {'epoch': {'flag_epoch': True, 'datetime': epoch_config}})
        wind_dict, output, directory = make_table_wind(rows, epoch=epoch_table, epoch_dict=epoch_dict)
        self.directories = getattr(self, 'directories', []) + [directory]
        return output

    def test_a_matching_epoch_is_quiet(self):
        output = self.build('2024-01-01T00:00:00.000Z', '2024-01-01T00:00:00Z')
        self.assertNotIn('Warning', output)

    def test_a_difference_within_the_tolerance_is_quiet(self):
        output = self.build('2024-01-01T00:00:00.000Z', '2024-01-01T00:30:00Z')
        self.assertNotIn('Warning', output)

    def test_a_different_epoch_warns(self):
        output = self.build('2024-01-01T00:00:00.000Z', '2024-07-01T00:00:00Z')
        self.assertIn('Warning', output)
        self.assertIn('different time', output)

    def test_a_table_without_an_epoch_is_quiet(self):
        output = self.build(None, '2024-01-01T00:00:00Z')
        self.assertNotIn('Warning', output)

    def test_no_epoch_in_the_config_is_quiet(self):
        output = self.build('2024-01-01T00:00:00.000Z', None)
        self.assertNotIn('Warning', output)


class TestTheShippedWindTables(unittest.TestCase):
    """同梱のテーブルが読めること。`database/wind/` はマスターなので default で引く。"""

    def read(self, filename, extrapolation='zero'):
        config = {'wind': {'flag_wind': True, 'kind_wind_model': 'fileread',
                           'directory_path_specify': 'default',
                           'filename_wind': filename,
                           'kind_extrapolation': extrapolation}}
        with quiet():
            return wind.initial_settings_wind(config)

    def test_the_ncep_sample_is_a_three_dimensional_field(self):
        wind_dict = self.read('wind_ncep_20240101_pacific.txt')
        self.assertEqual(wind_dict[wind.KEY_AXIS], [1, 2, 3])
        self.assertEqual(len(wind_dict[wind.KEY_LONGITUDE]), 9)
        self.assertEqual(len(wind_dict[wind.KEY_LATITUDE]), 9)
        self.assertEqual(len(wind_dict[wind.KEY_ALTITUDE]), 17)
        self.assertLess(wind_dict[wind.KEY_ALTITUDE][-1], 32.0)
        # 実データなので風速は数十 m/s の桁におさまる
        self.assertLess(np.abs(wind_dict[wind.KEY_WIND]).max(), 100.0)
        # 鉛直風はモデル化していない
        np.testing.assert_array_equal(wind_dict[wind.KEY_WIND][:,:,:,:,2], 0.0)

    def test_the_ncep_sample_carries_its_epoch(self):
        wind_dict = self.read('wind_ncep_20240101_pacific.txt')
        self.assertIsNotNone(wind_dict[wind.KEY_EPOCH])
        self.assertIn('2024-01-01', wind_dict[wind.KEY_EPOCH])

    def test_the_profile_sample_is_one_dimensional(self):
        wind_dict = self.read('wind_profile_sample.txt')
        self.assertEqual(wind_dict[wind.KEY_AXIS], [3])
        self.assertEqual(len(wind_dict[wind.KEY_LONGITUDE]), 1)
        # 対流圏界面付近にジェットの最大がある
        east = wind_dict[wind.KEY_WIND][0,0,0,:,0]
        altitude = wind_dict[wind.KEY_ALTITUDE]
        self.assertAlmostEqual(altitude[int(np.argmax(east))], 10.0, places=6)


class TestTheGeneratorScript(unittest.TestCase):
    """
    生成スクリプトのオフラインで動く経路。NCEP を取ってくる経路は
    ネットワークが要るのでここでは走らせない（CI で落ちるため）。
    """

    def setUp(self):
        import tempfile

        self.directory = tempfile.mkdtemp()

    def tearDown(self):
        import shutil

        shutil.rmtree(self.directory, ignore_errors=True)

    def script(self):
        return os.path.join(ROOT_DIR, 'database', 'wind', 'generate_wind_table.py')

    def test_the_profile_source_writes_a_table_the_solver_can_read(self):
        import subprocess
        import sys

        source = os.path.join(self.directory, 'profile.txt')
        with open(source, 'w') as f:
            f.write('# altitude[km] east[m/s] north[m/s]\n')
            f.write('  0.0   5.0  0.0\n')
            f.write(' 10.0  40.0  2.0\n')
            f.write(' 30.0 -10.0  0.0\n')

        output = os.path.join(self.directory, 'windmodel.txt')
        result = subprocess.run([sys.executable, self.script(), '--source', 'profile',
                                 '-i', source, '--centre', '140.0', '35.0',
                                 '--datetime', '2024-01-01T00:00:00Z', '-o', output],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.exists(output))

        with quiet():
            wind_dict = wind.initial_settings_wind(table_config(self.directory, 'windmodel.txt'))
        self.assertEqual(wind_dict[wind.KEY_AXIS], [3])
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 10.0), [40.0, 2.0, 0.0], atol=1.e-9)
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 5.0), [22.5, 1.0, 0.0], atol=1.e-9)
        self.assertIn('2024-01-01', wind_dict[wind.KEY_EPOCH])

    def test_the_calm_source_writes_a_readable_table(self):
        import subprocess
        import sys

        output = os.path.join(self.directory, 'windmodel.txt')
        result = subprocess.run([sys.executable, self.script(), '--source', 'calm',
                                 '--altitude-max', '50.0', '-o', output],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with quiet():
            wind_dict = wind.initial_settings_wind(table_config(self.directory, 'windmodel.txt'))
        np.testing.assert_array_equal(wind_dict[wind.KEY_WIND], 0.0)

    def test_the_ncep_source_needs_a_datetime(self):
        import subprocess
        import sys

        result = subprocess.run([sys.executable, self.script(), '--source', 'ncep',
                                 '-o', os.path.join(self.directory, 'windmodel.txt')],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--datetime', result.stdout + result.stderr)


class TestTheSolverWithAWindTable(unittest.TestCase):

    def setUp(self):
        import tempfile

        self.directory = tempfile.mkdtemp()

    def tearDown(self):
        import shutil

        shutil.rmtree(self.directory, ignore_errors=True)

    def run_with_table(self, rows, extrapolation='zero', time_max=2700.0, timestep=1.0):
        write_wind_table(os.path.join(self.directory, 'windmodel.txt'), rows)
        config = copy.deepcopy(load_config(CONFIG_REENTRY))
        config['computational_setup']['time_elapsed_maximum'] = time_max
        config['time_integration']['timestep_constant'] = timestep
        config['wind'] = table_config(self.directory, 'windmodel.txt', extrapolation)['wind']

        with quiet():
            atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
            aerodynamic_dict = satellite.initial_settings_satellite(config)
            wind_dict = wind.initial_settings_wind(config)
            orb = orbital()
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                orb.initial_settings(config)
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                solver.solve_equation_motion(config, iteration, time_elapsed,
                                             coordinate_dict, velocity_dict, trajectory_dict,
                                             atmosphere_dict, aerodynamic_dict, None, wind_dict)
        return coordinate_dict, velocity_dict

    def test_a_profile_table_of_zero_wind_reproduces_the_run_without_wind(self):
        rows = grid_rows([0.0], [0.0], [0.0, 200.0], lambda lo, la, al: [0.0, 0.0, 0.0])
        coordinate_table, velocity_table = self.run_with_table(rows)
        config = load_config(CONFIG_REENTRY)
        _, coordinate_calm, velocity_calm = run_solver(config, 2700.0, 1.0)
        np.testing.assert_array_equal(np.array(coordinate_table['cartesian']),
                                      np.array(coordinate_calm['cartesian']))
        np.testing.assert_array_equal(np.array(velocity_table['cartesian']),
                                      np.array(velocity_calm['cartesian']))

    def test_a_uniform_profile_table_matches_the_constant_model(self):
        # 同じ風を 2 通りの経路で与えて突き合わせる
        rows = grid_rows([0.0], [0.0], [0.0, 200.0], lambda lo, la, al: [30.0, -10.0, 0.0])
        coordinate_table, velocity_table = self.run_with_table(rows)
        config = load_config(CONFIG_REENTRY)
        _, coordinate_constant, velocity_constant = run_solver(config, 2700.0, 1.0,
                                                               wind_velocity=[30.0, -10.0, 0.0])
        np.testing.assert_allclose(np.array(coordinate_table['cartesian']),
                                   np.array(coordinate_constant['cartesian']), rtol=1.e-12, atol=1.e-6)

    def test_a_table_that_stops_below_the_trajectory_only_acts_where_it_has_data(self):
        #
        # 高度 20 km までしか持たないテーブル。zero 外挿ならそれより上では風なし。
        # 20 km 以下だけで流されるので、全高度に風を置いた場合より効きが小さくなる
        #
        rows_low = grid_rows([0.0], [0.0], [0.0, 20.0], lambda lo, la, al: [40.0, 0.0, 0.0])
        coordinate_low, _ = self.run_with_table(rows_low)
        rows_all = grid_rows([0.0], [0.0], [0.0, 200.0], lambda lo, la, al: [40.0, 0.0, 0.0])
        coordinate_all, _ = self.run_with_table(rows_all)
        config = load_config(CONFIG_REENTRY)
        _, coordinate_calm, _ = run_solver(config, 2700.0, 1.0)

        drift_low = coordinate_low['geodetic'][-1][0] - coordinate_calm['geodetic'][-1][0]
        drift_all = coordinate_all['geodetic'][-1][0] - coordinate_calm['geodetic'][-1][0]
        self.assertGreater(drift_low, 0.0)
        self.assertGreater(drift_all, drift_low)

    def test_the_wind_follows_the_position_along_the_trajectory(self):
        #
        # 緯度に比例する風を置くと、南へ下る軌道に沿って風も変わる。
        # 出力の風の列を、位置から独立に引き直して突き合わせる
        #
        rows = grid_rows([-120.0, 0.0, 120.0], [-60.0, 0.0, 60.0], [0.0, 200.0],
                         lambda lo, la, al: [la, 0.0, 0.0])
        coordinate_dict, _ = self.run_with_table(rows, time_max=600.0, timestep=1.0)

        config = copy.deepcopy(load_config(CONFIG_REENTRY))
        config['wind'] = table_config(self.directory, 'windmodel.txt')['wind']
        with quiet():
            wind_dict = wind.initial_settings_wind(config)

        latitude_list = []
        east_list = []
        for index in range(0, len(coordinate_dict['geodetic']), 100):
            geodetic = coordinate_dict['geodetic'][index]
            latitude_list.append(geodetic[1]*180.0/np.pi)
            east_list.append(wind.get_wind_local(geodetic, wind_dict)[0])
        # 風の東成分は緯度そのもの（緯度は測地、風の格子も測地の緯度で引いている）
        np.testing.assert_allclose(east_list, latitude_list, atol=1.e-9)
        self.assertGreater(max(latitude_list) - min(latitude_list), 5.0)


class TestTheTimeAxisOfAWindTable(unittest.TestCase):
    """
    7 列のテーブルは先頭が時刻（テーブルのエポックからの秒）。
    時刻軸は 4 番目の次元というだけで、縮退の扱いも範囲外の扱いも他の軸と同じ。
    """

    def tearDown(self):
        import shutil

        for directory in getattr(self, 'directories', []):
            shutil.rmtree(directory, ignore_errors=True)

    def build(self, rows, **keyword):
        wind_dict, output, directory = make_table_wind(rows, **keyword)
        self.directories = getattr(self, 'directories', []) + [directory]
        return wind_dict, output

    def series_rows(self):
        # 東成分が時刻に比例して増える表（0 s で 0、3600 s で 36）
        return grid_rows_time([0.0, 1800.0, 3600.0], [0.0], [0.0], [0.0, 100.0],
                              lambda t, lo, la, al: [t/100.0, 0.0, 0.0])

    def test_the_time_axis_is_read(self):
        wind_dict, _ = self.build(self.series_rows(),
                                  epoch='2024-01-01T00:00:00Z',
                                  epoch_dict=epoch_of('2024-01-01T00:00:00Z'))
        np.testing.assert_allclose(wind_dict[wind.KEY_TIME], [0.0, 1800.0, 3600.0])
        self.assertEqual(wind_dict[wind.KEY_AXIS], [0, 3])

    def test_the_nodes_of_the_time_axis_are_reproduced(self):
        wind_dict, _ = self.build(self.series_rows(),
                                  epoch='2024-01-01T00:00:00Z',
                                  epoch_dict=epoch_of('2024-01-01T00:00:00Z'))
        for value_time in (0.0, 1800.0, 3600.0):
            np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 50.0, value_time)[0],
                                       value_time/100.0, atol=1.e-9)

    def test_the_interpolation_in_time_is_linear(self):
        wind_dict, _ = self.build(self.series_rows(),
                                  epoch='2024-01-01T00:00:00Z',
                                  epoch_dict=epoch_of('2024-01-01T00:00:00Z'))
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 50.0, 900.0)[0], 9.0, atol=1.e-9)
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 50.0, 2700.0)[0], 27.0, atol=1.e-9)

    def test_the_epoch_of_the_run_shifts_the_time_axis(self):
        # config のエポックがテーブルのエポックより 30 分後なら、経過 0 s は表の 1800 s
        wind_dict, output = self.build(self.series_rows(),
                                       epoch='2024-01-01T00:00:00Z',
                                       epoch_dict=epoch_of('2024-01-01T00:30:00Z'))
        self.assertAlmostEqual(wind_dict[wind.KEY_OFFSET], 1800.0, places=6)
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 50.0, 0.0)[0], 18.0, atol=1.e-9)
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 50.0, 900.0)[0], 27.0, atol=1.e-9)

    def test_the_time_axis_is_clamped_outside_its_range(self):
        wind_dict, _ = self.build(self.series_rows(),
                                  epoch='2024-01-01T00:00:00Z',
                                  epoch_dict=epoch_of('2024-01-01T00:00:00Z'))
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 50.0, 100000.0)[0], 36.0, atol=1.e-9)
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 50.0, -100.0)[0], 0.0, atol=1.e-9)

    def test_a_single_snapshot_ignores_the_elapsed_time(self):
        rows = grid_rows([0.0], [0.0], [0.0, 100.0], lambda lo, la, al: [12.0, 0.0, 0.0])
        wind_dict, _ = self.build(rows)
        for value_time in (0.0, 1000.0, 1.e6):
            np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 50.0, value_time)[0], 12.0, atol=0.0)

    def test_a_time_axis_without_an_epoch_in_the_table_stops_the_run(self):
        with self.assertRaises(SystemExit):
            self.build(self.series_rows(), epoch_dict=epoch_of('2024-01-01T00:00:00Z'))

    def test_a_time_axis_without_an_epoch_in_the_config_stops_the_run(self):
        with self.assertRaises(SystemExit):
            self.build(self.series_rows(), epoch='2024-01-01T00:00:00Z')

    def test_an_epoch_outside_the_time_range_warns(self):
        _, output = self.build(self.series_rows(),
                               epoch='2024-01-01T00:00:00Z',
                               epoch_dict=epoch_of('2024-01-02T00:00:00Z'))
        self.assertIn('outside the time range', output)

    def test_mixing_six_and_seven_column_rows_stops_the_run(self):
        rows = self.series_rows()
        rows[0] = rows[0][1:]
        with self.assertRaises(SystemExit):
            self.build(rows)

    def test_a_four_dimensional_table_reproduces_its_nodes(self):
        time = [0.0, 3600.0]
        longitude = [-10.0, 10.0]
        latitude = [0.0, 20.0]
        altitude = [0.0, 50.0]

        def field(t, lo, la, al):
            return [t/3600.0, lo + la, al]

        wind_dict, _ = self.build(grid_rows_time(time, longitude, latitude, altitude, field),
                                  epoch='2024-01-01T00:00:00Z',
                                  epoch_dict=epoch_of('2024-01-01T00:00:00Z'))
        self.assertEqual(wind_dict[wind.KEY_AXIS], [0, 1, 2, 3])
        for value_time in time:
            for value_lon in longitude:
                for value_lat in latitude:
                    for value_alt in altitude:
                        np.testing.assert_allclose(
                            local(wind_dict, value_lon, value_lat, value_alt, value_time),
                            field(value_time, value_lon, value_lat, value_alt), atol=1.e-9)


class TestTheSolverWithATimeDependentTable(unittest.TestCase):

    def setUp(self):
        import tempfile

        self.directory = tempfile.mkdtemp()

    def tearDown(self):
        import shutil

        shutil.rmtree(self.directory, ignore_errors=True)

    def run_case(self, rows, epoch=None, time_max=600.0, timestep=1.0, output=False):
        write_wind_table(os.path.join(self.directory, 'windmodel.txt'), rows, epoch=epoch)
        config = copy.deepcopy(load_config(CONFIG_REENTRY))
        config['computational_setup']['time_elapsed_maximum'] = time_max
        config['time_integration']['timestep_constant'] = timestep
        config['wind'] = table_config(self.directory, 'windmodel.txt')['wind']
        config['epoch'] = {'flag_epoch': True, 'datetime': '2024-01-01T00:00:00Z'}
        config['post_process']['directory_output'] = self.directory
        config['post_process']['kml']['flag_output'] = False
        config['post_process']['tecplot']['flag_output'] = bool(output)
        config['post_process']['tecplot']['frequency_output'] = 100

        with quiet():
            epoch_dict = epoch_of('2024-01-01T00:00:00Z')
            atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
            aerodynamic_dict = satellite.initial_settings_satellite(config)
            wind_dict = wind.initial_settings_wind(config, epoch_dict)
            orb = orbital()
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                orb.initial_settings(config)
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                solver.solve_equation_motion(config, iteration, time_elapsed,
                                             coordinate_dict, velocity_dict, trajectory_dict,
                                             atmosphere_dict, aerodynamic_dict, None, wind_dict)
            if output:
                orb.output_tecplot(config, iteration, time_elapsed, coordinate_dict, velocity_dict,
                                   trajectory_dict, None, epoch_dict, wind_dict)
        return coordinate_dict, velocity_dict

    def rising_wind(self):
        # 東成分が 0 から 60 m/s まで 600 s で増える
        return grid_rows_time([0.0, 600.0], [0.0], [0.0], [0.0, 200.0],
                              lambda t, lo, la, al: [60.0*t/600.0, 0.0, 0.0])

    def test_a_time_dependent_table_differs_from_its_first_snapshot(self):
        _, _ = self.run_case(self.rising_wind(), epoch='2024-01-01T00:00:00Z')
        coordinate_rising, _ = self.run_case(self.rising_wind(), epoch='2024-01-01T00:00:00Z')
        rows_static = grid_rows([0.0], [0.0], [0.0, 200.0], lambda lo, la, al: [0.0, 0.0, 0.0])
        coordinate_static, _ = self.run_case(rows_static)
        self.assertNotAlmostEqual(coordinate_rising['geodetic'][-1][0],
                                  coordinate_static['geodetic'][-1][0], places=9)

    def test_the_drift_lies_between_the_two_snapshots(self):
        #
        # 0 -> 60 m/s と増える風で流される量は、ずっと 0 m/s の場合と
        # ずっと 60 m/s の場合の間になければならない
        #
        coordinate_rising, _ = self.run_case(self.rising_wind(), epoch='2024-01-01T00:00:00Z')
        rows_zero = grid_rows([0.0], [0.0], [0.0, 200.0], lambda lo, la, al: [0.0, 0.0, 0.0])
        rows_full = grid_rows([0.0], [0.0], [0.0, 200.0], lambda lo, la, al: [60.0, 0.0, 0.0])
        coordinate_zero, _ = self.run_case(rows_zero)
        coordinate_full, _ = self.run_case(rows_full)
        longitude_rising = coordinate_rising['geodetic'][-1][0]
        self.assertGreater(longitude_rising, coordinate_zero['geodetic'][-1][0])
        self.assertLess(longitude_rising, coordinate_full['geodetic'][-1][0])

    def test_the_runge_kutta_stages_use_their_own_time(self):
        #
        # 時刻に依る力を RK4 で正しく積むには、各段が自分の時刻の風を見なければならない
        # （段の時刻は t, t+dt/2, t+dt/2, t+dt）。
        #
        # 1 ステップだけ（dt = 窓の長さ）進めると、これが露わになる。0 から W へ
        # 線形に増える風なら、RK4 の重み (1/6, 2/6, 2/6, 1/6) と段の時刻 (0, 1/2, 1/2, 1)
        # から加速度の実効値は「風 W/2 のときの加速度」に等しく、**速度の増分が一致する**。
        # 全段を 1 段目の時刻で評価すると実効の風は 0、つまり無風と同じ答えになる。
        #
        # 位置で比べてはいけない。加速度が時刻に比例するときの変位は A T^2/6 で、
        # 一定加速度 A/2 の変位 A T^2/4 の 2/3 にしかならない（実測でも 33 % ずれる）。
        # 一致するのは速度の増分だけである。
        #
        window = 60.0
        rows_ramp = grid_rows_time([0.0, window], [0.0], [0.0], [0.0, 200.0],
                                   lambda t, lo, la, al: [200.0*t/window, 0.0, 0.0])
        rows_zero = grid_rows([0.0], [0.0], [0.0, 200.0], lambda lo, la, al: [0.0, 0.0, 0.0])
        rows_half = grid_rows([0.0], [0.0], [0.0, 200.0], lambda lo, la, al: [100.0, 0.0, 0.0])

        _, velocity_ramp = self.run_case(rows_ramp, epoch='2024-01-01T00:00:00Z',
                                         time_max=window, timestep=window)
        _, velocity_zero = self.run_case(rows_zero, time_max=window, timestep=window)
        _, velocity_half = self.run_case(rows_half, time_max=window, timestep=window)

        final_ramp = np.array(velocity_ramp['cartesian'][-1])
        final_zero = np.array(velocity_zero['cartesian'][-1])
        final_half = np.array(velocity_half['cartesian'][-1])

        # 無風とは違う（全段を 1 段目の時刻で見ていれば一致してしまう）
        self.assertGreater(np.linalg.norm(final_ramp - final_zero), 1.e-9)
        # 風 W/2 のときの速度増分と一致する。抗力が |v_air| v_air なので厳密ではなく、
        # 1 ステップが 60 s ＝ 450 km 分あって密度も動くため、実測の残差は 1.3 %。
        # 段の時刻を 1 段目で固定すると比は 1.0 になるので、5 % で十分に分かれる
        difference = np.linalg.norm(final_ramp - final_half)
        scale = np.linalg.norm(final_half - final_zero)
        self.assertGreater(scale, 0.0)
        self.assertLess(difference/scale, 0.05)

    def test_the_wind_written_to_the_output_follows_the_time_of_each_row(self):
        # 出力の風の列を、その行の時刻から独立に引き直して突き合わせる
        self.run_case(self.rising_wind(), epoch='2024-01-01T00:00:00Z', output=True)
        variables, rows = read_tecplot(os.path.join(self.directory, 'tecplot.dat'))
        time_column = rows[:, variables.index('Time[s]')]
        east_column = rows[:, variables.index('WindE[m/s]')]
        np.testing.assert_allclose(east_column, 60.0*time_column/600.0, atol=1.e-9)
        self.assertGreater(east_column.max() - east_column.min(), 50.0)


class TestMergingTwoWindTables(unittest.TestCase):
    """
    NCEP（下層）と HWM14（上層）を遷移層で繋ぐ経路。
    境目でそのまま切り替えると風速が飛ぶので、重みを線形に移す。
    """

    def setUp(self):
        import tempfile

        self.directory = tempfile.mkdtemp()

    def tearDown(self):
        import shutil

        shutil.rmtree(self.directory, ignore_errors=True)

    def script(self):
        return os.path.join(ROOT_DIR, 'database', 'wind', 'generate_wind_table.py')

    def write(self, name, rows, epoch=None):
        path = os.path.join(self.directory, name)
        write_wind_table(path, rows, epoch=epoch)
        return path

    def merge(self, low, high, transition=('20', '30'), name='merged.txt'):
        import subprocess
        import sys

        output = os.path.join(self.directory, name)
        result = subprocess.run([sys.executable, self.script(), '--source', 'merge',
                                 '-i', low, high, '--transition', transition[0], transition[1],
                                 '-o', output], capture_output=True, text=True)
        return result, output

    def test_the_transition_layer_blends_the_two_tables(self):
        #
        # 節点を遷移層の内側（22, 24, 28 km）に置く。節点が遷移層の外にしか無いと
        # 合成の重みが結果に現れず、階段状に切り替えても同じ答えになってしまう。
        # 見る点も 22 と 28 の非対称な 2 つにする（中点だけだと重みの反転を見逃す）。
        #
        altitude = [0.0, 10.0, 20.0, 22.0, 24.0, 28.0, 30.0]
        low = self.write('low.txt', grid_rows([0.0], [0.0], altitude,
                                              lambda lo, la, al: [10.0, 0.0, 0.0]),
                         epoch='2024-01-01T00:00:00Z')
        high = self.write('high.txt', grid_rows([0.0], [0.0], altitude + [100.0],
                                                lambda lo, la, al: [-50.0, 0.0, 0.0]),
                          epoch='2024-01-01T00:00:00Z')
        result, output = self.merge(low, high)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

        with quiet():
            wind_dict = wind.initial_settings_wind(table_config(self.directory, 'merged.txt'))
        # 遷移層の下は下層の値、上は上層の値、途中は重み (z-20)/10 の線形結合
        for altitude_query, weight in ((10.0, 0.0), (20.0, 0.0), (22.0, 0.2), (24.0, 0.4),
                                       (28.0, 0.8), (30.0, 1.0), (100.0, 1.0)):
            expected = (1.0 - weight)*10.0 + weight*(-50.0)
            np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, altitude_query)[0], expected,
                                       atol=1.e-5, err_msg='altitude {:g} km'.format(altitude_query))

    def test_the_horizontal_nodes_stay_aligned_through_the_merge(self):
        # 格子点ごとに違う値を入れて、合成で並びが崩れないことを見る
        altitude = [0.0, 20.0, 25.0, 30.0]
        low = self.write('low.txt', grid_rows([-10.0, 10.0], [0.0, 20.0], altitude,
                                              lambda lo, la, al: [lo, la, 0.0]))
        high = self.write('high.txt', grid_rows([-10.0, 10.0], [0.0, 20.0], altitude,
                                                lambda lo, la, al: [10.0*lo, 10.0*la, 0.0]))
        result, _ = self.merge(low, high)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with quiet():
            wind_dict = wind.initial_settings_wind(table_config(self.directory, 'merged.txt'))
        for longitude in (-10.0, 10.0):
            for latitude in (0.0, 20.0):
                np.testing.assert_allclose(local(wind_dict, longitude, latitude, 0.0),
                                           [longitude, latitude, 0.0], atol=1.e-5)
                np.testing.assert_allclose(local(wind_dict, longitude, latitude, 30.0),
                                           [10.0*longitude, 10.0*latitude, 0.0], atol=1.e-5)

    def test_the_merged_table_spans_both_altitude_ranges(self):
        low = self.write('low.txt', grid_rows([0.0], [0.0], [0.0, 31.0],
                                              lambda lo, la, al: [5.0, 0.0, 0.0]))
        high = self.write('high.txt', grid_rows([0.0], [0.0], [20.0, 500.0],
                                                lambda lo, la, al: [-5.0, 0.0, 0.0]))
        result, output = self.merge(low, high)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with quiet():
            wind_dict = wind.initial_settings_wind(table_config(self.directory, 'merged.txt'))
        self.assertAlmostEqual(wind_dict[wind.KEY_ALTITUDE][0], 0.0, places=6)
        self.assertAlmostEqual(wind_dict[wind.KEY_ALTITUDE][-1], 500.0, places=6)

    def test_tables_on_different_horizontal_grids_are_rejected(self):
        low = self.write('low.txt', grid_rows([0.0, 10.0], [0.0], [0.0, 30.0],
                                              lambda lo, la, al: [1.0, 0.0, 0.0]))
        high = self.write('high.txt', grid_rows([0.0], [0.0], [20.0, 100.0],
                                                lambda lo, la, al: [1.0, 0.0, 0.0]))
        result, _ = self.merge(low, high)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('longitude axis', result.stdout + result.stderr)

    def test_tables_with_the_same_number_of_nodes_but_different_values_are_rejected(self):
        # 節点数が同じだと長さの比較では見つからない。値まで見ていないと通ってしまう
        low = self.write('low.txt', grid_rows([0.0, 10.0], [0.0], [0.0, 30.0],
                                              lambda lo, la, al: [1.0, 0.0, 0.0]))
        high = self.write('high.txt', grid_rows([0.0, 20.0], [0.0], [20.0, 100.0],
                                                lambda lo, la, al: [1.0, 0.0, 0.0]))
        result, _ = self.merge(low, high)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('longitude axis', result.stdout + result.stderr)

    def test_a_transition_layer_outside_a_table_is_rejected(self):
        low = self.write('low.txt', grid_rows([0.0], [0.0], [0.0, 10.0],
                                              lambda lo, la, al: [1.0, 0.0, 0.0]))
        high = self.write('high.txt', grid_rows([0.0], [0.0], [20.0, 100.0],
                                                lambda lo, la, al: [1.0, 0.0, 0.0]))
        result, _ = self.merge(low, high)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('transition layer', result.stdout + result.stderr)

    def test_merging_needs_two_tables(self):
        import subprocess
        import sys

        low = self.write('low.txt', grid_rows([0.0], [0.0], [0.0, 30.0],
                                              lambda lo, la, al: [1.0, 0.0, 0.0]))
        result = subprocess.run([sys.executable, self.script(), '--source', 'merge',
                                 '-i', low, '-o', os.path.join(self.directory, 'merged.txt')],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('two input tables', result.stdout + result.stderr)

    def test_a_time_axis_survives_the_merge(self):
        low = self.write('low.txt', grid_rows_time([0.0, 3600.0], [0.0], [0.0], [0.0, 30.0],
                                                   lambda t, lo, la, al: [t/360.0, 0.0, 0.0]),
                         epoch='2024-01-01T00:00:00Z')
        high = self.write('high.txt', grid_rows_time([0.0, 3600.0], [0.0], [0.0], [20.0, 100.0],
                                                     lambda t, lo, la, al: [-t/360.0, 0.0, 0.0]),
                          epoch='2024-01-01T00:00:00Z')
        result, _ = self.merge(low, high)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with quiet():
            wind_dict = wind.initial_settings_wind(table_config(self.directory, 'merged.txt'),
                                                   epoch_of('2024-01-01T00:00:00Z'))
        np.testing.assert_allclose(wind_dict[wind.KEY_TIME], [0.0, 3600.0])
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 10.0, 3600.0)[0], 10.0, atol=1.e-5)
        np.testing.assert_allclose(local(wind_dict, 0.0, 0.0, 100.0, 3600.0)[0], -10.0, atol=1.e-5)


class TestTheHwm14Source(unittest.TestCase):
    """
    HWM14 は Fortran と NRL の係数バイナリなので依存にしない。
    使うときにだけ確かめ、無ければ導入方法を示して止まる。
    """

    def script(self):
        return os.path.join(ROOT_DIR, 'database', 'wind', 'generate_wind_table.py')

    def run_script(self, *argument):
        import subprocess
        import sys

        return subprocess.run([sys.executable, self.script()] + list(argument),
                              capture_output=True, text=True)

    def test_it_needs_a_datetime(self):
        result = self.run_script('--source', 'hwm14', '-o', os.devnull)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--datetime', result.stdout + result.stderr)

    def test_a_missing_hwm14_is_reported_with_the_way_to_get_it(self):
        try:
            import hwm14  # noqa: F401
            self.skipTest('HWM14 is installed in this environment')
        except ImportError:
            pass
        try:
            from pyhwm2014 import hwm14  # noqa: F401
            self.skipTest('HWM14 is installed in this environment')
        except ImportError:
            pass

        result = self.run_script('--source', 'hwm14', '--datetime', '2024-01-01T00:00:00Z',
                                 '-o', os.devnull)
        self.assertNotEqual(result.returncode, 0)
        output = result.stdout + result.stderr
        self.assertIn('HWM14 is not available', output)
        self.assertIn('f2py', output)
        self.assertIn('hwm123114.bin', output)


class TestMultilinearMatchesScipy(unittest.TestCase):
    """
    手書きの多次元線形補間が RegularGridInterpolator とビット一致すること。

    RGI をやめたのは速度のため（3 次元で 1 点 39 us -> 4 us。風のテーブルを読む
    チュートリアルでは計算時間の半分がそこだった）。速度のためだけの置き換えなので、
    検査することは「値が動かないこと」だけになる。空力表の
    test_atmosphere.TestBilinearMatchesScipy と同じ趣旨で、次元数が変わるぶんだけ広い。

    頂点を回る順序・重みの掛ける順序・足し込む順序を RGI の _evaluate_linear に
    合わせてあるので一致する（まとめ方を変えると最下位桁が動き、チュートリアルと
    validation の出力がバイト一致しなくなる ── それが高速化の合格条件だった）。

    節点そのものも見る（格子の内側では区間の選び方によらず同じ値になるはずで、
    そこがずれるのは重みの作り方が違うとき）。

    変異テスト（2026-09-16）: 頂点を回る順序を裏返す・重みの積を右から結合する・
    軸の対応をずらす、はいずれも落ちる。**区間の選び方を bisect_right に変えても
    落ちない**が、これは期待どおりで、節点では一方の重みが 1、他方が 0 になって
    同じ値に行き着く（空力表の evaluate_bilinear と同じ話）。
    """

    def setUp(self):
        self.rng = np.random.default_rng(20260916)

    def reference(self, node, values):
        return scipy.interpolate.RegularGridInterpolator(
            tuple(np.asarray(axis) for axis in node), values,
            method='linear', bounds_error=False, fill_value=None)

    def check(self, node, values, query):
        interpolator = self.reference(node, values)
        for point in query:
            np.testing.assert_array_equal(
                wind.evaluate_multilinear(node, values, point),
                interpolator(np.array([point]))[0])

    def test_it_matches_on_synthetic_grids(self):
        # 時刻・経度・緯度・高度の 4 軸まで。軸が 1 本に縮退した表（鉛直プロファイル）
        # から全部そろった表まで、set_interpolator が作りうる形をすべて通す
        for size in ((37,), (11, 37), (9, 7, 37), (3, 9, 7, 37)):
            node = [sorted(float(value) for value in self.rng.uniform(-180.0, 180.0, number))
                    for number in size]
            values = self.rng.normal(size=tuple(size) + (3,))

            query = [[float(self.rng.uniform(axis[0], axis[-1])) for axis in node]
                     for _ in range(300)]
            # 節点そのもの（軸ごとに別の節点を選ぶ組み合わせも混ぜる）
            query += [[float(axis[self.rng.integers(len(axis))]) for axis in node]
                      for _ in range(300)]

            with self.subTest(size=size):
                self.check(node, values, query)

    def test_it_matches_on_the_shipped_table(self):
        config = load_config(os.path.join(ROOT_DIR, 'tutorial', 'work_reentry_wind_table',
                                          'config.yml'))
        with quiet():
            wind_dict = wind.initial_settings_wind(config)

        interpolator = wind_dict[wind.KEY_INTERP]
        node = interpolator[wind.KEY_NODE]
        values = interpolator[wind.KEY_VALUE]
        self.assertEqual(len(node), 3, '経度 x 緯度 x 高度の表であること')

        query = [[float(self.rng.uniform(axis[0], axis[-1])) for axis in node]
                 for _ in range(300)]
        query += [[float(axis[self.rng.integers(len(axis))]) for axis in node]
                  for _ in range(300)]

        self.check(node, values, query)


if __name__ == '__main__':
    unittest.main()
