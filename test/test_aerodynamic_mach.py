#!/usr/bin/env python3
"""
3 自由度の CD を Mach 数でも変える空力表（"MACH <値>" でブロックを分けた表）のテスト。

Mach 軸の無い表は従来どおり Kn だけで CD を引き、出力も 1 バイトも変わらない（それは
参照出力との突き合わせ test_regression_3dof.py が見ている）。ここで確かめるのは:

  1. 表の読み取り。昇順への並べ替え、AOA との混在・Mach 1 点だけ・Mach の前のデータ行・
     重複した Mach は止める
  2. (Mach, Kn) の双一次。節点では表の値、Mach 方向に線形、両軸とも端でクランプ
  3. 音速 a = sqrt(gamma R_u T / M)。国際標準大気の海面の 340.294 m/s と照らす。
     比熱比・分子量が無ければ止め、6 自由度に Mach の表を渡しても止める
  4. **恒等式**: どの Mach でも同じ CD(Kn) を並べた表は、Kn だけの表と同じ軌道を与える。
     Mach 数の作り方（速さ・温度・音速の係数）を間違えても、この恒等式は落ちない。
     そこは 5 が見る
  5. 出力の Mach 列が、出力の速さと温度から作った V/a に一致し、かつ CD が本当に
     Mach で変わっている（Mach で CD を変えた表では軌道が動く。動く向きも見る）
  6. RK4 の各段が自分の速度で Mach 数を作っている（4 次の収束で見る）
"""

import copy
import os
import shutil
import tempfile
import unittest

import numpy as np

from context import ROOT_DIR, load_config, quiet

import atmosphere.atmosphere as atmosphere
import satellite.satellite as satellite
import solver.solver as solver
import wind.wind as wind
from orbital.orbital import orbital


CONFIG_REENTRY = os.path.join(ROOT_DIR, 'tutorial', 'work_reentry', 'config.yml')

# 空気（国際標準大気の海面で音速 340.294 m/s）
GAMMA_AIR  = 1.4
WEIGHT_AIR = 28.9644

ROW_TAIL = '0 0 0 0 0  0 0 0 0 0 0  0'   # CFy..CMz (5), SDV_CFx..SDV_CMz (6), Altitude


def write_table(path, blocks, header='MACH'):
    """blocks は [(見出しの値, [(Kn, CD), ...]), ...]。"""
    with open(path, 'w') as f:
        f.write('Test table\n')
        for value, rows in blocks:
            if value is not None:
                f.write('{} {}\n'.format(header, value))
            for knudsen, cd in rows:
                f.write('{!r} {!r} {}\n'.format(float(knudsen), float(cd), ROW_TAIL))


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


class TemporaryDirectoryMixin:

    def setUp(self):
        self.directory = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.directory, ignore_errors=True)

    def path(self, name='table.txt'):
        return os.path.join(self.directory, name)


class TestTheMachTableIsRead(TemporaryDirectoryMixin, unittest.TestCase):

    KNUDSEN = [1.e-4, 1.e-2, 1.0]

    def test_the_blocks_are_sorted_by_the_mach_number(self):
        write_table(self.path(), [(10.0, [(k, 1.7) for k in self.KNUDSEN]),
                                  (2.0,  [(k, 1.5) for k in self.KNUDSEN]),
                                  (5.0,  [(k, 1.6) for k in self.KNUDSEN])])
        angle, block, _, mach = satellite.parse_aerodynamic_file(self.path())
        np.testing.assert_array_equal(mach, [2.0, 5.0, 10.0])
        np.testing.assert_array_equal(angle, [0.0])
        self.assertEqual([b[0, 1] for b in block], [1.5, 1.6, 1.7])

    def test_a_table_without_mach_lines_has_no_mach_axis(self):
        write_table(self.path(), [(None, [(k, 1.5) for k in self.KNUDSEN])])
        _, _, _, mach = satellite.parse_aerodynamic_file(self.path())
        self.assertIsNone(mach)

    def test_mixing_aoa_and_mach_stops(self):
        with open(self.path(), 'w') as f:
            f.write('AOA 0\n1e-4 1.5 {}\nMACH 2\n1e-4 1.5 {}\n'.format(ROW_TAIL, ROW_TAIL))
        with quiet(), self.assertRaises(SystemExit):
            satellite.parse_aerodynamic_file(self.path())

    def test_data_before_the_first_mach_line_stops(self):
        with open(self.path(), 'w') as f:
            f.write('1e-4 1.5 {}\nMACH 2\n1e-4 1.5 {}\nMACH 3\n1e-4 1.5 {}\n'.format(ROW_TAIL, ROW_TAIL, ROW_TAIL))
        with quiet(), self.assertRaises(SystemExit):
            satellite.parse_aerodynamic_file(self.path())

    def test_a_single_mach_block_stops(self):
        write_table(self.path(), [(2.0, [(k, 1.5) for k in self.KNUDSEN])])
        with quiet(), self.assertRaises(SystemExit):
            satellite.parse_aerodynamic_file(self.path())

    def test_a_duplicated_mach_number_stops(self):
        write_table(self.path(), [(2.0, [(k, 1.5) for k in self.KNUDSEN]),
                                  (2.0, [(k, 1.6) for k in self.KNUDSEN])])
        with quiet(), self.assertRaises(SystemExit):
            satellite.parse_aerodynamic_file(self.path())

    def test_an_empty_mach_block_stops(self):
        with open(self.path(), 'w') as f:
            f.write('MACH 2\nMACH 3\n1e-4 1.5 {}\n'.format(ROW_TAIL))
        with quiet(), self.assertRaises(SystemExit):
            satellite.parse_aerodynamic_file(self.path())

    def test_blocks_with_different_knudsen_numbers_stop(self):
        write_table(self.path(), [(2.0, [(1.e-4, 1.5), (1.0, 1.9)]),
                                  (5.0, [(1.e-3, 1.6), (1.0, 1.9)])])
        config = {'satellite': {'directory_path_specify': 'manual', 'directory_aerodynamic': self.directory,
                                'filename_aerodynamic': 'table.txt', 'characteristic_length': 1.0}}
        with quiet(), self.assertRaises(SystemExit):
            satellite.read_aerodynamic_file(config)


class TestTheMachDrag(TemporaryDirectoryMixin, unittest.TestCase):

    MACH    = [2.0, 5.0, 10.0]
    KNUDSEN = [1.e-4, 1.e-2, 1.0]
    # CD[i, j] = Mach i, Kn j。どちらの向きにも変わる値にして、軸の取り違えを落とす
    CD = [[1.50, 1.55, 1.90],
          [1.60, 1.66, 1.95],
          [1.70, 1.74, 2.00]]

    def setUp(self):
        super().setUp()
        write_table(self.path(), [(m, list(zip(self.KNUDSEN, row))) for m, row in zip(self.MACH, self.CD)])
        config = {'satellite': {'directory_path_specify': 'manual', 'directory_aerodynamic': self.directory,
                                'filename_aerodynamic': 'table.txt', 'characteristic_length': 1.0,
                                'kind_aerodynamic_model': 'fileread'},
                  'atmosphere': {'specific_heat_ratio': GAMMA_AIR, 'molecular_weight': WEIGHT_AIR}}
        with quiet():
            self.aerodynamic_dict = satellite.initial_settings_satellite(config)

    def cd(self, knudsen, mach):
        return satellite.get_aerodynamic_coefficient_mach(knudsen, mach, self.aerodynamic_dict)

    def test_the_nodes_give_the_table(self):
        for i, mach in enumerate(self.MACH):
            for j, knudsen in enumerate(self.KNUDSEN):
                self.assertEqual(self.cd(knudsen, mach), self.CD[i][j])

    def test_it_is_linear_in_the_mach_number(self):
        for fraction in (0.25, 0.5, 0.8):
            mach = 5.0 + fraction*(10.0 - 5.0)
            expected = self.CD[1][1] + fraction*(self.CD[2][1] - self.CD[1][1])
            self.assertAlmostEqual(self.cd(1.e-2, mach), expected, places=14)

    def test_it_is_linear_in_the_knudsen_number_as_the_old_table(self):
        # Kn 方向は Kn だけの CD（np.interp、Kn について線形）と同じ扱い
        knudsen = 0.3
        expected = np.interp(knudsen, self.KNUDSEN, self.CD[0])
        self.assertAlmostEqual(self.cd(knudsen, 2.0), expected, places=14)

    def test_both_axes_are_clamped(self):
        self.assertEqual(self.cd(1.e-2, 0.5), self.CD[0][1])
        self.assertEqual(self.cd(1.e-2, 40.0), self.CD[2][1])
        self.assertEqual(self.cd(1.e-8, 5.0), self.CD[1][0])
        self.assertEqual(self.cd(1.e+3, 5.0), self.CD[1][2])
        self.assertEqual(self.cd(1.e+3, 40.0), self.CD[2][2])

    def test_the_speed_of_sound_is_that_of_the_standard_atmosphere(self):
        # 国際標準大気の海面: 288.15 K で 340.294 m/s（gamma 1.4、M 28.9644 g/mol）。
        # ISA は気体定数 8.31432 で作られていて、CODATA 2018 の 8.314462618 と 1.7e-5 違う。
        # 音速ではその半分の 9e-6 がずれとして残る
        mach = satellite.get_mach_number(340.294, 288.15, self.aerodynamic_dict)
        self.assertAlmostEqual(mach, 1.0, delta=1.e-5)
        mach_isa = mach*np.sqrt(satellite.CONSTANT_GAS_UNIVERSAL/8.31432)
        self.assertAlmostEqual(mach_isa, 1.0, delta=2.e-6)

    def test_the_mach_number_goes_as_one_over_the_square_root_of_temperature(self):
        ratio = satellite.get_mach_number(1000.0, 200.0, self.aerodynamic_dict) \
              / satellite.get_mach_number(1000.0, 800.0, self.aerodynamic_dict)
        self.assertAlmostEqual(ratio, 2.0, places=14)

    def test_the_old_keys_are_not_left_for_the_6dof_functions(self):
        # 6 自由度の関数が Mach 軸を黙って捨てた値を返せないように、格子を置いていない
        self.assertIsNone(self.aerodynamic_dict[satellite.KEY_CD_MEAN])
        self.assertIsNone(self.aerodynamic_dict[satellite.KEY_CF])


class TestTheMachTableNeedsItsSettings(TemporaryDirectoryMixin, unittest.TestCase):

    def config(self, atmosphere_section, flag_attitude=False):
        write_table(self.path(), [(2.0, [(1.e-4, 1.5), (1.0, 1.9)]),
                                  (5.0, [(1.e-4, 1.6), (1.0, 1.9)])])
        return {'satellite': {'directory_path_specify': 'manual', 'directory_aerodynamic': self.directory,
                              'filename_aerodynamic': 'table.txt', 'characteristic_length': 1.0,
                              'kind_aerodynamic_model': 'fileread'},
                'atmosphere': atmosphere_section,
                'attitude': {'flag_attitude': flag_attitude}}

    def test_a_missing_specific_heat_ratio_stops(self):
        with quiet(), self.assertRaises(SystemExit):
            satellite.initial_settings_satellite(self.config({'molecular_weight': WEIGHT_AIR}))

    def test_a_missing_molecular_weight_stops(self):
        with quiet(), self.assertRaises(SystemExit):
            satellite.initial_settings_satellite(self.config({'specific_heat_ratio': GAMMA_AIR}))

    def test_an_unphysical_specific_heat_ratio_stops(self):
        with quiet(), self.assertRaises(SystemExit):
            satellite.initial_settings_satellite(self.config({'specific_heat_ratio': 1.0,
                                                              'molecular_weight': WEIGHT_AIR}))

    def test_the_6dof_refuses_a_mach_table(self):
        with quiet(), self.assertRaises(SystemExit):
            satellite.initial_settings_satellite(self.config({'specific_heat_ratio': GAMMA_AIR,
                                                              'molecular_weight': WEIGHT_AIR},
                                                             flag_attitude=True))

    def test_a_table_without_mach_lines_does_not_ask_for_them(self):
        write_table(self.path(), [(None, [(1.e-4, 1.5), (1.0, 1.9)])])
        config = self.config({})
        write_table(self.path(), [(None, [(1.e-4, 1.5), (1.0, 1.9)])])
        with quiet():
            aerodynamic_dict = satellite.initial_settings_satellite(config)
        self.assertIsNone(aerodynamic_dict[satellite.KEY_MACH])


class TestTheSolverUsesTheAirSpeed(TemporaryDirectoryMixin, unittest.TestCase):
    """solver.get_drag_coefficient は風があれば対気速度で Mach 数を作る。"""

    def test_the_air_relative_velocity_is_used(self):
        write_table(self.path(), [(1.0, [(1.e-4, 1.0), (1.0, 1.0)]),
                                  (3.0, [(1.e-4, 3.0), (1.0, 3.0)])])
        config = {'satellite': {'directory_path_specify': 'manual', 'directory_aerodynamic': self.directory,
                                'filename_aerodynamic': 'table.txt', 'characteristic_length': 1.0,
                                'kind_aerodynamic_model': 'fileread'},
                  'atmosphere': {'specific_heat_ratio': GAMMA_AIR, 'molecular_weight': WEIGHT_AIR}}
        with quiet():
            aerodynamic_dict = satellite.initial_settings_satellite(config)
        temperature = 288.15
        sound = 340.294
        velocity = np.array([2.0*sound, 0.0, 0.0])
        # 風が無ければ ECEF 速度の Mach 2 で CD 2
        cd_calm = solver.get_drag_coefficient(1.e-2, temperature, velocity, None, True,
                                              None, None, aerodynamic_dict)
        self.assertAlmostEqual(cd_calm, 2.0, places=4)
        # 追い風で対気速度が Mach 1 まで落ちれば CD 1
        velocity_air = velocity - np.array([sound, 0.0, 0.0])
        cd_wind = solver.get_drag_coefficient(1.e-2, temperature, velocity, velocity_air, True,
                                              None, None, aerodynamic_dict)
        self.assertAlmostEqual(cd_wind, 1.0, places=4)


def run_reentry(table_path, time_max=600.0, output_directory=None, timestep=None):
    """tutorial/work_reentry を空力表だけ差し替えて 3 自由度で回す。"""
    config = copy.deepcopy(load_config(CONFIG_REENTRY))
    config['computational_setup']['time_elapsed_maximum'] = time_max
    if timestep is not None:
        config['time_integration']['timestep_constant'] = timestep
    if table_path is not None:
        config['satellite']['directory_aerodynamic'] = os.path.dirname(table_path)
        config['satellite']['filename_aerodynamic'] = os.path.basename(table_path)
    config['atmosphere']['specific_heat_ratio'] = GAMMA_AIR
    config['atmosphere']['molecular_weight'] = WEIGHT_AIR
    if output_directory is not None:
        config['post_process']['directory_output'] = output_directory
        config['post_process']['kml']['flag_output'] = False
        config['post_process']['tecplot']['flag_output'] = True
        config['post_process']['tecplot']['frequency_output'] = 1

    with quiet():
        atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
        aerodynamic_dict = satellite.initial_settings_satellite(config)
        wind_dict = wind.initial_settings_wind(config)
        orb = orbital()
        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = orb.initial_settings(config)
        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            solver.solve_equation_motion(config, iteration, time_elapsed,
                                         coordinate_dict, velocity_dict, trajectory_dict,
                                         atmosphere_dict, aerodynamic_dict, None, wind_dict)
        if output_directory is not None:
            orb.output_tecplot(config, iteration, time_elapsed, coordinate_dict, velocity_dict,
                               trajectory_dict, None, None, wind_dict, 0.0, aerodynamic_dict)
    return config, np.array(coordinate_dict['cartesian']), np.array(velocity_dict['cartesian'])


def read_reentry_table():
    """tutorial/work_reentry の Kn だけの表 (Kn, CD)。"""
    config = load_config(CONFIG_REENTRY)
    with quiet():
        aerodynamic_dict = satellite.initial_settings_satellite(config)
    return aerodynamic_dict[satellite.KEY_KN], aerodynamic_dict[satellite.KEY_CD_MEAN]


class TestAMachIndependentTableGivesTheSameTrajectory(TemporaryDirectoryMixin, unittest.TestCase):
    """
    恒等式: どの Mach でも同じ CD(Kn) を並べた Mach の表は、Kn だけの表と同じ軌道を与える。

    ビット一致はしない（双一次は Mach 方向の重み w と 1-w で同じ値を足し直すので、
    最下位桁が動きうる）。600 s の再突入で位置の相対差 1e-12 以内を見る。
    """

    def test_the_trajectories_agree(self):
        knudsen, cd = read_reentry_table()
        write_table(self.path(), [(mach, list(zip(knudsen, cd))) for mach in (0.5, 3.0, 12.0, 40.0)])
        _, coordinate_old, velocity_old = run_reentry(None)
        _, coordinate_new, velocity_new = run_reentry(self.path())
        self.assertEqual(coordinate_old.shape, coordinate_new.shape)
        scale = np.max(np.abs(coordinate_old))
        self.assertLess(np.max(np.abs(coordinate_new - coordinate_old))/scale, 1.e-12)
        scale = np.max(np.abs(velocity_old))
        self.assertLess(np.max(np.abs(velocity_new - velocity_old))/scale, 1.e-12)


class TestEachStageUsesItsOwnMachNumber(TemporaryDirectoryMixin, unittest.TestCase):
    """
    RK4 の各段は、その段の仮想速度で Mach 数を作る。1 段目の速度を使い回しても軌道は
    ほとんど変わらないが、刻みを半分にしたときの差の縮み方が 16 倍（4 次）から崩れる
    （実測: 正しい実装で 16.0 / 16.3、使い回すと 8.0 / 1.4）。

    CD が Mach 数に線形で、Kn にも Mach の範囲にも折れ目の無い表を使う。区分線形の
    折れ目をまたぐと、それだけで次数が落ちるため。
    """

    def test_the_scheme_stays_fourth_order(self):
        write_table(self.path(), [(0.0,  [(1.e-6, 1.0), (1.e+6, 1.0)]),
                                  (60.0, [(1.e-6, 4.0), (1.e+6, 4.0)])])
        position = {}
        for timestep in (2.0, 1.0, 0.5):
            _, coordinate, _ = run_reentry(self.path(), time_max=400.0, timestep=timestep)
            position[timestep] = coordinate[-1]
        ratio = np.linalg.norm(position[2.0] - position[1.0]) / np.linalg.norm(position[1.0] - position[0.5])
        self.assertGreater(ratio, 12.0)
        self.assertLess(ratio, 20.0)


class TestTheMachColumnAndTheMachDrag(unittest.TestCase):
    """Mach で CD を変えた表: 出力の Mach 列と、CD が Mach で本当に効いていること。"""

    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.mkdtemp()
        knudsen, cd = read_reentry_table()
        # 超音速より下では CD を 2 倍にする。極超音速（Mach 5 以上）は元の表のまま
        blocks = [(mach, list(zip(knudsen, cd*factor))) for mach, factor in ((1.0, 2.0), (2.0, 2.0),
                                                                             (5.0, 1.0), (40.0, 1.0))]
        cls.table = os.path.join(cls.directory, 'table_mach.txt')
        write_table(cls.table, blocks)
        cls.table_flat = os.path.join(cls.directory, 'table_flat.txt')
        write_table(cls.table_flat, [(mach, list(zip(knudsen, cd))) for mach in (1.0, 40.0)])

        cls.output = os.path.join(cls.directory, 'output')
        os.makedirs(cls.output)
        cls.config, cls.coordinate, cls.velocity = run_reentry(cls.table, time_max=2000.0,
                                                               output_directory=cls.output)
        _, cls.coordinate_flat, cls.velocity_flat = run_reentry(cls.table_flat, time_max=2000.0)
        cls.variables, cls.rows = read_tecplot(os.path.join(cls.output, 'tecplot.dat'))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.directory, ignore_errors=True)

    def column(self, name):
        return self.rows[:, self.variables.index(name)]

    def test_the_mach_column_is_written_after_the_knudsen_number(self):
        self.assertEqual(self.variables.index('Mach'), self.variables.index('Kn') + 1)

    def test_the_mach_column_is_the_speed_over_the_speed_of_sound(self):
        sound = np.sqrt(GAMMA_AIR*satellite.CONSTANT_GAS_UNIVERSAL/(WEIGHT_AIR*1.e-3)*self.column('Temp[K]'))
        np.testing.assert_allclose(self.column('Mach'), self.column('VelplAbs[m/s]')/sound, rtol=1.e-12)

    def test_the_flight_crosses_the_mach_range_of_the_table(self):
        mach = self.column('Mach')
        self.assertGreater(mach.max(), 5.0)
        self.assertLess(mach.min(), 2.0)

    def test_the_hypersonic_part_is_not_touched(self):
        # Mach 5 より上では両方の表の CD が同じなので、そこを抜けるまでは同じ軌道
        mach = self.column('Mach')
        index = int(np.argmax(mach < 5.0)) - 1
        self.assertGreater(index, 10)
        scale = np.max(np.abs(self.coordinate_flat))
        self.assertLess(np.max(np.abs(self.coordinate[:index] - self.coordinate_flat[:index]))/scale, 1.e-12)

    def test_a_larger_supersonic_drag_slows_the_vehicle_down_sooner(self):
        # 超音速で CD を 2 倍にしたほうが、同じ時刻の速さは小さい
        speed      = np.linalg.norm(self.velocity, axis=1)
        speed_flat = np.linalg.norm(self.velocity_flat, axis=1)
        mach = self.column('Mach')
        index = int(np.argmax(mach < 1.5))
        count = min(len(speed), len(speed_flat))
        self.assertTrue(np.all(speed[index:count] < speed_flat[index:count]))

    def test_a_table_without_mach_lines_writes_no_mach_column(self):
        directory = tempfile.mkdtemp()
        try:
            run_reentry(None, time_max=10.0, output_directory=directory)
            variables, _ = read_tecplot(os.path.join(directory, 'tecplot.dat'))
        finally:
            shutil.rmtree(directory, ignore_errors=True)
        self.assertNotIn('Mach', variables)


if __name__ == '__main__':
    unittest.main()
