#!/usr/bin/env python3
"""
リスタート（computational_setup.flag_initial: False）のテスト。

リスタートは v2.5.1 まで「明示的に止まる」だけだった（CODE_REVIEW A-4）。
実装したのは**最終点からの再開**で、リスタートファイルの履歴は復元しない:

  - ファイルは初期計算からの全履歴を持っているが、大気量（密度・温度・Kn）は
    書かれておらず、restart_process.frequency_output で間引いたファイルでは
    行と反復回数も対応しない。再開後の出力は再開点から始まる
  - 反復回数は 0 に戻し、**経過時間だけを継ぐ**。出力の時刻列は再開時刻から連続し、
    風のテーブルも同じ絶対時刻で引かれる

受け入れ条件は「**通しで計算した後半と、リスタートで計算した続きがビット一致すること**」。
これが成り立つのは、リスタートファイルが ECEF 直交系の状態をそのまま持っており、
再開時に測地系を経由しない（往復変換で最下位桁が動かない）ためで、
ここを壊すと一致しなくなる。

ヘッダは行番号ではなく内容で探す。姿勢・エポック・出力間隔の有無で行数が変わり、
以前の実装は 3 行目を決め打ちしていて 6 自由度のファイルで ValueError になっていた。
"""

import copy
import os
import shutil
import tempfile
import unittest

import numpy as np

from context import ROOT_DIR, load_config, quiet

import atmosphere.atmosphere as atmosphere
import epoch.epoch as epoch_module
import satellite.satellite as satellite
import solver.solver as solver
from orbital.orbital import orbital


def case_config(case, time_elapsed_maximum):
    """チュートリアルのケースを短く切った config を返す。"""
    config = load_config(os.path.join(ROOT_DIR, 'tutorial', case, 'config.yml'))
    config = copy.deepcopy(config)
    config['computational_setup']['time_elapsed_maximum'] = time_elapsed_maximum
    return config


def run_case(config, directory):
    """1 ケース走らせ、Tecplot と restart を directory に書く。"""
    config = copy.deepcopy(config)
    config['post_process']['directory_output'] = os.path.join(directory, 'output_result')
    config['restart_process']['directory_output'] = os.path.join(directory, 'output_restart')

    with quiet():
        atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
        aerodynamic_dict = satellite.initial_settings_satellite(config)
        epoch_dict = epoch_module.initial_settings_epoch(config)

        orb = orbital()
        orb.make_directory_output(config)

        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            orb.initial_settings(config, epoch_dict)
        time_elapsed_initial = time_elapsed

        attitude_dict = orb.initial_settings_attitude(config, coordinate_dict, velocity_dict)

        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            solver.solve_equation_motion(config, iteration, time_elapsed,
                                         coordinate_dict, velocity_dict, trajectory_dict,
                                         atmosphere_dict, aerodynamic_dict, attitude_dict)

        orb.output_tecplot(config, iteration, time_elapsed,
                           coordinate_dict, velocity_dict, trajectory_dict, attitude_dict,
                           epoch_dict, None, time_elapsed_initial)
        orb.output_restart(config, iteration, time_elapsed,
                           coordinate_dict['cartesian'], velocity_dict['cartesian'],
                           attitude_dict, epoch_dict)

    return orb, iteration, coordinate_dict, velocity_dict, attitude_dict


def read_data_lines(path):
    """Tecplot ファイルの数値行を文字列のまま返す（バイト一致を見るため）。"""
    lines = []
    with open(path) as f:
        for line in f:
            if line.startswith(('#', 'Variables', 'zone')) or not line.split():
                continue
            lines.append(line.rstrip('\n'))
    return lines


def write_restart_file(path, body, iteration=10, time_elapsed=10.0,
                       attitude=True, epoch=None, frequency=None):
    """
    ヘッダの並びを変えられるリスタートファイルを書く。

    output_restart が書くのと同じ並びで、姿勢・エポック・出力間隔の行は
    あったり無かったりする（それが行番号で探せない理由）。
    """
    with open(path, 'w') as f:
        f.write('# Restart data (ECEF, cartesian system)\n')
        if epoch is not None:
            f.write('# Epoch (UTC): ' + epoch + '\n')
        if attitude:
            f.write('# X, Y, Z, U, V, W, q0, q1, q2, q3, P, Q, R\n')
            f.write('# --Quaternion: ECEF to body. Angular velocity: body axes,'
                    ' relative to the inertial frame\n')
        if frequency is not None:
            f.write('# Output frequency: ' + str(frequency) + '\n')
        f.write('# Iteration, Elapsed time\n')
        f.write('# ' + str(iteration) + ' ' + str(time_elapsed) + '\n')
        for row in body:
            f.write(' '.join(str(value) for value in row) + '\n')
    return path


class TestRestartFilename(unittest.TestCase):
    """書いた名前と読む名前が一致すること。"""

    def setUp(self):
        self.orb = orbital()
        self.config = {'restart_process': {'directory_output': 'output_restart',
                                           'file_restart': 'restart.dat',
                                           'flag_time_series': False,
                                           'digid_step': 4,
                                           'restart_step': 500}}

    def test_a_single_file_is_overwritten(self):
        self.assertEqual(self.orb.get_restart_filename(self.config, 20),
                         'output_restart/restart.dat')

    def test_a_time_series_embeds_the_step(self):
        self.config['restart_process']['flag_time_series'] = True
        # 書く側（反復回数）と読む側（restart_step）が同じ規則であること。
        # かつては書き出しが restart_0500.dat、読み取りが restart_s0500.dat を探していた
        self.assertEqual(self.orb.get_restart_filename(self.config, 500),
                         'output_restart/restart_0500.dat')
        self.assertEqual(
            self.orb.get_restart_filename(self.config, 500),
            self.orb.get_restart_filename(self.config,
                                          self.config['restart_process']['restart_step']))


class TestRestartHeaderIsFoundByKey(unittest.TestCase):
    """ヘッダの行数が変わっても反復回数と経過時間を拾えること。"""

    def setUp(self):
        self.directory = tempfile.mkdtemp(prefix='tacode_restart_')
        self.addCleanup(shutil.rmtree, self.directory, True)
        self.orb = orbital()
        self.config = {'restart_process': {'directory_output': self.directory,
                                           'file_restart': 'restart.dat',
                                           'flag_time_series': False,
                                           'digid_step': 4,
                                           'restart_step': 0}}
        self.path = os.path.join(self.directory, 'restart.dat')

    def read(self, **kwargs):
        with quiet():
            return self.orb.read_restart(self.config)

    def test_a_3dof_file(self):
        write_restart_file(self.path, [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
                           iteration=7, time_elapsed=3.5, attitude=False)
        state = self.read()
        self.assertEqual(state[self.orb.KEY_RESTART_ITERATION], 7)
        self.assertEqual(state[self.orb.KEY_RESTART_TIME], 3.5)
        np.testing.assert_array_equal(state[self.orb.KEY_RESTART_COORDINATE], [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(state[self.orb.KEY_RESTART_VELOCITY], [4.0, 5.0, 6.0])
        self.assertIsNone(state[self.orb.KEY_RESTART_QUATERNION])

    def test_a_6dof_file_with_every_optional_header_line(self):
        # 姿勢・エポック・出力間隔が全部あると、反復回数の行は 6 行目に来る。
        # 3 行目を決め打ちしていた実装はここで ValueError になっていた
        row = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3]
        write_restart_file(self.path, [row], iteration=20000, time_elapsed=1000.0,
                           attitude=True, epoch='2015-01-01T01:30:00+00:00', frequency=20)
        state = self.read()
        self.assertEqual(state[self.orb.KEY_RESTART_ITERATION], 20000)
        self.assertEqual(state[self.orb.KEY_RESTART_TIME], 1000.0)
        np.testing.assert_array_equal(state[self.orb.KEY_RESTART_QUATERNION], [1.0, 0.0, 0.0, 0.0])
        np.testing.assert_array_equal(state[self.orb.KEY_RESTART_OMEGA], [0.1, 0.2, 0.3])
        self.assertEqual(state[self.orb.KEY_RESTART_EPOCH], '2015-01-01T01:30:00+00:00')

    def test_the_last_line_is_the_resume_point(self):
        body = [[float(n), 0.0, 0.0, 0.0, 0.0, 0.0] for n in range(0, 5)]
        write_restart_file(self.path, body, iteration=4, time_elapsed=4.0, attitude=False)
        state = self.read()
        self.assertEqual(state[self.orb.KEY_RESTART_COORDINATE][0], 4.0)

    def test_a_missing_file_stops(self):
        with self.assertRaises(SystemExit) as raised:
            with quiet():
                self.orb.read_restart(self.config)
        self.assertEqual(raised.exception.code, 1)

    def test_a_file_without_the_iteration_line_stops(self):
        with open(self.path, 'w') as f:
            f.write('# Restart data (ECEF, cartesian system)\n')
            f.write('1.0 2.0 3.0 4.0 5.0 6.0\n')
        with self.assertRaises(SystemExit) as raised:
            self.read()
        self.assertEqual(raised.exception.code, 1)

    def test_a_file_without_data_stops(self):
        write_restart_file(self.path, [], attitude=False)
        with self.assertRaises(SystemExit) as raised:
            self.read()
        self.assertEqual(raised.exception.code, 1)

    def test_an_unexpected_number_of_columns_stops(self):
        write_restart_file(self.path, [[1.0, 2.0, 3.0, 4.0]], attitude=False)
        with self.assertRaises(SystemExit) as raised:
            self.read()
        self.assertEqual(raised.exception.code, 1)

    def test_a_quaternion_off_the_unit_sphere_stops(self):
        # 正規化して黙って進めない。発散した計算の残骸を「続き」にしないため
        row = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        write_restart_file(self.path, [row])
        with self.assertRaises(SystemExit) as raised:
            self.read()
        self.assertEqual(raised.exception.code, 1)


class TestRestartEpoch(unittest.TestCase):
    """経過時間はエポックからの秒なので、エポックの食い違いは止める。"""

    def setUp(self):
        self.directory = tempfile.mkdtemp(prefix='tacode_restart_')
        self.addCleanup(shutil.rmtree, self.directory, True)
        self.orb = orbital()
        self.config = {'restart_process': {'directory_output': self.directory,
                                           'file_restart': 'restart.dat',
                                           'flag_time_series': False,
                                           'digid_step': 4,
                                           'restart_step': 0}}
        self.path = os.path.join(self.directory, 'restart.dat')

    def epoch_dict(self, string_datetime):
        with quiet():
            return epoch_module.initial_settings_epoch(
                {'epoch': {'flag_epoch': True, 'datetime': string_datetime}})

    def test_the_same_epoch_passes(self):
        epoch_dict = self.epoch_dict('2015-01-01T01:30:00Z')
        write_restart_file(self.path, [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], attitude=False,
                           epoch=epoch_module.get_string(epoch_dict, 0.0))
        with quiet():
            state = self.orb.read_restart(self.config, epoch_dict)
        self.assertIsNotNone(state[self.orb.KEY_RESTART_EPOCH])

    def test_a_different_epoch_stops(self):
        write_restart_file(self.path, [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], attitude=False,
                           epoch='2015-01-01T01:30:00+00:00')
        with self.assertRaises(SystemExit) as raised:
            with quiet():
                self.orb.read_restart(self.config, self.epoch_dict('2016-06-01T00:00:00Z'))
        self.assertEqual(raised.exception.code, 1)

    def test_a_missing_epoch_only_warns(self):
        write_restart_file(self.path, [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], attitude=False)
        with quiet() as output:
            self.orb.read_restart(self.config, self.epoch_dict('2015-01-01T01:30:00Z'))
        self.assertIn('Warning', output.getvalue())


class TestRestartContinuesTheTrajectory(unittest.TestCase):
    """通しで計算した後半と、リスタートで計算した続きが一致すること。"""

    def continue_and_compare(self, case, time_restart, time_end):
        directory = tempfile.mkdtemp(prefix='tacode_restart_')
        self.addCleanup(shutil.rmtree, directory, True)

        directory_short = os.path.join(directory, 'short')
        directory_long = os.path.join(directory, 'long')
        directory_continued = os.path.join(directory, 'continued')

        run_case(case_config(case, time_restart), directory_short)
        run_case(case_config(case, time_end), directory_long)

        config_continued = case_config(case, time_end)
        config_continued['computational_setup']['flag_initial'] = False
        # 短い計算が書いた restart をそのまま入力にする
        os.makedirs(directory_continued)
        shutil.copytree(os.path.join(directory_short, 'output_restart'),
                        os.path.join(directory_continued, 'output_restart'))
        run_case(config_continued, directory_continued)

        lines_long = read_data_lines(os.path.join(directory_long, 'output_result', 'tecplot.dat'))
        lines_continued = read_data_lines(
            os.path.join(directory_continued, 'output_result', 'tecplot.dat'))

        # 時刻列は再開時刻から連続する（0 に戻らない）
        self.assertEqual(float(lines_continued[0].split()[0]), time_restart)

        tail = [line for line in lines_long if float(line.split()[0]) >= time_restart]
        self.assertEqual(len(tail), len(lines_continued))
        self.assertGreater(len(tail), 1)
        for index in range(0, len(tail)):
            self.assertEqual(tail[index], lines_continued[index],
                             '%s: 行 %d が通し計算と一致しない' % (case, index))

    def test_three_degrees_of_freedom(self):
        self.continue_and_compare('work', 20.0, 40.0)

    def test_six_degrees_of_freedom(self):
        # 姿勢（クォータニオンと角速度）も再開点から続くこと
        self.continue_and_compare('work_reentry_6dof', 2.0, 4.0)


class TestRestartOfTheAttitude(unittest.TestCase):
    """姿勢の有無の食い違い。"""

    def setUp(self):
        self.directory = tempfile.mkdtemp(prefix='tacode_restart_')
        self.addCleanup(shutil.rmtree, self.directory, True)

    def test_a_3dof_file_cannot_restart_a_6dof_run(self):
        config = case_config('work_reentry_6dof', 1.0)
        config['attitude']['flag_attitude'] = False
        run_case(config, self.directory)

        config_continued = case_config('work_reentry_6dof', 2.0)
        config_continued['computational_setup']['flag_initial'] = False
        config_continued['restart_process']['directory_output'] = \
            os.path.join(self.directory, 'output_restart')

        with self.assertRaises(SystemExit) as raised:
            with quiet():
                orb = orbital()
                iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                    orb.initial_settings(config_continued)
                orb.initial_settings_attitude(config_continued, coordinate_dict, velocity_dict)
        self.assertEqual(raised.exception.code, 1)

    def test_a_6dof_file_read_by_a_3dof_run_warns(self):
        run_case(case_config('work_reentry_6dof', 1.0), self.directory)

        config_continued = case_config('work_reentry_6dof', 2.0)
        config_continued['computational_setup']['flag_initial'] = False
        config_continued['attitude']['flag_attitude'] = False
        config_continued['restart_process']['directory_output'] = \
            os.path.join(self.directory, 'output_restart')

        with quiet() as output:
            orb = orbital()
            orb.initial_settings(config_continued)
        self.assertIn('Warning', output.getvalue())


if __name__ == '__main__':
    unittest.main()
