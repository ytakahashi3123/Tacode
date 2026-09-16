#!/usr/bin/env python3
"""
src_helper の後処理ツールの設定ファイル（config_helper.yml）のテスト。

長い引数の列を毎回打たずに済ませるための仕組みなので、確かめるのは

  1. 優先順位が **コマンドライン > YAML > 既定値** であること
  2. 綴りを間違えたキーで**止まる**こと（黙って無視されると、設定したつもりで
     効いていない、という一番困る壊れ方になる）
  3. --save-config で書き出したものが、そのまま読み直せること
  4. 3 つのツールがどれも同じ流儀で読むこと（セクション名だけが違う）
"""

import argparse
import os
import subprocess
import sys
import tempfile
import unittest

import yaml

from context import ROOT_DIR

GENERAL_DIR = os.path.join(ROOT_DIR, 'src_helper', 'general')
if GENERAL_DIR not in sys.path:
    sys.path.insert(0, GENERAL_DIR)

import helper_config  # noqa: E402

SECTION = 'sample_tool'

SCRIPT = {'animate_trajectory': os.path.join(ROOT_DIR, 'src_helper', 'animate_trajectory',
                                             'animate_trajectory.py'),
          'montecarlo_dispersion': os.path.join(ROOT_DIR, 'src_helper', 'montecarlo_dispersion',
                                                'montecarlo_dispersion.py'),
          'montecarlo_animation': os.path.join(ROOT_DIR, 'src_helper', 'montecarlo_animation',
                                               'montecarlo_animation.py')}

# 読み込むだけで matplotlib を要求するのはこのツールだけ。他の 2 つは
# 入っていなければ絵を飛ばして正常終了する。
NEEDS_MATPLOTLIB = ('animate_trajectory',)

try:
    import matplotlib  # noqa: F401
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def build_parser():
    """ツールを 1 つ模した、小さな parser。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', nargs='?', default=None)
    parser.add_argument('-o', '--output', default='out.mp4')
    parser.add_argument('--frames', type=int, default=180)
    parser.add_argument('--view', default='flat')
    parser.add_argument('--mark', action='append', default=None)
    helper_config.add_argument(parser)
    return parser


def get_setting(argv, content=None, directory=None):
    """argv と設定ファイルの内容を与えて、まとまった設定を返す。"""
    path = os.path.join(directory, helper_config.FILENAME_DEFAULT) if directory else None
    if content is not None:
        with open(path, 'w') as f:
            yaml.safe_dump(content, f)

    argument = list(argv)
    if path is not None:
        argument = ['-file', path] + argument

    saved = sys.argv
    try:
        sys.argv = ['tool'] + argument
        return helper_config.get_setting(build_parser(), SECTION)
    finally:
        sys.argv = saved


class TestTheOrderOfPrecedence(unittest.TestCase):

    def test_the_defaults_are_used_without_a_file(self):
        argument = get_setting([])
        self.assertEqual(argument.frames, 180)
        self.assertEqual(argument.output, 'out.mp4')
        self.assertIsNone(argument.directory)

    def test_a_missing_file_is_not_an_error(self):
        # 設定ファイルは必須ではない（コマンドラインだけで従来どおり使える）
        saved = sys.argv
        try:
            sys.argv = ['tool', '-file', '/no/such/config.yml', '--frames', '7']
            argument = helper_config.get_setting(build_parser(), SECTION)
        finally:
            sys.argv = saved
        self.assertEqual(argument.frames, 7)

    def test_the_file_is_read(self):
        with tempfile.TemporaryDirectory() as directory:
            argument = get_setting([], {SECTION: {'frames': 30, 'view': '3d',
                                                  'directory': 'work'}}, directory)
            self.assertEqual(argument.frames, 30)
            self.assertEqual(argument.view, '3d')
            self.assertEqual(argument.directory, 'work')

    def test_the_command_line_wins(self):
        with tempfile.TemporaryDirectory() as directory:
            argument = get_setting(['--frames', '5', 'work_cli'],
                                   {SECTION: {'frames': 30, 'directory': 'work_file'}},
                                   directory)
            self.assertEqual(argument.frames, 5)
            self.assertEqual(argument.directory, 'work_cli')

    def test_the_command_line_wins_even_with_the_default_value(self):
        #
        # 既定値と同じ値をコマンドラインで明示したときも、設定ファイルではなく
        # コマンドラインが勝つこと。「既定値と違うかどうか」で判定していた頃は、
        # --frames 180 と書いても設定ファイルの 30 になっていた。
        #
        with tempfile.TemporaryDirectory() as directory:
            argument = get_setting(['--frames', '180', '--view', 'flat'],
                                   {SECTION: {'frames': 30, 'view': '3d'}},
                                   directory)
            self.assertEqual(argument.frames, 180)
            self.assertEqual(argument.view, 'flat')

    def test_another_section_is_ignored(self):
        with tempfile.TemporaryDirectory() as directory:
            argument = get_setting([], {'another_tool': {'frames': 30},
                                        SECTION: {'view': '3d'}}, directory)
            self.assertEqual(argument.frames, 180)
            self.assertEqual(argument.view, '3d')

    def test_a_list_is_read_as_a_list(self):
        # --mark のような繰り返せる引数は、YAML ではリストで書く
        with tempfile.TemporaryDirectory() as directory:
            argument = get_setting([], {SECTION: {'mark': ['a=1.dat', 'b=2.dat']}}, directory)
            self.assertEqual(argument.mark, ['a=1.dat', 'b=2.dat'])


class TestTheFileIsChecked(unittest.TestCase):

    def test_an_unknown_key_stops_the_run(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(SystemExit):
                get_setting([], {SECTION: {'frame': 30}}, directory)

    def test_a_section_which_is_not_a_mapping_stops_the_run(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(SystemExit):
                get_setting([], {SECTION: [1, 2, 3]}, directory)

    def test_a_file_which_is_not_a_mapping_stops_the_run(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, helper_config.FILENAME_DEFAULT)
            with open(path, 'w') as f:
                f.write('- just\n- a list\n')
            with self.assertRaises(SystemExit):
                get_setting([], None, directory)

    def test_an_empty_file_is_taken_as_no_setting(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, helper_config.FILENAME_DEFAULT)
            with open(path, 'w') as f:
                f.write('# nothing here\n')
            self.assertEqual(get_setting([], None, directory).frames, 180)

    def test_a_missing_required_value_stops_the_run(self):
        argument = get_setting([])
        with self.assertRaises(SystemExit):
            helper_config.require(argument, 'directory', SECTION)

    def test_a_given_required_value_is_returned(self):
        argument = get_setting(['work'])
        self.assertEqual(helper_config.require(argument, 'directory', SECTION), 'work')


class TestWritingTheSettings(unittest.TestCase):

    def test_what_is_written_is_read_back(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'written.yml')
            argument = get_setting(['work', '--frames', '12', '--view', '3d',
                                    '--save-config', path])
            self.assertTrue(helper_config.save_file(argument, SECTION))

            with open(path) as f:
                content = yaml.safe_load(f)
            self.assertEqual(content[SECTION]['frames'], 12)
            self.assertEqual(content[SECTION]['view'], '3d')
            self.assertEqual(content[SECTION]['directory'], 'work')
            # 設定の与え方そのものは書き出さない
            self.assertNotIn('file', content[SECTION])
            self.assertNotIn('save_config', content[SECTION])

            # 書いたものをそのまま読み直せること
            saved = sys.argv
            try:
                sys.argv = ['tool', '-file', path]
                again = helper_config.get_setting(build_parser(), SECTION)
            finally:
                sys.argv = saved
            self.assertEqual(again.frames, 12)
            self.assertEqual(again.view, '3d')

    def test_the_other_sections_are_kept(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'written.yml')
            with open(path, 'w') as f:
                yaml.safe_dump({'another_tool': {'frames': 3}}, f)

            argument = get_setting(['--save-config', path])
            helper_config.save_file(argument, SECTION)

            with open(path) as f:
                content = yaml.safe_load(f)
            self.assertEqual(content['another_tool'], {'frames': 3})
            self.assertIn(SECTION, content)

    def test_nothing_is_written_without_the_option(self):
        argument = get_setting([])
        self.assertFalse(helper_config.save_file(argument, SECTION))


class TestTheToolsReadTheSameFile(unittest.TestCase):
    """3 つのツールが同じ流儀で設定を読むこと。セクション名だけが違う。"""

    def run_script(self, name, arguments):
        return subprocess.run([sys.executable, SCRIPT[name]] + arguments,
                              capture_output=True, text=True)

    def check_round_trip(self, name):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'config_helper.yml')
            completed = self.run_script(name, ['--save-config', path])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)

            with open(path) as f:
                content = yaml.safe_load(f)
            self.assertEqual(list(content.keys()), [name])

            # 書き出したものを読ませても通ること（未知のキーが無い）
            completed = self.run_script(name, ['-file', path, '--save-config', path])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_animate_trajectory(self):
        self.check_round_trip('animate_trajectory')

    def test_montecarlo_dispersion(self):
        self.check_round_trip('montecarlo_dispersion')

    def test_montecarlo_animation(self):
        self.check_round_trip('montecarlo_animation')

    def test_a_tool_without_its_input_says_where_to_put_it(self):
        completed = self.run_script('montecarlo_dispersion', [])
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn('montecarlo_dispersion.directory', completed.stdout)

    def test_an_unknown_key_stops_a_tool(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'config_helper.yml')
            with open(path, 'w') as f:
                yaml.safe_dump({'montecarlo_dispersion': {'directry': 'work'}}, f)
            completed = self.run_script('montecarlo_dispersion', ['-file', path])
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn('Unknown setting', completed.stdout)


class TestTheTutorialSettings(unittest.TestCase):
    """同梱の config_helper.yml が、実際にそのツールで読めること。"""

    def test_the_shipped_files_are_read_by_their_tools(self):
        path_list = [os.path.join(ROOT_DIR, 'tutorial', 'work_montecarlo_wind',
                                  'config_helper.yml'),
                     os.path.join(ROOT_DIR, 'tutorial', 'work_reentry_6dof',
                                  'config_helper.yml')]
        for path in path_list:
            self.assertTrue(os.path.exists(path), path)
            with open(path) as f:
                content = yaml.safe_load(f)
            for section in content:
                self.assertIn(section, SCRIPT, path)
                if section in NEEDS_MATPLOTLIB and not HAS_MATPLOTLIB:
                    continue
                with tempfile.TemporaryDirectory() as directory:
                    written = os.path.join(directory, 'config_helper.yml')
                    completed = subprocess.run(
                        [sys.executable, SCRIPT[section], '-file', path,
                         '--save-config', written],
                        capture_output=True, text=True, cwd=os.path.dirname(path))
                    self.assertEqual(completed.returncode, 0,
                                     path + ': ' + completed.stdout + completed.stderr)


if __name__ == '__main__':

    unittest.main()
