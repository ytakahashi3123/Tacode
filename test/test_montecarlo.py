#!/usr/bin/env python3
"""
モンテカルロ（src/tacode-montecarlo.py）と、その結果をまとめるヘルパーのテスト。

モンテカルロはテンプレートの制御ファイルを **テキストとして** 書き換える。
そこで確かめるのは 3 点:

  1. 書き換えが狙った 1 か所だけに効くこと。キー名はセクションをまたいで重複するので
     （initial_settings.velocity と wind.velocity）、セクションで絞れないと
     風のつもりで初期速度を動かしてしまう。値の同じ行が並ぶとき（[0.0, 0.0, 0.0]）に
     巻き添えで一緒に書き換わらないことも見る
  2. 風のチュートリアル（tutorial/work_montecarlo_wind と tutorial_template_wind）が
     互いに整合していること。ばらつかせる基準値はモンテカルロ側の config から読まれ、
     書き換えられるのはテンプレート側の行なので、両者がずれていると黙って別の値で走る
  3. 短くしたケースを実際に 2 つ走らせ、ケースごとに違う風で違う軌道が出ること
  4. 走り終わったケースを 1 つの Tecplot ファイルにまとめる postprocess

src_helper/montecarlo_dispersion のほうは、終端点の分解（東・北）と、
テンプレートとの差分からばらついた入力を拾う経路を見る。
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import yaml

try:
    import matplotlib  # noqa: F401
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from context import ROOT_DIR, load_config

from montecarlo.montecarlo import montecarlo

HELPER_DIR = os.path.join(ROOT_DIR, 'src_helper', 'montecarlo_dispersion')
SCRIPT = os.path.join(HELPER_DIR, 'montecarlo_dispersion.py')
if HELPER_DIR not in sys.path:
    sys.path.insert(0, HELPER_DIR)

import montecarlo_dispersion  # noqa: E402

ANIMATION_DIR = os.path.join(ROOT_DIR, 'src_helper', 'montecarlo_animation')
SCRIPT_ANIMATION = os.path.join(ANIMATION_DIR, 'montecarlo_animation.py')
if ANIMATION_DIR not in sys.path:
    sys.path.insert(0, ANIMATION_DIR)

import montecarlo_animation  # noqa: E402

CONFIG_MONTECARLO_WIND = os.path.join(ROOT_DIR, 'tutorial', 'work_montecarlo_wind', 'config.yml')
TEMPLATE_WIND = os.path.join(ROOT_DIR, 'tutorial_template_wind')
CONFIG_TEMPLATE_WIND = os.path.join(TEMPLATE_WIND, 'config.yml')


def rewrite(text, key, values, section=None):
    """制御ファイルの断片を書き換え、書き換え後のテキストを返す。"""
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, 'config.yml')
        with open(path, 'w') as f:
            f.write(text)
        montecarlo().rewrite_control(path, key, len(values), values, section)
        with open(path) as f:
            return f.read()


TEXT_TWO_SECTIONS = """initial_settings:
  velocity:
    - 7450
    - 0
    - 0

wind:
  # [East, North, Up], m/s
  velocity:
    - 0.0
    - 0.0
    - 0.0
"""


TEXT_SECTION_WITHOUT_THE_KEY = """computational_setup:
  timestep_constant: 0.5
  time_elapsed_maximum: 5000.0

initial_settings:
  velocity:
    - 7450
    - 0
    - 0
"""


def rewrite_expecting_a_stop(text, key, values, section=None):
    """書き換えが停止することを確かめ、停止したあとのファイルの中身を返す。"""
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, 'config.yml')
        with open(path, 'w') as f:
            f.write(text)
        try:
            montecarlo().rewrite_control(path, key, len(values), values, section)
        except SystemExit:
            with open(path) as f:
                return f.read()
        raise AssertionError('rewrite_control did not stop')


class TestRewritingTheControlFile(unittest.TestCase):

    def test_the_section_decides_which_key_is_rewritten(self):
        text = rewrite(TEXT_TWO_SECTIONS, 'velocity', ['11.0', '22.0', '33.0'], 'wind')
        config = yaml.safe_load(text)
        self.assertEqual(config['wind']['velocity'], [11.0, 22.0, 33.0])
        self.assertEqual(config['initial_settings']['velocity'], [7450, 0, 0])

    def test_the_first_section_is_rewritten_as_well(self):
        text = rewrite(TEXT_TWO_SECTIONS, 'velocity', ['1.0', '2.0', '3.0'], 'initial_settings')
        config = yaml.safe_load(text)
        self.assertEqual(config['initial_settings']['velocity'], [1.0, 2.0, 3.0])
        self.assertEqual(config['wind']['velocity'], [0.0, 0.0, 0.0])

    def test_without_a_section_the_first_key_wins(self):
        # 従来の呼び方（セクションを渡さない）。ファイルの先頭から探す
        text = rewrite(TEXT_TWO_SECTIONS, 'velocity', ['1.0', '2.0', '3.0'])
        config = yaml.safe_load(text)
        self.assertEqual(config['initial_settings']['velocity'], [1.0, 2.0, 3.0])
        self.assertEqual(config['wind']['velocity'], [0.0, 0.0, 0.0])

    def test_equal_values_are_not_rewritten_together(self):
        # [0.0, 0.0, 0.0] の 3 行は文字列としては同じ。行を指定して置換しないと
        # 3 行とも同じ値になる（文字列置換で実際に起きていた壊れ方）
        text = rewrite(TEXT_TWO_SECTIONS, 'velocity', ['1.0', '2.0', '3.0'], 'wind')
        self.assertEqual(yaml.safe_load(text)['wind']['velocity'], [1.0, 2.0, 3.0])

    def test_only_the_target_lines_change(self):
        text = rewrite(TEXT_TWO_SECTIONS, 'velocity', ['1.0', '2.0', '3.0'], 'wind')
        before = TEXT_TWO_SECTIONS.splitlines()
        after = text.splitlines()
        self.assertEqual(len(before), len(after))
        differing = [i for i in range(0, len(before)) if before[i] != after[i]]
        self.assertEqual(len(differing), 3)

    def test_the_indent_and_the_comment_are_kept(self):
        text = """wind:
  velocity:
      - 0.0   # East
      - 0.0   # North
      - 0.0   # Up
"""
        result = rewrite(text, 'velocity', ['1.0', '2.0', '3.0'], 'wind')
        # 先頭のインデントと行末のコメントは残る（値のあとの空白は 1 個に詰まる）
        self.assertIn('      - 1.0 # East', result)
        self.assertIn('      - 3.0 # Up', result)
        self.assertEqual(yaml.safe_load(result)['wind']['velocity'], [1.0, 2.0, 3.0])

    def test_a_comment_mentioning_the_key_is_not_taken_for_it(self):
        text = """wind:
  # the initial velocity uses the same frame
  velocity:
    - 0.0
    - 0.0
"""
        result = rewrite(text, 'velocity', ['1.0', '2.0'], 'wind')
        self.assertEqual(yaml.safe_load(result)['wind']['velocity'], [1.0, 2.0])

    def test_a_one_element_list_is_rewritten(self):
        # テーブルの風を振るための wind.velocity_factor は 1 要素のリスト。
        # density_factor と同じ形なので、driver から見れば要素数が 1 なだけ
        text = """wind:
  velocity_factor:
    - 1.0
"""
        result = rewrite(text, 'velocity_factor', ['0.8'], 'wind')
        self.assertEqual(yaml.safe_load(result)['wind']['velocity_factor'], [0.8])

    def test_a_missing_section_stops_the_run(self):
        with self.assertRaises(SystemExit):
            rewrite(TEXT_TWO_SECTIONS, 'velocity', ['1.0'], 'atmosphere')

    def test_a_key_of_another_section_is_not_rewritten(self):
        # セクションの終わりで探索を打ち切らないと、指定したセクションにキーが無いとき
        # ファイル末尾まで走り、あとに並ぶ別のセクションの同名キーを書き換えてしまう
        # （computational_setup の velocity を指定して initial_settings.velocity が
        # 動いた、という壊れ方が実際にあった）
        text = rewrite_expecting_a_stop(TEXT_SECTION_WITHOUT_THE_KEY, 'velocity',
                                        ['1.0', '2.0', '3.0'], 'computational_setup')
        self.assertEqual(text, TEXT_SECTION_WITHOUT_THE_KEY)
        self.assertEqual(yaml.safe_load(text)['initial_settings']['velocity'], [7450, 0, 0])

    def test_a_comment_between_the_sections_does_not_close_the_section(self):
        # 区切りはインデントの無い "名前:" 行。空行とコメントはセクションの内側にも
        # 現れるので、そこで打ち切ってはならない
        text = """wind:
  flag_wind: True

# [East, North, Up], m/s
# （桁を詰めて書かれたコメントでも、セクションはまだ終わっていない）

  velocity:
    - 0.0
    - 0.0
    - 0.0

atmosphere:
  velocity:
    - 9.9
"""
        result = rewrite(text, 'velocity', ['1.0', '2.0', '3.0'], 'wind')
        config = yaml.safe_load(result)
        self.assertEqual(config['wind']['velocity'], [1.0, 2.0, 3.0])
        self.assertEqual(config['atmosphere']['velocity'], [9.9])

    def test_a_missing_key_stops_the_run(self):
        with self.assertRaises(SystemExit):
            rewrite(TEXT_TWO_SECTIONS, 'density_factor', ['1.0'], 'wind')

    def test_a_key_which_is_not_a_list_stops_the_run(self):
        text = """wind:
  flag_wind: True
"""
        with self.assertRaises(SystemExit):
            rewrite(text, 'flag_wind', ['False'], 'wind')


class TestTheWindMonteCarloTutorial(unittest.TestCase):

    def test_the_template_and_the_case_agree(self):
        # ばらつかせる基準値はモンテカルロ側の config から読まれ、書き換えられるのは
        # テンプレート側の行。両者がずれていると黙って別の値で走る
        config = load_config(CONFIG_MONTECARLO_WIND)
        template = load_config(CONFIG_TEMPLATE_WIND)

        self.assertTrue(config['wind']['flag_wind'])
        self.assertTrue(template['wind']['flag_wind'])
        self.assertEqual(config['wind']['velocity'], template['wind']['velocity'])
        self.assertEqual(config['initial_settings']['velocity'],
                         template['initial_settings']['velocity'])
        self.assertEqual(config['initial_settings']['coordinate'],
                         template['initial_settings']['coordinate'])

    def test_the_template_path_points_at_the_wind_template(self):
        config = load_config(CONFIG_MONTECARLO_WIND)
        self.assertEqual(config['montecarlo']['template_path_specify'], 'manual')
        path = os.path.normpath(os.path.join(os.path.dirname(CONFIG_MONTECARLO_WIND),
                                             config['montecarlo']['template_path']))
        self.assertEqual(path, os.path.normpath(TEMPLATE_WIND))

    def test_every_target_variable_exists_in_both_files(self):
        config = load_config(CONFIG_MONTECARLO_WIND)
        with open(CONFIG_TEMPLATE_WIND) as f:
            text_template = f.read()

        for name, section, dispersion in config['montecarlo']['target_variable']:
            self.assertIn(section, config)
            self.assertIn(name, config[section])
            self.assertIsInstance(config[section][name], list)
            self.assertGreater(dispersion, 0.0)
            # テンプレート側でも、そのセクションのそのキーが書き換えられること
            values = ['{:.1f}'.format(float(i)) for i in range(0, len(config[section][name]))]
            rewritten = yaml.safe_load(rewrite(text_template, name, values, section))
            self.assertEqual(rewritten[section][name], [float(i) for i in range(0, len(values))])

    def test_the_template_carries_the_tables_it_needs(self):
        template = load_config(CONFIG_TEMPLATE_WIND)
        for section, key, name in (('atmosphere', 'directory_atmosphere', 'filename_atmosphere'),
                                   ('satellite', 'directory_aerodynamic', 'filename_aerodynamic')):
            self.assertEqual(template[section]['directory_path_specify'], 'manual')
            # load_config が config の位置基準の絶対パスに直している
            self.assertTrue(os.path.exists(os.path.join(template[section][key],
                                                        template[section][name])))

    def test_two_shortened_cases_run_and_differ(self):
        # 制御ファイル・テンプレート・driver・出力の結線を丸ごと通す。
        # 1000 s に切り詰めても着地はしないが、風の効きは軌道に出る
        with tempfile.TemporaryDirectory() as directory:
            template = os.path.join(directory, 'template')
            shutil.copytree(TEMPLATE_WIND, template)
            shorten_the_run(os.path.join(template, 'config.yml'), 1000.0)

            with open(CONFIG_MONTECARLO_WIND) as f:
                config = yaml.safe_load(f)
            config['montecarlo']['number_iteration'] = 2
            config['montecarlo']['maximum_number_execution'] = 2
            config['montecarlo']['template_path_specify'] = 'manual'
            config['montecarlo']['template_path'] = template
            config['montecarlo']['work_dir'] = 'work'
            path_config = os.path.join(directory, 'config.yml')
            with open(path_config, 'w') as f:
                yaml.safe_dump(config, f)

            environment = dict(os.environ)
            environment['TACODE_PYTHON'] = sys.executable
            subprocess.run([sys.executable, os.path.join(ROOT_DIR, 'src', 'tacode-montecarlo.py'),
                            '-file', 'config.yml'],
                           cwd=directory, env=environment, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            wind_list = []
            position_list = []
            for name in ('case0001', 'case0002'):
                case = os.path.join(directory, 'work', name)
                with open(os.path.join(case, 'config.yml')) as f:
                    case_config = yaml.safe_load(f)
                wind_list.append(case_config['wind']['velocity'])
                state = montecarlo_dispersion.get_terminal_state(
                    os.path.join(case, 'output_result', 'tecplot.dat'))
                position_list.append([state['X'], state['Y'], state['Z']])
                # 風の列が出ており、その行の風が config のものと合う
                self.assertAlmostEqual(state['WindE'], case_config['wind']['velocity'][0], places=6)
                self.assertAlmostEqual(state['WindN'], case_config['wind']['velocity'][1], places=6)
                # ばらつかせていない初期条件は基準のまま
                self.assertEqual(case_config['initial_settings']['velocity'],
                                 config['initial_settings']['velocity'])

            self.assertNotEqual(wind_list[0], wind_list[1])
            self.assertGreater(np.linalg.norm(np.array(position_list[0])-np.array(position_list[1])),
                               1.e-6)

            # driver が最後に全ケースをまとめること（1 ケース = 1 ゾーン）
            gathered = os.path.join(directory, config['montecarlo']['result_dir'],
                                    config['montecarlo']['filename_tecplot'])
            self.assertTrue(os.path.exists(gathered))
            with open(gathered) as f:
                text = f.read()
            self.assertEqual(text.count('zone t="case'), 2)
            # 変数行はケースの出力から引き継ぐので、風の列がそのまま出ている
            self.assertEqual(text.count('Variables ='), 1)
            self.assertIn('WindE[m/s]', text)


class TestFindingTheCaseDirectories(unittest.TestCase):
    """
    `montecarlo_dispersion.find_case_directory` が拾うのはケースだけであること。

    ケースの名前は montecarlo.get_case_directory が作る「<case_dir> + 4 桁」。
    「テンプレート以外の全部」で拾うと、作業ディレクトリに置いた結果ディレクトリまで
    ケースとして数え、統計の母数が変わる。
    """

    def test_only_the_numbered_directories_are_taken(self):
        with tempfile.TemporaryDirectory() as directory:
            for name in ('case0001', 'case0002', 'case_template', 'result_tacode', 'figure'):
                os.makedirs(os.path.join(directory, name))
            with open(os.path.join(directory, 'case0003'), 'w') as f:
                f.write('a file, not a directory\n')

            case_list = [os.path.basename(path)
                         for path in montecarlo_dispersion.find_case_directory(directory)]
            self.assertEqual(case_list, ['case0001', 'case0002'])


class TestDetectingCasesWhichFail(unittest.TestCase):
    """
    ケースが落ちたことを親が検知すること。

    子の終了コードを見ずに待つだけだと、100 ケースのうち 3 ケースが落ちても
    親は exit 0 を返し、統計は残った 97 ケースから静かに作られる。
    """

    def build(self, directory, flag_allow_failure=False):
        driver = montecarlo()
        driver.work_dir = directory
        driver.case_dir = 'case'
        driver.filename_trajectory_tacode = os.path.join('output_result', 'tecplot.dat')
        driver.flag_allow_failure = flag_allow_failure
        return driver

    def test_a_non_zero_exit_code_is_recorded(self):
        driver = self.build('work')
        process = subprocess.Popen([sys.executable, '-c', 'raise SystemExit(3)'])
        driver.process_list = [['work/case0001', process]]
        driver.wait_tacode()
        self.assertEqual(driver.case_failed, [['work/case0001', 3]])
        self.assertEqual(driver.process_list, [])

    def test_a_case_which_succeeds_is_not_recorded(self):
        driver = self.build('work')
        process = subprocess.Popen([sys.executable, '-c', 'pass'])
        driver.process_list = [['work/case0001', process]]
        driver.wait_tacode()
        self.assertEqual(driver.case_failed, [])

    def test_a_case_without_a_result_file_is_found(self):
        # 子が 0 を返しても出力を書いていないことがある（そのまま統計から消える）
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))
            case = write_case(directory, 'case0002', (6378.137, -1.0, 0.0), (15.0, 10.0, 0.0))
            os.remove(os.path.join(case, 'output_result', 'tecplot.dat'))

            driver = self.build(directory)
            driver.check_case_result()
            self.assertEqual(driver.case_no_result, [os.path.join(directory, 'case0002')])

    def test_the_report_stops_the_run(self):
        driver = self.build('work')
        driver.case_failed = [['work/case0001', 1]]
        with self.assertRaises(SystemExit) as raised:
            driver.report_failure()
        self.assertEqual(raised.exception.code, 1)

    def test_the_report_is_silent_when_every_case_completes(self):
        self.assertEqual(self.build('work').report_failure(), 0)

    def test_the_flag_lets_the_run_go_on(self):
        driver = self.build('work', flag_allow_failure=True)
        driver.case_failed = [['work/case0001', 1]]
        driver.case_no_result = ['work/case0002']
        self.assertEqual(driver.report_failure(), 2)

    def test_a_case_which_fails_and_writes_nothing_counts_once(self):
        # 落ちたケースは結果ファイルも書いていないので、両方の記録に挙がる
        driver = self.build('work', flag_allow_failure=True)
        driver.case_failed = [['work/case0001', 1]]
        driver.case_no_result = ['work/case0001']
        self.assertEqual(driver.report_failure(), 1)

    def test_a_case_is_not_counted_twice(self):
        driver = self.build('work')
        driver.add_case_no_result('work/case0001')
        driver.add_case_no_result('work/case0001')
        self.assertEqual(driver.case_no_result, ['work/case0001'])

    def test_the_driver_returns_a_non_zero_exit_code(self):
        # 結線を丸ごと通す: テンプレートの設定を壊し、driver が失敗を返すこと。
        # ケースは初期化で落ちるので、計算そのものは走らない
        with tempfile.TemporaryDirectory() as directory:
            template = os.path.join(directory, 'template')
            shutil.copytree(TEMPLATE_WIND, template)
            path_template = os.path.join(template, 'config.yml')
            with open(path_template) as f:
                text = f.read()
            with open(path_template, 'w') as f:
                f.write(text.replace('kind_atmosphere_model: fileread',
                                     'kind_atmosphere_model: bogus'))

            with open(CONFIG_MONTECARLO_WIND) as f:
                config = yaml.safe_load(f)
            config['montecarlo']['number_iteration'] = 2
            config['montecarlo']['maximum_number_execution'] = 2
            config['montecarlo']['template_path_specify'] = 'manual'
            config['montecarlo']['template_path'] = template
            config['montecarlo']['work_dir'] = 'work'
            config['montecarlo']['result_dir'] = 'result'
            with open(os.path.join(directory, 'config.yml'), 'w') as f:
                yaml.safe_dump(config, f)

            environment = dict(os.environ)
            environment['TACODE_PYTHON'] = sys.executable
            done = subprocess.run([sys.executable,
                                   os.path.join(ROOT_DIR, 'src', 'tacode-montecarlo.py'),
                                   '-file', 'config.yml'],
                                  cwd=directory, env=environment,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn('Cases that did not complete:', done.stdout.decode())


def shorten_the_run(path, time_elapsed_maximum):
    """テンプレートの計算時間を短くする（テストを軽くするため）。"""
    with open(path) as f:
        config = yaml.safe_load(f)
    config['computational_setup']['time_elapsed_maximum'] = time_elapsed_maximum
    with open(path, 'w') as f:
        yaml.safe_dump(config, f)
    return


TECPLOT_HEADER = ('# Tecplot data: Tacode\n'
                  'Variables = Time[s],X[km],Y[km],Z[km],Long[deg.],Lati[deg.],Alti[km]\n'
                  'zone t=time i= 2 f=point\n')


def write_case(directory, name, position, wind_velocity):
    """終端点だけが意味を持つ、最小のケースディレクトリを作る。"""
    case = os.path.join(directory, name)
    os.makedirs(os.path.join(case, 'output_result'))
    with open(os.path.join(case, 'output_result', 'tecplot.dat'), 'w') as f:
        f.write(TECPLOT_HEADER)
        f.write('0.0 6378.137000000 0.000000000 0.000000000 0.0 0.0 100.0\n')
        f.write('1.0 {:.9f} {:.9f} {:.9f} 0.0 0.0 0.0\n'.format(*position))
    with open(os.path.join(case, 'config.yml'), 'w') as f:
        yaml.safe_dump({'wind': {'velocity': list(wind_velocity)}}, f)
    return case


class TestTheDispersionHelper(unittest.TestCase):

    def test_the_offsets_are_east_and_north_of_the_origin(self):
        # 赤道上、経度 0 の真上を基準にとると、ECEF の +y が東、+z が北になる
        radius = 6378.137
        state_list = [{'X': radius, 'Y': 0.0, 'Z': 0.0},
                      {'X': radius, 'Y': 2.0, 'Z': 0.0},
                      {'X': radius, 'Y': 0.0, 'Z': 3.0}]
        origin, east, north, up = montecarlo_dispersion.get_dispersion(
            state_list, np.array([radius, 0.0, 0.0]))
        np.testing.assert_allclose(east, [0.0, 2.0, 0.0], atol=1.e-9)
        np.testing.assert_allclose(north, [0.0, 0.0, 3.0], atol=1.e-9)
        np.testing.assert_allclose(up, [0.0, 0.0, 0.0], atol=1.e-9)

    def test_the_default_origin_is_the_mean(self):
        radius = 6378.137
        state_list = [{'X': radius, 'Y': -1.0, 'Z': 0.0},
                      {'X': radius, 'Y': +1.0, 'Z': 0.0}]
        origin, east, north, up = montecarlo_dispersion.get_dispersion(state_list)
        self.assertAlmostEqual(float(np.mean(east)), 0.0, places=9)
        np.testing.assert_allclose(east, [-1.0, +1.0], atol=1.e-6)

    def test_the_dispersed_input_is_found_from_the_template(self):
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case_template', (6378.137, 0.0, 0.0), (20.0, 10.0, 0.0))
            case_list = [write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0)),
                         write_case(directory, 'case0002', (6378.137, -1.0, 0.0), (15.0, 10.0, 0.0))]

            template = montecarlo_dispersion.find_template_directory(directory)
            self.assertEqual(os.path.basename(template), 'case_template')
            self.assertEqual(montecarlo_dispersion.find_case_directory(directory), case_list)

            name_list, value_dict = montecarlo_dispersion.get_dispersed_variable(
                case_list, template, 'config.yml')
            # 動いたのは東向き成分だけ。北と上はテンプレートと同じなので出てこない
            self.assertEqual(name_list, ['wind.velocity[0]'])
            self.assertEqual(value_dict[case_list[0]]['wind.velocity[0]'], 25.0)
            self.assertEqual(value_dict[case_list[1]]['wind.velocity[0]'], 15.0)

    def test_the_command_line_writes_a_table_and_a_csv(self):
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case_template', (6378.137, 0.0, 0.0), (20.0, 10.0, 0.0))
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))
            write_case(directory, 'case0002', (6378.137, -1.0, 0.0), (15.0, 10.0, 0.0))
            path_csv = os.path.join(directory, 'dispersion.csv')

            result = subprocess.run([sys.executable, SCRIPT, directory, '-o', path_csv],
                                    check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            text = result.stdout.decode()
            self.assertIn('case0001', text)
            self.assertIn('wind.velocity[0]', text)
            self.assertIn('Number of cases:  2', text)

            with open(path_csv) as f:
                lines = f.read().splitlines()
            self.assertEqual(len(lines), 3)
            self.assertIn('East[km]', lines[0])

    def test_the_covariance_ellipse_follows_the_scatter(self):
        # 東西に広く南北に狭い雲。長軸は東西を向き、半径は各方向の標準偏差になる
        angle = np.linspace(0.0, 2.0*np.pi, 360, endpoint=False)
        east = 10.0*np.cos(angle)
        north = 2.0*np.sin(angle)
        semi_major, semi_minor, tilt = montecarlo_dispersion.get_covariance_ellipse(east, north)
        self.assertAlmostEqual(semi_major, float(np.std(east, ddof=1)), places=6)
        self.assertAlmostEqual(semi_minor, float(np.std(north, ddof=1)), places=6)
        self.assertAlmostEqual(abs(tilt) % 180.0, 0.0, places=6)

    def test_the_ellipse_turns_with_the_scatter(self):
        # 45 度に伸びた雲なら長軸も 45 度
        value = np.linspace(-1.0, 1.0, 101)
        semi_major, semi_minor, tilt = montecarlo_dispersion.get_covariance_ellipse(value, value)
        self.assertAlmostEqual(tilt % 180.0, 45.0, places=6)
        self.assertLess(semi_minor, 1.e-9)

    def test_a_long_path_is_shortened_for_the_plot(self):
        self.assertEqual(montecarlo_dispersion.shorten_path('/a/b/c/d/e.dat'),
                         os.path.join('...', 'c', 'd', 'e.dat'))
        self.assertEqual(montecarlo_dispersion.shorten_path('d/e.dat'),
                         os.path.join('d', 'e.dat'))

    def test_a_mark_is_placed_but_left_out_of_the_statistics(self):
        # --mark は「モンテカルロの外で走らせた計算」を図に重ねるためのもの。
        # 統計（ケース数、平均、標準偏差）には入らない
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_mark:
            write_case(directory, 'case_template', (6378.137, 0.0, 0.0), (20.0, 10.0, 0.0))
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))
            write_case(directory, 'case0002', (6378.137, 3.0, 0.0), (15.0, 10.0, 0.0))
            mark = os.path.join(write_case(directory_mark, 'table', (6378.137, -5.0, 0.0),
                                           (0.0, 0.0, 0.0)),
                                'output_result', 'tecplot.dat')

            result = subprocess.run([sys.executable, SCRIPT, directory,
                                     '--mark', 'table='+mark],
                                    check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            text = result.stdout.decode()
            # ケースの平均 (+2 km) を原点にとるので、比較点は -7 km
            self.assertIn('Mark: table  East -7.000 km,  North +0.000 km', text)
            self.assertIn('Number of cases:  2', text)

    def test_a_missing_mark_file_stops_the_run(self):
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case_template', (6378.137, 0.0, 0.0), (20.0, 10.0, 0.0))
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))
            result = subprocess.run([sys.executable, SCRIPT, directory,
                                     '--mark', 'x=/no/such/file.dat'],
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertNotEqual(result.returncode, 0)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_the_plot_is_written(self):
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case_template', (6378.137, 0.0, 0.0), (20.0, 10.0, 0.0))
            for i in range(0, 5):
                write_case(directory, 'case{:04d}'.format(i+1),
                           (6378.137, 1.0*i, 0.5*i), (20.0+i, 10.0, 0.0))
            path_plot = os.path.join(directory, 'dispersion.png')

            subprocess.run([sys.executable, SCRIPT, directory, '--plot', path_plot],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            self.assertTrue(os.path.exists(path_plot))
            self.assertGreater(os.path.getsize(path_plot), 1000)

    def test_a_reference_case_moves_the_origin(self):
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            write_case(directory, 'case_template', (6378.137, 0.0, 0.0), (20.0, 10.0, 0.0))
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))
            write_case(directory, 'case0002', (6378.137, 3.0, 0.0), (15.0, 10.0, 0.0))
            # 基準は走査するディレクトリの外に置く（ケースとして数えられないように）
            reference = os.path.join(write_case(directory_reference, 'calm', (6378.137, 0.0, 0.0),
                                                (0.0, 0.0, 0.0)),
                                     'output_result', 'tecplot.dat')

            result = subprocess.run([sys.executable, SCRIPT, directory, '--reference', reference],
                                    check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            text = result.stdout.decode()
            # 基準からの平均のずれ（東 (1+3)/2 = 2 km）が出る。ケースの平均を原点に
            # とったのでは 0 になってしまう量で、風の平均的な効きはこちらで読む
            self.assertIn('East +2.000 km', text)


class TestGatheringTheResults(unittest.TestCase):
    """
    postprocess は各ケースの Tecplot 出力を 1 ファイルにまとめる（1 ケース = 1 ゾーン）。

    列の中身には立ち入らず、変数行はケースの出力からそのまま引き継ぐ
    （風や姿勢で列が増えるので、形式をここにもう 1 つ持たない）。
    """

    def build(self, directory, flag_tecplot=True):
        """postprocess に必要な属性だけを持たせた montecarlo を返す。"""
        driver = montecarlo()
        driver.work_dir = directory
        driver.case_dir = 'case'
        driver.filename_trajectory_tacode = os.path.join('output_result', 'tecplot.dat')
        driver.result_dir = os.path.join(directory, 'result')
        driver.flag_tecplot = flag_tecplot
        driver.filename_tecplot = 'gathered.dat'
        os.makedirs(driver.result_dir)
        return driver

    def test_every_case_becomes_a_zone(self):
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))
            write_case(directory, 'case0002', (6378.137, -1.0, 0.0), (15.0, 10.0, 0.0))

            path = self.build(directory).postprocess()
            with open(path) as f:
                lines = f.read().splitlines()

            # 変数行は 1 度だけ。ケースの出力のものがそのまま出る
            variables = [line for line in lines if line.lower().startswith('variables')]
            self.assertEqual(len(variables), 1)
            self.assertIn('Variables = Time[s],X[km]', variables[0])

            zone = [line for line in lines if line.startswith('zone')]
            self.assertEqual(zone, ['zone t="case0001" i= 2 f=point',
                                    'zone t="case0002" i= 2 f=point'])
            # 宣言した点数と実際の行数が合っていること（Tecplot が読めなくなる）
            data = [line for line in lines
                    if not line.startswith(('zone', '#')) and not line.lower().startswith('variables')]
            self.assertEqual(len(data), 4)
            self.assertIn('1.0 6378.137000000 1.000000000 0.000000000 0.0 0.0 0.0', data)

    def test_each_zone_declares_its_own_number_of_points(self):
        # ケースごとに行数は違う（着地する時刻が違う）。ゾーンの i= と実際の行数が
        # ずれると Tecplot が読めなくなるので、ケースごとに数え直す
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))
            case = write_case(directory, 'case0002', (6378.137, -1.0, 0.0), (15.0, 10.0, 0.0))
            with open(os.path.join(case, 'output_result', 'tecplot.dat'), 'a') as f:
                f.write('2.0 6378.137000000 -1.000000000 0.000000000 0.0 0.0 0.0\n')

            path = self.build(directory).postprocess()
            with open(path) as f:
                lines = f.read().splitlines()

            zone = [line for line in lines if line.startswith('zone')]
            self.assertEqual(zone, ['zone t="case0001" i= 2 f=point',
                                    'zone t="case0002" i= 3 f=point'])

    def test_the_template_is_left_out(self):
        # テンプレートは連番でないので拾われない（拾うとケース数が 1 つ増える）
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case_template', (6378.137, 0.0, 0.0), (20.0, 10.0, 0.0))
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))

            path = self.build(directory).postprocess()
            with open(path) as f:
                text = f.read()
            self.assertEqual(text.count('zone t="'), 1)
            self.assertIn('zone t="case0001"', text)

    def test_a_case_without_a_result_is_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))
            case = write_case(directory, 'case0002', (6378.137, -1.0, 0.0), (15.0, 10.0, 0.0))
            os.remove(os.path.join(case, 'output_result', 'tecplot.dat'))

            path = self.build(directory).postprocess()
            with open(path) as f:
                text = f.read()
            self.assertEqual(text.count('zone t="'), 1)

    def test_cases_with_different_columns_stop_the_run(self):
        # 風を一部のケースだけで入れた、というような取り違え。黙って混ぜない
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))
            case = write_case(directory, 'case0002', (6378.137, -1.0, 0.0), (15.0, 10.0, 0.0))
            path_case = os.path.join(case, 'output_result', 'tecplot.dat')
            with open(path_case) as f:
                text = f.read()
            with open(path_case, 'w') as f:
                f.write(text.replace('Alti[km]', 'Alti[km],WindE[m/s]'))

            with self.assertRaises(SystemExit):
                self.build(directory).postprocess()

    def test_nothing_is_written_when_the_flag_is_off(self):
        with tempfile.TemporaryDirectory() as directory:
            write_case(directory, 'case0001', (6378.137, 1.0, 0.0), (25.0, 10.0, 0.0))
            driver = self.build(directory, flag_tecplot=False)
            self.assertIsNone(driver.postprocess())
            self.assertEqual(os.listdir(driver.result_dir), [])

    def test_an_empty_working_directory_is_reported_and_not_written(self):
        with tempfile.TemporaryDirectory() as directory:
            driver = self.build(directory)
            self.assertIsNone(driver.postprocess())
            self.assertEqual(os.listdir(driver.result_dir), [])


class TestTheDispersionAnimation(unittest.TestCase):
    """
    src_helper/montecarlo_animation は、ケースの軌跡と分散円を動かして見せる。
    幾何（基準ケースから見た東・北）は描かなくても検査できるので分けてある。
    """

    def test_the_offsets_are_east_and_north_of_the_reference_at_the_same_time(self):
        # 赤道上、経度 0 の真上から見ると ECEF の +y が東、+z が北。
        # 基準は「同時刻の基準ケースの位置」なので、時々刻々のずれが出る
        radius = 6378.137
        position = np.array([[[radius, 0.0, 0.0], [radius, 2.0, 3.0]]])
        position_reference = np.array([[radius, 0.0, 0.0], [radius, 0.0, 0.0]])

        east, north = montecarlo_animation.get_offset(position, position_reference)
        np.testing.assert_allclose(east, [[0.0, 2.0]], atol=1.e-9)
        np.testing.assert_allclose(north, [[0.0, 3.0]], atol=1.e-9)

    def test_the_window_holds_the_origin_and_every_point(self):
        east = np.array([[0.0, 20.0], [0.0, 30.0]])
        north = np.array([[0.0, 5.0], [0.0, -4.0]])
        window_east, window_north = montecarlo_animation.get_window(east, north)

        self.assertLess(window_east[0], 0.0)
        self.assertGreater(window_east[1], 30.0)
        self.assertLess(window_north[0], -4.0)
        self.assertGreater(window_north[1], 5.0)
        # 縦横の縮尺を揃えて描くので、窓も正方形にしておく
        self.assertAlmostEqual(window_east[1]-window_east[0],
                               window_north[1]-window_north[0], places=9)

    def test_the_altitude_is_stretched_about_the_surface(self):
        # globe ビューの高度の誇張。地表（高度 0）は動かず、高度だけが factor 倍になる
        radius = montecarlo_animation.RADIUS_PLANET
        position = np.array([[radius, 0.0, 0.0], [radius + 100.0, 0.0, 0.0]])

        np.testing.assert_allclose(montecarlo_animation.stretch_altitude(position, 1.0),
                                   position, rtol=0.0, atol=0.0)
        stretched = montecarlo_animation.stretch_altitude(position, 10.0)
        self.assertAlmostEqual(float(np.linalg.norm(stretched[0])), radius, places=6)
        self.assertAlmostEqual(float(np.linalg.norm(stretched[1])), radius + 1000.0, places=6)
        # 向きは変えない
        np.testing.assert_allclose(stretched[1]/np.linalg.norm(stretched[1]),
                                   position[1]/np.linalg.norm(position[1]), atol=1.e-12)

    def test_the_longitude_is_unwrapped_and_labelled_back(self):
        # +-180 度をまたぐ軌道でケースの平均を取るため、経度は連続化して持つ
        self.assertAlmostEqual(montecarlo_animation.normalize_longitude(240.0), -120.0, places=9)
        self.assertAlmostEqual(montecarlo_animation.normalize_longitude(-190.0), 170.0, places=9)

    def run_script(self, arguments):
        return subprocess.run([sys.executable, SCRIPT_ANIMATION] + arguments,
                              capture_output=True, text=True)

    def build_run(self, directory, directory_reference):
        """最小のモンテカルロ結果（3 ケース）と、基準ケースを作る。"""
        write_case(directory, 'case_template', (6378.137, 0.0, 0.0), (20.0, 10.0, 0.0))
        for i in range(0, 3):
            write_case(directory, 'case{:04d}'.format(i+1),
                       (6378.137, 1.0*i, 0.5*i), (20.0+i, 10.0, 0.0))
        # 基準は走査するディレクトリの外に置く（ケースとして数えられないように）
        return os.path.join(write_case(directory_reference, 'calm', (6378.137, 0.0, 0.0),
                                       (0.0, 0.0, 0.0)),
                            'output_result', 'tecplot.dat')

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_a_snapshot_of_the_last_frame_is_written(self):
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            output = os.path.join(directory, 'snapshot.png')
            completed = self.run_script([directory, '--reference', reference,
                                         '--snapshot', '-1', '--dpi', '50', '-o', output])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertTrue(os.path.exists(output))
            self.assertGreater(os.path.getsize(output), 1000)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_an_html_animation_is_written_as_a_single_file(self):
        # ffmpeg の要らない出力。フレーム数だけ埋め込まれること
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            output = os.path.join(directory, 'anim.html')
            completed = self.run_script([directory, '--reference', reference,
                                         '--frames', '3', '--dpi', '40', '-o', output])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            with open(output) as f:
                page = f.read()
            self.assertEqual(page.count('data:image/png;base64'), 3)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_a_marked_run_is_drawn_but_left_out_of_the_statistics(self):
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_mark, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            mark = os.path.join(write_case(directory_mark, 'table', (6378.137, -5.0, 0.0),
                                           (0.0, 0.0, 0.0)),
                                'output_result', 'tecplot.dat')
            output = os.path.join(directory, 'snapshot.png')
            completed = self.run_script([directory, '--reference', reference,
                                         '--mark', 'table='+mark, '--snapshot', '-1',
                                         '--dpi', '50', '-o', output])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertIn('Cases read:  3', completed.stdout)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_the_three_dimensional_view_is_drawn(self):
        # --view 3d は軌跡そのもの（経度・緯度・高度）を箱に描き、分散円を床に落とす
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            output = os.path.join(directory, 'box.png')
            completed = self.run_script([directory, '--reference', reference, '--view', '3d',
                                         '--snapshot', '-1', '--dpi', '50', '-o', output])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertGreater(os.path.getsize(output), 1000)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_the_relative_view_is_drawn(self):
        # --view 3d-relative は同じ箱を「基準ケースから見た東・北」で描く。
        # 散らばりが育つ様子はこちらでしか読めない
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            output = os.path.join(directory, 'relative.png')
            completed = self.run_script([directory, '--reference', reference,
                                         '--view', '3d-relative', '--snapshot', '-1',
                                         '--dpi', '50', '-o', output])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertGreater(os.path.getsize(output), 1000)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_the_altitude_limit_clips_the_three_dimensional_view(self):
        # 箱の上を越える点は描かない（matplotlib の 3 次元は線を箱で切ってくれない）ので、
        # 上限を与えても最後まで書けること
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            output = os.path.join(directory, 'box.png')
            completed = self.run_script([directory, '--reference', reference, '--view', '3d',
                                         '--altitude-max', '50', '--snapshot', '-1',
                                         '--dpi', '50', '-o', output])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertGreater(os.path.getsize(output), 1000)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_the_globe_view_is_drawn(self):
        # --view globe は基準からの差ではなく ECEF の軌道そのものを、地球ごと描く
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            output = os.path.join(directory, 'globe.png')
            completed = self.run_script([directory, '--reference', reference,
                                         '--view', 'globe', '--exaggerate', '20',
                                         '--snapshot', '-1', '--dpi', '50', '-o', output])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertGreater(os.path.getsize(output), 1000)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_the_follow_view_is_drawn(self):
        # --view follow はカメラが基準ケースを追い、窓の中に「基準からのずれ」を描く
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            output = os.path.join(directory, 'follow.png')
            completed = self.run_script([directory, '--reference', reference,
                                         '--view', 'follow', '--snapshot', '-1',
                                         '--dpi', '50', '-o', output])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertGreater(os.path.getsize(output), 1000)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_the_follow_view_takes_a_fixed_window(self):
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            output = os.path.join(directory, 'follow_fixed.png')
            completed = self.run_script([directory, '--reference', reference,
                                         '--view', 'follow', '--window', '25',
                                         '--snapshot', '-1', '--dpi', '50', '-o', output])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertGreater(os.path.getsize(output), 1000)

    def test_the_follow_window_is_centred_between_the_reference_and_the_cases(self):
        # 基準ケースに窓を載せると、風下に寄ったケース群で窓の半分が空く
        position_reference = np.array([[7000.0, 0.0, 0.0], [7000.0, 0.0, 0.0]])
        position = np.array([[[7000.0, 10.0, 0.0], [7000.0, 30.0, 0.0]],
                             [[7000.0, 10.0, 0.0], [7000.0, 30.0, 0.0]]])
        centre, half = montecarlo_animation.get_follow_window(position, position_reference,
                                                              np.zeros(len(position_reference)))

        np.testing.assert_allclose(centre[0], [7000.0, 5.0, 0.0])
        np.testing.assert_allclose(centre[1], [7000.0, 15.0, 0.0])
        # 半幅は中心からいちばん遠い点（この場合はケースと基準が同距離）に余白を掛けたもの
        self.assertAlmostEqual(half[1], montecarlo_animation.MARGIN_FOLLOW*15.0, places=9)

    def test_the_follow_window_holds_every_case_and_the_reference(self):
        # 基準ケースは別に見ていない。中心が中点なので、基準ケースの距離は
        # 最遠のケースの距離を超えられない（この性質に頼っている）
        position_reference = np.array([[7000.0, 0.0, 0.0]])
        position = np.array([[[7000.0, 12.0, 3.0]], [[7000.0, -4.0, -9.0]]])
        centre, half = montecarlo_animation.get_follow_window(position, position_reference,
                                                              np.zeros(len(position_reference)))

        for point in (position[0, 0], position[1, 0], position_reference[0]):
            self.assertLessEqual(float(np.linalg.norm(point - centre[0])), float(half[0]) + 1.e-9)

    def test_the_reference_is_inside_the_window_whatever_the_cases(self):
        # 上の性質を乱数で確かめる（成り立たなければ基準ケースが窓から出る）
        generator = np.random.default_rng(20260911)
        for _ in range(0, 200):
            position_reference = generator.normal(scale=50.0, size=(1, 3))
            position = generator.normal(scale=50.0, size=(generator.integers(1, 8), 1, 3))
            centre, half = montecarlo_animation.get_follow_window(
                position, position_reference, np.zeros(len(position_reference)))
            self.assertLessEqual(float(np.linalg.norm(position_reference[0] - centre[0])),
                                 float(half[0]) + 1.e-9)

    def test_the_follow_window_never_narrows(self):
        # 散らばりが一時的に縮むたびに寄っては引いてを繰り返すと、何が動いているのか読めない
        position_reference = np.zeros((3, 3))
        position = np.array([[[0.0, 40.0, 0.0], [0.0, 4.0, 0.0], [0.0, 8.0, 0.0]]])
        centre, half = montecarlo_animation.get_follow_window(position, position_reference,
                                                              np.zeros(len(position_reference)))

        self.assertTrue(np.all(np.diff(half) >= 0.0))

    def test_the_follow_window_has_a_floor_and_can_be_fixed(self):
        # 突入直後は散らばりが 0 なので、下限が無いと窓が潰れる
        position_reference = np.zeros((2, 3))
        position = np.zeros((3, 2, 3))
        centre, half = montecarlo_animation.get_follow_window(position, position_reference,
                                                              np.zeros(len(position_reference)))
        np.testing.assert_allclose(half, montecarlo_animation.WINDOW_FOLLOW_MINIMUM)

        centre, half = montecarlo_animation.get_follow_window(
            position, position_reference, np.zeros(len(position_reference)), 12.5)
        np.testing.assert_allclose(half, 12.5)

    def test_the_window_keeps_the_ground_in_sight(self):
        # 散らばりだけで決めると、突入直後は窓が数 km なのに地表は 150 km 下にあり、
        # 何も無い空間に点が 1 つという絵になる
        position_reference = np.array([[6378.137 + 150.0, 0.0, 0.0],
                                       [6378.137 + 10.0, 0.0, 0.0]])
        altitude_reference = np.array([150.0, 10.0])
        position = position_reference[np.newaxis, :, :].copy()

        centre, half = montecarlo_animation.get_follow_window(position, position_reference,
                                                              altitude_reference,
                                                              0.0, 1.15)
        # 地表（高度 0）が窓の中にあること
        self.assertGreaterEqual(float(half[0]), 150.0)
        self.assertGreaterEqual(float(half[1]), 10.0)
        # 高度に追わせているので、降りるにつれて窓は閉じる
        self.assertLess(float(half[1]), float(half[0]))

        # 0 を与えると地表は考えない（散らばりだけ）
        centre, half = montecarlo_animation.get_follow_window(position, position_reference,
                                                              altitude_reference, 0.0, 0.0)
        self.assertLess(float(half[0]), 150.0)

    def test_a_fixed_window_ignores_the_ground_and_the_dispersion(self):
        position_reference = np.array([[6378.137 + 150.0, 0.0, 0.0]])
        position = np.array([[[6378.137 + 150.0, 40.0, 0.0]]])
        centre, half = montecarlo_animation.get_follow_window(position, position_reference,
                                                              np.array([150.0]), 8.0, 1.15)
        np.testing.assert_allclose(half, 8.0)

    def test_the_scale_bar_carries_ticks_at_both_ends(self):
        # 目盛りが無いと、ただの線が何を意味するのか分からない
        origin = np.zeros(3)
        unit_along = np.array([1.0, 0.0, 0.0])
        unit_up = np.array([0.0, 0.0, 1.0])
        bar = montecarlo_animation.get_scale_bar(origin, 10.0, unit_along, unit_up,
                                                 ratio_tick=0.1)

        # NaN で 3 つの区間に切れていること（目盛り・本体・目盛り）
        self.assertEqual(int(np.count_nonzero(np.isnan(bar[:, 0]))), 2)
        finite = bar[np.isfinite(bar[:, 0])]
        # 本体は長さ 10 km、目盛りは上下に 1 km
        self.assertAlmostEqual(float(np.max(finite[:, 0]) - np.min(finite[:, 0])), 10.0, places=9)
        self.assertAlmostEqual(float(np.max(finite[:, 2]) - np.min(finite[:, 2])), 2.0, places=9)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_the_follow_view_draws_a_scale_bar_when_it_is_asked_for(self):
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            output = os.path.join(directory, 'follow_bar.png')
            completed = self.run_script([directory, '--reference', reference,
                                         '--view', 'follow', '--scale-bar',
                                         '--snapshot', '-1', '--dpi', '50', '-o', output])
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertGreater(os.path.getsize(output), 1000)

    def test_a_point_outside_the_window_is_dropped(self):
        # matplotlib の 3 次元は線を箱で切ってくれないので、外れる点は NaN にして落とす
        centre = np.array([0.0, 0.0, 0.0])
        line = np.array([[1.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 0.0, -2.0]])
        clipped = montecarlo_animation.clip_to_window(line, centre, 3.0)

        self.assertFalse(np.any(np.isnan(clipped[0])))
        self.assertTrue(np.all(np.isnan(clipped[1])))
        self.assertFalse(np.any(np.isnan(clipped[2])))
        # 元の配列は書き換えない（同じ軌跡を先頭の点にも使う）
        self.assertFalse(np.any(np.isnan(line)))

    def test_the_scale_bar_is_a_round_number_within_the_window(self):
        for half, expected in ((0.6, 0.5), (2.0, 2.0), (7.0, 5.0), (30.0, 20.0), (99.0, 50.0)):
            self.assertAlmostEqual(montecarlo_animation.get_scale_length(half), expected,
                                   places=9)
        self.assertEqual(montecarlo_animation.get_scale_length(0.0), 0.0)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_the_planet_radius_agrees_with_the_other_tool(self):
        # 地球半径は animate_trajectory にもある。二重に持っている値なので、ずれたら困る
        sys.path.insert(0, os.path.join(ROOT_DIR, 'src_helper', 'animate_trajectory'))
        import animate_trajectory
        self.assertEqual(montecarlo_animation.RADIUS_PLANET, animate_trajectory.RADIUS_PLANET)

    @unittest.skipUnless(HAS_MATPLOTLIB, 'matplotlib is not installed')
    def test_an_unknown_output_format_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as directory_reference:
            reference = self.build_run(directory, directory_reference)
            completed = self.run_script([directory, '--reference', reference, '--frames', '2',
                                         '-o', os.path.join(directory, 'a.xyz')])
            self.assertNotEqual(completed.returncode, 0)

    def test_a_directory_without_a_case_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            completed = self.run_script([directory, '-o',
                                         os.path.join(directory, 'a.gif')])
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn('No case', completed.stdout)


if __name__ == '__main__':

    unittest.main()
