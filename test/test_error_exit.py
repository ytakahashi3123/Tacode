#!/usr/bin/env python3
"""
エラー経路の終了コードのテスト。

設定ミスを見つけたとき、ソルバーは "Program stopped." と出して止まる。**そのとき
終了コードが 0 でないこと**をここで固定する。

組み込みの `exit()` は `SystemExit(None)` を投げるので、終了コードは 0 になる。
つまり「止まったのに成功扱い」で、呼び出し側（シェルスクリプト、モンテカルロの
driver、CI）は失敗を検知できない。実際に wind / epoch / satellite / montecarlo の
23 か所がこの形で、風のモデル名を打ち間違えた計算が exit 0 を返していた。

そこで 2 つの角度から押さえる:

  1. **構造**: src/ に裸の `exit()` が無いこと（あれば 0 を返す経路が戻ったということ）
  2. **挙動**: 代表的な設定ミスで終了コードが 1 になること。ソルバーを実際に
     子プロセスとして起動する経路も 1 つ通す
"""

import glob
import os
import re
import subprocess
import sys
import unittest

from context import ROOT_DIR, SRC_DIR, load_config, quiet

import epoch.epoch as epoch
import wind.wind as wind


# 行として裸の exit() を呼んでいるもの（コメントや文字列の中は見ない）
PATTERN_BARE_EXIT = re.compile(r'^\s*exit\s*\(\s*\)\s*$')


class TestNoBareExitInTheSource(unittest.TestCase):

    def test_the_source_does_not_call_the_builtin_exit(self):
        offending = []
        for path in sorted(glob.glob(os.path.join(SRC_DIR, '**', '*.py'), recursive=True)):
            with open(path) as f:
                for number, line in enumerate(f, start=1):
                    if PATTERN_BARE_EXIT.match(line):
                        offending.append('{:s}:{:d}'.format(os.path.relpath(path, ROOT_DIR), number))

        self.assertEqual(offending, [],
                         'bare exit() returns the exit code 0; use sys.exit(1) in an error path')


class TestTheExitCodeOfAnErrorPath(unittest.TestCase):

    def test_an_unknown_wind_model_exits_with_one(self):
        config = load_config()
        config['wind']['flag_wind'] = True
        config['wind']['kind_wind_model'] = 'bogus'
        with quiet():
            with self.assertRaises(SystemExit) as raised:
                wind.initial_settings_wind(config)
        self.assertEqual(raised.exception.code, 1)

    def test_a_wind_velocity_of_the_wrong_length_exits_with_one(self):
        config = load_config()
        config['wind']['flag_wind'] = True
        config['wind']['kind_wind_model'] = 'constant'
        config['wind']['velocity'] = [1.0, 2.0]
        with quiet():
            with self.assertRaises(SystemExit) as raised:
                wind.initial_settings_wind(config)
        self.assertEqual(raised.exception.code, 1)

    def test_a_malformed_epoch_exits_with_one(self):
        config = load_config()
        config['epoch']['flag_epoch'] = True
        config['epoch']['datetime'] = 'yesterday afternoon'
        with quiet():
            with self.assertRaises(SystemExit) as raised:
                epoch.initial_settings_epoch(config)
        self.assertEqual(raised.exception.code, 1)

    def test_the_solver_returns_a_non_zero_exit_code(self):
        # 子プロセスとして起動する経路（シェルスクリプトとモンテカルロの driver が
        # 見ているのはこれ）。設定を読んだ時点で落ちるので計算は走らない
        import shutil
        import tempfile

        import yaml

        case = os.path.join(ROOT_DIR, 'tutorial', 'work')
        with tempfile.TemporaryDirectory() as directory:
            with open(os.path.join(case, 'config.yml')) as f:
                config = yaml.safe_load(f)
            config['wind']['flag_wind'] = True
            config['wind']['kind_wind_model'] = 'bogus'
            with open(os.path.join(directory, 'config.yml'), 'w') as f:
                yaml.safe_dump(config, f)
            # ケースは自分の database/ を参照するので、それも一緒に運ぶ
            shutil.copytree(os.path.join(case, 'database'), os.path.join(directory, 'database'))

            done = subprocess.run([sys.executable, os.path.join(SRC_DIR, 'tacode.py')],
                                  cwd=directory, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        self.assertEqual(done.returncode, 1)
        self.assertIn('Program stopped.', done.stdout.decode())


if __name__ == '__main__':
    unittest.main()
