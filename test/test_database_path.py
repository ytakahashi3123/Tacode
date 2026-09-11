#!/usr/bin/env python3
"""
データベース（大気・空力・風）のパス解決と、ケースが持つテーブルの素性のテスト。

パス解決は `general.get_database_directory` に一本化してある。以前は 3 モジュールが
別々に `if auto/default ... elif manual ... else` を書いており、しかも else の既定が
モジュールごとに違っていた（片方は `database/aerodynamice` という綴り誤り）。
綴りを間違えた `directory_path_specify` が黙って既定へ落ちると、ケースが自分の
`database/` に置いたテーブルではなくリポジトリのマスターを読む。

ケース側の `database/` はマスターの複製という約束なので、その一致もここで見る。
"""

import glob
import hashlib
import os
import tempfile
import unittest

from context import ROOT_DIR, load_config, quiet

import atmosphere.atmosphere as atmosphere
import satellite.satellite as satellite
import wind.wind as wind
from general.general import DIRECTORY_DATABASE, get_database_directory


# ケースが複製して持つテーブルのうち、マスターと違っていてよいもの（今は無し）
CASE_TABLE_EXCEPTION = ()


class TestPathResolution(unittest.TestCase):

    def test_default_and_auto_point_at_the_master(self):
        for specify in ('default', 'auto'):
            for name in ('atmosphere', 'aerodynamic', 'wind'):
                section = {'directory_path_specify': specify}
                directory = get_database_directory(section, name, name, 'directory_' + name)
                self.assertEqual(os.path.normpath(directory),
                                 os.path.join(DIRECTORY_DATABASE, name))

    def test_a_missing_section_falls_back_to_the_master(self):
        # 省略可能なセクション（wind）では section そのものが None で来る
        directory = get_database_directory(None, 'wind', 'wind', 'directory_wind', 'database/wind')
        self.assertEqual(os.path.normpath(directory), os.path.join(DIRECTORY_DATABASE, 'wind'))

    def test_the_master_directory_exists(self):
        for name in ('atmosphere', 'aerodynamic', 'wind'):
            self.assertTrue(os.path.isdir(os.path.join(DIRECTORY_DATABASE, name)), name)

    def test_manual_is_taken_as_given(self):
        section = {'directory_path_specify': 'manual', 'directory_atmosphere': 'database/atmosphere'}
        directory = get_database_directory(section, 'atmosphere', 'atmosphere', 'directory_atmosphere')
        self.assertEqual(directory, 'database/atmosphere')

    def test_a_misspelt_specification_stops_the_run(self):
        for specify in ('Manual', 'MANUAL', 'defalut', 'file', ''):
            section = {'directory_path_specify': specify, 'directory_atmosphere': 'database/atmosphere'}
            with self.assertRaises(SystemExit) as raised:
                with quiet():
                    get_database_directory(section, 'atmosphere', 'atmosphere', 'directory_atmosphere')
            self.assertEqual(raised.exception.code, 1)

    def test_manual_without_the_directory_key_stops_the_run(self):
        section = {'directory_path_specify': 'manual'}
        with self.assertRaises(SystemExit) as raised:
            with quiet():
                get_database_directory(section, 'atmosphere', 'atmosphere', 'directory_atmosphere')
        self.assertEqual(raised.exception.code, 1)

    def test_the_wind_keeps_its_default_directory(self):
        # 風だけは manual のときの既定値を持つ（README に書いてある挙動）
        section = {'directory_path_specify': 'manual'}
        directory = get_database_directory(section, 'wind', 'wind', 'directory_wind', 'database/wind')
        self.assertEqual(directory, 'database/wind')


class TestPathResolutionThroughTheModules(unittest.TestCase):
    """3 つのモジュールが同じ解決を通っていること。"""

    def test_the_atmosphere_stops_on_a_misspelt_specification(self):
        config = load_config()
        config['atmosphere']['directory_path_specify'] = 'Manual'
        with self.assertRaises(SystemExit):
            with quiet():
                atmosphere.read_atmosphere_file(config)

    def test_the_aerodynamics_stops_on_a_misspelt_specification(self):
        config = load_config()
        config['satellite']['directory_path_specify'] = 'Manual'
        with self.assertRaises(SystemExit):
            with quiet():
                satellite.read_aerodynamic_file(config)

    def test_the_wind_stops_on_a_misspelt_specification(self):
        config = load_config()
        config['wind']['directory_path_specify'] = 'Manual'
        config['wind']['kind_wind_model'] = 'fileread'
        with self.assertRaises(SystemExit):
            with quiet():
                wind.read_wind_file(config)

    def test_default_reads_the_master_table(self):
        # ケースの database/ を通さずにマスターを読めること
        config = load_config()
        config['atmosphere']['directory_path_specify'] = 'default'
        config['atmosphere']['filename_atmosphere'] = 'atmospheremodel.txt'
        with quiet():
            atmosphere_dict = atmosphere.read_atmosphere_file(config)
        self.assertGreater(atmosphere_dict[atmosphere.KEY_DATA], 0)


class TestCaseTablesMatchTheMaster(unittest.TestCase):
    """
    各ケースの `database/` はマスターの複製であること。

    ケースごとに複製する方針そのものは意図したもの（ケースが自分の使うテーブルを
    持ち運ぶ）だが、片方だけ直すと静かにずれる。
    """

    @staticmethod
    def _digest(path):
        with open(path, 'rb') as f:
            return hashlib.md5(f.read()).hexdigest()

    def test_every_case_table_is_identical_to_the_master(self):
        pattern = os.path.join(ROOT_DIR, '*', 'database', '*', '*.txt')
        pattern_case = os.path.join(ROOT_DIR, '*', '*', 'database', '*', '*.txt')
        path_list = sorted(glob.glob(pattern) + glob.glob(pattern_case))
        path_list = [path for path in path_list
                     if os.path.relpath(path, ROOT_DIR).split(os.sep)[0] != 'database']
        self.assertGreater(len(path_list), 0, 'no case table was found')

        for path in path_list:
            relative = os.path.relpath(path, ROOT_DIR)
            if relative in CASE_TABLE_EXCEPTION:
                continue
            name_database, name_file = relative.split(os.sep)[-2:]
            path_master = os.path.join(DIRECTORY_DATABASE, name_database, name_file)
            self.assertTrue(os.path.exists(path_master),
                            relative + ' has no master in database/' + name_database)
            self.assertEqual(self._digest(path), self._digest(path_master),
                             relative + ' differs from the master copy')


if __name__ == '__main__':
    unittest.main()
