#!/usr/bin/env python3
"""
エポック（絶対時刻）のテスト。

`epoch` は風（NCEP の時間内挿、HWM14 の通日と世界時）と三体摂動の共通の土台で、
ここを取り違えると全部が静かにずれる。ユリウス日は astropy と突き合わせる
（あれば。無ければ skip する。`src_helper` の matplotlib と同じ任意依存の扱い）。

既定が無効であること、無効なら出力が従来と 1 バイトも変わらないことも見る。
"""

import datetime
import unittest

from context import quiet

import epoch.epoch as epoch


def make_epoch(value):
    with quiet():
        return epoch.initial_settings_epoch({'epoch': {'flag_epoch': True, 'datetime': value}})


class TestTheEpochIsOffByDefault(unittest.TestCase):

    def test_a_config_without_the_section_gives_none(self):
        self.assertIsNone(epoch.initial_settings_epoch({}))

    def test_an_explicit_false_gives_none(self):
        config = {'epoch': {'flag_epoch': False, 'datetime': '2026-03-21T03:00:00Z'}}
        self.assertIsNone(epoch.initial_settings_epoch(config))

    def test_a_none_section_gives_none(self):
        self.assertIsNone(epoch.initial_settings_epoch({'epoch': None}))

    def test_the_tutorial_configs_leave_it_off(self):
        # 参照出力を壊さないための条件。ここが True になった時点で回帰テストが落ちる
        import glob
        import os

        from context import ROOT_DIR, load_config

        paths = sorted(glob.glob(os.path.join(ROOT_DIR, 'tutorial', '*', 'config.yml')))
        paths.append(os.path.join(ROOT_DIR, 'tutorial', 'template', 'config.yml'))
        paths.append(os.path.join(ROOT_DIR, 'src', 'config.yml'))
        self.assertTrue(len(paths) >= 6)
        for path in paths:
            config = load_config(path)
            self.assertIn('epoch', config, path)
            self.assertFalse(config['epoch']['flag_epoch'], path)
            self.assertIsNone(epoch.initial_settings_epoch(config), path)


class TestParsingTheEpoch(unittest.TestCase):

    def test_an_iso_string_with_a_z_is_read_as_utc(self):
        # fromisoformat は Python 3.11 未満で 'Z' を読めないので自前で処理している
        epoch_dict = make_epoch('2026-03-21T03:00:00Z')
        self.assertEqual(epoch.get_string(epoch_dict, 0.0), '2026-03-21T03:00:00.000Z')

    def test_a_string_without_a_time_zone_is_taken_as_utc(self):
        epoch_dict = make_epoch('2026-03-21T03:00:00')
        self.assertEqual(epoch.get_string(epoch_dict, 0.0), '2026-03-21T03:00:00.000Z')

    def test_an_offset_is_converted_to_utc(self):
        epoch_dict = make_epoch('2026-03-21T12:00:00+09:00')
        self.assertEqual(epoch.get_string(epoch_dict, 0.0), '2026-03-21T03:00:00.000Z')

    def test_a_datetime_object_from_yaml_is_accepted(self):
        # PyYAML は引用符の無いタイムスタンプを datetime に解決してしまう
        import yaml

        config = yaml.safe_load('epoch:\n  flag_epoch: True\n  datetime: 2026-03-21T03:00:00Z\n')
        self.assertIsInstance(config['epoch']['datetime'], datetime.datetime)
        with quiet():
            epoch_dict = epoch.initial_settings_epoch(config)
        self.assertEqual(epoch.get_string(epoch_dict, 0.0), '2026-03-21T03:00:00.000Z')

    def test_a_yaml_timestamp_with_an_offset_is_accepted(self):
        import yaml

        config = yaml.safe_load('epoch:\n  flag_epoch: True\n  datetime: 2026-03-21T12:00:00+09:00\n')
        with quiet():
            epoch_dict = epoch.initial_settings_epoch(config)
        self.assertEqual(epoch.get_string(epoch_dict, 0.0), '2026-03-21T03:00:00.000Z')

    def test_a_bare_date_is_accepted_as_midnight(self):
        import yaml

        config = yaml.safe_load('epoch:\n  flag_epoch: True\n  datetime: 2026-03-21\n')
        self.assertIsInstance(config['epoch']['datetime'], datetime.date)
        with quiet():
            epoch_dict = epoch.initial_settings_epoch(config)
        self.assertEqual(epoch.get_string(epoch_dict, 0.0), '2026-03-21T00:00:00.000Z')

    def test_a_malformed_string_stops_the_run(self):
        with self.assertRaises(SystemExit):
            with quiet():
                epoch.initial_settings_epoch({'epoch': {'flag_epoch': True, 'datetime': 'yesterday'}})

    def test_a_missing_datetime_stops_the_run(self):
        with self.assertRaises(SystemExit):
            with quiet():
                epoch.initial_settings_epoch({'epoch': {'flag_epoch': True}})


class TestAdvancingTheEpoch(unittest.TestCase):

    def test_the_elapsed_time_is_added(self):
        epoch_dict = make_epoch('2026-03-21T03:00:00Z')
        self.assertEqual(epoch.get_string(epoch_dict, 3600.0), '2026-03-21T04:00:00.000Z')

    def test_a_fractional_second_survives(self):
        # 姿勢の dt = 0.05 s を刻むので、秒未満が丸められては困る
        epoch_dict = make_epoch('2026-03-21T03:00:00Z')
        self.assertEqual(epoch.get_string(epoch_dict, 0.05), '2026-03-21T03:00:00.050Z')

    def test_the_day_rolls_over(self):
        epoch_dict = make_epoch('2026-12-31T23:59:00Z')
        self.assertEqual(epoch.get_string(epoch_dict, 120.0), '2027-01-01T00:01:00.000Z')
        self.assertEqual(epoch.get_day_of_year(epoch_dict, 120.0), 1)


class TestTheDerivedQuantities(unittest.TestCase):

    def test_the_day_of_year_counts_from_one(self):
        self.assertEqual(epoch.get_day_of_year(make_epoch('2026-01-01T00:00:00Z'), 0.0), 1)

    def test_the_day_of_year_accounts_for_a_leap_year(self):
        # 2024 は閏年、2026 はそうでない。HWM14 に渡す通日がここでずれると季節が 1 日動く
        self.assertEqual(epoch.get_day_of_year(make_epoch('2024-03-01T00:00:00Z'), 0.0), 61)
        self.assertEqual(epoch.get_day_of_year(make_epoch('2026-03-01T00:00:00Z'), 0.0), 60)

    def test_the_second_of_day_is_universal_time(self):
        epoch_dict = make_epoch('2026-03-21T12:34:56Z')
        self.assertAlmostEqual(epoch.get_second_of_day(epoch_dict, 0.0), 12*3600 + 34*60 + 56, places=9)

    def test_the_second_of_day_is_taken_after_the_offset_is_removed(self):
        epoch_dict = make_epoch('2026-03-21T12:00:00+09:00')
        self.assertAlmostEqual(epoch.get_second_of_day(epoch_dict, 0.0), 3*3600, places=9)

    def test_the_julian_date_of_the_reference_epoch(self):
        # J2000.0 = 2000-01-01T12:00:00 TT。UTC で見た 2000-01-01T12:00:00 の JD は 2451545.0
        epoch_dict = make_epoch('2000-01-01T12:00:00Z')
        self.assertAlmostEqual(epoch.get_julian_date(epoch_dict, 0.0), 2451545.0, places=9)

    def test_the_julian_date_advances_by_one_per_day(self):
        epoch_dict = make_epoch('2026-03-21T03:00:00Z')
        difference = epoch.get_julian_date(epoch_dict, 86400.0) - epoch.get_julian_date(epoch_dict, 0.0)
        self.assertAlmostEqual(difference, 1.0, places=12)

    def test_the_julian_date_matches_astropy(self):
        try:
            from astropy.time import Time
        except ImportError:
            self.skipTest('astropy is not installed')
        for string in ('2026-03-21T03:00:00Z', '1999-12-31T23:59:59Z', '2024-02-29T06:30:00Z'):
            epoch_dict = make_epoch(string)
            reference = Time(string.replace('Z', ''), scale='utc').jd
            self.assertAlmostEqual(epoch.get_julian_date(epoch_dict, 0.0), reference, places=9, msg=string)


if __name__ == '__main__':
    unittest.main()
