#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/10
#
# 絶対時刻（エポック）。
#
# 従来このコードは 0 起点の経過秒（time_elapsed）しか持たない。しかし
#
#   - NCEP の風は 6 時間格子なので、時間内挿に暦日時が要る
#   - HWM14 は通日（day of year）と世界時（UT）を引数に取る
#   - 月・太陽の暦（三体摂動）も同じく暦日時が要る
#
# といずれも絶対時刻を要求するので、その共通の土台をここに置く。
#
# 既定は無効。config に epoch セクションが無い（あるいは flag_epoch: False）なら
# initial_settings_epoch は None を返し、出力を含めて従来と完全に同じ挙動になる。
# 姿勢（attitude_dict）と同じ設計で、None かどうかだけで分岐する。
#
# 時刻はすべて UTC で扱う。**閏秒は考慮しない**（UTC を等間隔の秒として扱い、
# 内部では datetime + timedelta で進める）。この用途では効かない
# （風の 6 時間格子に対して数十秒、月の位置に対して 1 秒あたり 1 km 未満）。
# 秒の定義まで問題になる計算をするなら TAI/TT との差を入れる必要がある。

import sys as sys
import datetime as datetime
from general.general import get_setting

# Dict key
KEY_DATETIME = 'datetime'

SECOND_PER_DAY = 86400.0

# ユリウス日の起点: 1858-11-17T00:00:00 UTC = MJD 0 = JD 2400000.5
JULIAN_DATE_MJD_ZERO = 2400000.5
DATETIME_MJD_ZERO    = datetime.datetime(1858, 11, 17, tzinfo=datetime.timezone.utc)


def initial_settings_epoch(config):
  #
  # config の epoch セクションを読む。
  # 戻り値は epoch_dict、既定（セクション無し / flag_epoch: False）では None。
  #
  section = get_setting(config, 'epoch', None)
  if not bool( get_setting(section, 'flag_epoch', False) ) :
    return None

  datetime_setting = get_setting(section, 'datetime', None)
  if datetime_setting is None :
    print('Epoch is enabled but epoch.datetime is not given.')
    print('--Give the UTC date and time in ISO 8601, e.g. 2026-03-21T03:00:00Z')
    print('Program stopped.')
    exit()

  epoch_dict = {KEY_DATETIME: parse_datetime(datetime_setting)}

  print('Setting the epoch...')
  print('--Epoch (UTC):', get_string(epoch_dict, 0.0))
  print('--Julian date:', '{:.6f}'.format(get_julian_date(epoch_dict, 0.0)))

  return epoch_dict


def parse_datetime(datetime_setting):
  #
  # config の値を UTC の datetime にする。
  #
  # PyYAML は引用符の無い ISO 8601 のタイムスタンプを datetime.datetime に、
  # 日付だけなら datetime.date に解決してしまうので、文字列とあわせて 3 通りを受ける。
  # PyYAML が返す datetime はオフセットを適用済みの naive な UTC なので、
  # タイムゾーンが無いものは UTC とみなしてよい（引用符つきの文字列で
  # タイムゾーンを省いた場合も、UTC で与えたものとして扱う）。
  #
  if isinstance(datetime_setting, datetime.datetime) :
    time_epoch = datetime_setting
  elif isinstance(datetime_setting, datetime.date) :
    time_epoch = datetime.datetime(datetime_setting.year, datetime_setting.month, datetime_setting.day)
  else :
    time_epoch = parse_datetime_string( str(datetime_setting) )

  if time_epoch.tzinfo is None :
    time_epoch = time_epoch.replace(tzinfo=datetime.timezone.utc)

  return time_epoch.astimezone(datetime.timezone.utc)


def parse_datetime_string(string_datetime):
  #
  # ISO 8601 の文字列を datetime にする。
  # 末尾の 'Z' は datetime.fromisoformat が Python 3.11 未満で読めないので自前で置き換える
  # （このコードは 3.9 から動かす）。
  #
  string_tmp = string_datetime.strip()
  if string_tmp.endswith('Z') or string_tmp.endswith('z') :
    string_tmp = string_tmp[:-1] + '+00:00'

  try:
    time_epoch = datetime.datetime.fromisoformat(string_tmp)
  except ValueError:
    print('Epoch is not in ISO 8601 format:', string_datetime)
    print('--Expected something like 2026-03-21T03:00:00Z')
    print('Program stopped.')
    exit()

  return time_epoch


def get_datetime(epoch_dict, time_elapsed):
  # エポックから time_elapsed 秒後の UTC
  return epoch_dict[KEY_DATETIME] + datetime.timedelta(seconds=float(time_elapsed))


def get_string(epoch_dict, time_elapsed):
  # 出力に添える ISO 8601 の文字列（UTC、ミリ秒まで）
  time_tmp = get_datetime(epoch_dict, time_elapsed)
  return time_tmp.strftime('%Y-%m-%dT%H:%M:%S') + '.{:03d}'.format(time_tmp.microsecond//1000) + 'Z'


def get_julian_date(epoch_dict, time_elapsed):
  # ユリウス日。暦（月・太陽）と、NCEP の時間軸の突き合わせに使う
  difference = get_datetime(epoch_dict, time_elapsed) - DATETIME_MJD_ZERO
  return JULIAN_DATE_MJD_ZERO + difference.total_seconds()/SECOND_PER_DAY


def get_day_of_year(epoch_dict, time_elapsed):
  # 通日（1 月 1 日が 1）。HWM14 と NRLMSISE-00 がこの形で受ける
  return get_datetime(epoch_dict, time_elapsed).timetuple().tm_yday


def get_second_of_day(epoch_dict, time_elapsed):
  # 世界時の当日通算秒。HWM14 が地方時を作るのに使う
  time_tmp = get_datetime(epoch_dict, time_elapsed)
  return float(time_tmp.hour*3600 + time_tmp.minute*60 + time_tmp.second) + time_tmp.microsecond*1.e-6
