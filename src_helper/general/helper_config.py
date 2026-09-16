#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2026/09/10
#
# src_helper の後処理ツールの設定を YAML で与えるための共通部分。
#
# ソルバーが config.yml を読むのと同じ流儀にしてある: 既定では**カレント
# ディレクトリの config_helper.yml** を読み、-file で別のファイルを指せる。
# 1 つのファイルにツールごとのセクションを並べる:
#
#   animate_trajectory:
#     filename: output_result/tecplot.dat
#     output: attitude.mp4
#
#   montecarlo_animation:
#     directory: work_montecarlo_wind
#     view: 3d
#
# 優先順位は **コマンドライン > YAML > 既定値**。長い引数の列を毎回打たなくて
# 済むようにするのが目的で、その場限りの上書きだけコマンドラインでやる。
#
# キーは argparse の dest 名（--altitude-max なら altitude_max）。綴りを間違えた
# キーは黙って無視せずに停止する（設定したつもりで効いていない、が一番困るため）。
#
# --save-config で、いま効いている設定をそのまま YAML に書き出せる。既存の
# ファイルがあれば自分のセクションだけ差し替える（**コメントは残らない**ので、
# 手で書いたファイルに上書きするときは注意）。

import argparse
import os
import sys

import yaml

# 既定の設定ファイル名（カレントディレクトリから読む。ソルバーの config.yml と同じ考え方）
FILENAME_DEFAULT = 'config_helper.yml'

# 設定ファイルにも --save-config にも現れないキー（ツールの動きではなく、設定の与え方そのもの）
KEY_INTERNAL = ('file', 'save_config')


def add_argument(parser):
  #
  # 設定ファイルまわりの引数を足す。各ツールの argument() から呼ぶ。
  #
  parser.add_argument('-file', '--file', dest='file', default=FILENAME_DEFAULT,
                      help='settings file, read from the current directory '
                           '(default: %(default)s). Missing is not an error: the command '
                           'line and the defaults are then used as they are')
  parser.add_argument('--save-config', dest='save_config', nargs='?', const='',
                      default=None, metavar='PATH',
                      help='write the settings in effect to this file (default: the file '
                           'given by -file) and exit. Other sections of an existing file '
                           'are kept, but its comments are not')

  return parser


def read_file(filename):
  #
  # 設定ファイル全体を読む。無ければ空の辞書（設定ファイルは必須ではない）。
  #
  if not os.path.exists(filename) :
    return {}

  with open(filename) as f:
    content = yaml.safe_load(f)

  if content is None :
    return {}
  if not isinstance(content, dict) :
    print('The settings file must be a mapping of sections:', filename)
    sys.exit(1)

  return content


def get_given(parser):
  #
  # コマンドラインで**実際に与えられた**ものだけを返す。
  #
  # 既定値を一時的に argparse.SUPPRESS にして同じ引数を読み直すと、与えられなかった
  # 項目は属性そのものが現れない。既定値との比較では「既定値と同じ値を明示した」
  # 場合を拾えないので、こちらで判定する。parser は呼び出し元のものなので、
  # 既定値は必ず戻す。
  #
  default_list = [action.default for action in parser._actions]
  try:
    for action in parser._actions:
      action.default = argparse.SUPPRESS
    given = vars(parser.parse_args())
  finally:
    for action, value in zip(parser._actions, default_list):
      action.default = value

  return given


def get_setting(parser, section):
  #
  # コマンドライン・設定ファイル・既定値をこの順の優先度でまとめ、
  # argparse.Namespace にして返す。
  #
  # 「コマンドラインで与えられたか」は、既定値を伏せた読み直しで判定する
  # （get_given）。既定値と同じ値を明示したときも、設定ファイルではなく
  # コマンドラインが勝つ。
  #
  argument = parser.parse_args()
  default  = vars(parser.parse_args([]))

  given = {key: value for key, value in get_given(parser).items()
           if key not in KEY_INTERNAL}

  content = read_file(argument.file)
  setting = content.get(section, {})
  if setting is None :
    setting = {}
  if not isinstance(setting, dict) :
    print('Section', section, 'of the settings file must be a mapping:', argument.file)
    sys.exit(1)

  unknown = [key for key in setting if key not in default or key in KEY_INTERNAL]
  if len(unknown) > 0 :
    print('Unknown setting in section', section, 'of', argument.file, ':', ', '.join(sorted(unknown)))
    print('--Available:', ', '.join(sorted(key for key in default if key not in KEY_INTERNAL)))
    sys.exit(1)

  value = dict(default)
  value.update(setting)
  value.update(given)
  value['file']        = argument.file
  value['save_config'] = argument.save_config

  return argparse.Namespace(**value)


def save_file(argument, section):
  #
  # --save-config。いま効いている設定を書き出して True を返す（呼んだ側はそこで終える）。
  # 指定が無ければ False。
  #
  if argument.save_config is None :
    return False

  filename = argument.save_config if argument.save_config != '' else argument.file

  content = read_file(filename)
  content[section] = {key: value for key, value in vars(argument).items()
                      if key not in KEY_INTERNAL}

  with open(filename, 'w') as f:
    f.write('# Settings for the post-processing tools in src_helper/\n')
    f.write('# --One section per tool; the keys are the long options without the dashes.\n')
    f.write('# --The command line overrides what is written here.\n')
    yaml.safe_dump(content, f, default_flow_style=False, sort_keys=True,
                   allow_unicode=True)

  print('Written: ', filename)

  return True


def require(argument, name, section):
  #
  # 設定ファイルにもコマンドラインにも無かった必須の項目を報告して止める。
  #
  value = getattr(argument, name, None)
  if value is not None :
    return value

  print('Give', name, 'on the command line, or as', section+'.'+name,
        'in the settings file:', argument.file)
  sys.exit(1)
