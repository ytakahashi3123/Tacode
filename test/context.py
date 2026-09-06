#!/usr/bin/env python3
"""
テスト共通のセットアップ。

src/ を import path に載せ、テスト用の config を組み立てるヘルパーを提供する。
各テストは `from context import ...` で使う。
"""

import contextlib
import copy
import io
import os
import sys

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.normpath(os.path.join(TEST_DIR, '..'))
SRC_DIR = os.path.join(ROOT_DIR, 'src')

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

import yaml  # noqa: E402


@contextlib.contextmanager
def quiet():
    """Tacode の print を飲み込む。テスト出力を読めるようにするため。"""
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        yield buffer


def load_config(path=None):
    """リポジトリの config を読む。既定はチュートリアルケース。"""
    if path is None:
        path = os.path.join(ROOT_DIR, 'tutorial', 'work', 'config.yml')
    with open(path) as f:
        config = yaml.safe_load(f)
    return absolutize_database_path(config, os.path.dirname(os.path.abspath(path)))


def absolutize_database_path(config, case_directory):
    """
    manual 指定のデータベースパスを、config が置かれたディレクトリ基準の絶対パスにする。

    各ケースは自分の作業ディレクトリに database/ を持っており、ソルバーはそれを
    カレントディレクトリ基準で解決する。テストは cd せずに走るので、ここで
    config の場所に読み替えないとリポジトリ直下の database/ を掴んでしまう。
    """
    for section, key in (('satellite', 'directory_aerodynamic'),
                         ('atmosphere', 'directory_atmosphere')):
        block = config.get(section)
        if not block or block.get('directory_path_specify') != 'manual':
            continue
        directory = block[key]
        if not os.path.isabs(directory):
            block[key] = os.path.normpath(os.path.join(case_directory, directory))
    return config


def two_body_config():
    """
    純粋な 2 体問題にした config を返す。

    J 項・自転（コリオリ／遠心力）・空力抗力をすべて 0 にするので、
    残る力は -GM/r^2 のみ。ケプラー運動の保存量を検査するのに使う。
    """
    config = copy.deepcopy(load_config())

    for key in config['planet']['potential_factor']:
        config['planet']['potential_factor'][key] = 0.0
    config['planet']['rotation_rate'] = 0.0

    # 大気を定数モデルにして密度 0 = 抗力なし
    config['atmosphere']['kind_atmosphere_model'] = 'constant'
    config['atmosphere']['density'] = 0.0
    config['atmosphere']['temperature'] = 300.0
    config['atmosphere']['knudsen'] = 1.0
    config['satellite']['kind_aerodynamic_model'] = 'constant'
    config['satellite']['drag_coefficient'] = 0.0

    return config


def gravitational_parameter(config):
    """GM [m^3/s^2]"""
    return config['planet']['gravitational_constant'] * config['planet']['mass']
