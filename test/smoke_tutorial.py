#!/usr/bin/env python3
"""
チュートリアルケースを実際に走らせて、出力ファイルを検査する。

単体テスト（test_*.py）は関数と不変条件しか見ないので、
「tacode.py を起動して 3 種類のファイルが正しく書けるか」はここで確かめる。
6 自由度（姿勢）ケースも短くして走らせ、姿勢の列が付くことを確かめる。
KML を書くので simplekml が必要。

    python3 test/smoke_tutorial.py

作業用ディレクトリは一時領域に作るので、リポジトリの中は汚さない。
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.normpath(os.path.join(TEST_DIR, '..'))
CONFIG = os.path.join(ROOT_DIR, 'tutorial', 'work', 'config.yml')
CONFIG_6DOF = os.path.join(ROOT_DIR, 'tutorial', 'work_reentry_6dof', 'config.yml')
TACODE = os.path.join(ROOT_DIR, 'src', 'tacode.py')

# 6 自由度ケースは短く切り詰めて走らせる
TIME_MAX_6DOF = 20.0
TIMESTEP_6DOF = 0.05
EXPECTED_ROWS_6DOF = int(round(TIME_MAX_6DOF/TIMESTEP_6DOF)) + 1
EXPECTED_COLUMN_6DOF = 27
EXPECTED_ALPHA_6DOF = 20.0

KML_NS = {'k': 'http://www.opengis.net/kml/2.2'}

# チュートリアルの設定: time_elapsed_maximum 5000 s / dt 1 s / KML は 10 ステップおき
EXPECTED_ROWS = 5001
EXPECTED_LAST_TIME = 5000.0
EXPECTED_RESTART_HEADER = '# 5000 5000.0'
EXPECTED_KML_POINTS = 501

failures = []


def check(condition, message):
    print('  %-4s %s' % ('OK' if condition else 'FAIL', message))
    if not condition:
        failures.append(message)


def read_tecplot(path):
    header = None
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('zone'):
                header = line
                continue
            words = line.split()
            if len(words) >= 14 and not line.startswith(('#', 'Variables')):
                rows.append([float(w) for w in words])
    return header, rows


def run_attitude_case():
    import yaml

    print()
    print('Running the 6-DOF (attitude) case')

    with open(CONFIG_6DOF) as f:
        config = yaml.safe_load(f)
    config['computational_setup']['time_elapsed_maximum'] = TIME_MAX_6DOF
    config['time_integration']['timestep_constant'] = TIMESTEP_6DOF
    config['post_process']['tecplot']['frequency_output'] = 1
    config['restart_process']['frequency_output'] = 1

    workdir = tempfile.mkdtemp(prefix='tacode_smoke_6dof_')
    try:
        with open(os.path.join(workdir, 'config.yml'), 'w') as f:
            yaml.safe_dump(config, f)

        completed = subprocess.run([sys.executable, TACODE],
                                   cwd=workdir, capture_output=True, text=True)
        if completed.returncode != 0:
            print(completed.stdout[-2000:])
            print(completed.stderr[-2000:], file=sys.stderr)
            check(False, 'tacode.py が終了コード %d で落ちた' % completed.returncode)
            return
        check(True, 'tacode.py が正常終了した')

        tecplot = os.path.join(workdir, 'output_result', 'tecplot.dat')
        restart = os.path.join(workdir, 'output_restart', 'restart.dat')
        check(os.path.exists(tecplot), '出力が生成された: tecplot.dat')
        if not os.path.exists(tecplot):
            return

        with open(tecplot) as f:
            lines = f.read().splitlines()
        variables = [line for line in lines if line.startswith('Variables')][0]
        check(variables.rstrip().endswith('AoAtotal[deg.]'),
              'Variables 行に姿勢の列が付いている')

        header, rows = read_tecplot(tecplot)
        check(len(rows) == EXPECTED_ROWS_6DOF,
              'tecplot の行数が %d（実際 %d）' % (EXPECTED_ROWS_6DOF, len(rows)))
        check(all(len(row) == EXPECTED_COLUMN_6DOF for row in rows),
              'すべての行が %d 列' % EXPECTED_COLUMN_6DOF)

        declared = header.split('i=')[1].split()[0] if header else ''
        check(declared.isdigit() and int(declared) == len(rows),
              'zone ヘッダの点数と行数が一致')

        norm = [abs(sum(value*value for value in row[14:18])**0.5 - 1.0) for row in rows]
        check(max(norm) < 1.e-9, 'クォータニオンが単位長を保っている（最大ずれ %.2e）' % max(norm))

        check(abs(rows[0][24] - EXPECTED_ALPHA_6DOF) < 1.e-6,
              '初期迎角が %.1f 度（実際 %.6f）' % (EXPECTED_ALPHA_6DOF, rows[0][24]))
        check(max(abs(row[24]) for row in rows) < EXPECTED_ALPHA_6DOF + 0.5,
              '迎角が初期振幅を超えて発散しない')
        check(max(abs(row[25]) for row in rows) < 0.05,
              '面内に置いた初期姿勢で横滑り角が立たない')

        with open(restart) as f:
            restart_lines = [line for line in f.read().splitlines() if not line.startswith('#')]
        check(all(len(line.split()) == 13 for line in restart_lines),
              'restart が位置・速度・姿勢の 13 列')
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def main():
    if not os.path.exists(CONFIG):
        print('config が見つからない: %s' % CONFIG, file=sys.stderr)
        return 1

    workdir = tempfile.mkdtemp(prefix='tacode_smoke_')
    try:
        shutil.copy(CONFIG, os.path.join(workdir, 'config.yml'))

        print('Running the tutorial case in %s' % workdir)
        completed = subprocess.run([sys.executable, TACODE],
                                   cwd=workdir, capture_output=True, text=True)

        if completed.returncode != 0:
            print(completed.stdout[-2000:])
            print(completed.stderr[-2000:], file=sys.stderr)
            print('  FAIL tacode.py が終了コード %d で落ちた' % completed.returncode)
            return 1
        print('  OK   tacode.py が正常終了した')

        tecplot = os.path.join(workdir, 'output_result', 'tecplot.dat')
        kml = os.path.join(workdir, 'output_result', 'geodetic.kml')
        restart = os.path.join(workdir, 'output_restart', 'restart.dat')

        for path in (tecplot, kml, restart):
            check(os.path.exists(path), '出力が生成された: %s' % os.path.basename(path))
        if failures:
            return 1

        # --- Tecplot -----------------------------------------------------
        header, rows = read_tecplot(tecplot)
        check(len(rows) == EXPECTED_ROWS,
              'tecplot の行数が %d（実際 %d）' % (EXPECTED_ROWS, len(rows)))
        check(rows[0][0] == 0.0, 'tecplot の先頭が t=0.0')
        check(rows[-1][0] == EXPECTED_LAST_TIME,
              'tecplot の末尾が t=%.1f（実際 %.1f）' % (EXPECTED_LAST_TIME, rows[-1][0]))

        declared = header.split('i=')[1].split()[0] if header else ''
        check(declared.isdigit(), 'zone ヘッダの点数が整数（実際 "%s"）' % declared)
        check(declared.isdigit() and int(declared) == len(rows),
              'zone ヘッダの点数と行数が一致')

        altitudes = [r[6] for r in rows]
        check(all(a > 0.0 for a in altitudes), '全ステップで高度が正')
        densities = [r[11] for r in rows]
        check(all(d > 0.0 for d in densities), '全ステップで密度が正')

        # --- restart -----------------------------------------------------
        with open(restart) as f:
            restart_lines = f.read().splitlines()
        check(restart_lines[2] == EXPECTED_RESTART_HEADER,
              'restart ヘッダが "%s"（実際 "%s"）' % (EXPECTED_RESTART_HEADER, restart_lines[2]))
        check(len(restart_lines) == 3 + EXPECTED_ROWS,
              'restart のデータ行数が %d（実際 %d）' % (EXPECTED_ROWS, len(restart_lines) - 3))

        # --- KML ---------------------------------------------------------
        root = ET.parse(kml).getroot()
        linestrings = root.findall('.//k:LineString', KML_NS)
        check(len(linestrings) == 1,
              'KML の LineString が 1 本（実際 %d 本）' % len(linestrings))

        if linestrings:
            coords = linestrings[0].find('k:coordinates', KML_NS).text.split()
            check(len(coords) == EXPECTED_KML_POINTS,
                  'LineString の点数が %d（実際 %d）' % (EXPECTED_KML_POINTS, len(coords)))
            check(len(coords) >= 2, '線分が描画される点数がある')

        colors = [c.text for c in root.findall('.//k:color', KML_NS)]
        check(bool(colors) and all(re.fullmatch(r'[0-9a-fA-F]{8}', c or '') for c in colors),
              'KML の色が aabbggrr の 16 進（実際 %s）' % colors)

        run_attitude_case()

        print()
        if failures:
            print('%d 件失敗' % len(failures))
            return 1
        print('スモークテスト成功')
        return 0
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main())
