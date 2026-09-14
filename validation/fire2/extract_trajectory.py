#!/usr/bin/env python3
#
# Project Fire flight II の再突入軌道（NASA TN D-3569 表 V）を読み取る。
#
# 表 V は 1617.75 s から 1929.75 s まで 0.5 秒刻み 625 点の**数値表**である。
# ただし手に入るのは 1966 年の走査版で、埋め込まれている OCR が壊れている:
# 桁が落ち、`E` が `EUR` に、`.` が `-` になり、数字 1 個ずつに砕けた行もある。
# **黙って読み飛ばすと、実際より密な表を返したことになる**ので、この読み取り器は
# 読めなかった値を数えて報告し、その場所には NaN を書く。
#
# 配置は文字の並び順ではなく `pdftotext -bbox-layout` の**座標**から組み直す。
# 行の間隔が開くと pdftotext は 1 レコードを列方向に歩くので、テキスト出力の
# 並び順は表の並び順ではない。1 レコードは 7 列 × 3 行:
#
#   1 行目 | t        緯度     経度     高度     速度     経路角   方位角
#   2 行目 | t - t_0  動圧     動圧     気圧     気圧     密度     密度
#          |          (N/m2)  (lb/ft2) (N/m2)  (lb/ft2)  (SI)    (slug/ft3)
#   3 行目 |          Reynolds Mach     加速度   温度
#
# **3 列は SI と英単位で二度印字されている。**これは同じ数の独立した 2 回目の
# 読み取りなので、片方が壊れていればもう片方から厳密に復元でき、両方生きていれば
# 相互検証になる。片割れの無い列（緯度・経度・高度・速度・経路角・方位角・温度・
# 加速度・Mach）は前後のレコードで押さえる。0.5 秒刻みの滑らかな軌道なので、
# 桁の落ちた値は近傍への当てはめから外れる。
#
# **内挿で埋めることはしない。**読めなかった値は NaN のまま出力し、報告に載せる。
#
# 使い方:
#   python3 extract_trajectory.py                    # reference/ を作り直す
#   python3 extract_trajectory.py --report-damage    # 壊れた値を 1 個ずつ並べる

import argparse
import math
import os
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))

# **参照 PDF はリポジトリの外**にある。読み取り済みの reference/ を置いてあるので
# 通常はこのスクリプトを走らせる必要は無い。
FILE_REPORT_DEFAULT = os.path.join(
  DIRECTORY_SCRIPT, '../../../references/pdf',
  'NASA-TN-D-3569_Lewis-Scallion_1966_Fire-II-Flight-Parameters.pdf')
FILE_OUTPUT_DEFAULT = os.path.join(DIRECTORY_SCRIPT, 'reference/fire2_trajectory.dat')

TIME_FIRST   = 1617.75    # 最初のレコード。高度 121 920 m の突入点
TIME_STEP    = 0.5
RECORD_COUNT = 625        # 1929.75 s まで

# 7 列の格子の中心（PDF の点）。ページごとに数点ずれるので許容幅を持たせる
COLUMN_CENTER    = [138.0, 195.0, 272.0, 351.0, 430.0, 506.0, 586.0]
COLUMN_TOLERANCE = 22.0

# これより近い断片は同じ数の一部。砕けた数の内部の隙間は 12 点を超えず、
# 列と列の隙間は 15 点を下回らない。**列の先頭に来た断片は隙間に関わらず別の数**
FRAGMENT_GAP_MAXIMUM = 13.0

# 指数部は必ず独立した語で、その手前だけ隙間が広いことがある。ここまでは繋ぐ
# （それでも列と列の最小の隙間より狭い）
EXPONENT_GAP_MAXIMUM = 14.5

# レコードの 3 行は 7.5 点間隔、次のレコードは 27.5 点下。
#
# **固定のしきい値で行を分けてはいけない。**走査が傾いているページがあり、
# そこでは同じ行でも右端の列が左端より 4〜6 点上に来る。行の間隔 7.5 点と
# 同じ桁なので、しきい値で切ると行が 1 つずれ、密度が経路角の列に、
# 温度が気圧の列に入る（実際に起きた）。代わりに**レコードごとに縦位置を
# まとめてから**、まとまりの並び順で 3 行に割り当てる。
ROW_SPACING     = 7.5     # 行の間隔（点）
ROW_GAP_MINIMUM = 1.5     # これ以上離れた縦位置は別の行
ROW_ABOVE_LIMIT = 6.5     # 時刻より上に来てよい量（傾きの分）
ROW_BELOW_LIMIT = 24.0    # これより下はどのレコードにも属さない

# 二度印字されている列の単位換算
PASCAL_PER_PSF   = 47.880259       # N/m2 / (lb/ft2)
KGM3_PER_SLUGFT3 = 515.378818      # kg/m3 / (slug/ft3)

# 両単位で読めた値が食い違いと見なされる相対差。表は 8 桁で、報告書は換算に
# 使った係数を書いていないので、その分の余裕を持たせてある
TOLERANCE_UNIT = 1.0e-4

# 数字と紛らわしい文字。**形の上でその数字にしかなり得ないものだけ**を入れる
DIGIT_LOOKALIKE = {'O': '0', 'o': '0', 'Q': '0', 'D': '0',
                   'l': '1', 'I': '1', 'i': '1', '|': '1', 'L': '1', '!': '1',
                   'Z': '2', 'z': '2',
                   'S': '5', 's': '5',
                   'G': '6', 'b': '6',
                   'T': '7', 'r': '7',
                   'B': '8',
                   'g': '9', 'q': '9'}

# 指数の目印
EXPONENT_LOOKALIKE = ['E', '\u20ac', '\u00a3', '\u00a9']

# 表の値はすべてこの形をしている: 符号・`0.`・**8 桁**・`E`・符号つき 2 桁の指数。
# 無傷で出てくる 8382 個のうち 8381 個がちょうど 8 桁なので、**8 桁を要求すること
# 自体が桁落ちの検出になる**（1 桁落ちた値は、緩めれば 10 分の 1 の数として通る）。
# match ではなく search を使うのは、数の手前に残った走査のゴミを落とすため
PATTERN_VALUE = re.compile(r'(-?)0\.(\d{8})E([+-]?)(\d{1,2})$')
PATTERN_TIME  = re.compile(r'^(1[6-9]\d\d)\.(\d\d)$')

# 出力する量。(名前, 行, 列, 片割れの行, 片割れの列, SI への係数)
FIELDS = [('latitude',          0, 1, None, None, None),
          ('longitude',         0, 2, None, None, None),
          ('altitude',          0, 3, None, None, None),
          ('velocity',          0, 4, None, None, None),
          ('flightpath',        0, 5, None, None, None),
          ('heading',           0, 6, None, None, None),
          ('dynamic_pressure',  1, 1,    1,    2, PASCAL_PER_PSF),
          ('pressure',          1, 3,    1,    4, PASCAL_PER_PSF),
          ('density',           1, 5,    1,    6, KGM3_PER_SLUGFT3),
          ('reynolds',          2, 1, None, None, None),
          ('mach',              2, 2, None, None, None),
          ('acceleration',      2, 3, None, None, None),
          ('temperature',       2, 4, None, None, None)]

FIELDS_OUTPUT = ['latitude', 'longitude', 'altitude', 'velocity', 'flightpath',
                 'heading', 'dynamic_pressure', 'pressure', 'density',
                 'temperature', 'mach', 'acceleration', 'reynolds']

# 近傍に当てた 2 次式から、当てはめの残差の OUTLIER_SIGMA 倍**かつ**下の下限を
# 超えて外れた値は捨てる。下限は、ばらつきの無い列が自分の丸め誤差で引っかからない
# ようにするためのもの
OUTLIER_SIGMA  = 12.0
OUTLIER_WINDOW = 4        # 当てはめに使う前後のレコード数
OUTLIER_PASSES = 4        # 1 個消すと隣が見えるので繰り返す
OUTLIER_FLOOR  = {'latitude': 1.0e-3, 'longitude': 1.0e-3, 'altitude': 5.0,
                  'velocity': 1.0, 'flightpath': 1.0e-2, 'heading': 1.0e-2,
                  'temperature': 2.0, 'mach': 5.0e-2, 'acceleration': 5.0e-3,
                  'reynolds': 2.0e-3, 'dynamic_pressure': 2.0e-3,
                  'pressure': 2.0e-3, 'density': 2.0e-3}

# 大気の列は突入の間に 10 桁動き、値ではなく**対数が**滑らかなので、そちらで当てる。
# 上の下限はこの 4 列では相対値になる
OUTLIER_LOGARITHMIC = ('density', 'pressure', 'dynamic_pressure', 'reynolds')

# 指数が落ちると値は桁ごと飛び、それが数レコード続くと上の当てはめは効かない
# （窓の中に壊れた値が何個も入り、残差が膨らんで何も外れて見えなくなる）。
# そこで**符号の変わらない列は先に近傍の中央値と比べる**。少数派の壊れた値は
# 中央値を動かせないので、桁の飛んだ値はここで落ちる
MAGNITUDE_SCREEN = ('altitude', 'velocity', 'temperature', 'mach', 'reynolds',
                    'density', 'pressure', 'dynamic_pressure', 'latitude',
                    'longitude', 'flightpath', 'heading')
MAGNITUDE_FACTOR = 3.0

# 中央値を取る窓は当てはめの窓より広くする。列の割り当てがずれると 5 レコードほど
# 続けて壊れることがあり、前後 4 点では壊れた側が多数派になって中央値が動く
MAGNITUDE_WINDOW = 10

# 符号を直すのは、窓の中の値のこの割合以上が逆符号で、かつ大きさの比が
# この範囲に入るとき（末端で値が 4 個未満しか残っていなければ全員一致を要求する）
SIGN_FRACTION  = 0.85
SIGN_MAGNITUDE = (0.5, 2.0)

# 両単位が食い違ったとき、近傍が推す側を採るのは、他方がこの倍率以上遠いときだけ
DISPUTE_MARGIN = 10.0

XHTML = '{http://www.w3.org/1999/xhtml}'


def run_pdftotext(file_report):
  # 報告書の座標つきテキスト（XML）を得る
  if not os.path.exists(file_report) :
    print('No such report:', file_report)
    print('--The scanned report lives outside the repository. Give its path with')
    print('  --report, or keep the reference/ file already in the case.')
    print('Program stopped.')
    sys.exit(1)

  handle, path = tempfile.mkstemp(suffix='.xml')
  os.close(handle)
  try:
    subprocess.run(['pdftotext', '-bbox-layout', file_report, path], check=True)
    with open(path, encoding='utf-8', errors='replace') as stream:
      return stream.read()
  except (OSError, subprocess.CalledProcessError) as instance:
    print('Cannot read the report with pdftotext:', instance)
    print('--Install poppler-utils.')
    print('Program stopped.')
    sys.exit(1)
  finally:
    os.unlink(path)


def read_pages(text_xml):
  # ページごとの語の一覧 (xmin, xmax, y, 文字列)
  root      = ET.fromstring(text_xml)
  page_list = []
  for page in root.iter(XHTML + 'page'):
    word_list = []
    for word in page.iter(XHTML + 'word'):
      if word.text is None :
        continue
      word_list.append((float(word.get('xMin')), float(word.get('xMax')),
                        float(word.get('yMin')), word.text.strip()))
    page_list.append(word_list)
  return page_list


def repair(text):
  #
  # 文字列を表の値に直す。直せなければ None を返す。
  # **値が一意に決まる置き換えしかしない**。その数字にしかなり得ない文字を
  # 数字にし、紛れ込んだ空白を取るだけで、それでも型に合わない語は捨てる
  # （呼び出し側が片割れの単位を見に行くか、NaN を書く）。
  #
  text = text.replace(' ', '')
  for marker in EXPONENT_LOOKALIKE:
    text = text.replace(marker, 'E')
  if not text :
    return None

  sign = ''
  if text[0] in '-\u2014\u2013' :
    sign, text = '-', text[1:]
  if not text :
    return None

  # 数の手前に走査が残したゴミ（"I0.4617" の I、"'0.1135" の '）。
  # 仮数は必ず `0.` で始まるので、その前にあるものは値の一部ではない
  opening = text.find('0.')
  if opening > 0 and not text[:opening].isdigit() :
    text = text[opening:]

  head = DIGIT_LOOKALIKE.get(text[0], text[0])
  if head != '0' :
    return None
  rest = text[1:]
  if not rest :
    return None
  # 2 文字目は必ず小数点。ここに来た `-` や `,` は負号ではなく小数点の読み違い
  if rest[0] in '-,\u00b7:;' :
    rest = '.' + rest[1:]
  if rest[0] != '.' :
    return None

  body  = rest[1:]
  split = body.find('E')
  if split < 0 :
    return None
  mantissa = ''.join(DIGIT_LOOKALIKE.get(c, c) for c in body[:split])
  exponent = ''.join(DIGIT_LOOKALIKE.get(c, c) for c in body[split + 1:])
  if exponent[:1] in ('\u2014', '\u2013') :
    exponent = '-' + exponent[1:]

  match = PATTERN_VALUE.search('%s0.%sE%s' % (sign, mantissa, exponent))
  if match is None :
    return None
  return float('%s0.%sE%s%s' % (match.group(1), match.group(2),
                                match.group(3) or '+', match.group(4)))


def repair_time(text):
  # 時刻の列は指数表記ではなく普通の小数
  text = ''.join(DIGIT_LOOKALIKE.get(c, c) if c != '.' else c
                 for c in text.replace(' ', ''))
  if PATTERN_TIME.match(text) is None :
    return None
  return float(text)


def column_of(xmin):
  # 値が始まる列の番号。どの列にも属さなければ None
  column   = None
  distance = COLUMN_TOLERANCE
  for index, center in enumerate(COLUMN_CENTER):
    if abs(xmin - center) < distance :
      column, distance = index, abs(xmin - center)
  return column


def starts_a_column(xmin, xmin_current):
  # 途中の断片が、いま組み立てている列とは別の列の先頭かどうか
  column = column_of(xmin)
  return column is not None and column != column_of(xmin_current)


def is_exponent(head, tail, gap):
  #
  # 断片が直前の数の指数部かどうか。指数は必ず独立した語で、走査が仮数の中より
  # 広い隙間を空けることがある。`E` の直後にしか来ず 1〜2 桁にしかならないので、
  # （別の列の先頭でないことは既に除いてある）その隙間は跨いでよい
  #
  if gap > EXPONENT_GAP_MAXIMUM :
    return False
  if not head or head[-1] not in EXPONENT_LOOKALIKE :
    return False
  digit = ''.join(DIGIT_LOOKALIKE.get(c, c) for c in tail.lstrip('+-'))
  return 1 <= len(digit) <= 2 and digit.isdigit()


def merge_fragments(word_list):
  # 1 行の語を、隣り合う断片を繋いで (xmin, y, 文字列) の値にまとめる
  value_list = []
  current    = None
  for xmin, xmax, y, text in word_list:
    joins = (current is not None
             and not starts_a_column(xmin, current[0])
             and (xmin - current[2] < FRAGMENT_GAP_MAXIMUM
                  or is_exponent(current[1], text, xmin - current[2])))
    if joins :
      current = (current[0], current[1] + text, xmax, current[3])
    else:
      if current is not None :
        value_list.append((current[0], current[3], current[1]))
      current = (xmin, text, xmax, y)
  if current is not None :
    value_list.append((current[0], current[3], current[1]))
  return value_list


def page_values(word_list):
  # 1 ページの値を (xmin, y, 文字列) で返す
  line = {}
  for word in word_list:
    line.setdefault(round(word[2]/2.0), []).append(word)
  value_list = []
  for key in sorted(line):
    value_list.extend(merge_fragments(sorted(line[key])))
  value_list.sort(key=lambda value: (value[1], value[0]))
  return value_list


def collect_records(page_list):
  #
  # レコード番号ごとに {(行, 列): [文字列]} を返す。
  #
  # レコードは時刻の列で見つける。表は TIME_FIRST から TIME_STEP 刻みなので、
  # 0 列目でその時刻に読める値がレコードの頭になる。時刻そのものが壊れた
  # レコードは**推測せずに落とす**ので、鍵は時刻ではなく番号にしてある。
  #
  # 他の値は、自分の**すぐ上にある時刻**のレコードに付ける。並び順で切ると
  # 1 行目の大半を落とす（走査が傾いていて、時刻と同じ行の値が時刻より
  # わずかに上に来ると、1 つ前のレコードに入ってしまう）。
  #
  record = {}
  for word_list in page_list:
    value_list = page_values(word_list)

    anchor_list = []
    for xmin, y, text in value_list:
      if column_of(xmin) != 0 :
        continue
      time = repair_time(text)
      if time is None :
        continue
      number = (time - TIME_FIRST)/TIME_STEP
      if abs(number - round(number)) > 1.0e-6 :
        continue
      number = int(round(number))
      if 0 <= number < RECORD_COUNT :
        anchor_list.append((y, number))
    anchor_list.sort()
    if not anchor_list :
      continue

    pending = {}
    for xmin, y, text in value_list:
      column = column_of(xmin)
      if column is None :
        continue
      chosen = None
      for y_anchor, number in anchor_list:
        if y_anchor > y + ROW_ABOVE_LIMIT :
          break
        chosen = (y_anchor, number)
      if chosen is None :
        continue
      offset = y - chosen[0]
      if offset > ROW_BELOW_LIMIT :
        continue
      pending.setdefault(chosen[1], []).append((offset, column, text))

    for number in pending:
      cells = record.setdefault(number, {})
      for row, column, text in assign_rows(pending[number]):
        cells.setdefault((row, column), []).append(text)
  return record


def assign_rows(entry_list):
  #
  # 1 レコードの値を 3 行に割り当てる。縦位置をまとまりに分け、**並び順**で
  # 行を決める（傾いたページでも、行の中の散らばりより行の間隔のほうが広い）。
  # まとまりが 2 つしか無いときだけ、離れ方を見て 1 行目と 3 行目の組かどうかを
  # 判断する。3 つより多ければ、近いもの同士から順に併合する
  #
  entry_list = sorted(entry_list)
  group      = []
  for entry in entry_list:
    if group and entry[0] - group[-1][-1][0] <= ROW_GAP_MINIMUM :
      group[-1].append(entry)
    else:
      group.append([entry])

  while len(group) > 3:
    distance = [group[index + 1][0][0] - group[index][-1][0]
                for index in range(len(group) - 1)]
    index = distance.index(min(distance))
    group[index] = group[index] + group[index + 1]
    del group[index + 1]

  if len(group) == 2 :
    centre = [sum(item[0] for item in one)/len(one) for one in group]
    row_of = [0, 2] if centre[1] - centre[0] > 1.5*ROW_SPACING else [0, 1]
  else:
    row_of = list(range(len(group)))

  result = []
  for index, one in enumerate(group):
    for offset, column, text in one:
      result.append((row_of[index], column, text))
  return result


class Damage(object):
  # 読めなかった値と、その扱い

  def __init__(self):
    self.entry_list = []
    self.count      = {}

  def add(self, time, name, kind, detail):
    self.entry_list.append((time, name, kind, detail))
    self.count[kind] = self.count.get(kind, 0) + 1

  def report(self, flag_verbose):
    print('  damaged or missing values: {:d}'.format(len(self.entry_list)))
    for kind in sorted(self.count):
      print('    {:<28s} {:d}'.format(kind, self.count[kind]))
    if flag_verbose :
      for time, name, kind, detail in self.entry_list:
        print('    t = {:8.2f}  {:<18s} {:<28s} {:s}'.format(time, name, kind, detail))


def first_readable(candidate_list):
  #
  # 1 つの升に複数の文字列が入ることがある（数の脇に走査のゴミが残った場合）。
  # 型に合う最初のものが数そのもの
  #
  if not candidate_list :
    return None, None
  for text in candidate_list:
    value = repair(text)
    if value is not None :
      return text, value
  return candidate_list[0], None


def read_field(cells, name, row, column, row_twin, column_twin, factor,
               time, damage, disputed):
  # 1 つの量を SI で返す。壊れていれば片割れの単位から復元する
  primary, value = first_readable(cells.get((row, column)))

  twin = None
  if row_twin is not None :
    dummy, converted = first_readable(cells.get((row_twin, column_twin)))
    if converted is not None :
      twin = converted*factor

  if value is not None and twin is not None :
    scale = max(abs(value), abs(twin))
    if scale > 0.0 and abs(value - twin)/scale > TOLERANCE_UNIT :
      damage.add(time, name, 'units disagree', '{:.8g} vs {:.8g}'.format(value, twin))
      disputed[name] = (value, twin)
      return float('nan')
    return value
  if value is not None :
    return value
  if twin is not None :
    damage.add(time, name, 'read from the other unit', repr(primary))
    return twin
  damage.add(time, name, 'unreadable',
             repr(primary) if primary is not None else 'absent')
  return float('nan')


def solve3(matrix, right):
  # 3 元 1 次連立を Gauss 消去で解く。解けなければ None
  row_list = [list(matrix[i]) + [right[i]] for i in range(3)]
  for step in range(3):
    pivot = max(range(step, 3), key=lambda i: abs(row_list[i][step]))
    if abs(row_list[pivot][step]) < 1.0e-30 :
      return None
    row_list[step], row_list[pivot] = row_list[pivot], row_list[step]
    for other in range(step + 1, 3):
      factor = row_list[other][step]/row_list[step][step]
      for column in range(step, 4):
        row_list[other][column] -= factor*row_list[step][column]
  answer = [0.0]*3
  for step in (2, 1, 0):
    total = row_list[step][3]
    for column in range(step + 1, 3):
      total -= row_list[step][column]*answer[column]
    answer[step] = total/row_list[step][step]
  return answer


def quadratic_prediction(window):
  # (ずれ, 値) に 2 次式を当て、ずれ 0 での値と残差を返す
  sums   = [0.0]*5
  moment = [0.0]*3
  for offset, value in window:
    power = 1.0
    for order in range(5):
      sums[order] += power
      power *= offset
    moment[0] += value
    moment[1] += value*offset
    moment[2] += value*offset*offset
  coefficient = solve3([[sums[0], sums[1], sums[2]],
                        [sums[1], sums[2], sums[3]],
                        [sums[2], sums[3], sums[4]]], moment)
  if coefficient is None :
    return window[0][1], float('inf')

  residual_list = []
  for offset, value in window:
    model = coefficient[0] + coefficient[1]*offset + coefficient[2]*offset*offset
    residual_list.append(abs(value - model))
  residual_list.sort()
  # 二乗平均ではなく**絶対残差の中央値**にしてある。走査は数レコード続けて
  # 壊すことがあり、二乗で測ると、暴きたい値そのものに残差を持ち上げられて
  # 何も引っかからなくなる。1.4826 は標準偏差と同じ尺度に載せるため
  return coefficient[0], residual_list[len(residual_list)//2]*1.4826


def transform(name, value):
  # 列が滑らかになる空間へ移す
  if name not in OUTLIER_LOGARITHMIC :
    return value
  if value <= 0.0 :
    return None
  return math.log(value)


def untransform(name, value):
  return math.exp(value) if name in OUTLIER_LOGARITHMIC else value


def neighbourhood(column, index, name):
  # 当てはめに使う前後の値
  window = []
  for offset in range(-OUTLIER_WINDOW, OUTLIER_WINDOW + 1):
    other = index + offset
    if offset == 0 or other < 0 or other >= len(column) :
      continue
    if math.isnan(column[other]) :
      continue
    mapped = transform(name, column[other])
    if mapped is not None :
      window.append((float(offset), mapped))
  return window


def screen_signs(table, damage):
  #
  # 負号が落ちた値を直す。走査は `-` を取りこぼすことがあり、経路角の
  # -35.5 度が +35.5 度になって出てくる（大きさは合っているので、大きさで測る
  # ふるいは素通りする）。**前後の値がすべて逆符号のときだけ**直す。これは
  # 値の推測ではなく、`O` を `0` と読むのと同じ、落ちた文字の復元である。
  # 符号の落ちた値が近くに 2 つあると互いに邪魔をするので、全員一致ではなく
  # SIGN_FRACTION の多数決にし、さらに 2 巡する
  #
  for dummy in range(2):
    screen_signs_once(table, damage)


def screen_signs_once(table, damage):
  for name in MAGNITUDE_SCREEN:
    column = [row[name] for row in table]
    for index in range(len(column)):
      value = column[index]
      if math.isnan(value) or value == 0.0 :
        continue
      neighbour = []
      for offset in range(-MAGNITUDE_WINDOW, MAGNITUDE_WINDOW + 1):
        other = index + offset
        if offset == 0 or other < 0 or other >= len(column) :
          continue
        if not math.isnan(column[other]) and column[other] != 0.0 :
          neighbour.append(column[other])
      if not neighbour :
        continue
      opposite = sum(1 for item in neighbour if item*value < 0.0)
      magnitude = sorted(abs(item) for item in neighbour)
      middle    = magnitude[len(magnitude)//2]
      # 逆符号が大多数**で、かつ大きさが近傍と揃っている**ときだけ直す。
      # 少数の値しか残っていない末端では全員一致を要求する
      if opposite < (len(neighbour) if len(neighbour) < 4
                     else SIGN_FRACTION*len(neighbour)) :
        continue
      if not SIGN_MAGNITUDE[0] < abs(value)/middle < SIGN_MAGNITUDE[1] :
        continue
      damage.add(table[index]['time'], name, 'sign restored',
                 '{:.8g} to {:.8g}'.format(value, -value))
      table[index][name] = -value
      column[index]      = -value


def screen_magnitudes(table, damage):
  # 近傍の中央値から桁ごと離れた値を落とす
  for name in MAGNITUDE_SCREEN:
    column = [row[name] for row in table]
    for index in range(len(column)):
      value = column[index]
      if math.isnan(value) or value == 0.0 :
        continue
      neighbour = []
      for offset in range(-MAGNITUDE_WINDOW, MAGNITUDE_WINDOW + 1):
        other = index + offset
        if offset == 0 or other < 0 or other >= len(column) :
          continue
        if not math.isnan(column[other]) and column[other] != 0.0 :
          neighbour.append(abs(column[other]))
      if len(neighbour) < 8 :
        continue
      neighbour.sort()
      middle = neighbour[len(neighbour)//2]
      ratio  = abs(value)/middle
      if ratio > MAGNITUDE_FACTOR or ratio < 1.0/MAGNITUDE_FACTOR :
        damage.add(table[index]['time'], name, 'decades from the neighbours',
                   '{:.8g} against {:.8g}'.format(value, middle))
        table[index][name] = float('nan')


def resolve_disputes(table, damage):
  #
  # 両単位が食い違った値を決める。どちらかの桁が間違っているので、
  # 外れ値の判定と同じく前後のレコードに決めさせる。**片方が明らかに遠いときだけ**
  # 採り、見分けが付かない組は推測せずに NaN のままにする
  #
  for name in FIELDS_OUTPUT:
    column = [row[name] for row in table]
    for index, row in enumerate(table):
      pair = row['disputed'].get(name)
      if pair is None :
        continue
      window = neighbourhood(column, index, name)
      if len(window) < 5 :
        continue
      predicted = quadratic_prediction(window)[0]
      scored    = []
      for candidate in pair:
        mapped = transform(name, candidate)
        if mapped is not None :
          scored.append((abs(mapped - predicted), candidate))
      if len(scored) < 2 :
        continue
      scored.sort()
      if scored[0][0]*DISPUTE_MARGIN > scored[1][0] :
        continue
      row[name] = scored[0][1]
      damage.add(row['time'], name, 'unit columns settled',
                 '{:.8g} over {:.8g}'.format(scored[0][1], scored[1][1]))


def flag_outliers(table, damage):
  # 近傍の滑らかさを壊す値を NaN にする。1 個消すと隣が見えるので繰り返す
  for dummy in range(OUTLIER_PASSES):
    if not flag_outliers_once(table, damage) :
      break


def flag_outliers_once(table, damage):
  flag_removed = False
  for name in FIELDS_OUTPUT:
    floor = OUTLIER_FLOOR.get(name)
    if floor is None :
      continue
    column = [row[name] for row in table]
    for index in range(len(column)):
      if math.isnan(column[index]) :
        continue
      window = neighbourhood(column, index, name)
      if len(window) < 5 :
        continue
      predicted, residual = quadratic_prediction(window)
      mapped = transform(name, column[index])
      error  = float('inf') if mapped is None else abs(mapped - predicted)
      if error > floor and error > OUTLIER_SIGMA*max(residual, 1.0e-12) :
        damage.add(table[index]['time'], name, 'breaks the neighbours',
                   '{:.8g}, neighbours give {:.8g}'.format(
                     column[index], untransform(name, predicted)))
        table[index][name] = float('nan')
        column[index]      = float('nan')
        flag_removed       = True
  return flag_removed


def check_identities(table):
  #
  # 表が満たすはずの関係。**これは診断であって選別ではない**。
  # p = rho R T は、ゾンデの大気を使った区間で R = 287.05 から 10 % ずれ、
  # 組成の変わる 95 km より上ではさらにずれる。走査ではなく表の性質なので、
  # 報告はするが、これで値を捨てることはしない
  #
  result  = {}
  dynamic, ideal, mach = [], [], []
  for row in table:
    density, velocity = row['density'], row['velocity']
    if not math.isnan(density) and not math.isnan(velocity) \
       and not math.isnan(row['dynamic_pressure']) :
      model = 0.5*density*velocity*velocity
      if model > 0.0 :
        dynamic.append(abs(row['dynamic_pressure'] - model)/model)
    if not math.isnan(density) and not math.isnan(row['temperature']) \
       and not math.isnan(row['pressure']) and row['pressure'] > 0.0 :
      ideal.append(abs(row['pressure'] - density*287.05*row['temperature'])/row['pressure'])
    if not math.isnan(row['temperature']) and not math.isnan(velocity) \
       and not math.isnan(row['mach']) and row['mach'] > 0.0 :
      model = velocity/math.sqrt(1.4*287.05*row['temperature'])
      mach.append(abs(row['mach'] - model)/row['mach'])
  for name, value_list in (('q = rho V^2 / 2', dynamic),
                           ('p = rho R T', ideal),
                           ('M = V / a', mach)):
    if value_list :
      value_list.sort()
      result[name] = (len(value_list), value_list[len(value_list)//2], value_list[-1])
  return result


def count_missing(table):
  total   = len(table)*len(FIELDS_OUTPUT)
  missing = sum(1 for row in table for name in FIELDS_OUTPUT if math.isnan(row[name]))
  return missing, total


def write_table(file_output, table):
  directory = os.path.dirname(file_output)
  if directory and not os.path.isdir(directory) :
    os.makedirs(directory)

  missing, total = count_missing(table)
  header = (
    '# Project Fire flight II reentry trajectory\n'
    '#\n'
    '# NASA TN D-3569 (Lewis and Scallion, 1966), table V:\n'
    '#   "REENTRY TRAJECTORY PARAMETERS OBTAINED FROM COMPUTER SIMULATION".\n'
    '#\n'
    "# This is NASA's own three-degree-of-freedom simulation of the reentry,\n"
    '# started from the Ascension Island TPQ-18 radar at t = 1608 s and flown\n'
    '# through an atmosphere measured by two Nike-Apache sounding rockets.\n'
    '# It is NOT a direct measurement. The report shows it merging with the\n'
    '# radar track on both sides of the blackout and reaching the measured\n'
    '# impact point within 500 m.\n'
    '#\n'
    '# Extracted from the scan by extract_trajectory.py.\n'
    '# {:d} of {:d} values could not be read and are written as nan.\n'
    '#\n'
    '# time             s, elapsed flight time (entry interface at {:.2f} s)\n'
    '# latitude         deg., geodetic, north positive\n'
    '# longitude        deg., east positive\n'
    '# altitude         m\n'
    '# velocity         m/s, relative to the rotating Earth\n'
    '# flightpath       deg. from the local geodetic horizon, negative downwards\n'
    '# heading          deg. clockwise from true north\n'
    '# dynamic_pressure N/m2\n'
    '# pressure         N/m2\n'
    '# density          kg/m3\n'
    '# temperature      K\n'
    '# mach             -\n'
    '# acceleration     g units, negative is deceleration\n'
    '# reynolds         -\n'.format(missing, total, TIME_FIRST))

  with open(file_output, 'w') as stream:
    stream.write(header)
    stream.write('#' + ' '.join(['{:>15s}'.format('time')]
                                + ['{:>15s}'.format(name) for name in FIELDS_OUTPUT]) + '\n')
    for row in table:
      stream.write(' ' + ' '.join(['{:15.6f}'.format(row['time'])]
                                  + ['{:15.8g}'.format(row[name]) for name in FIELDS_OUTPUT])
                   + '\n')


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--report', type=str, default=FILE_REPORT_DEFAULT,
                      help='Scanned NASA TN D-3569')
  parser.add_argument('--output', type=str, default=FILE_OUTPUT_DEFAULT,
                      help='Where to write the extracted table')
  parser.add_argument('--report-damage', action='store_true',
                      help='List every value that could not be read')
  args = parser.parse_args()

  print('Reading ' + os.path.normpath(args.report))
  record = collect_records(read_pages(run_pdftotext(args.report)))
  print('  {:d} of {:d} records found'.format(len(record), RECORD_COUNT))
  missing = [index for index in range(RECORD_COUNT) if index not in record]
  if missing :
    print('  records with no readable time column: '
          + ', '.join('{:.2f}'.format(TIME_FIRST + TIME_STEP*index) for index in missing))

  damage = Damage()
  table  = []
  for index in sorted(record):
    time = TIME_FIRST + TIME_STEP*index
    row  = {'time': time, 'disputed': {}}
    for name, row_table, column, row_twin, column_twin, factor in FIELDS:
      row[name] = read_field(record[index], name, row_table, column, row_twin,
                             column_twin, factor, time, damage, row['disputed'])
    table.append(row)

  screen_signs(table, damage)
  screen_magnitudes(table, damage)
  resolve_disputes(table, damage)
  flag_outliers(table, damage)
  damage.report(args.report_damage)

  print('  values kept, of {:d} records:'.format(len(table)))
  for name in FIELDS_OUTPUT:
    good = sum(1 for row in table if not math.isnan(row[name]))
    print('    {:<18s} {:5d}  ({:5.1f} %)'.format(name, good, 100.0*good/len(table)))

  print('  identities (reported, not enforced -- see check_identities):')
  for name, (count, median, worst) in sorted(check_identities(table).items()):
    print('    {:<18s} {:5d} records, median {:8.2e}, worst {:8.2e}'.format(
          name, count, median, worst))

  write_table(args.output, table)
  print('Wrote {:s} ({:d} records)'.format(os.path.relpath(args.output), len(table)))
  if not args.report_damage and damage.entry_list :
    print('--Run again with --report-damage to list the damaged values')

  return


if __name__ == '__main__':
  main()
