#!/usr/bin/env python3

import sys as sys
import numpy as np
import os as os
import scipy.interpolate   # set_interpolator のコメントに残した
                          # RegularGridInterpolator を戻すときに要る
import attitude.attitude as attitude
from general.general import get_database_directory

# Dict key
KEY_LENGTH   = 'characteristic_length'
KEY_AREA     = 'characteristic_area'
KEY_MASS     = 'mass'
KEY_DRAG     = 'drag_coefficient'

KEY_KN      = 'Knudsen_number'
KEY_CD_MEAN = 'CD_mean'
KEY_ALT     = 'Altitude'
KEY_INTERP  = 'Interpolator'

# 6 自由度計算で使う迎角依存の係数表
KEY_AOA     = 'Angle_of_attack'
KEY_CF      = 'Force_coefficient'
KEY_CM      = 'Moment_coefficient'
KEY_COEF    = 'Coefficient'
KEY_FILE    = 'Filename'
KEY_LENGTH_REF = 'Length_reference'   # 表の Knudsen 数を作った代表長さ（書いていなければ None）

# 空力データベースの列数（Kn, CFx..CMz, SDV_CFx..SDV_CMz, Altitude）
NUM_COLUMN_AERODYNAMIC = 14

# 軸対称性の判定: 面内成分に対する比の閾値と、丸め誤差を拾わないための下限
TOLERANCE_AXISYMMETRY = 1.e-3
FLOOR_AXISYMMETRY     = 1.e-12

# 係数表の見出しに書ける「この表の Knudsen 数を作った代表長さ」の目印と、
# config の characteristic_length と食い違っているとみなす相対差
MARKER_LENGTH_REFERENCE   = 'Reference length'
TOLERANCE_LENGTH_REFERENCE = 1.e-6


def initial_settings_satellite(config):

  # kind_aerodynamic_model: constant は係数表を使わない（solver は config の
  # drag_coefficient を読む）。それでも表を読みに行くと、使いもしないファイルが
  # 無いだけで計算が止まる
  if config['satellite']['kind_aerodynamic_model'] == 'constant' :
    print('Aerodynamic model: constant. The coefficient table is not read.')
    return {}

  aerodynamic_dict = read_aerodynamic_file(config)

  aerodynamic_dict = set_interpolator(aerodynamic_dict)

  # 軸対称性の検査は、姿勢を解いていて、かつ係数を表から引くときだけ意味を持つ
  # （3 自由度では表から CFx しか使わず、constant の空力モデルでは表そのものを使わない）
  section = attitude.get_setting(config, 'attitude', None)
  if bool( attitude.get_setting(section, 'flag_attitude', False) ) and \
     config['satellite']['kind_aerodynamic_model'] == 'fileread' :
    warn_if_not_axisymmetric(aerodynamic_dict)

  return aerodynamic_dict


def check_axisymmetry(aerodynamic_dict):
  #
  # 係数表が軸対称の機体のものかどうかを見る。
  #
  # 6 自由度では表を全迎角 alpha_total だけで引き、attitude.matrix_aerodynamic_roll で
  # 係数ベクトルを実際の横流れ面へ回す。この扱いが厳密なのは軸対称の機体だけで、
  # そのとき表は全迎角・全 Knudsen 数で CFy = CMx = CMz = 0 になる（横流れ面の中では
  # 横力もロール・ヨーのモーメントも立たない）。表にこれらが残っていると、ロール行列は
  # それを面内の量として黙って誤った向きへ回す。表に迎角以外の姿勢変数が無い以上、
  # 正しい向きを復元する手段は無いので、検出して知らせるしかない。
  #
  # 判定は面内成分に対する比で行う（係数の絶対値は代表長さ・代表面積の取り方で変わる）。
  # CFy は max(|CFx|, |CFz|) に対して、CMx と CMz は max(|CMy|) に対して見る。
  # 代表スケールが 0 に潰れている表もあるので、下限 FLOOR_AXISYMMETRY を併せて使う
  # （解析モデルの表に残る 1e-17 級の丸め誤差を拾わないため）。
  #
  # 戻り値は逸脱した成分の一覧（軸対称なら空）。
  #
  angle_of_attack    = aerodynamic_dict[KEY_AOA]
  knudsen            = aerodynamic_dict[KEY_KN]
  coefficient_force  = aerodynamic_dict[KEY_CF]
  coefficient_moment = aerodynamic_dict[KEY_CM]

  scale_force  = np.abs( coefficient_force[:,:,[0,2]] ).max()
  scale_moment = np.abs( coefficient_moment[:,:,1] ).max()

  component_asymmetric = ( ('CFy', coefficient_force[:,:,1],  scale_force ),
                           ('CMx', coefficient_moment[:,:,0], scale_moment),
                           ('CMz', coefficient_moment[:,:,2], scale_moment) )

  violation = []
  for name, component, scale in component_asymmetric:
    magnitude = np.abs(component).max()
    if magnitude <= max( TOLERANCE_AXISYMMETRY*scale, FLOOR_AXISYMMETRY ) :
      continue
    index_aoa, index_knudsen = np.unravel_index( np.argmax(np.abs(component)), component.shape )
    violation.append( {'name'           : name,
                       'magnitude'      : magnitude,
                       'ratio'          : magnitude/scale if scale > 0.0 else float('inf'),
                       'angle_of_attack': angle_of_attack[index_aoa],
                       'knudsen'        : knudsen[index_knudsen]} )

  return violation


def warn_if_not_axisymmetric(aerodynamic_dict):
  # 軸対称でない係数表を 6 自由度で使おうとしていることを知らせる。
  # 停止はしない（表そのものは読めているし、迎角面内の力とモーメントは使えるため）。

  violation = check_axisymmetry(aerodynamic_dict)
  if len(violation) == 0 :
    return False

  print('Warning: the aerodynamic table is not that of an axisymmetric body.')
  print('--File:', aerodynamic_dict[KEY_FILE])
  for item in violation:
    print('--{} reaches {:.4e} ({:.2e} of the in-plane coefficients) at AOA {:g} deg., Kn {:.4e}'.format(
          item['name'], item['magnitude'], item['ratio'], item['angle_of_attack'], item['knudsen']))
  print('--The 6-DOF computation looks the table up with the total angle of attack alone and')
  print('--rotates the coefficients into the actual sideslip plane. That is exact only for an')
  print('--axisymmetric body, whose table has CFy = CMx = CMz = 0 everywhere. The components')
  print('--above are rotated as if they lay in that plane, so the side force and the rolling')
  print('--and yawing moments will point in the wrong direction.')
  print('--Use an axisymmetric table, or extend the table and the rotation to the sideslip')
  print('--angle as well.')

  return True


def set_interpolator(aerodynamic_dict):
  # 空力データベースの補間器を初期化時に 1 度だけ構築する
  # （大気側と同じく、毎ステップの再構築を避ける）
  #
  # 3 自由度の CD（Knudsen 数の 1 次元線形）はここに置かない。
  # np.interp が補間器を作らずに同じ値を返すので、構築するものが無い

  print('Setting aerodynamic interpolator...')

  interpolator = {}

  # 迎角依存の係数表があるときだけ (AOA, Kn) の格子を用意する。
  # 迎角が 1 つしかない表（従来の AOA 0 のみのファイル）では作らず、
  # 6 自由度計算では表の値を迎角によらず使う（復元モーメントは立たない）。
  #
  # ここに置くのは力とモーメントを 1 本にまとめた (AOA, Kn, 6) の配列で、
  # 引くのは get_aerodynamic_coefficient_attitude の evaluate_bilinear。
  # scipy の RegularGridInterpolator は使わない: 1 点を引くのに 22 us かかり、
  # 6 自由度計算の 1 割を占めていた（手書きの双一次は 5 us）。RGI は任意次元の
  # 超立方体を itertools で回す造りなので、2 次元 1 点では手間が計算を覆う。
  #
  # **値は RGI とビット単位で一致する。** evaluate_bilinear は RGI の
  # _evaluate_linear と同じ順序で足し込んである（実テーブル 3 種・節点を含む
  # 9936 点で確認済み。test_atmosphere.TestBilinearMatchesScipy が検査する）。
  #
  # 戻すときはここを次の 2 行に替え、evaluate_bilinear の呼び出しを
  # interpolator[KEY_COEF](np.array([[aoa, kn]]))[0] に戻せばよい:
  #
  #   interpolator[KEY_COEF] = scipy.interpolate.RegularGridInterpolator(
  #                              (angle_of_attack, aerodynamic_dict[KEY_KN]), coefficient,
  #                              method='linear', bounds_error=False, fill_value=None)
  #
  angle_of_attack = aerodynamic_dict[KEY_AOA]
  if len(angle_of_attack) > 1 :
    interpolator[KEY_COEF] = np.concatenate([aerodynamic_dict[KEY_CF], aerodynamic_dict[KEY_CM]], axis=2)

  aerodynamic_dict[KEY_INTERP] = interpolator

  return aerodynamic_dict


def read_aerodynamic_file(config):

  directory_path = get_database_directory(config['satellite'], 'satellite', 'aerodynamic', 'directory_aerodynamic')

  filename_tmp = directory_path + '/' + config['satellite']['filename_aerodynamic']
  print('Reading aerodynamic model...:', filename_tmp)

  angle_of_attack, data_block, length_reference = parse_aerodynamic_file(filename_tmp)

  # 各ブロックは同じ Knudsen 数の並びでなければ (AOA, Kn) の格子にならない
  kn_aero = data_block[0][:,0]
  for index in range(1, len(data_block)):
    if data_block[index].shape != data_block[0].shape or \
       not np.array_equal(data_block[index][:,0], kn_aero) :
      print('The aerodynamic table has inconsistent Knudsen numbers between the AOA blocks.')
      print('--Every AOA block must list the same Knudsen numbers in the same order.')
      print('--File:', filename_tmp)
      print('Program stopped.')
      sys.exit(1)
  if np.any(np.diff(kn_aero) <= 0.0) :
    print('The Knudsen numbers in the aerodynamic table are not in ascending order.')
    print('--File:', filename_tmp)
    print('Program stopped.')
    sys.exit(1)

  # 迎角ごと・Knudsen 数ごとの係数（力・モーメント）
  coefficient_force  = np.array([block[:,1:4] for block in data_block])
  coefficient_moment = np.array([block[:,4:7] for block in data_block])
  alt_aero = data_block[0][:,13]

  # 3 自由度計算で使う CD は迎角 0 に最も近いブロックの CFx とする
  # （従来の AOA 0 のみのファイルではそのブロックそのもの）。
  #
  # **同じ表を 2 通りに使っていることに注意**（CODE_REVIEW B-8）。3 自由度は
  # ここで取り出す 1 本の CD(Kn) だけを使い（get_aerodynamic_coefficient）、
  # 6 自由度は (迎角, Kn) の格子として表全体を引く
  # （get_aerodynamic_coefficient_attitude）。迎角 0 では両者はビット一致する。
  index_zero = int(np.argmin(np.abs(angle_of_attack)))
  cfx_mean   = coefficient_force[index_zero,:,0]

  print('--Angles of attack in the table (deg.):', ', '.join(['{:g}'.format(value) for value in angle_of_attack]))

  warn_if_length_differs(config, filename_tmp, length_reference)

  aerodynamic_dict = {KEY_KN:kn_aero, KEY_CD_MEAN:cfx_mean, KEY_ALT:alt_aero,
                      KEY_AOA:angle_of_attack, KEY_CF:coefficient_force, KEY_CM:coefficient_moment,
                      KEY_FILE:filename_tmp, KEY_LENGTH_REF:length_reference}

  return aerodynamic_dict


def parse_length_reference(line, filename_tmp):
  #
  # 見出しの "# Reference length: <値> m" を読む。その行でなければ None を返す。
  #
  # 目印はあるのに数値が読めない行は止める（書いたつもりの値が黙って無視され、
  # 検査も黙って行われなくなるため）。
  #
  text = line.lstrip('#').strip()
  if not text.lower().startswith( MARKER_LENGTH_REFERENCE.lower() ) :
    return None

  field = text.split(':', 1)
  word  = field[1].split() if len(field) > 1 else []
  try:
    return float(word[0])
  except (IndexError, ValueError):
    print('The "{}" line of the aerodynamic table has no number.'.format(MARKER_LENGTH_REFERENCE))
    print('--Line:', line.strip())
    print('--File:', filename_tmp)
    print('--Write it as "# {}: 0.8 m", or leave the line out.'.format(MARKER_LENGTH_REFERENCE))
    print('Program stopped.')
    sys.exit(1)


def warn_if_length_differs(config, filename_tmp, length_reference):
  #
  # 表の Knudsen 数を作った代表長さと、この計算の代表長さが食い違っていないか
  # （CODE_REVIEW B-7）。
  #
  # Kn = lambda/L なので、**表の Kn 軸は表を作ったときの L に紐づいている**。別の L で
  # 引くと、同じ高度に対して表の別の場所を読むことになり、CD が黙ってずれる。
  # 代表長さはモーメント係数の規格化にも使われるので、6 自由度ではモーメントも動く。
  #
  # **停止はしない。**手元にある表を別の機体に当てるのは、近似と割り切れば実際に行う
  # ことで（同梱のチュートリアルと Apollo の 3 自由度ケースがまさにそれをしている）、
  # 軸対称でない表を使うとき（warn_if_not_axisymmetric）と同じ扱いにしてある。
  #
  # 代表長さを書いていない表は検査しない（出所の分からない古いファイルが読めなくなる）。
  #
  if length_reference is None :
    return False

  length = float( config['satellite'][KEY_LENGTH] )
  if abs( length - length_reference ) <= TOLERANCE_LENGTH_REFERENCE*max( abs(length_reference), abs(length) ) :
    return False

  print('Caution: the Knudsen axis of the aerodynamic table was built for a different body.')
  print('--File :', filename_tmp)
  print('--Table: reference length {:g} m'.format(length_reference))
  print('--Case : satellite.characteristic_length = {:g} m'.format(length))
  print('--At a given altitude this run enters the table at {:.4g} times the Knudsen number'.format(
        length_reference/length))
  print('--the table associates with it, so the drag is read off the wrong part of the curve.')
  print('--The 6-DOF moments are scaled by the same length as well.')
  print('--Use the length the table was built with, or a table built for this body.')

  return True


def parse_aerodynamic_file(filename_tmp):
  #
  # 空力データベースを読む。迎角ごとのブロックに対応する。
  #
  #   <見出し行（数値でない行は読み飛ばす）>
  #   # Reference length: <値> m      <- 任意。あれば代表長さの食い違いを検査する
  #   AOA 0
  #   <Kn, CFx, CFy, CFz, CMx, CMy, CMz, SDV_CFx ... SDV_CMz, Altitude>
  #   ...
  #   AOA 10
  #   ...
  #
  # "AOA" 行が 1 つも無いファイル（迎角 0 のみの従来形式）もそのまま読める。
  #
  with open(filename_tmp) as f:
    lines = f.readlines()

  angle_of_attack  = []
  data_block       = []
  data_current     = None
  length_reference = None

  for line in lines:
    words = line.split()
    if len(words) == 0 :
      continue
    if words[0].startswith('#') :
      value = parse_length_reference(line, filename_tmp)
      if value is not None :
        length_reference = value
      continue

    # 迎角の見出し行
    if words[0].upper() == 'AOA' :
      if len(words) < 2 :
        print('An "AOA" line in the aerodynamic table has no value.')
        print('--File:', filename_tmp)
        print('Program stopped.')
        sys.exit(1)
      angle_of_attack.append( float(words[1]) )
      data_current = []
      data_block.append( data_current )
      continue

    # データ行は「列数がちょうど揃った数値の行」だけを採る。
    # それ以外（見出しなど）は無害に読み飛ばす。
    if len(words) != NUM_COLUMN_AERODYNAMIC :
      continue
    try:
      values = [float(word) for word in words]
    except ValueError:
      continue

    if data_current is None :
      # "AOA" 行が無いまま数値が現れた場合は迎角 0 のブロックとみなす
      angle_of_attack.append( 0.0 )
      data_current = []
      data_block.append( data_current )
    data_current.append( values )

  if len(data_block) == 0 or len(data_block[0]) == 0 :
    print('No data row was found in the aerodynamic table.')
    print('--File:', filename_tmp)
    print('Program stopped.')
    sys.exit(1)

  # 迎角の昇順に並べ替える
  angle_of_attack = np.array(angle_of_attack)
  order = np.argsort(angle_of_attack, kind='stable')
  angle_of_attack = angle_of_attack[order]
  data_block = [np.array(data_block[index]) for index in order]

  if len(angle_of_attack) > 1 and np.any(np.diff(angle_of_attack) <= 0.0) :
    print('The aerodynamic table has duplicated angles of attack.')
    print('--File:', filename_tmp)
    print('Program stopped.')
    sys.exit(1)

  return angle_of_attack, data_block, length_reference


def get_aerodynamic_coefficient(knudsen, knudsen_aerodynamic, cdmean_aerodynamic):

  # Knudsen 数から 3 自由度の CD（迎角 0 ブロックの CFx）を線形補間する。
  # 表の外側は端の値でクランプする（np.interp の既定の挙動そのもの。
  # 6 自由度側の get_aerodynamic_coefficient_attitude と同じ扱い）。
  #
  # かつては scipy.interpolate.interp1d(kind='linear') を初期化時に作っていたが、
  # interp1d は SciPy 1.10 以降 legacy で、np.interp が同じ値をビット単位で返す。

  return np.interp(knudsen, knudsen_aerodynamic, cdmean_aerodynamic)




def get_aerodynamic_coefficient_attitude(knudsen, angle_attack_total, aerodynamic_dict):
  #
  # 全迎角（deg）と Knudsen 数から、係数表の面（機体軸 x-z 面）における
  # 力・モーメントの係数ベクトルを得る。
  #
  #   coefficient_force  = [CFx, CFy, CFz]
  #   coefficient_moment = [CMx, CMy, CMz]
  #
  # 表の外側は端の値でクランプする（1 次元の CD と同じ扱い）。
  #
  # **3 自由度とは同じ表の使い方が違う**（CODE_REVIEW B-8）。3 自由度が使うのは
  # read_aerodynamic_file が取り出した「迎角 0 ブロックの CFx」1 本だけで、
  # ここは (迎角, Kn) の格子として表全体を引く。迎角 0 では両者はビット一致する。
  #
  knudsen_table = aerodynamic_dict[KEY_KN]
  aoa_table     = aerodynamic_dict[KEY_AOA]

  knudsen_clamped = min( max(knudsen, knudsen_table[0]), knudsen_table[-1] )

  if len(aoa_table) == 1 :
    # 迎角 1 点だけの表。迎角によらず同じ値を使う
    coefficient_force  = np.array([np.interp(knudsen_clamped, knudsen_table, aerodynamic_dict[KEY_CF][0,:,m]) for m in range(0,3)])
    coefficient_moment = np.array([np.interp(knudsen_clamped, knudsen_table, aerodynamic_dict[KEY_CM][0,:,m]) for m in range(0,3)])
    return coefficient_force, coefficient_moment

  aoa_clamped = min( max(angle_attack_total, aoa_table[0]), aoa_table[-1] )

  coefficient = evaluate_bilinear(aoa_table, knudsen_table,
                                  aerodynamic_dict[KEY_INTERP][KEY_COEF],
                                  aoa_clamped, knudsen_clamped)

  return coefficient[0:3], coefficient[3:6]


def evaluate_bilinear(grid_first, grid_second, values, point_first, point_second):
  #
  # 2 次元の線形補間。格子は昇順、点は格子の内側にあること（呼び出し側でクランプ済み）。
  #
  # scipy.interpolate.RegularGridInterpolator の代わりに手で書いてある。理由と
  # 戻し方は set_interpolator のコメントを見ること。**足し込む順序を RGI の
  # _evaluate_linear に合わせてあるので、値はビット単位で一致する**:
  #
  #   value = 0 + v[i,j]*((1-y0)*(1-y1)) + v[i,j+1]*((1-y0)*y1)
  #             + v[i+1,j]*(y0*(1-y1))   + v[i+1,j+1]*(y0*y1)
  #
  # 区間の決め方（searchsorted の既定 side='left' から 1 を引き、両端で丸める）も
  # RGI に合わせてある。格子の内側では side の取り方で値は変わらないが
  # （節点では一方が重み 1、他方が重み 0 になって同じ値に行き着く）、
  # 読み比べられるように揃えておく。
  #
  index_first  = min( max(int(np.searchsorted(grid_first,  point_first )) - 1, 0), len(grid_first)  - 2 )
  index_second = min( max(int(np.searchsorted(grid_second, point_second)) - 1, 0), len(grid_second) - 2 )

  distance_first  = (point_first  - grid_first[index_first]  )/(grid_first[index_first+1]   - grid_first[index_first]  )
  distance_second = (point_second - grid_second[index_second])/(grid_second[index_second+1] - grid_second[index_second])

  value = 0.0
  for index_i, weight_i in ((index_first, 1.0 - distance_first), (index_first + 1, distance_first)):
    for index_j, weight_j in ((index_second, 1.0 - distance_second), (index_second + 1, distance_second)):
      value = value + values[index_i, index_j]*(weight_i*weight_j)

  return value
