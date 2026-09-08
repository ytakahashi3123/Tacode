#!/usr/bin/env python3

import numpy as np
import os as os
import scipy.interpolate
import attitude.attitude as attitude

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

# 空力データベースの列数（Kn, CFx..CMz, SDV_CFx..SDV_CMz, Altitude）
NUM_COLUMN_AERODYNAMIC = 14

# 軸対称性の判定: 面内成分に対する比の閾値と、丸め誤差を拾わないための下限
TOLERANCE_AXISYMMETRY = 1.e-3
FLOOR_AXISYMMETRY     = 1.e-12


def initial_settings_satellite(config):

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

  print('Setting aerodynamic interpolator...')

  interpolator = {}
  interpolator[KEY_CD_MEAN] = scipy.interpolate.interp1d(aerodynamic_dict[KEY_KN], aerodynamic_dict[KEY_CD_MEAN], kind="linear")

  # 迎角依存の係数表があるときだけ (AOA, Kn) の 2 次元補間器を作る。
  # 迎角が 1 つしかない表（従来の AOA 0 のみのファイル）では作らず、
  # 6 自由度計算では表の値を迎角によらず使う（復元モーメントは立たない）。
  angle_of_attack = aerodynamic_dict[KEY_AOA]
  if len(angle_of_attack) > 1 :
    coefficient = np.concatenate([aerodynamic_dict[KEY_CF], aerodynamic_dict[KEY_CM]], axis=2)
    interpolator[KEY_COEF] = scipy.interpolate.RegularGridInterpolator(
                               (angle_of_attack, aerodynamic_dict[KEY_KN]), coefficient,
                               method='linear', bounds_error=False, fill_value=None)

  aerodynamic_dict[KEY_INTERP] = interpolator

  return aerodynamic_dict


def set_satellite_property(config):

  mass_satellite      = config['satellite']['mass']
  dragcoef_satellite  = config['satellite']['drag_coefficient']
  area_satellite      = config['satellite']['characteristic_aree']
  length_satellite    = config['satellite']['characteristic_length']

  return


def read_aerodynamic_file(config):

  script_directory = os.path.dirname(os.path.realpath(__file__))
  if config['satellite']['directory_path_specify'] == 'auto' or config['satellite']['directory_path_specify'] == 'default':
    directory_path = script_directory + '/../../database/aerodynamic' 
  elif config['satellite']['directory_path_specify'] == 'manual':
    directory_path = config['satellite']['directory_aerodynamic']
  else :
    directory_path = script_directory + '/../../database/aerodynamice' 

  filename_tmp = directory_path + '/' + config['satellite']['filename_aerodynamic']
  print('Reading aerodynamic model...:', filename_tmp)

  angle_of_attack, data_block = parse_aerodynamic_file(filename_tmp)

  # 各ブロックは同じ Knudsen 数の並びでなければ (AOA, Kn) の格子にならない
  kn_aero = data_block[0][:,0]
  for index in range(1, len(data_block)):
    if data_block[index].shape != data_block[0].shape or \
       not np.array_equal(data_block[index][:,0], kn_aero) :
      print('The aerodynamic table has inconsistent Knudsen numbers between the AOA blocks.')
      print('--Every AOA block must list the same Knudsen numbers in the same order.')
      print('--File:', filename_tmp)
      print('Program stopped.')
      exit()
  if np.any(np.diff(kn_aero) <= 0.0) :
    print('The Knudsen numbers in the aerodynamic table are not in ascending order.')
    print('--File:', filename_tmp)
    print('Program stopped.')
    exit()

  # 迎角ごと・Knudsen 数ごとの係数（力・モーメント）
  coefficient_force  = np.array([block[:,1:4] for block in data_block])
  coefficient_moment = np.array([block[:,4:7] for block in data_block])
  alt_aero = data_block[0][:,13]

  # 3 自由度計算で使う CD は迎角 0 に最も近いブロックの CFx とする
  # （従来の AOA 0 のみのファイルではそのブロックそのもの）
  index_zero = int(np.argmin(np.abs(angle_of_attack)))
  cfx_mean   = coefficient_force[index_zero,:,0]

  print('--Angles of attack in the table (deg.):', ', '.join(['{:g}'.format(value) for value in angle_of_attack]))

  aerodynamic_dict = {KEY_KN:kn_aero, KEY_CD_MEAN:cfx_mean, KEY_ALT:alt_aero,
                      KEY_AOA:angle_of_attack, KEY_CF:coefficient_force, KEY_CM:coefficient_moment,
                      KEY_FILE:filename_tmp}

  return aerodynamic_dict


def parse_aerodynamic_file(filename_tmp):
  #
  # 空力データベースを読む。迎角ごとのブロックに対応する。
  #
  #   <見出し行（数値でない行は読み飛ばす）>
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

  angle_of_attack = []
  data_block      = []
  data_current    = None

  for line in lines:
    words = line.split()
    if len(words) == 0 :
      continue
    if words[0].startswith('#') :
      continue

    # 迎角の見出し行
    if words[0].upper() == 'AOA' :
      if len(words) < 2 :
        print('An "AOA" line in the aerodynamic table has no value.')
        print('--File:', filename_tmp)
        print('Program stopped.')
        exit()
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
    exit()

  # 迎角の昇順に並べ替える
  angle_of_attack = np.array(angle_of_attack)
  order = np.argsort(angle_of_attack, kind='stable')
  angle_of_attack = angle_of_attack[order]
  data_block = [np.array(data_block[index]) for index in order]

  if len(angle_of_attack) > 1 and np.any(np.diff(angle_of_attack) <= 0.0) :
    print('The aerodynamic table has duplicated angles of attack.')
    print('--File:', filename_tmp)
    print('Program stopped.')
    exit()

  return angle_of_attack, data_block


def get_aerodynamic_coefficient(knudsen, knudsen_aerodynamic, cdmean_aerodynamic, interpolator_aero):

  # Interpolate aerodynamic data from knudsen number of satellite
  # 補間器は initial_settings_satellite で構築済みのものを評価するだけ

  if knudsen < knudsen_aerodynamic[0] :
    cdmean = cdmean_aerodynamic[0]
  elif knudsen > knudsen_aerodynamic[-1] :
    cdmean = cdmean_aerodynamic[-1]
  else :
    cdmean = interpolator_aero[KEY_CD_MEAN](knudsen)

  return cdmean




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
  knudsen_table = aerodynamic_dict[KEY_KN]
  aoa_table     = aerodynamic_dict[KEY_AOA]

  knudsen_clamped = min( max(knudsen, knudsen_table[0]), knudsen_table[-1] )

  if len(aoa_table) == 1 :
    # 迎角 1 点だけの表。迎角によらず同じ値を使う
    coefficient_force  = np.array([np.interp(knudsen_clamped, knudsen_table, aerodynamic_dict[KEY_CF][0,:,m]) for m in range(0,3)])
    coefficient_moment = np.array([np.interp(knudsen_clamped, knudsen_table, aerodynamic_dict[KEY_CM][0,:,m]) for m in range(0,3)])
    return coefficient_force, coefficient_moment

  aoa_clamped = min( max(angle_attack_total, aoa_table[0]), aoa_table[-1] )

  coefficient = aerodynamic_dict[KEY_INTERP][KEY_COEF]( np.array([[aoa_clamped, knudsen_clamped]]) )[0]

  return coefficient[0:3], coefficient[3:6]
