#!/usr/bin/env python3

import sys as sys
import numpy as np
import os as os
import glob as glob
import shutil as shutil
import random as random
from orbital.orbital import orbital
from general.general import get_setting


class montecarlo(orbital):

  def __init__(self):

    print("Constructing class: montecarlo")

    # 失敗したケースの記録
    # --親が失敗を検知できるようにするためのもの。終了コードを見ずに待つだけだと、
    #   ケースが落ちても親は次へ進み、統計は残ったケースから静かに作られる。
    # --__init__で持つのは、initial_settingsを通さずにpostprocessだけを呼ぶ
    #   使い方（テスト）でも属性が在ることを保証するため。
    self.case_failed        = []     # [ケースディレクトリ, 終了コード]
    self.case_no_result     = []     # 結果ファイルが無い（または空の）ケース
    self.flag_allow_failure = False

    # 分散を振る乱数。**シードを与えると同じ config が同じケース群を作る**
    # （報告書の図を作り直せるようにするため）。既定は None で、OS のエントロピーから
    # 種を取る従来どおりの挙動。initial_settings で config の値に置き換わる
    self.random_seed      = None
    self.random_generator = random.Random()

    return


  def initial_settings(self, config):

    self.work_dir      = config['montecarlo']['work_dir']
    self.case_dir      = config['montecarlo']['case_dir']

    path_specify = config['montecarlo']['template_path_specify']
    default_path = '/../../tutorial/template' 
    manual_path  = config['montecarlo']['template_path']
    self.template_path = self.get_directory_path(path_specify, default_path, manual_path)

    # Make directory
    super().make_directory_rm(self.work_dir)
  
    # Copy template for tacode run
    self.work_dir_template = self.work_dir+'/'+self.case_dir+'_template'
    shutil.copytree(self.template_path, self.work_dir_template)

    self.filename_control_tacode    = config['montecarlo']['filename_control']
    self.filename_trajectory_tacode = config['montecarlo']['filename_trajectory']

    self.cmd_tacode = config['montecarlo']['cmd_shell']
    self.root_dir   = os.getcwd()
    self.cmd_home = os.path.dirname(os.path.realpath(__file__)) + '/..'

    # 同時に走らせるケース数（'auto' ならこのマシンのコア数）
    self.maximum_number_execution = self.get_maximum_number_execution(config)

    # Counter
    self.iter = 1

    self.process_list = []

    # Result file
    # --全ケースの計算が終わったあと、postprocessが1つのTecplotファイルにまとめる。
    #   古いconfigにこれらのキーが無くても動くようにget_settingで読む
    section = config['montecarlo']
    self.flag_allow_failure = bool( get_setting(section, 'flag_allow_failure', False) )

    self.set_random_generator(section)
    self.result_dir       = get_setting(section, 'result_dir', 'result_tacode')
    self.flag_tecplot     = bool( get_setting(section, 'flag_tecplot', False) )
    self.filename_tecplot = get_setting(section, 'filename_tecplot', 'tecplot_montecarlo.dat')

    # Make directory
    if self.flag_tecplot :
      super().make_directory(self.result_dir)

    return


  def set_random_generator(self, section):
    #
    # 分散を振る乱数を用意する。
    #
    # **montecarlo.random_seed を与えると、同じ config が同じケース群を作る。**
    # 既定は None で、OS のエントロピーから種を取る従来どおりの挙動（毎回違う）。
    # 報告書に載せた図を作り直すには、そのときのシードが要る。
    #
    self.random_seed      = get_setting(section, 'random_seed', None)
    self.random_generator = random.Random(self.random_seed)
    if self.random_seed is not None :
      print('--Random seed:', self.random_seed)

    return


  def get_dispersed_value(self, value_default, dispersion):
    #
    # 基準値に分散を掛ける。**分散は乗算（相対）**なので、基準値 0 の成分は動かない
    # （加算の分散は当面不要というユーザー判断、2026-09-10）。
    #
    return [ value*(1.0 + (self.random_generator.random() - 0.50)*dispersion)
             for value in value_default ]


  def is_toplevel(self, line):
    #
    # YAMLのトップレベルのキー（インデントの無い"名前:"行）かどうか。
    # 空行とコメント行はセクションの内側にも現れるので、区切りとは見ない。
    #
    line_strip = line.strip()
    if line_strip == '' or line_strip.startswith('#') :
      return False

    return not line[0].isspace()


  def rewrite_control(self,filename,txt_indentified,ele_indentified,txt_replaced,txt_root=None):
    #
    # txt_rootのセクションからtxt_indentifiedのキーを探し、その直後に並ぶ
    # (ele_indentified)個のリスト要素（"- 値"の行）の値をtxt_replacedで置換する。
    #
    # セクションで絞るのは、同じキー名が複数のセクションに現れるため
    # （initial_settings.velocityとwind.velocity）。txt_rootを省略すると
    # ファイルの先頭から探す（従来の挙動）。
    #
    # 置換は行番号を指定して行う。値の同じ行が複数あるとき（風の[0.0, 0.0, 0.0]など）、
    # 文字列置換だと無関係な行まで一緒に書き換わってしまうため。

    # Reading control file
    with open(filename) as f:
      lines = f.readlines()

    # リストとして取得
    lines_strip = [line.strip() for line in lines]

    # 探索の範囲を決める（セクション名はトップレベルなのでインデントが無い）
    #
    # 終わりで閉じるのが要点である。閉じないと、指定したセクションにキーが無いとき
    # 探索がファイル末尾まで走り、別のセクションの同名キーを書き換えてしまう
    # （wind.velocityとinitial_settings.velocityのような組）。
    i_start = 0
    i_end   = len(lines)
    if txt_root is not None:
      i_root = [i for i, line in enumerate(lines) if line.startswith(txt_root+':')]
      if len(i_root) == 0:
        print('Section is not found in the control file:', txt_root, ',File:', filename)
        sys.exit(1)
      i_start = i_root[0]
      # 次のトップレベルのキー（インデントの無い行。空行とコメントは除く）でセクションが終わる
      i_next = [i for i in range(i_start+1, len(lines)) if self.is_toplevel(lines[i])]
      if len(i_next) > 0:
        i_end = i_next[0]

    # キーの行を特定する
    i_key = [i for i in range(i_start, i_end) if lines_strip[i].startswith(txt_indentified+':')]
    if len(i_key) == 0:
      print('Variable is not found in the control file:', txt_indentified, ',Section:', txt_root, ',File:', filename)
      sys.exit(1)
    i_key = i_key[0]

    # キーの直後に並ぶリスト要素を置換する
    for m in range(0,ele_indentified):
      i_line = i_key+m+1
      words  = lines_strip[i_line].split() if i_line < len(lines) else []
      # Replace (words[0]に該当する'-'は置換しない、その次のwords[1]を置換する)
      if len(words) < 2 or words[0] != '-':
        print('Variable is not a list of', ele_indentified, 'elements:', txt_indentified, ',Section:', txt_root, ',File:', filename)
        sys.exit(1)
      words[1] = txt_replaced[m]
      # インデントを考慮して新しい行を構築する
      indent       = lines[i_line][:len(lines[i_line])-len(lines[i_line].lstrip())]
      lines[i_line] = indent + ' '.join(words) + '\n'

    # Update the file
    with open(filename, mode="w") as f:
      f.writelines(lines)

    return


  def get_maximum_number_execution(self, config):
    #
    # 同時に走らせるケースの数。
    #
    # **'auto' ならこのマシンの論理コア数**にする。台数ぶんの数字を tutorial の
    # config に焼き込むと、書いた人の機械でしか合わないため。ケース数より多く
    # 走らせることはできないので、そこで頭を打つ。
    #
    # 走るのは互いに独立した子プロセスなので、上げると同時に使うメモリも比例して
    # 増える（1 ケースおよそ 90 MB）。
    #
    setting          = config['montecarlo']['maximum_number_execution']
    number_iteration = config['montecarlo']['number_iteration']

    if isinstance(setting, str) :
      if setting.strip().lower() != 'auto' :
        print('montecarlo.maximum_number_execution must be a positive integer or "auto":', setting)
        print('Program stopped.')
        sys.exit(1)
      number = os.cpu_count() or 1
      print('--Cases at the same time:', number, '(auto)')
    else :
      number = int(setting)

    if number < 1 :
      print('montecarlo.maximum_number_execution must be a positive integer or "auto":', setting)
      print('Program stopped.')
      sys.exit(1)

    return min(number, number_iteration)


  def run_tacode(self,config):
    import subprocess
    # Tacodeの実行

    num_iteration = config['montecarlo']['number_iteration']
    maximum_number_execution = self.maximum_number_execution
    
    # 計算ディレクトリに移動、実行、元ディレクトリに戻る
    # --戻すのはtry/finallyの中。途中で例外が出たまま戻らないと、以降のケースが
    #   前のケースのディレクトリの中に作られる
    os.chdir( self.work_dir_case )
    try:
      # Get relative path
      current_path  = os.getcwd()
      relative_path = os.path.relpath(self.cmd_home, current_path)

      # Run Tacode（子プロセスのカレントディレクトリは起動時に決まるので、
      # このあとで親が戻っても影響しない）
      process = subprocess.Popen([self.cmd_tacode, relative_path])
      self.process_list.append( [self.work_dir_case, process] )
    finally:
      os.chdir( self.root_dir )

    if self.iter%maximum_number_execution == 0 or self.iter == num_iteration:
      self.wait_tacode()

    return


  def wait_tacode(self):
    #
    # 走らせたケースの終了を待ち、終了コードを確認する。
    #
    # 待つだけで戻り値を捨てると、ケースが落ちたことが誰にも伝わらない。
    # ここで数えておいて、最後にreport_failureが報告する。
    #
    for case_dir, process in self.process_list:
      returncode = process.wait()
      if returncode != 0 :
        print('--Caution: the case failed, exit code', returncode, ':', case_dir)
        self.case_failed.append( [case_dir, returncode] )

    self.process_list = []

    return


  def add_case_no_result(self, case_dir):
    #
    # 結果ファイルが無い（または空の）ケースを記録する。二重に数えない。
    #
    if case_dir not in self.case_no_result :
      self.case_no_result.append(case_dir)

    return


  def check_case_result(self):
    #
    # 全ケースに結果ファイルが在ることを確かめる。
    #
    # 終了コードだけでは足りない: 子が0を返しても出力を書いていないことがあり、
    # flag_tecplotがFalseだとpostprocessが走らないのでそこの検査も通らない。
    #
    for case in self.get_case_directory():
      filename_tmp = case+'/'+self.filename_trajectory_tacode
      if not os.path.exists(filename_tmp) :
        self.add_case_no_result(case)

    return


  def report_failure(self):
    #
    # 失敗したケースを報告する。既定では1つでも失敗していれば停止する（exit 1）。
    #
    # montecarlo.flag_allow_failureをTrueにしたときだけ、件数を報告して続ける。
    # 「100ケースのうち3ケースが黙って落ちて、統計は97ケースから作られていた」
    # という結果を残さないためのもの。
    #
    # 同じケースが「終了コードが 0 でない」と「結果ファイルが無い」の両方に
    # 挙がるので、ケースの数として数え直す
    case_all = [case_dir for case_dir, returncode in self.case_failed]
    for case_dir in self.case_no_result:
      if case_dir not in case_all :
        case_all.append(case_dir)

    number_failed = len(case_all)
    if number_failed == 0 :
      return 0

    print('Cases that did not complete: ', number_failed)
    for case_dir, returncode in self.case_failed:
      print('--Exit code', returncode, ':', case_dir)
    for case_dir in self.case_no_result:
      print('--No result file:', case_dir+'/'+self.filename_trajectory_tacode)

    if self.flag_allow_failure :
      print('--montecarlo.flag_allow_failure is True, so the run goes on with the cases that did complete.')
      return number_failed

    print('--Set montecarlo.flag_allow_failure: True to go on with the cases that did complete.')
    print('Program stopped.')
    sys.exit(1)


  def get_case_directory(self):
    #
    # 走らせたケースディレクトリを名前順に返す（case_dir + 4桁の連番）。
    # テンプレート（case_dir + '_template'）は連番でないので拾われない。
    #
    pattern = self.work_dir+'/'+self.case_dir+'[0-9]'*4

    return sorted( glob.glob(pattern) )


  def read_tecplot_case(self, filename):
    #
    # ケースのTecplot出力を「変数行」と「データ行」に分ける。
    #
    # 列の中身には立ち入らない。風や姿勢の有無で列が増えるうえ、形式を
    # ここにもう1つ持つと出力側と二重管理になるため（変数行はそのまま引き継ぐ）。
    #
    with open(filename) as f:
      lines = f.readlines()

    variables = None
    data      = []
    for line in lines:
      line_strip = line.strip()
      if line_strip == '' or line_strip.startswith('#') :
        continue
      if line_strip.lower().startswith('variables') :
        variables = line_strip
        continue
      if line_strip.lower().startswith('zone') :
        continue
      data.append(line_strip)

    return variables, data


  def postprocess(self):
    #
    # 全ケースのTecplot出力を1つのファイルにまとめる（1ケース＝1ゾーン）。
    #
    # まとめるだけで、統計（着地点のばらつき）は取らない。そちらは
    # src_helper/montecarlo_dispersion/の仕事で、あちらはケースの制御ファイルとの
    # 差分からばらついた入力も拾う。ここはTecplotで全ケースを一度に見るためのもの。
    #
    if not self.flag_tecplot :
      return None

    case_list = self.get_case_directory()
    if len(case_list) == 0 :
      print('--No case directory to gather:', self.work_dir)
      return None

    # ケースは1つずつ読んで書き出す（全ケースを抱えると100ケースで数百MBになる）
    filename_out  = self.result_dir+'/'+self.filename_tecplot
    variables_all = None
    number_zone   = 0
    file = None
    for case in case_list:
      filename_tmp = case+'/'+self.filename_trajectory_tacode
      if not os.path.exists(filename_tmp) :
        print('--Caution: the case has no result file, skipped:', filename_tmp)
        self.add_case_no_result(case)
        continue

      variables, data = self.read_tecplot_case(filename_tmp)
      if len(data) == 0 :
        print('--Caution: the case has no data line, skipped:', filename_tmp)
        self.add_case_no_result(case)
        continue

      if file is None :
        # 最初に中身のあったケースの変数行を、そのまま全体の変数行にする
        variables_all = variables
        print('--Writing the gathered Tecplot file... :', filename_out)
        file = open(filename_out, mode='w')
        file.write('# Tecplot data: Tacode-MonteCarlo'+self.newline_code)
        if variables_all is not None :
          file.write(variables_all+self.newline_code)
      elif variables != variables_all :
        # 列がケースごとに違うと1つのファイルにまとめられない
        # （風や姿勢を一部のケースだけで入れた、など）
        file.close()
        print('The cases do not share the same variables:', filename_tmp)
        print('Program stopped.')
        sys.exit(1)

      # ゾーンの点数はケースごとに数え直す（着地する時刻が違うので行数が揃わない）
      file.write('zone t="'+os.path.basename(case)+'" i= '+str(len(data))+' f=point'+self.newline_code)
      for line in data:
        file.write(line+self.newline_code)
      number_zone += 1

    if file is None :
      print('--No result file to gather in:', self.work_dir)
      return None

    file.close()
    print('--Done, ', number_zone, 'zones')

    return filename_out


  def f_tacode(self,config):
    # Tacodeのコントロールファイルを適切に修正して、tacodeを実行する。

    print('Iteration: ', self.iter)

    # Caseディレクトリの作成
    number_padded      = '{0:04d}'.format(self.iter)
    self.work_dir_case = self.work_dir+'/'+self.case_dir+number_padded
    print('--Case directory: ', self.work_dir_case)
    shutil.copytree(self.work_dir_template, self.work_dir_case)

    # コントロールファイルの書き換え 
    print('--Modification: control file')
    filename_ctl = self.work_dir_case+'/'+self.filename_control_tacode
    var_montecarlo = config['montecarlo']['target_variable']
    for n in range(0, len(var_montecarlo)):
      var_name_ctl   = var_montecarlo[n][0]
      var_root_ctl   = var_montecarlo[n][1]
      var_dispersion = var_montecarlo[n][2]

      var_default = config[var_root_ctl][var_name_ctl]

      txt_indentified = var_name_ctl
      ele_indentified = len(var_default)
      txt_replaced = [ str(value) for value in self.get_dispersed_value(var_default, var_dispersion) ]

      print('Variable:',var_name_ctl,'in',var_root_ctl, ',Default:',var_default, ',With dispersion:',txt_replaced)
      self.rewrite_control(filename_ctl, txt_indentified, ele_indentified, txt_replaced, var_root_ctl)

    # Tacodeの実行
    print('--Start Tacode')
    self.run_tacode(config)
    print('--End Tacode')

    # Trajectoryファイルの読み込みと誤差評価
    #print('--Postprocess based on results')
    #error = self.evaluate_error()
    #print('--Done, Postprocess based on results')

    # Update counter
    self.iter += 1

    return #error
  

  def montecarlo_routine(self,config):

    for n in range(0,config['montecarlo']['number_iteration']):
      self.f_tacode(config)

    # 結果のとりまとめ（最後のケースまで待ってから。run_tacodeが
    # self.iter == number_iterationで待つので、ここでは全ケースが終わっている）
    print('Postprocess: gathering the results')
    self.postprocess()

    # 落ちたケースの報告（既定では停止する）
    self.check_case_result()

    return self.report_failure()
