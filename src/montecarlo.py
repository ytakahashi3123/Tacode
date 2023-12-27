#!/usr/bin/env python3

import numpy as np
import os as os
import shutil as shutil
from orbital.orbital import orbital


class montecarlo(orbital):

  def __init__(self):

    print("Constructing class: montecarlo")

    return


  def initial_settings(self, config):

    self.work_dir      = config['montecarlo']['work_dir']
    self.case_dir      = config['montecarlo']['case_dir']

    path_specify = config['montecarlo']['template_path_specify']
    default_path = '/../../testcase_template' 
    manual_path  = config['montecarlo']['template_path']
    self.template_path = orbital.get_directory_path(path_specify, default_path, manual_path)

    # Make directory
    super().make_directory_rm(self.work_dir)
  
    # Copy template for tacode run
    #if ( os.path.exists(self.work_dir) ):
    #  os.rmdir(self.work_dir)
    self.work_dir_template = self.work_dir+'/'+self.case_dir+'_template'
    shutil.copytree(self.template_path, self.work_dir_template)

    self.filename_control_tacode    = config['montecarlo']['filename_control']
    self.filename_trajectory_tacode = config['montecarlo']['filename_trajectory']

    self.cmd_tacode = config['montecarlo']['cmd_shell']
    self.root_dir   = os.getcwd()
    self.cmd_home = os.path.dirname(os.path.realpath(__file__))

    # Tacodeコントロールファイルにおける置換関連
    self.txt_indentified = 'density_factor:' 
    self.ele_indentified = 1

    # Counter
    self.iter = 1

    # Result file
    self.result_dir       = config['montecarlo']['result_dir']
    self.flag_tecplot     = config['montecarlo']['flag_tecplot']
    self.header_tecplot   = config['montecarlo']['header_tecplot']
    self.filename_tecplot = config['montecarlo']['filename_tecplot']

    # Make directory
    #super().make_directory(self.result_dir)

    return


  def rewrite_control(self,filename,txt_indentified,ele_indentified,txt_replaced):
    # 
    # txt_indentifiedの文字列を含む行を抽出し、その(ele_indentified)列目要素を置換する。

    # Reading control file
    with open(filename) as f:
      lines = f.readlines()

    # リストとして取得 
    lines_strip = [line.strip() for line in lines]

    # 置換する行を特定する
    # i_line    = [i for i, line in enumerate(lines_strip) if txt_indentified in line]
    line_both        = [(i, line) for i, line in enumerate(lines_strip) if txt_indentified in line]
    i_line, str_line = list(zip(*line_both))
    # 抽出した行をスペース・タブで分割する。そのele_indentified列目を置換し、line_replacedというstr型に戻す。
    words = lines_strip[i_line[0]+1].split()
    #for n in range(0,len(ele_indentified)):
    #  words[n] = txt_replaced[n]
    words[self.ele_indentified] = txt_replaced
    line_replaced  = ' '.join(words)

    # lines_newはリストになることに注意。そのため、'',joinでstr型に戻す
    lines_new     = [item.replace( lines_strip[i_line[0]+1], line_replaced ) for item in lines]
    str_lines_new = ''.join(lines_new)

    # Update the file
    with open(filename, mode="w") as f:
      f.write(str_lines_new)
    
    return


  def run_tacode(self):
    import subprocess
    # Tacodeの実行
    
    # 計算ディレクトリに移動、実行、元ディレクトリに戻る
    os.chdir( self.work_dir_case )
    #subprocess.call('pwd')

    # Get relative path
    current_path=os.getcwd()
    relative_path = os.path.relpath(self.cmd_home, current_path)
    
    # Run Tacode
    try:
      subprocess.run([self.cmd_tacode, relative_path], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
      print(f"Error occurred. Error code: {e.returncode}")
      print("Error output:", e.stdout)

    os.chdir( self.root_dir )

    return


  def evaluate_error(self):

    # Trajectoryデータの読み込み
    filename_tmp = self.work_dir_case+'/'+self.filename_trajectory_tacode
    print('--Reading trajectory file (tacode)... :',filename_tmp)
    time_sec, longitude, latitude, altitude, velocity,  density, temperature, kn = super().read_inputdata_tacode(filename_tmp)
    velocity_long = velocity[0]
    velocity_lat  = velocity[1]
    velocity_alt  = velocity[2]
    velocity_mag  = velocity[3]

    return error_tmp


  def f_tacode(self):
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
    # --Atmospheric density factor
    txt_replaced = str(0.5)
    self.rewrite_control(filename_ctl,self.txt_indentified,self.ele_indentified,txt_replaced)

    # Tacodeの実行
    print('--Start Tacode')
    self.run_tacode()
    print('--End Tacode')

    # Trajectoryファイルの読み込みと誤差評価
    print('--Evaluating error between Tacode and GPR results')
   # error = self.evaluate_error()
    print('--Done, Error between Tacode and GPR: ')

    # カウンタの更新
    self.iter += 1

    return #error


def main():

  # 設定ファイルの読み込み
  file_control_default = orbital.file_control_default
  arg                  = orbital.argument(file_control_default)
  file_control         = arg.file
  config               = orbital.read_config_yaml(file_control)

  montecarlo.initial_settings(config)

  for n in range(0,10):
    montecarlo.f_tacode()

  return


if __name__ == '__main__':

  print('Initializing Tacode-MonteCarlo')

  # Call classes
  orbital = orbital()
  montecarlo = montecarlo()

  # Main
  main()

  print('Finalizing Tacode-MonteCarlo')

  exit()