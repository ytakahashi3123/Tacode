#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date; 2023/12/31

import numpy as np
import os as os
import shutil as shutil
from orbital.orbital import orbital

class montecarlo(orbital):

  def __init__(self):

    print("Constructing class: montecarlo")

    return


  def initialize(self, config):

    # GPR results
    #self.unit_covert_timeunit = unit_covert_timeunit
    #self.time_sec_opt  = time_opt
    #self.time_day_opt  = time_opt/unit_covert_timeunit
    #self.longitude_opt = geodetic_coord_opt[0]
    #self.latitude_opt  = geodetic_coord_opt[1]
    #self.altitude_opt  = geodetic_coord_opt[2]

    #self.i_target_opt         = super().getNearestIndex(self.time_sec_opt, config['tacode']['target_time']*self.unit_covert_timeunit)
    #self.target_time_opt      = self.time_day_opt[self.i_target_opt]
    #self.target_longitude_opt = self.longitude_opt[self.i_target_opt]
    #self.target_latitude_opt  = self.latitude_opt[self.i_target_opt]
    #self.target_altitude_opt  = self.altitude_opt[self.i_target_opt]

    #print('Time at start position:       ',self.target_time_opt, 'day')
    #print('Initial Geodetic coordinate, :',self.target_longitude_opt,'deg.,',self.target_latitude_opt,'deg.,',self.target_altitude_opt,'km')

    self.work_dir      = config['montecarlo']['work_dir']
    self.case_dir      = config['montecarlo']['case_dir']
    self.template_path = config['montecarlo']['template_path']

    # Make directory
    super().make_directory(self.work_dir)
  
    # Copy template for tacode run
    if ( os.path.exists(self.work_dir) ):
      os.rmdir(self.work_dir)
    self.work_dir_template = self.work_dir+'/'+self.case_dir+'_template'
    shutil.copytree(self.template_path, self.work_dir_template)

    self.filename_control_tacode    = config['montecarlo']['filename_control']
    self.filename_trajectory_tacode = config['montecarlo']['filename_trajectory']

    self.cmd_tacode = config['montecarlo']['cmd_tacode']
    self.root_dir   = os.getcwd()

    # Tacodeコントロールファイルにおける置換関連
    self.txt_coord_indentified = '# Initial coordinate' 
    self.ele_coord_indentified = [0,1,2]
    self.txt_veloc_indentified = '# Initial velocity' 
    self.ele_veloc_indentified = [0,1,2]

    # Counter
    self.iter = 1

    # Used in Tacode trajectory (:trajectory.dat) for comparison with GPR result
    #self.time_start      = config['tacode']['time_start']
    #self.time_end        = config['tacode']['time_end']
    #self.target_time_set = config['tacode']['target_time']

    # Result file
    self.result_dir       = config['montecarlo']['result_dir']
    self.flag_tecplot     = config['montecarlo']['flag_tecplot']
    self.header_tecplot   = config['montecarlo']['header_tecplot']
    self.filename_tecplot = config['montecarlo']['filename_tecplot']

    # Make directory
    super().make_directory(self.result_dir)

    return


  def rewrite_control(self,filename,txt_indentified,ele_indentified,txt_replaced):
    # 
    # txt_indentifiedの文字列を含む行を抽出し、その(ele_indentified+1)番目要素を置換する。

    # ファイル読み込み
    with open(filename) as f:
      lines = f.readlines()

    # リストとして取得 
    lines_strip = [line.strip() for line in lines]

    # 置換する行を特定する
    # i_line    = [i for i, line in enumerate(lines_strip) if txt_indentified in line]
    line_both        = [(i, line) for i, line in enumerate(lines_strip) if txt_indentified in line]
    i_line, str_line = list(zip(*line_both))
    # 抽出した行をスペース・タブで分割する。そのele_indentified+1番目を置換し、line_replacedというstr型に戻す。
    words = lines_strip[i_line[0]].split()
    for n in range(0,len(ele_indentified)):
      words[n] = txt_replaced[n]
    line_replaced  = ' '.join(words)

    # lines_newはリストになることに注意。そのため、'',joinでstr型に戻す
    lines_new     = [item.replace( lines_strip[i_line[0]], line_replaced ) for item in lines]
    str_lines_new = ''.join(lines_new)

    # 同じファイル名で保存
    with open(filename, mode="w") as f:
      f.write(str_lines_new)
    
    return


  def run_tacode(self):
    import subprocess
    # Tacodeの実行
    
    # 計算ディレクトリに移動、実行、元ディレクトリに戻る
    os.chdir( self.work_dir_case )
#    subprocess.call('pwd')
    subprocess.call( self.cmd_tacode )
    os.chdir( self.root_dir )

    return


  def evaluate_error(self):
    # Tacodeによるトラジェクトリ結果とGPR結果の誤差を評価する

    # Trajectoryデータの読み込み
    filename_tmp = self.work_dir_case+'/'+self.filename_trajectory_tacode
    print('--Reading trajectory file (tacode)... :',filename_tmp)
    time_sec, longitude, latitude, altitude, velocity,  density, temperature, kn = super().read_inputdata_tacode(filename_tmp)
    velocity_long = velocity[0]
    velocity_lat  = velocity[1]
    velocity_alt  = velocity[2]
    velocity_mag  = velocity[3]

    # 開始時刻をGPRデータと合わせる
    #time_day        = time_sec/orbital.unit_covert_day2sec
    #time_day_offset = time_day+self.target_time_set
    #time_sec_offset = time_sec+self.target_time_set*orbital.unit_covert_day2sec

    #i_start = super().getNearestIndex(time_day_offset, self.time_start)
    #i_end   = super().getNearestIndex(time_day_offset, self.time_end)

    # 誤差評価の計算
    #error_tmp = 0.0
    #count_tmp = 0
    #for n in range(i_start, i_end):
    #  count_tmp = count_tmp+1
    #  for m in range(0,len(self.time_day_opt)):
    #    if (self.time_day_opt[m] >= time_day_offset[n] ):
    #      m_opt = m
    #      break
    #  grad_fact = ( time_day_offset[n] - self.time_day_opt[m_opt-1] )/( self.time_day_opt[m_opt] - self.time_day_opt[m_opt-1] )
#   #   longitude_opt_cor = ( self.longitude_opt[m_opt] - self.longitude_opt[m_opt-1] )*grad_fact + self.longitude_opt[m_opt-1]
#   #   latitude_opt_cor  = ( self.latitude_opt[m_opt]  - self.latitude_opt[m_opt-1]  )*grad_fact + self.latitude_opt[m_opt-1]
    #  altitude_opt_cor  = ( self.altitude_opt[m_opt]  - self.altitude_opt[m_opt-1]  )*grad_fact + self.altitude_opt[m_opt-1]
    #  error_tmp = error_tmp + ( altitude[n] - altitude_opt_cor )**2 

#   #     error_tmp = error_tmp/float(count_tmp)
    #error_tmp = np.sqrt( error_tmp/float(count_tmp) )

    #if( self.flag_tecplot ):
    #  number_padded     = '{0:04d}'.format(self.iter)
    #  filename_tmp      = super().split_file(self.result_dir+'/'+self.filename_tecplot,'_case'+number_padded,'.')
    #  header_tmp        = self.header_tecplot
    #  print_message_tmp = '--Writing tecplot file... '
    #  delimiter_tmp     = '\t'
    #  comments_tmp      = ''
    #  output_tmp        = np.c_[time_day_offset, longitude, latitude, altitude, velocity_long, velocity_lat, velocity_alt, velocity_mag, density, temperature, kn, time_sec_offset]
    #  super().write_tecplotdata( filename_tmp, print_message_tmp, header_tmp, delimiter_tmp, comments_tmp, output_tmp )

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
    # --Initial Coordinate
    txt_replaced = [str(self.target_longitude_opt), str(self.target_latitude_opt), str(self.target_altitude_opt*orbital.unit_covert_km2m)]
    self.rewrite_control(filename_ctl,self.txt_coord_indentified,self.ele_coord_indentified,txt_replaced)
    # --Initial Velocity
    #u_lon, u_lat,u_alt  = x[:,0],x[:,1],x[:,2]
    #x_tmp        = [u_lon[0].tolist(),u_lat[0].tolist(),u_alt[0].tolist()]
    #txt_replaced = [str(n) for n in x_tmp]
    #print('--Initial velocity: ',x_tmp)
    self.rewrite_control(filename_ctl,self.txt_veloc_indentified,self.ele_veloc_indentified,txt_replaced)

    # Tacodeの実行
    print('--Start Tacode')
    self.run_tacode()
    print('--End Tacode')

    # Trajectoryファイルの読み込みと誤差評価
    print('--Evaluating error between Tacode and GPR results')
   # error = self.evaluate_error()
    print('--Done, Error between Tacode and GPR: ',error)

    # カウンタの更新
    self.iter += 1

    return #error


def main():

  # 設定ファイルの読み込み
  file_control_default = self.file_control_default
  arg                  = self.argument(file_control_default)
  file_control         = arg.file
  config               = self.read_config_yaml(file_control)

  self.initial_settings(config)

  for n in range(0,10):
    self.f_tacode()

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