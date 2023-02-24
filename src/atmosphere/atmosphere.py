#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/05/31

import numpy as np


# 単位変換
# ’read_atmosphere_file’でKey errorがあったときは下記の単位変換が正しいかチェックする
unit_convert ={'km':1.0, 'cm-3':1.e6, 'g/cm-3': 1.e3, 'K': 1.0}

# Dict key
KEY_DATA   = 'Number_Data'
KEY_ATM    = 'Number_Atmosphere'
KEY_KN     = 'Knudsen_number'

KEY_Height = 'Height'
KEY_N2     = 'N2'
KEY_O2     = 'O2'
KEY_N      = 'N'
KEY_O      = 'O'
KEY_Mass_density = 'Mass_density'
KEY_Temperature_neutral = 'Temperature_neutral'

LIST_MOLECULAR_KIND = [KEY_N2, KEY_O2, KEY_N, KEY_O]
DICT_DIAMETER_MOLECULAR = {KEY_N2: 3.75e-10, KEY_O2: 3.54e-10, KEY_N: 3.10e-10, KEY_O: 3.04e-10}


def initial_settings_atmosphere(config):

  atmosphere_dict = read_atmosphere_file(config)

  atmosphere_dict = set_knudsen_number(config, atmosphere_dict)

  return atmosphere_dict 


def read_atmosphere_file(config):

  filename_tmp = config['atmosphere']['directory_atmosphere'] + '/' + config['atmosphere']['filename_atmosphere']
  print('Reading atmosphere model...:', filename_tmp)

  # File open
  with open(filename_tmp) as f:
    lines = f.readlines()
  # リストとして取得 
  lines_strip = [line.strip() for line in lines]

  # 大気モデルパラメーターの取得
  atmosphere_name = []
  atmosphere_unit = []
  for i in range(0, len(lines_strip) ):
    if lines_strip[i] == 'Selected parameters are:' :
      i_count = i+1
      break
    if i == len(lines_strip)-1 :
      print('There is no specific word: "Selected parameters are:", Check atmosphere model file.')
      print('Program stopped.')
      exit()

  for i in range(i_count, len(lines_strip) ):
    # --ブランク行が来たらループから出る
    if lines_strip[i] == '':
      break
    # --空白で分割
    words = lines_strip[i].split()
    # --コンマ削除
    words = [n.replace(",","") for n in words]
    atmosphere_name.append( words[1] )
    atmosphere_unit.append( words[2] )
  

  num_atmosphere_var = len(atmosphere_name)
  
  print('Atmosphere name, unit')
  for i in range(0,num_atmosphere_var):
    print(i,atmosphere_name[i],',', atmosphere_unit[i])
  
  # Check unit convert  
  for i in range(0,num_atmosphere_var):
    try:
      unit_convert[atmosphere_unit[i]]
    except KeyError as instance:
      print(instance)
      print('Key error for unit convert. Check module: atmosphere.py, or atmosphere file')
      print('Program stopped.')
      exit()


  # 大気モデルデータの取得
  i_count = i_count + len(atmosphere_name) + 2
  num_array_tmp    = len(lines_strip) - i_count
  atmosphere_model = np.zeros(num_atmosphere_var*num_array_tmp).reshape(num_atmosphere_var,num_array_tmp) 
  for i in range(i_count, len(lines_strip) ):
    words = lines_strip[i].split()
    for j in range(0,num_atmosphere_var):
      atmosphere_model[j,i-i_count] = float( words[j] )*unit_convert[atmosphere_unit[j]]


  atmosphere_dict = {KEY_DATA:num_array_tmp, KEY_ATM:num_atmosphere_var}
  for i in range(0,num_atmosphere_var):
    atmosphere_dict[ atmosphere_name[i] ] = atmosphere_model[i]

  f.close()

  #print( atmosphere_dict['Temperature_neutral'])

  return atmosphere_dict


def set_knudsen_number(config, atmosphere_dict):

  print('Setting Knudsen number...')

  length   = config['satellite']['characteristic_length'] 
  num_data = atmosphere_dict[KEY_DATA]
  num_atm  = atmosphere_dict[KEY_ATM]

  d2_nd_total = np.zeros(num_data).reshape(num_data)
  for m in range(0, len(LIST_MOLECULAR_KIND) ):
    key_tmp      = LIST_MOLECULAR_KIND[m]
    try:
      numb_density = atmosphere_dict[key_tmp]
      diamter      = DICT_DIAMETER_MOLECULAR[key_tmp]
    except KeyError as instance:
      continue
    # Calculate sum( nd*diamter^2 )
    for n in range(0,num_data):
      d2_nd_total[n] = d2_nd_total[n] + numb_density[n]*diamter**2
  
  knudsen_number = 1.0/( np.sqrt(2.0)*np.pi*d2_nd_total*length)
  atmosphere_dict[KEY_KN] = knudsen_number

  return atmosphere_dict


def get_atmosphere_property(altitude, altitude_atm, density_atm, temperature_atm, knudsen_atm):
  
  import scipy
  from scipy import interpolate

  # Interpolate atmosphere data from altitude data of satellite

  if altitude < altitude_atm[0] :
    density     = density_atm[0] 
    temperature = temperature_atm[0]
    knudsen     = knudsen_atm[0]
  elif altitude > altitude_atm[-1] :
    density     = density_atm[-1] 
    temperature = temperature_atm[-1]
    knudsen     = knudsen_atm[-1]
  else :
    #f_linear = scipy.interpolate.interp1d(altitude_atm, density_atm)
    #f_quad   = scipy.interpolate.interp1d(altitude_atm, density_atm ,kind="quadratic")
    f_density     = scipy.interpolate.interp1d(altitude_atm, density_atm, kind="cubic")
    f_temperature = scipy.interpolate.interp1d(altitude_atm, temperature_atm, kind="cubic")
    f_knudsen     = scipy.interpolate.interp1d(altitude_atm, knudsen_atm, kind="cubic")
    density       = f_density(altitude)
    temperature   = f_temperature(altitude)
    knudsen       = f_knudsen(altitude)

  #x_new = np.linspace(0.5, 300, num=500)
  #print(x_new)
  #print( f_cubic(x_new) )

  return density, temperature, knudsen


