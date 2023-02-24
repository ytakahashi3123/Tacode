#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/05/23

import numpy as np

# Dict key
KEY_LENGTH   = 'characteristic_length'
KEY_AREA     = 'characteristic_area'
KEY_MASS     = 'mass'
KEY_DRAG     = 'drag_coefficient'

KEY_KN      = 'Knudsen_number'
KEY_CD_MEAN = 'CD_mean'
KEY_ALT     = 'Altitude'


def initial_settings_satellite(config):

  aerodynamic_dict = read_aerodynamic_file(config)

  return aerodynamic_dict


def set_satellite_property(config):

  mass_satellite      = config['satellite']['mass']
  dragcoef_satellite  = config['satellite']['drag_coefficient']
  area_satellite      = config['satellite']['characteristic_aree']
  length_satellite    = config['satellite']['characteristic_length']

  return


def read_aerodynamic_file(config):

  filename_tmp = config['satellite']['directory_aerodynamic'] + '/' + config['satellite']['filename_aerodynamic']
  print('Reading aerodynamic model...:', filename_tmp)

  # File open
  delimiter = None
  comments = '#'
  skiprows = 3
  data_input = np.loadtxt(filename_tmp, delimiter=delimiter, comments=comments, skiprows=skiprows)
  kn_aero  = data_input[:,0]
  cfx_mean = data_input[:,1]
  cfy_mean = data_input[:,2]
  cfz_mean = data_input[:,3]
  cmx_mean = data_input[:,4]
  cmy_mean = data_input[:,5]
  cmz_mean = data_input[:,6]
  cfx_std = data_input[:,7]
  cfy_std = data_input[:,8]
  cfz_std = data_input[:,9]
  cmx_std = data_input[:,10]
  cmy_std = data_input[:,11]
  cmz_std = data_input[:,12]
  alt_aero = data_input[:,13]

  aerodynamic_dict = {KEY_KN:kn_aero, KEY_CD_MEAN:cfx_mean, KEY_ALT:alt_aero}

  return aerodynamic_dict


def get_aerodynamic_coefficient(knudsen, knudsen_aerodynamic, cdmean_aerodynamic):

  import scipy
  from scipy import interpolate

  # Interpolate aerodynamic data from knudsen number of satellite

  if knudsen < knudsen_aerodynamic[0] :
    cdmean = cdmean_aerodynamic[0] 
  elif knudsen > knudsen_aerodynamic[-1] :
    cdmean = cdmean_aerodynamic[-1] 
  else :
    f_cdmean = scipy.interpolate.interp1d(knudsen_aerodynamic, cdmean_aerodynamic ,kind="linear")
    cdmean   = f_cdmean(knudsen)

  return cdmean


