#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/05/23

import numpy as np

tiny_value = 1.e-50

def convert_cartesian_geodetic(config, cartesian_coord):
  #
  # Cartesian Coordinate --> Longitude, Latitude, Altitude coordinate (WGS84)
  #
  # carcoord_x,y,z: meter
  # Longitude: Radius
  # Latitude: Radius
  # Altitude: meter

  radius_equat_planet = config['planet']['radius']
  ellipticity_planet  = config['planet']['ellipticity']

  carcoord_x = cartesian_coord[0]
  carcoord_y = cartesian_coord[1]
  carcoord_z = cartesian_coord[2]

  radius_proj = np.sqrt( carcoord_x**2 + carcoord_y**2 )
  B_geo_tmp   = np.sign( carcoord_z ) * radius_equat_planet * (1.0 - ellipticity_planet)
  E_geo_tmp   = ( (carcoord_z + B_geo_tmp)*B_geo_tmp/radius_equat_planet - radius_equat_planet )/radius_proj
  F_geo_tmp   = ( (carcoord_z - B_geo_tmp)*B_geo_tmp/radius_equat_planet + radius_equat_planet )/radius_proj

  P_geo_tmp   = 4.0*(E_geo_tmp*F_geo_tmp + 1.0)/3.0
  Q_geo_tmp   = 2.0*(E_geo_tmp**2 - F_geo_tmp**2)
  D_geo_tmp   = P_geo_tmp**3 + Q_geo_tmp**2

  if( D_geo_tmp >= 0.0 ):
  #if( 1 in np.sign(D_geo_tmp) ):
    S_geo_tmp = np.sqrt(D_geo_tmp) + Q_geo_tmp
    S_geo_tmp = np.sign(S_geo_tmp) * np.exp( np.log( np.abs(S_geo_tmp+tiny_value))/3.0 )
    V_geo_tmp = - ( 2.0*Q_geo_tmp + (P_geo_tmp/(S_geo_tmp+tiny_value) - S_geo_tmp)**3 )/(3.0*P_geo_tmp+tiny_value)
  else:
    V_geo_tmp = 2.0 * np.sqrt( -P_geo_tmp ) * np.cos( np.arccos(Q_geo_tmp/np.sqrt( - P_geo_tmp**3 )) / 3.0)

  G_geo_tmp = 0.5*( E_geo_tmp + np.sqrt(E_geo_tmp**2 + V_geo_tmp) )
  T_geo_tmp = np.sqrt( G_geo_tmp**2 + (F_geo_tmp - V_geo_tmp*G_geo_tmp)/(2.0*G_geo_tmp - E_geo_tmp )) - G_geo_tmp

# Longitude
  longitude = np.sign(carcoord_y) * np.arccos( carcoord_x / np.sqrt( carcoord_x**2 + carcoord_y**2 ) )
# Latitude
  latitude  = np.arctan( (1.0 - T_geo_tmp**2)*radius_equat_planet / (2.0 * B_geo_tmp * T_geo_tmp + tiny_value ))
# Altitude
  altitude  = (radius_proj - radius_equat_planet * T_geo_tmp) * np.cos(latitude) + (carcoord_z - B_geo_tmp) * np.sin(latitude)

  geodetic_coord = [longitude, latitude, altitude]

  return geodetic_coord


def convert_geodetic_cartesian(config, geodetic_coord):
  #
  # Longitude, Latitude, Altitude coordinate (WGS84) --> Cartesian Coordinate
  #
  # Longitude: Radius
  # Latitude: Radius
  # Altitude: meter
  # carcoord_x,y,z: meter
  
  radius_equat_planet = config['planet']['radius']
  ellipticity_planet  = config['planet']['ellipticity']

  longitude = geodetic_coord[0]
  latitude  = geodetic_coord[1]
  altitude  = geodetic_coord[2]

  eccentricity2  = ellipticity_planet*(2.0-ellipticity_planet)
  altitude_geoid = radius_equat_planet / np.sqrt( 1.0 - eccentricity2 * np.sin( latitude )**2 )

  carcoord_x = ( altitude_geoid + altitude ) * np.cos( latitude ) * np.cos( longitude )
  carcoord_y = ( altitude_geoid + altitude ) * np.cos( latitude ) * np.sin( longitude )
  carcoord_z = ( altitude_geoid * (1.0 - eccentricity2 ) + altitude ) * np.sin( latitude )

  cartesian_coord = [carcoord_x, carcoord_y, carcoord_z]

  return cartesian_coord


def set_angle_polar(config, coord):

  # 直交座標系（回転座標系）から極座標における経度(alpha)・緯度(beta)を計算する。WGSの経度・緯度とは厳密にことなるので注意
  #
  radius      = np.sqrt( coord[0]**2 + coord[1]**2 + coord[2]**2 ) 
  angle_beta  = np.arcsin( coord[2]/np.sqrt( coord[0]**2 + coord[1]**2 + coord[2]**2 ) )
  angle_alpha = np.arccos( coord[0]/np.sqrt( coord[0]**2 + coord[1]**2 )  ) * np.sign( coord[1] )

  polar_coord =[radius, angle_beta, angle_alpha]

  return polar_coord


def convert_carteasian_polar(config,vec_input,longitude,latitude):

  # 直交座標上の３次元ベクトルを極座標に変換
  vec_tmp = convert_coordinate_rxyz(vec_input, longitude,'z')
  vec_tmp = convert_coordinate_rxyz(vec_tmp  ,-latitude,'y')
  # この座標回転で1,2,3成分が[高度、経度、緯度]の順番で出るので[経度、緯度、高度]に補正する
  vec_res = exchanege_axis_zxy2xyz(vec_tmp)

  return vec_res


def convert_polar_carteasian(config,vec_input,longitude,latitude):
  
  # 極座標の３次元ベクトルを直交座標上のに変換
  vec_tmp = exchanege_axis_xyz2zxy(vec_input)
  vec_tmp = convert_coordinate_rxyz(vec_tmp,  latitude ,'y')
  vec_res = convert_coordinate_rxyz(vec_tmp ,-longitude,'z')

  return vec_res


def convert_coordinate_rxyz(vec,angle_rot,axis):
  #
  # Coordinate rotation
  # angle_rot: rad
  # axis: x or y or z
  #
  vec_x = vec[0]
  vec_y = vec[1]
  vec_z = vec[2]

  if (axis == 'x') :
    vec_conv_x = vec_x
    vec_conv_y = vec_y*np.cos(angle_rot) + vec_z*np.sin(angle_rot)
    vec_conv_z =-vec_y*np.sin(angle_rot) + vec_z*np.cos(angle_rot)
  elif (axis == 'y') :
    vec_conv_x = vec_x*np.cos(angle_rot) - vec_z*np.sin(angle_rot)
    vec_conv_y = vec_y
    vec_conv_z = vec_x*np.sin(angle_rot) + vec_z*np.cos(angle_rot)
  elif (axis == 'z') :
    vec_conv_x =  vec_x*np.cos(angle_rot) + vec_y*np.sin(angle_rot)
    vec_conv_y = -vec_x*np.sin(angle_rot) + vec_y*np.cos(angle_rot)
    vec_conv_z =  vec_z
  else :
    print( 'Axis is not correct.' )
    print( 'Program stopped. ')
    exit()
  
  return [vec_conv_x,vec_conv_y,vec_conv_z]


def exchanege_axis_xyz2zxy(vec_input):

  vec_conv_1 = vec_input[2]
  vec_conv_2 = vec_input[0]
  vec_conv_3 = vec_input[1]

  return [vec_conv_1,vec_conv_2,vec_conv_3]


def exchanege_axis_zxy2xyz(vec_input):

  vec_conv_1 = vec_input[1]
  vec_conv_2 = vec_input[2]
  vec_conv_3 = vec_input[0]

  return [vec_conv_1,vec_conv_2,vec_conv_3]