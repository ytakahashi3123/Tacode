#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/05/23

import sys as sys
import numpy as np
from orbital.orbital import orbital
from general.general import get_setting


def get_bank_angle(config, time_elapsed):
  #
  # バンク角 [deg.]。既定は satellite.bank_angle（定数）。
  #
  # **satellite.bank_angle_table があればそちらを使う。** 表は [時刻 s, バンク角 deg.] の
  # 並びで、時刻は計算の経過時間。時刻について線形内挿し、両端の外は端の値で止める。
  # 実機のロール変調を入力にするためのもので（validation/apollo10）、誘導則そのものは
  # 持たない。
  #
  section = config['satellite']
  table = get_setting(section, 'bank_angle_table', None)
  if table is None :
    return float( get_setting(section, 'bank_angle', 0.0) )

  node = np.array(table, dtype=float)
  if node.ndim != 2 or node.shape[1] != 2 or node.shape[0] < 2 :
    print('satellite.bank_angle_table must be a list of at least two [time, angle] pairs.')
    print('Program stopped.')
    sys.exit(1)
  if np.any( np.diff(node[:,0]) <= 0.0 ) :
    print('The time of satellite.bank_angle_table must increase.')
    print('Program stopped.')
    sys.exit(1)

  return float( np.interp(time_elapsed, node[:,0], node[:,1]) )


def get_lift_direction(coordinate, velocity_aero, angle_bank):
  #
  # 揚力の単位ベクトル（ECEF 成分）。
  #
  # 揚力は対気速度に直交する。バンク角 0 で「速度と鉛直上向きが張る面の中で上向き側」、
  # バンク角は**速度ベクトルまわり**に測り、**正で進行方向の右側**へ倒す
  # （bank 90 deg. で水平右、180 deg. で下向き）。鉛直は地心の上向き r/|r| を使う
  # （初期速度・風と同じ地心ローカル系の規約。README の Reference frame の節）。
  #
  # 速度が鉛直と平行なときは「面」が定まらない。そのとき揚力は 0 とする
  # （真上・真下に飛んでいる瞬間だけで、軌道計算では起きない）。
  #
  direction_velocity = velocity_aero/np.linalg.norm(velocity_aero)
  direction_up       = coordinate/np.linalg.norm(coordinate)

  # 速度に直交する成分（バンク 0 の向き）
  component_up = direction_up - np.dot(direction_up, direction_velocity)*direction_velocity
  magnitude_up = np.linalg.norm(component_up)
  if magnitude_up <= 0.0 :
    return np.zeros(3)
  component_up = component_up/magnitude_up

  # 進行方向の右側（東向きに飛んでいれば南）
  component_right = np.cross(direction_velocity, component_up)

  return np.cos(angle_bank)*component_up + np.sin(angle_bank)*component_right


def force_initialsettings(config):

  force = np.zeros(5*3).reshape(5,3)

  return force


def force_routine(config, coordinate, velocity, mass_satellite, area_satellite, cdmean_aerodynamic, density_factor, density, force, force_aerodynamic=None, velocity_air=None, angle_bank=None):
  
  potential_factor     = config['planet']['potential_factor']
  radius_equat_planet  = config['planet']['radius']
  gravity_const_planet = config['planet']['gravitational_constant']
  mass_planet          = config['planet']['mass']
  rotation_rate_planet = config['planet']['rotation_rate']


  # Gravitational force
  radius_pmass    = np.sqrt( coordinate[0]**2 + coordinate[1]**2 + coordinate[2]**2 )
  radius_by_coord = radius_equat_planet/radius_pmass                   # a_e/r
  gme_by_radius2  = gravity_const_planet*mass_planet/(radius_pmass**2) # G*M_e/r2

  sin_beta = coordinate[2]/radius_pmass
  cos_beta = np.sqrt( coordinate[0]**2 + coordinate[1]**2 )/radius_pmass

  # arctan2 を使う（経度 180 度での誤算出を避ける）
  angle_long = np.arctan2( coordinate[1], coordinate[0] )
  sin_labd   = np.sin( angle_long )
  cos_labd   = np.cos( angle_long )

  J2  = potential_factor['J2']
  J22 = potential_factor['J22']
  J3  = potential_factor['J3']
  J4  = potential_factor['J4']
  labd22           = potential_factor['Lambda22']*orbital.deg2rad
  sin_2labd_labd22 = np.sin( 2.0*(angle_long + labd22 ) )
  cos_2labd_labd22 = np.cos( 2.0*(angle_long + labd22 ) )

  force_g_r = gme_by_radius2 * (- 1.0                                                                           \
                                + 1.5    *radius_by_coord**2*J2  *( 3.0*sin_beta**2 -  1.0 )                    \
                                + 9.0    *radius_by_coord**2*J22 *(     cos_beta**2         )*cos_2labd_labd22  \
                                + 2.0    *radius_by_coord**3*J3  *( 5.0*sin_beta**3 -  3.0*sin_beta )           \
                                + 5.0/8.0*radius_by_coord**4*J4  *(35.0*sin_beta**4 - 30.0*sin_beta**2 + 3.0 )  \
                                )
  force_g_a = gme_by_radius2 * (  6.0    *radius_by_coord**2*J22 *(     cos_beta            )*sin_2labd_labd22 )
  force_g_b = gme_by_radius2 * (                                                                                \
                                - 1.0    *radius_by_coord**2*J2  *( 3.0*sin_beta*cos_beta  )                    \
                                + 6.0    *radius_by_coord**2*J22 *(     sin_beta*cos_beta  )*cos_2labd_labd22   \
                                - 0.5    *radius_by_coord**3*J3  *(15.0*sin_beta**2 -  3.0 )*cos_beta           \
                                - 0.5    *radius_by_coord**4*J4  *(35.0*sin_beta**3 - 15.0*sin_beta )*cos_beta  \
                               ) 

  force[1,0] =  ( force_g_r*cos_beta - force_g_b*sin_beta )*cos_labd - force_g_a*sin_labd
  force[1,1] =  ( force_g_r*cos_beta - force_g_b*sin_beta )*sin_labd + force_g_a*cos_labd
  force[1,2] =    force_g_r*sin_beta + force_g_b*cos_beta


  # Coriolis
  force[2,0]  = 2.0*rotation_rate_planet*velocity[1]
  force[2,1]  =-2.0*rotation_rate_planet*velocity[0]
  force[2,2]  = 0.0


  # Centrifugal
  force[3,0] = rotation_rate_planet**2*coordinate[0]
  force[3,1] = rotation_rate_planet**2*coordinate[1]
  force[3,2] = 0.0


  # Aerodynamic (Fx = 1/2 rho U^2 * Ux/U)
  # --6 自由度計算では姿勢に応じた空力加速度（ECEF 成分）が呼び出し側で求まっているので、
  #   それを与える。3 自由度では従来どおり速度方向の抗力のみ。
  if force_aerodynamic is None :
    # 風があるときは対気速度で抗力を作る。上の重力・コリオリ力・遠心力は ECEF 速度のまま
    # （velocity_air が None なら velocity をそのまま使うので、風を切れば従来とビット単位で同じ）
    velocity_aero  = velocity if velocity_air is None else velocity_air
    velocity_mag   = np.sqrt(velocity_aero[0]**2+velocity_aero[1]**2+velocity_aero[2]**2)
    fact_aero      = 0.50*density_factor*density*velocity_mag*area_satellite*cdmean_aerodynamic/mass_satellite
    #"density factor" added by Tomoki Sakai 2023/2/3
    force[4,0] = -fact_aero * velocity_aero[0]
    force[4,1] = -fact_aero * velocity_aero[1]
    force[4,2] = -fact_aero * velocity_aero[2]

    # 揚力。**既定は satellite.lift_coefficient = 0 で、そのときここは一切通らない**
    # （＝揚力を入れなければ従来の出力とビット単位で同じ）。
    # 6 自由度では姿勢から力が決まるので、この枝には来ない。
    coefficient_lift = float( get_setting(config['satellite'], 'lift_coefficient', 0.0) )
    if coefficient_lift != 0.0 :
      # バンク角は呼び出し側（solver）が段の時刻で引いて渡す。None のときは
      # 時刻に依らない設定として読む
      angle_tmp  = get_bank_angle(config, 0.0) if angle_bank is None else angle_bank
      fact_lift  = 0.50*density_factor*density*velocity_mag**2*area_satellite*coefficient_lift/mass_satellite
      direction_lift = get_lift_direction(coordinate, velocity_aero, angle_tmp*orbital.deg2rad)
      force[4,0] = force[4,0] + fact_lift*direction_lift[0]
      force[4,1] = force[4,1] + fact_lift*direction_lift[1]
      force[4,2] = force[4,2] + fact_lift*direction_lift[2]
  else :
    force[4,0] = force_aerodynamic[0]
    force[4,1] = force_aerodynamic[1]
    force[4,2] = force_aerodynamic[2]


  #print(force[1,:],force[2,:],force[3,:])

  # Total
  force[0,:] = force[1,:] + force[2,:] + force[3,:] + force[4,:]


  return force

