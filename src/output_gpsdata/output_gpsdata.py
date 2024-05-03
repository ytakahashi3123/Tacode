#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/05/31

import numpy as np

gpx_standard_year   = 2015
gpx_standard_month  = 1
gpx_standard_day    = 1
gpx_standard_hour   = 0
gpx_standard_minute = 0
gpx_standard_second = 0

gpx_case_name = 'Satellite_trajectory'

gpx_xml = '<?xml version="1.0" encoding="UTF-8"?>'
gpx_header = '<gpx xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns="http://www.topografix.com/GPX/1/0" xsi:schemaLocation="http://www.topografix.com/GPX/1/1 http://www.topografix.com/GPX/1/1/gpx.xsd">'
gpx_footer = '</gpx>'

gpx_trk_s = '<trk>'
gpx_trk_e = '</trk>'
gpx_name_s = '<name>'
gpx_name_e = '</name>'
gpx_number_s = '<number>'
gpx_number_e = '</number>'
gpx_trkseg_s = '<trkseg>'
gpx_trkseg_e = '</trkseg>'

gpx_trkpt_1 = '<trkpt lat="'
gpx_trkpt_2 = '" lon="'
gpx_trkpt_3 = '"><ele>'
gpx_trkpt_4 = '</ele><time>'
gpx_trkpt_5 = '</time></trkpt>'

newline_code ='\n'

rad2deg = 180.0/np.pi


def output_routine(config, iteration, coordinate_dict, velocity_dict):

  if config['post_process']['flag_output_gpx'] :
    output_gpx(config, iteration, coordinate_dict, velocity_dict)

  if config['post_process']['kml']['flag_output'] :
    output_kml(config, iteration, coordinate_dict, velocity_dict)

  return


def output_kml(config, iteration, coordinate_dict, velocity_dict):

  import simplekml as simplekml

  coordinate_geod = coordinate_dict['geodetic']

  config_kml           = config['post_process']['kml']
  filename_tmp         = config['post_process']['directory_output'] + '/' + config_kml['filename_output']
  linestyle_color_kml  = config_kml['linestyle_color']
  linestyle_width_kml  = config_kml['linestyle_width']
  extrude_kml          = config_kml['extrude']
  frequency_output_kml = config_kml['frequency_output']

  dt = config['time_integration']['timestep_constant']


  print('Writing KML file... ', filename_tmp)

  gps_linestrings = []
  for n in range(0,iteration):
    if n%frequency_output_kml == 0 :
      time_tmp = float(n)*dt
      lat_tmp = coordinate_geod[n][1]*rad2deg
      lon_tmp = coordinate_geod[n][0]*rad2deg
      alt_tmp = int(coordinate_geod[n][2])
      gps_linestrings.append( [str(time_tmp), str(lon_tmp), str(lat_tmp), str(alt_tmp), str(linestyle_color_kml)] )

  kml = simplekml.Kml()

  for linestring in gps_linestrings:
#    ls = kml.newlinestring(name=unicode(linestring[0], 'utf-8'))
    ls = kml.newlinestring(name='Time_'+linestring[0])
    ls.style.linestyle.color = linestring[4]
    ls.style.linestyle.width = linestyle_width_kml
    ls.extrude = extrude_kml
    ls.altitudemode = simplekml.AltitudeMode.absolute
    ls.coords = [(float(linestring[1]), float(linestring[2]), float(linestring[3]))]

  kml.save(filename_tmp)

  return


def output_gpx(config, iteration, coordinate_dict, velocity_dict):

  import gpxpy as gpxpy

  if config['post_process']['flag_output_gpx'] :

    coordinate_geod = coordinate_dict['geodetic']

    filename_tmp = config['post_process']['directory_output'] + '/' + config['post_process']['filename_output_gpx']

    gpx_file = open(filename_tmp, 'w')

    gpx = gpxpy.gpx.GPX()

    # Create first track in our GPX:
    gpx_track = gpxpy.gpx.GPXTrack()
    gpx.tracks.append(gpx_track)

    # Time
    #time_set = datetime(2018,1,1,0,0)   # (1)initial date
    #gps_delta = timedelta(seconds=1)    # (2)gps period: 1s

    # Create first segment in our GPX track:
    gpx_segment = gpxpy.gpx.GPXTrackSegment()
    gpx_track.segments.append(gpx_segment)

    # Create points:
    for n in range(0,iteration):
      lat_tmp = coordinate_geod[n][1]*rad2deg
      lon_tmp = coordinate_geod[n][0]*rad2deg
      alt_tmp = int(coordinate_geod[n][2])
      gpx_segment.points.append(gpxpy.gpx.GPXTrackPoint(lat_tmp, lon_tmp, elevation=alt_tmp)) 

    #gpx_segment.points.append(gpxpy.gpx.GPXTrackPoint(2.1234, 5.1234, elevation=1234))
    #gpx_segment.points.append(gpxpy.gpx.GPXTrackPoint(2.1235, 5.1235, elevation=1235))
    #gpx_segment.points.append(gpxpy.gpx.GPXTrackPoint(2.1236, 5.1236, elevation=1236))

    # You can add routes and waypoints, too...
    #print('Created GPX:', gpx.to_xml())
    gpx_file.write(gpx.to_xml())
    gpx_file.close()

  return


def output_gpx_native(config, iteration, coordinate_dict, velocity_dict):

  if config['post_process']['flag_output_gpx'] :

    coordinate_geod = coordinate_dict['geodetic']

    filename_tmp = config['post_process']['directory_output'] + '/' + config['post_process']['filename_output_gpx']
    print( 'Writing GPX file... ', filename_tmp )
    file = open(filename_tmp, "w")

    # --Header
    file.write( gpx_xml    + newline_code )
    file.write( gpx_header + newline_code )
    file.write( gpx_trk_s  + newline_code )
    file.write( gpx_name_s + gpx_case_name+gpx_name_e + newline_code )
    file.write( gpx_number_s + '1' + gpx_number_e + newline_code )
    file.write( gpx_trkseg_s  + newline_code )

    for n in range(0,iteration):
      lat_tmp = coordinate_geod[n][1]*rad2deg
      lon_tmp = coordinate_geod[n][0]*rad2deg
      alt_tmp = coordinate_geod[n][2]
      gpx_tmp = gpx_trkpt_1 + str(lat_tmp) + gpx_trkpt_2 + str(lon_tmp) +  gpx_trkpt_3 + str(alt_tmp) + gpx_trkpt_4 + '2015-01-01T00:00:' + str(n) + 'Z' + gpx_trkpt_5
      file.write(gpx_tmp + newline_code)

    file.write( gpx_trkseg_e + newline_code )
    file.write( gpx_trk_e + newline_code )
    file.write( gpx_footer )

  return 




