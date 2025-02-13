import requests
import numpy as np
import math
from io import BytesIO
import geo
from instruments import Camera, Radar, Position

cameras = {'3-7': Camera.from_filename('cam1_calibration.yml',
                                       Position(-1.439252, 51.146668, 0.1),
                                       rotation=10, scale_factor=3040/500,
                                       camera_type='allsky'),
           '3-8': Camera.from_filename('cam1_calibration.yml',
                                       Position(-1.438511, 51.149064, 0.1),
                                       rotation=10, scale_factor=3040/500,
                                       camera_type='allsky'),
           '5-1': Camera.from_filename('cam1_calibration.yml',
                                       Position(-1.363346, 50.952369, 0.1),
                                       rotation=10, scale_factor=3040/500,
                                       camera_type='allsky'),
           '5-2': Camera.from_filename('cam1_calibration.yml',
                                       Position(-1.141453, 51.147445, 0.1),
                                       rotation=10, scale_factor=3040/500,
                                       camera_type='allsky'),
           '5-3': Camera.from_filename('cam1_calibration.yml',
                                       Position(-1.443658, 51.320375, 0.1),
                                       rotation=10, scale_factor=3040/500,
                                       camera_type='allsky'),
           '5-4': Camera.from_filename('cam1_calibration.yml',
                                       Position(-1.10155, 51.29585, 0.1),
                                       rotation=10, scale_factor=3040/500,
                                       camera_type='allsky-planned'),
           '5-5': Camera.from_filename('cam1_calibration.yml',
                                       Position(-1.78568, 51.09967, 0.1),
                                       rotation=10, scale_factor=3040/500,
                                       camera_type='allsky-planned'),
           '5-6': Camera.from_filename('cam1_calibration.yml',
                                       Position(-1.77448, 51.34033, 0.1),
                                       rotation=10, scale_factor=3040/500,
                                       camera_type='allsky-planned'),
           }

radars = {'kepler': Radar(0.6, # beamwidth in degrees
                          Position(-1.439339, 51.144752, 0.1), rotation=0)}

def svg_path(points):
    return 'M'+' L'.join(['{:.1f} {:.1f}'.format(b[0], b[1]) for b in points])

def range_rings(camera, ranges):
    range_out = {}
    # Plot range rings
    alt = 10
    for rd in ranges:
        rl = []
        for az in range(0, 361, 10):
            elev, dist = np.rad2deg(np.arctan2(alt, rd)), np.sqrt(alt**2 +rd**2)
            rl.append(camera.iead_to_pix(*camera.gead_to_iead(elev, az, dist)))
        rl = np.array(rl)
        range_out[rd] = rl
    return range_out


def radar_beam(camera: Camera, radar: Radar, radar_elev, radar_azimuth):
    # Get radar direction
    radar_elev, radar_azimuth = 40, -180

    # Plot radar direction
    radar_beam_locs = radar.beam(radar_elev, radar_azimuth, [0, 8, 35])

    radar_out = {}

    radar_beam_px = np.array([camera.target_pix(pt) for pt in radar_beam_locs[1]])
    radar_out['radar_near'] = radar_beam_px

    radar_beam_px = np.array([camera.target_pix(pt) for pt in radar_beam_locs[2]])
    radar_out['radar_far'] = radar_beam_px

    radar_beam_px = np.array([camera.target_pix(pt) for pt in radar_beam_locs[:, 0]])
    radar_out['radar_ray'] = radar_beam_px

    return radar_out


def aircraft_trails(camera: Camera):
    # Get aircraft data
    response = requests.get('https://opendata.adsb.fi/api/v2/'+
                            f'lat/{camera.position.lat}/' +
                            f'lon/{camera.position.lon}/dist/100')
    aircraft = response.json()

    acft_out = {}

    for acft in aircraft['aircraft']:
        try:
            if acft['alt_geom']<26000:
                continue
            dist = geo.haversine(camera.position.lon, camera.position.lat, acft['lon'], acft['lat'])
            acft_out[acft['hex']] = {
                'acft_dist': dist,
                'alt': acft['alt_geom'],
                'flight': acft['flight'],
                't': acft['t']}

            if dist < 90:
                pl = camera.target_pix(Position(acft['lon'], acft['lat'], acft['alt_geom']/(3.33*1000)))
                acft_out[acft['hex']]['acft_px'] = pl

            # gs is in knots, not 100% about wind, but probably also knots
            wind_u, wind_v = (0.51444*acft['ws']*np.sin(np.deg2rad(acft['wd'])),
                              0.51444*acft['ws']*np.cos(np.deg2rad(acft['wd'])))
            acft_u, acft_v = (0.51444*acft['gs']*np.sin(np.deg2rad(acft['track'])),
                              0.51444*acft['gs']*np.cos(np.deg2rad(acft['track'])))

            #To get the track location, we wind back the aircraft position,
            # then integrate the wind forwards. A simple solution in then
            times = np.arange(-3600, 0.1, 10) # Seconds
            
            track_offset_km = np.array([
                acft_u-wind_u, acft_v-wind_v])[:, None] * times /1000
            track_pos = np.array(geo.xy_offset_to_ll(acft['lon'], acft['lat'], *track_offset_km))
            pl_track = np.array([camera.target_pix(Position(tp[0], tp[1], acft['alt_geom']/(3.33*1000)))
                                 for tp in track_pos.T])
            dist = geo.haversine(camera.position.lon, camera.position.lat, track_pos[0], track_pos[1])
            acft_out[acft['hex']]['track_px'] = pl_track.T

            spread_velocity = 3 # Spreading velocity in ms-1 (https://doi.org/10.1029/98JD02594) 
            spread_direction = acft['track']+90 # Simple for now - mostly governed by aircraft speed/direction
            spread_u, spread_v = spread_velocity*(np.array([np.sin(np.deg2rad(spread_direction)),
                                                            np.cos(np.deg2rad(spread_direction)),]))
            track_leftoffset_km = np.array([
                acft_u-wind_u+spread_u, acft_v-wind_v+spread_v])[:, None] * times /1000
            track_leftpos = np.array(geo.xy_offset_to_ll(acft['lon'], acft['lat'], *track_leftoffset_km))
            pl_track = np.array([camera.target_pix(Position(tp[0], tp[1], acft['alt_geom']/(3.33*1000)))
                                 for tp in track_leftpos.T])
            acft_out[acft['hex']]['track_left_px'] = pl_track.T

            track_rightoffset_km = np.array([
                acft_u-wind_u-spread_u, acft_v-wind_v-spread_v])[:, None] * times /1000
            track_rightpos = np.array(geo.xy_offset_to_ll(acft['lon'], acft['lat'], *track_rightoffset_km))
            pl_track = np.array([camera.target_pix(Position(tp[0], tp[1], acft['alt_geom']/(3.33*1000)))
                                 for tp in track_rightpos.T])
            acft_out[acft['hex']]['track_right_px'] = pl_track.T

            

        except KeyError:
            pass

    return acft_out
