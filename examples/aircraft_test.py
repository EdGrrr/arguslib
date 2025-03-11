import arguslib
import arguslib.aircraft
import arguslib.instruments
import arguslib.instruments.statuslib
from arguslib.instruments.instruments import Position
from arguslib.aircraft.fleet import Fleet
import datetime
import matplotlib.pyplot as plt
import numpy as np

def get_pcl_details(fleet, acft, dtime, seconds_since_formation):
    trail = fleet.aircraft[acft].get_trail(dtime, tlen=seconds_since_formation)
    metadata = fleet.aircraft[acft].get_data(dtime, ['alt_geom', 'oat'], tlen=seconds_since_formation)
    data = {v:metadata[v][0] for v in metadata.keys()}
    data['lon'] = trail[0][0]
    data['lat'] = trail[1][0]
    return data

def calc_ead(camera, acft_data):
    return camera.position.target_ead(Position(
        acft_data['lon'],
        acft_data['lat'],
        acft_data['alt_geom']/(3.33*1000)))

fleet = Fleet(variables=[
    'lon', 'lat',
    'alt_baro', 'alt_geom', 'geom_rate',
    'tas', 'gs', 'ws',
    'track', 'true_heading', 'wd',
    'oat'])


fleet_file_pattern = '/disk1/Data/ADS-B/COBALT/20250305_ADS-B'

#Track contrail segment location and radar scans for this aircraft
time_range = [datetime.datetime(2025, 3, 5, 12),
              datetime.datetime(2025, 3, 5, 14)]
track_step = 60
acft = '405f13'
contrail_form_time = datetime.datetime(2025, 3, 5, 12, 38, 5)

fleet.load_output(fleet_file_pattern)

camstr = '3-7' # As closest to Kepler
camera = arguslib.instruments.statuslib.cameras[camstr]


dtime = time_range[0]
output = []
while dtime< time_range[1]:
    if dtime > contrail_form_time+datetime.timedelta(seconds=60):
        acft_data = get_pcl_details(fleet, acft, dtime, (dtime-contrail_form_time).total_seconds())
        ead_data = calc_ead(camera, acft_data)
        acft_data['elevation'] = ead_data[0]
        azi = ead_data[1]
        if azi < 0: azi += 360
        acft_data['azimuth'] = azi
        acft_data['range'] = ead_data[2]
        output.append(acft_data)
        print((dtime-contrail_form_time).total_seconds())
    dtime += datetime.timedelta(seconds=track_step)

output = {v: [output[a][v] for a in range(len(output))] for v
          in output[0].keys()}

kepler_commands = np.array([
    [
        (datetime.datetime(2025, 3, 5, 12, 44, 0)-contrail_form_time).total_seconds(),
        266, 66, 51.14, -1.51],
    [
        (datetime.datetime(2025, 3, 5, 12, 49, 0)-contrail_form_time).total_seconds(),
        256, 39, 51.12, -1.62],
    [
        (datetime.datetime(2025, 3, 5, 12, 54, 00)-contrail_form_time).total_seconds(),
        254, 26.1, 51.09, -1.74],
    [
        (datetime.datetime(2025, 3, 5, 12, 59, 00)-contrail_form_time).total_seconds(),
        253, 19.36, 51.06, -1.86],
    [
        (datetime.datetime(2025, 3, 5, 13, 4, 0)-contrail_form_time).total_seconds(),
        252.5, 15.3, 51.04, -1.97]
])


plt.subplot(311)
plt.plot(output['lon'], output['lat'], label='trail')
plt.scatter(-1.44, 51.16, c='green', label='Kepler') # Chilbolton
plt.scatter(kepler_commands[:, 4], kepler_commands[:, 3], c='r', label='Kepler Commands')
plt.xlabel('lon')
plt.ylabel('lat')
plt.legend()

plt.subplot(312)
plt.plot(output['elevation'])
plt.scatter(kepler_commands[:, 0]/60, kepler_commands[:, 2], c='r')
plt.xlabel('Time since formation')
plt.ylabel('Elevation')

plt.subplot(313)
plt.plot(output['azimuth'])
plt.scatter(kepler_commands[:, 0]/60, kepler_commands[:, 1], c='r')
plt.xlabel('Time since formation')
plt.ylabel('Azimuth')


plt.gcf().suptitle(acft + ' ' + contrail_form_time.isoformat())
plt.show()

# acft_data = {}
# for n, name in enumerate(fleet.variables):
#     acft_data[name] = fleet.aircraft[acft].pos.positions[:, n]



# daysec = (dtime.hour*3600+dtime.minute*60+dtime.second)
# offset = daysec%15
# index = daysec//15
# length = (tlen/15)+2 # For before+after slots
# startind = int(max(0, index-length))
# acft_time = fleet.aircraft[acft].pos.times

# times = (daysec-acft_time)[startind:index]

# gs = self.positions[startind:index, self.variables.index('gs')]
# track = self.positions[startind:index, self.variables.index('track')]
# ws = self.positions[startind:index, self.variables.index('ws')]
# wd = self.positions[startind:index, self.variables.index('wd')]
# wind_u, wind_v = (0.51444*ws*np.sin(np.deg2rad(wd)),
#                   0.51444*ws*np.cos(np.deg2rad(wd)))
# if wind_filter>0:
#     def wind_conv_filter(wval, wind_filter):
#         filtered = (np.convolve(np.ones(int(wind_filter)), wval, mode='same') /
#                     np.convolve(np.ones(int(wind_filter)), np.ones(len(wval)), mode='same'))
#         if wind_filter>len(wval):
#             difflen = int(wind_filter)-len(wval)
#             filtered = filtered[difflen//2:(difflen//2+len(wval))]
#         return filtered
#     wind_u = wind_conv_filter(wind_u, wind_filter)
#     wind_v = wind_conv_filter(wind_v, wind_filter)

# acft_u, acft_v = (0.51444*gs*np.sin(np.deg2rad(track)),
#                   0.51444*gs*np.cos(np.deg2rad(track)))

# track_offset_km = np.array([
#     wind_u, wind_v]) * times /1000

# lon = self.positions[startind:index, self.variables.index('lon')]
# lat = self.positions[startind:index, self.variables.index('lat')]

# track_pos = np.array(geo.xy_offset_to_ll(lon, lat, *track_offset_km))
