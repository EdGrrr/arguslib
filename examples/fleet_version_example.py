'''This file is for testing the different 'fleet' classes'''
import datetime as dt
from datetime import time
import pandas as pd
import os
import re
from pathlib import Path
import random
import calendar
from collections import defaultdict
from arguslib import AutomaticADSBAircraftInterface, DirectUndistortedCamera
from arguslib.camera.camera_array import CameraArray
from arguslib.camera.undistorted_camera import UndistortedCamera
from arguslib.radar.radar import Radar
from arguslib.radar.radar_interface import RadarInterface
from arguslib.misc.times import get_camera_configs, get_available_scans_for_day

# We need the locator to find scan times
from arguslib.radar.locator import RadarData, initialise_locator
import matplotlib.pyplot as plt

# Compare the different fleet implementations
from arguslib.aircraft.fleet import Fleet, FleetOld


campaign_name = 'COBALT'
tlen = 60

camera_configs = get_camera_configs()
camera_config = camera_configs['Single Camera (3-7)']

radar = Radar.from_config(campaign_name)
cri = RadarInterface(radar, camera_config)
cai = AutomaticADSBAircraftInterface(
    UndistortedCamera.from_config("COBALT", "3-7"))
# Update the instrument that the AircraftInterface wraps.
# This is what you meant by "replace the aircraftinterface.camera property"
cai.camera = cri

# Setup some scan times
year, month, day = 2025, 5, 1
scan_date = dt.datetime(year, month, day)
start_time = time(6, 0)
end_time = time(18, 0)
scan_times = get_available_scans_for_day(
    campaign_name, "rhi",
    scan_date, start_time, end_time
)


# This is the 'classic' implementation, based on dicts/lists
fleet = FleetOld(variables=['lon', 'lat', 'alt_geom', 'oat', 'ws', 'wd'])
fleet.load_output('/net/hardin/disk1/Data/ADS-B/COBALT/20250501_ADS-B')


# The trackerlib-based implementation (which should be faster)
nfleet = Fleet(variables=['lon', 'lat', 'alt_geom', 'oat', 'ws', 'wd'])
# This implementation pre-calculates the wind filters, so it is
# specified where the aircraft data is loaded.
nfleet.load_output('/net/hardin/disk1/Data/ADS-B/COBALT/20250501_ADS-B', wind_filter=10)

# Make a set of intersection plots, showing the two fleet objects on a radar and camera.
for ind in range(10, 100, 10):
    scan_time = scan_times[ind]
    print("\n Scan time: "+f"{scan_time}" )

    # New Fleet object
    cai.fleet = nfleet
    ax = cai.show(
        scan_time,
        color_icao=True,
        trail_kwargs=dict(label_acft=True, advection_winds='aircraft', wind_filter=10),
        tlen=tlen * 60)
    ax[-1].legend(fontsize='x-small', ncol=3)
    plt.figure(1).savefig(f'./test_{scan_time.isoformat()}.jpg')
    plt.close('all')

    # Old/classic object
    cai.fleet = fleet
    ax = cai.show(
        scan_time,
        color_icao=True,
        trail_kwargs=dict(label_acft=True, advection_winds='aircraft', wind_filter=10),
        tlen=tlen * 60)
    ax[-1].legend(fontsize='x-small', ncol=3)
    plt.figure(1).savefig(f'./test_{scan_time.isoformat()}_old.jpg')
    plt.close('all')

