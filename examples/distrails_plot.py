import datetime as dt
from tqdm import tqdm
import pandas as pd
import numpy as np
import os
import re
from pathlib import Path
import random
from arguslib import AutomaticADSBAircraftInterface, DirectUndistortedCamera, AircraftInterface
from arguslib.aircraft.fleet import FleetOld
from arguslib.camera.camera_array import CameraArray
from arguslib.camera.undistorted_camera import UndistortedCamera
from arguslib.radar.radar import Radar
from arguslib.radar.radar_interface import RadarInterface
from arguslib.radar.locator import RadarData, initialise_locator
from arguslib.misc.times import (get_camera_configs, get_available_scans_for_day,
                                 find_active_year_months, find_available_dates_in_month)

import matplotlib.pyplot as plt
from glob import glob
import pyart


campaign_name = 'COBALT'
tlen = 60*60
intersect_chunk_size = 20  # s


camera_configs = get_camera_configs()
camera_config = camera_configs['Single Camera (3-7)']

radar = Radar.from_config(campaign_name)
cri = RadarInterface(radar, camera_config)
cam = UndistortedCamera.from_config("COBALT", "3-7")

time1 = dt.datetime.fromisoformat("2025-05-01T07:25:06")
time2 = dt.datetime.fromisoformat("2025-05-01T07:30:23")

ax = plt.subplot(211, projection="polar")
cam.show(time2, ax=ax)

ax = plt.subplot(413)
pradar = radar.data_loader.get_pyart_radar(time1)
display = pyart.graph.RadarDisplay(pradar)
display.plot('DBZ', ax=ax, vmax=0, vmin=-30, colorbar_flag=False, title_flag=False)
display.set_limits(xlim=(-7, -1), ylim=(6, 10))
ax.set_ylabel('Altitude')
ax.set_xlabel('')
ax.text(0.02, 0.02, time1, transform=ax.transAxes)

ax = plt.subplot(414)
pradar = radar.data_loader.get_pyart_radar(time2)
display = pyart.graph.RadarDisplay(pradar)
display.plot('DBZ', ax=ax, vmax=0, vmin=-30, colorbar_flag=False, title_flag=False)
display.set_limits(xlim=(-3, 3), ylim=(6, 10))
ax.set_ylabel('Altitude')
ax.text(0.02, 0.02, time2, transform=ax.transAxes)

plt.gcf().set_size_inches(4, 7)
plt.savefig('distrail_example.jpg', dpi=200, bbox_inches='tight')
plt.show()

