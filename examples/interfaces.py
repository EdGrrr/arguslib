"""This script contains increasingly complex examples of how to use the
interface classes to plot radar and aircraft data."""

# %%
import datetime
from pathlib import Path
from arguslib.aircraft import CameraAircraftInterface
from arguslib.instruments.radar import Radar
from arguslib.instruments.camera_array import CameraArray
from arguslib.instruments.camera import Camera
from arguslib.radar.camera_radar_interface import CameraRadarInterface


# %% Plotting the aircraft tracks on a single camera
cai = CameraAircraftInterface.from_campaign("COBALT", "3-7")
dt = datetime.datetime(2025, 3, 9)
adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
cai.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))
cai.show(dt.replace(hour=13, minute=5), tlen=10 * 60)


# %% Plotting a single camera with simultaneous radar data
dt = datetime.datetime(2025, 3, 25, 9, 50, 49)
cri = CameraRadarInterface.from_campaign("COBALT", "3-7")
cri.show(dt)


# %% Plotting multiple cameras simultaneously with radar data
dt = datetime.datetime(2025, 3, 25, 9, 50, 49)
multicam = CameraArray(
    [
        Camera.from_config("COBALT", "3-8"),
        Camera.from_config("COBALT", "5-1"),
        Camera.from_config("COBALT", "5-2"),
        Camera.from_config("COBALT", "5-3"),
        Camera.from_config("COBALT", "5-4"),
    ],
    layout_shape=(3, 3),
)
radar = Radar.from_config("COBALT")
cri = CameraRadarInterface(radar, multicam)
cri.show(dt)


# %% Plot multiple cameras with flight tracks
cai = CameraAircraftInterface(multicam)
dt = datetime.datetime(2025, 3, 9)
adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
cai.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))
cai.show(dt.replace(hour=13, minute=5), tlen=10 * 60, color_icao=True)


# %% Plot multiple cameras with flight tracks and radar data
dt = datetime.datetime(2025, 3, 25, 9, 50, 49)
cri = CameraRadarInterface(radar, cai)
cri.show(dt, kwargs_camera=dict(tlen=10 * 60, color_icao=True))


# %% Plot multiple cameras with flight tracks and radar data with flight tracks
dt = datetime.datetime(2025, 3, 25, 9, 50, 49)
cri = CameraRadarInterface(radar, multicam)
cai = CameraAircraftInterface(cri)
cai.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))
cai.show(dt, tlen=10 * 60, color_icao=True)


# %%
