"""This script contains increasingly complex examples of how to use the
interface classes to plot radar and aircraft data."""

# %%
import datetime
from pathlib import Path
from arguslib.aircraft import AircraftInterface
from arguslib.instruments.radar import Radar
from arguslib.instruments.camera_array import CameraArray
from arguslib.radar.radar_interface import RadarInterface

adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")


# %% Plotting the aircraft tracks on a single camera
cai = AircraftInterface.from_campaign("COBALT", "3-7")
dt = datetime.datetime(2025, 3, 9)
cai.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))
cai.show(dt.replace(hour=13, minute=5), tlen=10 * 60)


# %% Plotting a single camera with simultaneous radar data
dt = datetime.datetime(2025, 3, 25, 9, 50, 49)
cri = RadarInterface.from_campaign("COBALT", "3-7")
cri.show(dt)


# %% Plotting multiple cameras simultaneously with radar data
dt = datetime.datetime(2025, 3, 25, 9, 50, 49)
multicam = CameraArray.from_config("COBALTArray")
radar = Radar.from_config("COBALT")
cri = RadarInterface(radar, multicam)
cri.show(dt)


# %% Plot multiple cameras with flight tracks
cai = AircraftInterface(multicam)
dt = datetime.datetime(2025, 3, 9)
cai.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))
cai.show(dt.replace(hour=13, minute=5), tlen=10 * 60, color_icao=True)


# %% Plot multiple cameras with flight tracks and radar data
dt = datetime.datetime(2025, 3, 25, 9, 50, 49)
cri = RadarInterface(radar, cai)
cri.show(dt, kwargs_camera=dict(tlen=10 * 60, color_icao=True))


# %% Plot multiple cameras with flight tracks and radar data with flight tracks
dt = datetime.datetime(2025, 3, 25, 9, 50, 49)
cri = RadarInterface(radar, multicam)
cai = AircraftInterface(cri)
cai.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))
cai.show(dt, tlen=10 * 60, color_icao=True)

# %%
