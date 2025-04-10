# %%

from datetime import datetime
from argusflightserver.fleet import Aircraft
from arguslib.aircraft.aircraft_interface import AircraftInterface
from arguslib.instruments.calibration import PerspectiveProjection, focal_length_px
from arguslib.instruments.instruments import Position
from arguslib.instruments.camera import Camera
from arguslib.video.locator import CameraData
import matplotlib.pyplot as plt
from pathlib import Path

dt = datetime(2025, 4, 2, 9, 50)

cam = Camera.from_config("COBALT", "2-11")

cam.show(dt)


# %%

cai = AircraftInterface(cam)

adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
cai.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))

cai.show(dt, tlen=10 * 60, trail_kwargs={"plot_kwargs": {"max_range_km": 1000}})
