# %%

from datetime import datetime
import matplotlib.pyplot as plt
from pathlib import Path

from arguslib import AircraftInterface, Camera
from arguslib.misc.plotting import show_the_sun

dt = datetime(2025, 5, 1, 16, 10)
cam = Camera.from_config("COBALT", "2-11")

# This is the "arguslib" way to show the camera data. At the moment, this is producing something mostly-right, but with the *left and right* sides flipped!
ax = cam.show(dt)
show_the_sun(ax, cam, plotting_method='scatter')


# %%

cai = AircraftInterface(cam)

adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
cai.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))

cai.show(dt, tlen=30 * 60, trail_kwargs={"wind_filter":10, "plot_kwargs": {"max_range_km": 1000, 'lw': 0.5, 'alpha': 0.5}})

# # %%

# %%
