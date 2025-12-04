# %%

from datetime import datetime
import matplotlib.pyplot as plt
from pathlib import Path

from arguslib import AircraftInterface, Camera
from arguslib.misc.plotting import show_the_sun

cam = Camera.from_config("COBALT", "2-11")

for h in range(13, 18):
    dt = datetime(2025, 5, 1, h, 10)
    ax = cam.show(dt)
    show_the_sun(ax, cam, plotting_method='scatter')


# %%

cai = AircraftInterface(cam)

adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
cai.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))

cai.show(dt, tlen=30 * 60, trail_kwargs={"wind_filter":10, "plot_kwargs": {"max_range_km": 1000, 'lw': 0.5, 'alpha': 0.5}})

# # %%

# %%
