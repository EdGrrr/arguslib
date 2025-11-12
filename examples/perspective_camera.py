# %%

from datetime import datetime
from pathlib import Path

from arguslib import AircraftInterface, Camera

dt = datetime(2025, 4, 2, 9, 50)

cam = Camera.from_config("COBALT", "2-11")

cam.show(dt)


# %%

cai = AircraftInterface(cam)

adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
cai.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))

cai.show(dt, tlen=10 * 60, trail_kwargs={"plot_kwargs": {"max_range_km": 1000}})
