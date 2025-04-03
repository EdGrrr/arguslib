# %%
import datetime
from pathlib import Path
from arguslib.aircraft.aircraft_interface import AircraftInterface

int = AircraftInterface.from_campaign("COBALT", "3-7")
dt = datetime.datetime(2025, 3, 9)
adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
int.fleet.load_output(str(adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")))
int.show(dt.replace(hour=13, minute=5), tlen=10 * 60)

# %%
