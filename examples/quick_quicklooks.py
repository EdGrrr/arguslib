# %%

# %%
import datetime
from pathlib import Path

from arguslib import AircraftInterface, RadarInterface, Radar, UndistortedCamera

adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
outdir = Path(__file__).parent / "output" / "camera_radar_quicklooks"
outdir.mkdir(exist_ok=True, parents=True)

radar = Radar.from_config("COBALT")
cam = UndistortedCamera.from_config("COBALT", "3-7")
# cam = CameraArray.from_config("COBALTArray")
cri = RadarInterface(radar, cam)
cai = AircraftInterface(cri)

dt = datetime.datetime(2025, 5, 11, 11, 57, 41)

cai.load_flight_data(dt)

# %%
ax_cam, ax_radar = cai.show(
    datetime.datetime(2025, 5, 11, 6, 9),
    # dt,
    tlen=15 * 60,
    color_icao=True,
    trail_kwargs={
        "label_acft": True,
        # "icao_include": ["06a11e", "3949e4"],
        # "icao_include": ["4cadfe"],
        # "icao_include":[
        # "407f7b",
        # ],
        "plot_kwargs": {
            "alpha": 0.5,
        },
        # "advection_winds":'aircraft'
    },
    kwargs_camera={"brightness_adjust": 2.0},
)

ax_radar.legend(ncol=4, fontsize="small")

# %%
