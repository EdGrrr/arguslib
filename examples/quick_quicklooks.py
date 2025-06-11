# %%

# %%
import datetime
from pathlib import Path

from arguslib.instruments.camera import Camera

from arguslib.aircraft import AircraftInterface
from arguslib.misc.geo import ft_to_km
from arguslib.radar import RadarInterface
from arguslib.instruments.radar import Radar
from arguslib.instruments.camera_array import CameraArray

from arguslib.instruments.undistorted_camera import UndistortedCamera

adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
outdir = Path(__file__).parent / "output" / "camera_radar_quicklooks"
outdir.mkdir(exist_ok=True, parents=True)

radar = Radar.from_config("COBALT")
cam = UndistortedCamera.from_config("COBALT", "3-7")
# cam = CameraArray.from_config("COBALTArray")
cri = RadarInterface(radar, cam)
cai = AircraftInterface(cri) 

dt = datetime.datetime(2025, 5, 11, 11, 46, 57)

adsb_file = adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")
if cai.fleet.loaded_file != str(adsb_file):
    cai.fleet.load_output(str(adsb_file))


# %% replace cai.show with cri.show to just plot the radar and camera
cai.fleet.assign_era5_winds()
# %%
ax_cam, ax_radar = cai.show(
    dt,
    tlen=2*60 * 60,
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
            "radar_kwargs": {"plotting_method": "intersect_plot"},
        },
        # "advection_winds":'aircraft'
    },
    # kwargs_camera={"brightness_adjust": 2.0},
)

ax_radar.legend(ncol=4, fontsize="small")

# %%
from arguslib.misc.geo import ft_to_km
ft_to_km(cai.fleet.get_tracks(dt)['407f7b'][-1])
# %%
from arguslib.instruments import Position
radar.position.target_ead(Position(radar.position.lon, radar.position.lat, 10))
# %%
