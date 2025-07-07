
# %%

from arguslib.misc.met import download_era5_winds
import cv2
from tqdm import trange
from arguslib.aircraft import AircraftInterface
from pathlib import Path

from arguslib.instruments.direct_camera import DirectUndistortedCamera
import datetime
# %%
cam = DirectUndistortedCamera.from_config("COBALT", "3-7")
cam.show(datetime.datetime(2025,5,11,11,57,57))
cam.image


aci = AircraftInterface(cam)
# Dates with ERA5 data:
# dt_start = datetime.datetime(2025,3,8,10)
dt_start = datetime.datetime(2025,5,11,7,30)
# dt_start = datetime.datetime(2025,4,1,11,00)
# dt_start = datetime.datetime(2025,3,29,7,00)



adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
outdir = Path(__file__).parent / "output" / "videos"
outdir.mkdir(exist_ok=True, parents=True)

adsb_file = adsb_datadir / (dt_start.strftime("%Y%m%d") + "_ADS-B")
if aci.fleet.loaded_file != str(adsb_file):
    aci.fleet.load_output(str(adsb_file))

# %%
aci.fleet.assign_era5_winds()

# %%
out = cv2.VideoWriter(outdir / f"{dt_start.isoformat(timespec='minutes')}.mp4", cv2.VideoWriter_fourcc(*'avc1'), 4, (3040,3040))
for i in trange(90):
    aci.show(dt_start+datetime.timedelta(minutes=i), tlen=15*60)
    out.write(aci.camera.to_image_array()[:,:,::-1])
out.release()


# %%
aci.show(dt_start)
aci.camera.image

# %%
# download_era5_winds(dt_start)
# %%
