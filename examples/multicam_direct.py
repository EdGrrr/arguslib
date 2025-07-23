# %%
import datetime as dt


from arguslib.camera.camera_array import CameraArray
from arguslib.camera.undistorted_camera import UndistortedCamera
from arguslib.instruments.instruments import Position


multicam2 = CameraArray(
    [
        UndistortedCamera.from_config("COBALT", "3-7"), 
        UndistortedCamera.from_config("COBALT", "5-1"),
        UndistortedCamera.from_config("COBALT", "5-2"),
        UndistortedCamera.from_config("COBALT", "5-3"),
        UndistortedCamera.from_config("COBALT", "5-4"),
    ],
    (3,3)
)


# multicam.show(dt.datetime(2025,3,12,11,50))
ax = multicam2.show(dt.datetime(2025,3,12,11,50))

# construct a 10km grid between the min/max lon/lat.
max_lon = multicam2.cameras[2].position.lon
min_lon = multicam2.cameras[4].position.lon
max_lat = multicam2.cameras[3].position.lat
min_lat = multicam2.cameras[1].position.lat

# %%
from arguslib.radar.radar import Radar
import datetime as dt
radar = Radar.from_config("COBALT")
radar.show(dt.datetime(2025, 5, 31, 15,26,7), var='DBZ')
# %%
from arguslib.radar.radar_interface import RadarInterface
cri = RadarInterface(radar, multicam2)
cri.show(dt.datetime(2025, 5, 31, 15,26,7))
# %%
