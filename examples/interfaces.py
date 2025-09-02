"""This script contains increasingly complex examples of how to use the
interface classes to plot radar and aircraft data."""

# %%
import datetime
from arguslib.aircraft import AircraftInterface
from arguslib.radar.radar import Radar
from arguslib.camera.camera_array import CameraArray
from arguslib.camera.undistorted_camera import UndistortedCamera
from arguslib.radar.radar_interface import RadarInterface


# %% Plotting the aircraft tracks on a single camera
cai = AircraftInterface.from_campaign("COBALT", "3-7")
dt = datetime.datetime(2025, 5, 1, 7, 40)
cai.load_flight_data(dt)
cai.show(dt, tlen=10 * 60)

# %% Plotting a single camera with simultaneous radar data
dt = datetime.datetime(2025, 3, 25, 9, 50, 49)
ri_on_camera = RadarInterface.from_campaign("COBALT", "3-7")
ri_on_camera.show(dt)


# %% Plotting multiple cameras simultaneously with radar data
dt = datetime.datetime(2025, 3, 25, 9, 50, 49)

# "Undistorted" cameras straighten out any distortion from the images.
# Standard cameras also exist: `arguslib.camera.camera.Camera`.
# Camera arrays can also be loaded from config: CameraArray.from_config("COBALT"), but then contain standard (distorted cameras).
multicam = CameraArray(
    [
        UndistortedCamera.from_config("COBALT", cam_id)
        for cam_id in ["5-1", "5-2", "5-3", "5-4", "3-7"]
    ],
    (3, 3),
)

radar = Radar.from_config("COBALT")
ri_on_multicam = RadarInterface(radar, multicam)
ri_on_multicam.show(dt)


# %% Plot multiple cameras with flight tracks
ai_on_multicam = AircraftInterface(multicam)
dt = datetime.datetime(2025, 5, 1, 7, 40)
ai_on_multicam.load_flight_data(dt)
ai_on_multicam.show(dt, tlen=10 * 60, color_icao=True)


# %% Put a radar interface on top of the aircraft interface to show the radar beams and radar data next to it
dt = datetime.datetime(2025, 5, 1, 7, 25, 6)
ri_on_ai_on_multicam = RadarInterface(radar, ai_on_multicam)
ri_on_ai_on_multicam.show(dt, kwargs_camera=dict(tlen=30 * 60, color_icao=True))

# %% Or put the aircraft interface on top of the radar to show the tracks on the radar too
dt = datetime.datetime(2025, 5, 1, 6, 6, 18)
cri = RadarInterface(radar, multicam)
ai_on_ri_on_multicam = AircraftInterface(cri)
ai_on_ri_on_multicam.load_flight_data(dt)
ai_on_ri_on_multicam.show(
    dt,
    tlen=30 * 60,
    color_icao=True,
    trail_kwargs=dict(plot_kwargs=dict(plotting_method="intersect_plot")),
)
