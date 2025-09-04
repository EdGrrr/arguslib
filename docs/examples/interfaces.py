"""
Combining Interfaces
====================

This example demonstrates the power of `arguslib`'s interface system by
stacking multiple interfaces together. We will create a complex visualization
that includes:

1. A `CameraArray` to show multiple camera views simultaneously.
2. A `RadarInterface` to overlay radar data on top of the camera array.
3. An `AircraftInterface` to plot flight tracks over both the cameras and radar data.

This showcases how different data sources can be combined into a single,
comprehensive plot.
"""

# %%
import datetime
from arguslib.aircraft import AircraftInterface
from arguslib.radar.radar import Radar
from arguslib.camera.camera_array import CameraArray
from arguslib.camera.undistorted_camera import UndistortedCamera
from arguslib.radar.radar_interface import RadarInterface

# %%
# Stacking Interfaces for a Composite View
# ----------------------------------------
#
# We can combine multiple cameras with flight tracks and radar data using nested
# interfaces. Here, we build up the final visualization step-by-step.
#
# First, we define the core instruments: a `Radar` and a `CameraArray` composed
# of `UndistortedCamera` instances to provide a wide, distortion-free field of view.
dt = datetime.datetime(2025, 5, 1, 7, 30, 23)
radar = Radar.from_config("COBALT")
multicam = CameraArray.from_config("COBALTArray", camera_class=UndistortedCamera)
#
# Next, we wrap the `CameraArray` with a `RadarInterface` to overlay the radar
# data.
cri = RadarInterface(radar, multicam)
#
# %%
# Finally, we wrap the `RadarInterface` with an `AircraftInterface`. This allows
# us to load flight data and plot aircraft tracks on top of the combined
# radar-on-camera view.
#
# The `show()` method on the outermost interface renders the entire scene,
# showing the aircraft tracks and their intersection with the radar scan over
# the multi-camera background.
ai_on_ri_on_multicam = AircraftInterface(cri)
ai_on_ri_on_multicam.load_flight_data(dt)
ax = ai_on_ri_on_multicam.show(
    dt,
    tlen=45 * 60,
    color_icao=True,
)

# %%
