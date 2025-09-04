"""
Correcting Lens Distortion
==========================

This example demonstrates how to use `UndistortedCamera` to correct for the
barrel or pincushion distortion often present in wide-angle camera lenses.
"""

# %%
import datetime
from arguslib.camera.camera import Camera
from arguslib.camera.undistorted_camera import UndistortedCamera

dt = datetime.datetime(2025, 3, 25, 9, 0, 0)

# %%
# Standard Camera with Lens Distortion
# ------------------------------------
#
# When we plot an image from a standard `Camera` object, especially one with a
# wide-angle lens, we can often see the effects of lens distortion. Notice how
# straight lines appear to curve near the edges of the image.
cam = Camera.from_config("COBALT", "3-7")
ax = cam.show(dt)

# %%
# Correcting Distortion with `UndistortedCamera`
# ----------------------------------------------
#
# The `UndistortedCamera` class applies a correction based on the camera's
# calibration parameters, effectively "unwarping" the image. The resulting
# plot shows a much more natural, rectilinear view, as if seen through a
# distortion-free lens.
undistorted_cam = UndistortedCamera.from_config("COBALT", "3-7")
ax = undistorted_cam.show(dt)
