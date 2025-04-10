# %%

from datetime import datetime
from arguslib.instruments.calibration import PerspectiveProjection, focal_length_px
from arguslib.instruments.instruments import Position
from arguslib.instruments.camera import Camera
from arguslib.video.locator import CameraData
import matplotlib.pyplot as plt

cam = Camera.from_config("COBALT", "2-11")

cam.show(datetime(2025, 4, 8, 9, 50))
# %%

fx = 2.75  # focal length in mm
fy = 2.75  # focal length in mm

sensor_size_mm = 7.4  # diag sensor size in mm

image_size_px = (4608, 2304)  # image size in pixels (width, height)

fx_px = focal_length_px(fx, image_size_px[0], sensor_size_mm)
fy_px = focal_length_px(fy, image_size_px[1], sensor_size_mm)


cx = 2304  # x coordinate of the principal point in pixels
cy = 1296  # y coordinate of the principal point in pixels


proj = PerspectiveProjection([fx_px, fy_px], [cx, cy], None)


cam = Camera(proj, Position(-0.1791071, 51.499189, 0.1), rotation=[30, 250, 190])

cam.data_loader = CameraData("COBALT", "2-11")

fig, ax = plt.subplots()

ax = cam.show(datetime(2025, 4, 8, 9, 50), ax=ax, lr_flip=False)

# %%
fx_px, fy_px
# %%
