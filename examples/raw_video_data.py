# %%
import cv2
import matplotlib.pyplot as plt
from arguslib.camera import Video
from arguslib.camera.video import extract_timestamp, extract_exposure


vid = Video(
    "/disk1/Data/ARGUS/COBALT/3-8/videos/2025-03-09/argus-3-8_20250309_090030_A.mp4"
)

image = vid.get_frame(5)

timestamp = extract_timestamp(image)
exposure = extract_exposure(image)

print(f"Image timestamp: {timestamp}")
print(f"Image exposure: {exposure}us")

fig, ax = plt.subplots()
ax.imshow(image[::-1, :, ::-1])
ax.set_title(f"Timestamp: {timestamp}, Exposure: {exposure}us")
fig.show()
# %%
