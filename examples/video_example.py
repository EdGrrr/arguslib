import cv2
import numpy as np
import datetime

ts_factor = 4

def extract_timestamp(image, ts_factor=ts_factor):
    """Return a datetime object with the image timestamp to the nearest second."""
    image_ts_array = (image[0, 0 : (31 * ts_factor) : ts_factor, 0] > 128).astype("int")
    int_timestamp = int("".join(str(a) for a in image_ts_array), 2)
    return datetime.datetime.fromtimestamp(int_timestamp)


def extract_exposure(image, ts_factor=ts_factor):
    image_exp_array = (
        image[ts_factor, 0 : (31 * ts_factor) : ts_factor, 0] > 128
    ).astype("int")
    return int("".join(str(a) for a in image_exp_array), 2)


vidcap = cv2.VideoCapture(
    "/disk1/Data/ARGUS/COBALT/5-1/videos/2025-01-01/argus-5-1_20250101_090030_A.mp4"
)

success, image = vidcap.read()

timestamp = extract_timestamp(image)
exposure = extract_exposure(image)

print(f"Image timestamp: {timestamp}")
print(f"Image exposure: {exposure}us")
