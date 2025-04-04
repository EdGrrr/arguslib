import numpy as np
import datetime
import cv2
from zoneinfo import ZoneInfo

from pytz import utc

# Size of the timestamp label
ts_factor = 4


class Video:
    def __init__(self, filepath):
        self.filepath = filepath
        self.cap = cv2.VideoCapture(filepath)
        self.n_frames = int(self.cap.get(cv2.CAP_PROP_FRAME_COUNT))
        self.time_bounds = self.get_timestamps([0, self.n_frames - 1])

    def get_data_time(self, dt, return_timestamp=False):
        n = self.estimate_frame_number(dt)
        n = np.round(n).astype(int)
        frame = self.get_frame(n)
        if return_timestamp:
            return frame, extract_timestamp(frame)
        return frame

    def get_frame(self, n):
        if n < 0:
            n = self.n_frames + n
        self.cap.set(cv2.CAP_PROP_POS_FRAMES, n)
        ret, frame = self.cap.read()
        if not ret:
            raise ValueError(f"Could not read frame {n} from {self.filepath}")

        return frame

    def get_timestamps(self, ns=None):
        timestamps = []

        if ns is None:
            ns = range(self.n_frames)

        for i in ns:
            frame = self.get_frame(i)
            timestamps.append(extract_timestamp(frame))
        return tuple(timestamps)

    def estimate_frame_number(self, dt):
        if dt < self.time_bounds[0] or dt > self.time_bounds[1]:
            raise ValueError(
                f"Timestamp {dt} is not in the video time bounds for the video {self.filepath}"
            )
        return np.interp(
            dt.timestamp(),
            [self.time_bounds[0].timestamp(), self.time_bounds[1].timestamp()],
            [0, self.n_frames - 1],
        )

    def __len__(self):
        return self.n_frames

    def __repr__(self):
        return f"Video({self.filepath}, {self.n_frames} frames)"


def timestamp_image(tstamp, data, ts_factor=ts_factor, exposure=None):
    stamp = int(tstamp)
    ts_array = np.array(list(np.binary_repr(stamp, 31))).astype("int")
    ts_array = ts_array[None, :].repeat(ts_factor, axis=1).repeat(ts_factor, axis=0)
    data[:ts_factor, : (31 * ts_factor), :] = 255 * ts_array[:, :, None]
    if exposure is not None:
        stamp = int(exposure)
        ts_array = np.array(list(np.binary_repr(stamp, 31))).astype("int")
        ts_array = ts_array[None, :].repeat(ts_factor, axis=1).repeat(ts_factor, axis=0)
        data[ts_factor : (2 * ts_factor), : (31 * ts_factor), :] = (
            255 * ts_array[:, :, None]
        )
    return data


def create_timestamp():
    pass


def extract_timestamp(image, ts_factor=ts_factor):
    """Return a datetime object with the image timestamp to the nearest second."""
    image_ts_array = (image[0, 0 : (31 * ts_factor) : ts_factor, 0] > 128).astype("int")
    int_timestamp = int("".join(str(a) for a in image_ts_array), 2)
    local_time = datetime.datetime.fromtimestamp(
        int_timestamp, tz=ZoneInfo("Europe/London")
    )
    utc_time = local_time.astimezone(utc)
    return utc_time.replace(tzinfo=None)  # Make "timezone naive" to be compatible.


def extract_exposure(image, ts_factor=ts_factor):
    image_exp_array = (
        image[ts_factor, 0 : (31 * ts_factor) : ts_factor, 0] > 128
    ).astype("int")
    return int("".join(str(a) for a in image_exp_array), 2)


def round_seconds(dt):
    if dt.microsecond >= 500000:
        dt += datetime.timedelta(seconds=1)
    return dt.replace(microsecond=0)


def round_up_interval(dt, interval):
    dt = dt.replace(microsecond=0) + datetime.timedelta(seconds=1)
    seconds_since_hour = (
        dt - dt.replace(minute=0, second=0, microsecond=0)
    ).total_seconds()
    extra_seconds = interval - seconds_since_hour % interval
    return dt + datetime.timedelta(seconds=extra_seconds)
