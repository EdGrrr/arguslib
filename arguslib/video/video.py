import numpy as np
import datetime

# Size of the timestamp label
ts_factor = 4


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
    return datetime.datetime.fromtimestamp(int_timestamp)


def extract_exposure(image, ts_factor=ts_factor):
    image_exp_array = (
        image[ts_factor, 0 : (31 * ts_factor) : ts_factor, 0] > 128
    ).astype("int")
    return int("".join(str(a) for a in image_exp_array), 2)


def round_seconds(dt):
    if dt.microsecond >= 500000:
        dt += datetime.timedelta(seconds == 1)
    return dt.replace(microsecond=0)


def round_up_interval(dt, interval):
    dt = dt.replace(microsecond=0) + datetime.timedelta(seconds=1)
    seconds_since_hour = (
        dt - dt.replace(minute=0, second=0, microsecond=0)
    ).total_seconds()
    extra_seconds = interval - seconds_since_hour % interval
    return dt + datetime.timedelta(seconds=extra_seconds)
