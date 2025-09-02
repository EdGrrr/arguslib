# %%
"""Finds the closest image to the one required"""

from arguslib.misc.times import convert_to_london_naive
from csat2.locator import FileLocator
import os

from .video import Video
from zoneinfo import ZoneInfo

from pytz import utc

video_filename_format = "/disk1/Data/ARGUS/{campaign}/{camstr}/videos/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}_{hour:0>2}{min:0>2}{second:0>2}_A.mp4"

cal_filename_format = "/disk1/Data/ARGUS/{campaign}/{camstr}/cal/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}T{hour:0>2}{min:0>2}{second:0>2}_CAL{im_index}.mp4"

os.environ["OPENCV_FFMPEG_LOGLEVEL"] = "0"  # suppress opencv warnings


def initialise_locator():
    locator = FileLocator()
    locator.search_paths["ARGUS"] = {}
    locator.search_paths["ARGUS"]["video"] = [video_filename_format]
    locator.search_paths["ARGUS"]["cal"] = [cal_filename_format]
    return locator


class CameraData:
    def __init__(self, campaign, camstr):
        self.campaign = campaign
        self.camstr = camstr
        self.locator = initialise_locator()

        self.video = None
        self.current_video_path = None

        self.image = None
        self.current_image_time = None

    def get_data_time(self, dt, return_timestamp=False):
        filepath = self.get_video_file(dt)

        if filepath is None:
            raise FileNotFoundError(
                f"No camera {self.camstr} video file found for {dt}"
            )

        if filepath != self.current_video_path:
            self.video = Video(filepath)
            self.current_video_path = filepath

        try:
            self.image, self.current_image_time = self.video.get_data_time(
                dt, return_timestamp=True
            )
        except ValueError as e:
            # Likely a "not in time bounds" error - happens when we are "between two frames"
            # if we are at the start, return the first frame if its within a minute
            if (
                dt < self.video.time_bounds[0]
                and (self.video.time_bounds[0] - dt).seconds < 60
            ):
                self.image, self.current_image_time = self.video.get_data_time(
                    self.video.time_bounds[0], return_timestamp=True
                )
            # if we are at the end, return the last frame if its within a minute
            elif (
                dt > self.video.time_bounds[1]
                and (dt - self.video.time_bounds[1]).seconds < 60
            ):
                self.image, self.current_image_time = self.video.get_data_time(
                    self.video.time_bounds[1], return_timestamp=return_timestamp
                )
            else:
                raise e

        if return_timestamp:
            return self.image, self.current_image_time
        else:
            return self.image

    def get_video_file(self, dt):
        # dt is in utc.
        # but (TEST) the files are **named** with local time?
        # get dt object which is timezone naive but
        # dt = dt.replace(hour=dt.hour-1)

        files = self.locator.search(
            "ARGUS",
            "video",
            campaign=self.campaign,
            camstr=self.camstr,
            year=dt.year,
            mon=dt.month,
            day=dt.day,
            hour=dt.hour,
            min="**",
            second="**",
        )

        # sometimes there are completely empty files – annoying. these are abou 258bytes
        files = [f for f in files if os.path.getsize(f) > 1000]

        # there are some other corrupt files... will need to check if these are broken by trying to load them
        mins = [int(f.split("_")[-2][2:4]) for f in files]
        secs = [int(f.split("_")[-2][4:6]) for f in files]

        tot_secs = [
            (dt.minute * 60 + dt.second) - (m * 60 + s) for m, s in zip(mins, secs)
        ]

        files = [f for f, delta_t in zip(files, tot_secs) if delta_t >= -60]
        tot_secs = [delta_t for delta_t in tot_secs if delta_t >= -60]

        if len(files) == 0:
            return None
        elif len(files) == 1:
            if is_mp4_file_corrupt(files[0]):
                return None
            return files[0]
        else:
            closest = min(tot_secs)
            fpath = files[tot_secs.index(closest)]
            if is_mp4_file_corrupt(fpath):
                return None
            return fpath


def is_mp4_file_corrupt(filepath):
    """
    Check if the mp4 file is corrupt by trying to open it with cv2.
    If it fails, return True, otherwise return False.
    """
    import cv2

    cap = cv2.VideoCapture(filepath)
    if not cap.isOpened():
        return True
    else:
        return False


# %%
