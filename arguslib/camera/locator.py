# %%
"""Finds the closest image to the one required"""

import os
import datetime as dtmod
from csat2.locator import FileLocator
from .sources import VideoFrameSource

# The FileLocator format for the default MP4 video files
video_filename_format = "/disk1/Data/ARGUS/{campaign}/{camstr}/videos/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}_{hour:0>2}{min:0>2}{second:0>2}_A.mp4"

cal_filename_format = "/disk1/Data/ARGUS/{campaign}/{camstr}/cal/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}T{hour:0>2}{min:0>2}{second:0>2}_CAL{im_index}.mp4"

os.environ["OPENCV_FFMPEG_LOGLEVEL"] = "0"  # suppress opencv warnings


class CameraData:
    """
    Default data loader for ARGUS cameras.
    This loader finds and reads data from timestamped MP4 files.
    """

    def __init__(self, campaign, camstr):
        from csat2.locator import FileLocator  # keep existing behaviour

        self.campaign = campaign
        self.camstr = camstr
        self.locator = FileLocator()
        self.locator.search_paths["ARGUS"] = {}
        self.locator.search_paths["ARGUS"]["video"] = [video_filename_format]
        self.locator.search_paths["ARGUS"]["cal"] = [cal_filename_format]

        self._source = None
        self._source_key = None  # (type, path)

    def _select_source_for_dt(self, dt: dtmod.datetime):
        """Finds the correct MP4 VideoFrameSource for the given datetime."""
        mp4_path = self.get_video_file(dt)

        if mp4_path is None:
            raise FileNotFoundError(
                f"No camera {self.camstr} video file found for {dt}"
            )

        key = ("mp4", mp4_path)
        if self._source_key != key:
            self._source = VideoFrameSource(mp4_path)
            self._source_key = key
        return self._source

    def get_data_time(self, dt: dtmod.datetime, return_timestamp: bool = False):
        """
        Gets image data from the appropriate MP4 file at the nearest possible timestamp.
        """
        try:
            source = self._select_source_for_dt(dt)
            data, timestamp = source.get_data_time(dt, return_timestamp=True)
        except (FileNotFoundError, ValueError) as e:
            # Re-raise FileNotFoundError, or catch "not in time bounds" ValueError
            if isinstance(e, FileNotFoundError):
                raise e
            raise FileNotFoundError(
                f"No video data for {self.camstr} at {dt}: {e}"
            ) from e

        if abs(timestamp - dt).total_seconds() > 180:
            raise FileNotFoundError(
                f"No image within 3 minute tolerance for {self.camstr} at {dt} (closest at {timestamp}). The nearest timestamp is {timestamp}."
            )

        if not return_timestamp:
            return data
        return data, timestamp

    def get_video_file(self, dt):
        """Finds the local MP4 file path for a given datetime."""
        # dt is in utc.
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
        cap.release()  # Release the file handle
        return False


# %%
