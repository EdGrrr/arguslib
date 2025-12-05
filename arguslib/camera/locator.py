# %%
"""Finds the closest image to the one required"""

import os
import datetime as dtmod
from .sources import VideoFrameSource
from .video import get_video_time_bounds

# The FileLocator format for the default MP4 video files
video_filename_format = "/disk1/Data/ARGUS/{campaign}/{camstr}/videos/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}_{hour:0>2}{min:0>2}{second:0>2}_A.mp4"

cal_filename_format = "/disk1/Data/ARGUS/{campaign}/{camstr}/cal/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}T{hour:0>2}{min:0>2}{second:0>2}_CAL{im_index}.mp4"

os.environ["OPENCV_FFMPEG_LOGLEVEL"] = "0"  # suppress opencv warnings


class CameraData:
    """
    Default data loader for ARGUS cameras.
    This loader finds and reads data from timestamped MP4 files.
    """

    def __init__(self, campaign, camstr, invert_axes=[False, False], timestamp_timezone="UTC"):
        from csat2.locator import FileLocator  # keep existing behaviour

        self.campaign = campaign
        self.camstr = camstr
        self.locator = FileLocator()
        self.locator.search_paths["ARGUS"] = {}
        self.locator.search_paths["ARGUS"]["video"] = [video_filename_format]
        self.locator.search_paths["ARGUS"]["cal"] = [cal_filename_format]
        self.timestamp_timezone = timestamp_timezone

        self._source = None
        self._source_key = None  # (type, path)
        self._invert_axes = invert_axes
        
        self.image = None  # placeholder for the current image
        self.current_image_time = None  # placeholder for the current image timestamp

    def _select_source_for_dt(self, dt: dtmod.datetime):
        """Finds the correct MP4 VideoFrameSource for the given datetime."""
        mp4_path = self.get_video_file(dt)

        if mp4_path is None:
            raise FileNotFoundError(
                f"No camera {self.camstr} video file found for {dt}"
            )

        key = ("mp4", mp4_path)
        if self._source_key != key:
            self._source = VideoFrameSource(mp4_path, timestamp_timezone=self.timestamp_timezone)
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
            # Likely caused by being "between" video files
            if isinstance(e, FileNotFoundError):
                raise e
            raise FileNotFoundError(
                f"No video data for {self.camstr} at {dt}: {e}"
            ) from e

        if abs(timestamp - dt).total_seconds() > 180:
            raise FileNotFoundError(
                f"No image within 3 minute tolerance for {self.camstr} at {dt} (closest at {timestamp}). The nearest timestamp is {timestamp}."
            )
        
        if self._invert_axes[0]:
            data = data[::-1, :]
        if self._invert_axes[1]:
            data = data[:, ::-1]
            
        self.image = data
        self.current_image_time = timestamp

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

        # only consider files whose nominal time is before the requested time.
        # Allow a 3 minute tolerance, in case of missing images.
        files = [f for f, delta_t in zip(files, tot_secs) if delta_t >= -180]
        tot_secs = [delta_t for delta_t in tot_secs if delta_t >= -180]

        if len(files) == 0:
            return None
        elif len(files) == 1:
            if is_mp4_file_corrupt(files[0]):
                return None
            return files[0]
        else:
            # this nees to be a little bit smarter.
            # We want the frame closest to the requested time.
            # candidate files are the two closest before the requested time.
            closest = min(tot_secs)
            
            # check the time bounds of this file
            dt_start, dt_end = get_video_time_bounds(files[tot_secs.index(closest)], timestamp_timezone=self.timestamp_timezone)
            if dt_start <= dt <= dt_end:
                fpath = files[tot_secs.index(closest)]
                if is_mp4_file_corrupt(fpath):
                    return None
                return fpath
            else:
                next_closest = min(
                    [ts for ts in tot_secs if ts != closest], default=None
                )
                dt_start_prev, dt_end_prev = get_video_time_bounds(
                    files[tot_secs.index(next_closest)], timestamp_timezone=self.timestamp_timezone
                )
                if next_closest is not None and dt_start_prev <= dt <= dt_end_prev:
                    fpath = files[tot_secs.index(next_closest)]
                    if is_mp4_file_corrupt(fpath):
                        return None
                    return fpath
                elif dt_end_prev < dt < dt_start:
                    # in between two files, return the closest
                    if (dt - dt_end_prev < dt_start - dt) and not is_mp4_file_corrupt(
                        files[tot_secs.index(next_closest)]
                    ):
                        fpath = files[tot_secs.index(next_closest)]
                        return fpath
                    else:
                        fpath = files[tot_secs.index(closest)]
                        if is_mp4_file_corrupt(fpath):
                            return None
                        return fpath
                else:
                    # this shouldn't happen?
                    raise RuntimeError("Logic error in get_video_file")
    

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
