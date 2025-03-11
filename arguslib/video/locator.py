# %%
"""Finds the closest image to the one required"""

from csat2.locator import FileLocator
import os

from .video import Video

video_filename_format = "/disk1/Data/ARGUS/{campaign}/{camstr}/videos/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}_{hour:0>2}{min:0>2}{second:0>2}_A.mp4"

cal_filename_format = "/disk1/Data/ARGUS/{campaign}/{camstr}/cal/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}T{hour:0>2}{min:0>2}{second:0>2}_CAL{im_index}.mp4"


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

    def get_data_time(self, dt):
        filepath = self.get_video_file(dt)

        if filepath is None:
            raise FileNotFoundError(f"No video file found for {dt}")

        if filepath != self.current_video_path:
            self.video = Video(filepath)
            self.current_video_path = filepath

        try:
            return self.video.get_data_time(dt)
        except ValueError as e:
            # Likely a "not in time bounds" error - happens when we are "between two frames"
            # if we are at the start, return the first frame if its within a minute
            if (
                dt < self.video.time_bounds[0]
                and (self.video.time_bounds[0] - dt).seconds < 60
            ):
                return self.video.get_data_time(self.video.time_bounds[0])
            # if we are at the end, return the last frame if its within a minute
            elif (
                dt > self.video.time_bounds[1]
                and (dt - self.video.time_bounds[1]).seconds < 60
            ):
                return self.video.get_data_time(self.video.time_bounds[1])
            else:
                raise e

    def get_video_file(self, dt):
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

        if len(files) == 0:
            return None
        elif len(files) == 1:
            return files[0]
        else:
            mins = [int(f.split("_")[-2][2:4]) for f in files]
            secs = [int(f.split("_")[-2][4:6]) for f in files]
            tot_secs = [m * 60 + s for m, s in zip(mins, secs)]
            closest = min(tot_secs, key=lambda x: abs(x - dt.minute * 60 - dt.second))
            return files[tot_secs.index(closest)]


# %%
