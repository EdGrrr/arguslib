import xarray as xr
import datetime
import pyart
from csat2.locator import FileLocator

# TODO: could update these to discriminate the specific radar... Probably a job for when we put the search paths in the config.
vpt_filename_format = "/disk1/Data/{campaign}/radar/L1/{year}{mon:0>2}{day:0>2}/ncas-mobile-ka-band-radar-1_cao_{year}{mon:0>2}{day:0>2}-{hour:0>2}**_vpt_l1_v1.0.0.nc"
rhi_filename_format = "/disk1/Data/{campaign}/radar/L1/{year}{mon:0>2}{day:0>2}/ncas-mobile-ka-band-radar-1_cao_{year}{mon:0>2}{day:0>2}-{hour:0>2}**_rhi_l1_v1.0.0.nc"


def initialise_locator():
    locator = FileLocator()
    locator.search_paths["ARGUS"] = {}
    locator.search_paths["ARGUS"]["vpt"] = [vpt_filename_format]
    locator.search_paths["ARGUS"]["rhi"] = [rhi_filename_format]
    return locator


class RadarData:
    def __init__(self, campaign, scan_type="rhi"):
        self.campaign = campaign
        self.scan_type = scan_type
        self.locator = initialise_locator()

        self.file = None
        self.current_file_path = None

    def get_gridded_data_time(self, dt, var) -> xr.DataArray:
        if self.scan_type == "vpt":
            return self._get_gridded_vpt_data(dt, var)
        elif self.scan_type == "rhi":
            return self._get_gridded_rhi_data(dt, var)
        else:
            raise ValueError(f"Unknown scan type {self.scan_type}")

    def _get_gridded_rhi_data(self, dt, var) -> xr.DataArray:
        radar = self.get_pyart_radar(dt)

        if var not in radar.fields:
            raise ValueError(f"Variable {var} not found in radar file")

        grid = pyart.map.grid_rhi_sweeps(
            radar, roi_func="dist_beam", nb=0.6, bsp=0.3, min_radius=30
        )  # Assumptions made about beamwidth, spacing, and min_radius.
        # spacing an width informed by radar parameters, min_radius guessed to look ok.

        return grid[var]

    def _get_gridded_vpt_data(self, dt, var) -> xr.DataArray:
        raise NotImplementedError("VPT data gridding not yet implemented")

    def get_pyart_radar(self, dt) -> pyart.core.Radar:
        filepath = self.get_filepath(dt)

        if filepath is None:
            raise FileNotFoundError(f"No file found for {dt}")

        if filepath != self.current_file_path:
            self.radar_data = pyart.io.read(filepath)
            self.current_file_path = filepath

        return self.radar_data

    def get_filepath(self, dt):
        if self.scan_type == "vpt":
            hour = "**"
        else:
            hour = dt.hour

        files = self.locator.search(
            "ARGUS",
            self.scan_type,
            campaign=self.campaign,
            year=dt.year,
            mon=dt.month,
            day=dt.day,
            hour=hour,
            min="**",
            second="**",
        )

        if len(files) == 0:
            return None
        elif len(files) == 1:
            return files[0]
        else:
            mins = [int(f.split("_")[-4][-4:-2]) for f in files]
            secs = [int(f.split("_")[-4][-2:]) for f in files]
            tot_secs = [m * 60 + s for m, s in zip(mins, secs)]
            closest = min(tot_secs, key=lambda x: abs(x - dt.minute * 60 - dt.second))
            return files[tot_secs.index(closest)]

    def get_next_time(self, dt, max_gap_hrs=48):
        """
        Get the next available time for the radar data.
        """
        hour = "**"

        files = self.locator.search(
            "ARGUS",
            self.scan_type,
            campaign=self.campaign,
            year=dt.year,
            mon=dt.month,
            day=dt.day,
            hour=hour,
            min="**",
            second="**",
        )
        files = sorted(files)

        dates = []
        for f in files:
            dates.append(datetime.datetime.strptime(f.split("_")[-4], "%Y%m%d-%H%M%S"))

        next_available = None
        for i, date in enumerate(dates):
            if date > dt:  # dont want to include the current file
                # check if the gap is less than max_gap_hrs
                if (date - dt).total_seconds() / 3600 > max_gap_hrs:
                    return None
                next_available = date
                return next_available

        if next_available is None:
            # no files today.
            next_day = dt + datetime.timedelta(days=1)
            next_day = next_day.replace(hour=0, minute=0, second=0)
            # deplete the max_gap_hrs
            max_gap_hrs -= (next_day - dt).total_seconds() / 3600
            if max_gap_hrs < 0:
                return None

            return self.get_next_time(
                datetime.datetime(dt.year, dt.month, dt.day, dt.hour + 1), max_gap_hrs
            )
