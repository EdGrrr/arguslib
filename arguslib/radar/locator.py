import xarray as xr
import datetime
import pyart
from csat2.locator import FileLocator
import os
import re  # Added for robust parsing

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

        self.radar_data: pyart.core.Radar | None = None  # Type hint for clarity
        self.current_file_path: str | None = None

    def get_gridded_data_time(self, dt, var) -> xr.DataArray:
        if self.scan_type == "vpt":
            return self._get_gridded_vpt_data(dt, var)
        elif self.scan_type == "rhi":
            return self._get_gridded_rhi_data(dt, var)
        else:
            raise ValueError(f"Unknown scan type {self.scan_type}")

    def _get_gridded_rhi_data(self, dt, var) -> xr.DataArray:
        radar = self.get_pyart_radar(dt)  # This will now use the smarter loader

        if var not in radar.fields:
            raise ValueError(
                f"Variable {var} not found in radar file: {self.current_file_path}"
            )

        grid = pyart.map.grid_rhi_sweeps(
            radar, roi_func="dist_beam", nb=0.6, bsp=0.3, min_radius=30
        )  # Assumptions made about beamwidth, spacing, and min_radius.
        # spacing an width informed by radar parameters, min_radius guessed to look ok.

        return grid[var]

    def _get_gridded_vpt_data(self, dt, var) -> xr.DataArray:
        raise NotImplementedError("VPT data gridding not yet implemented")

    def get_pyart_radar(self, dt) -> pyart.core.Radar:
        # Check cache first: if currently loaded radar data covers dt
        if self.radar_data and self.current_file_path:
            try:
                if (
                    hasattr(self.radar_data, "time")
                    and "data" in self.radar_data.time
                    and self.radar_data.time["data"].size > 0
                ):

                    scan_start_time = pyart.util.datetime_from_radar(self.radar_data)
                    # Ensure time data is not empty for duration calculation
                    # Duration from the start of the volume to the time of the last ray recorded.
                    scan_duration_seconds = (
                        self.radar_data.time["data"][-1]
                        - self.radar_data.time["data"][0]
                    )
                    scan_end_time = scan_start_time + datetime.timedelta(
                        seconds=scan_duration_seconds
                    )

                    if scan_start_time <= dt <= scan_end_time:
                        return self.radar_data
                else:
                    # Cached data is invalid or incomplete, invalidate
                    self.radar_data = None
                    self.current_file_path = None
            except Exception:  # Broad exception if there's any issue with cached data
                self.radar_data = None  # Invalidate cache
                self.current_file_path = None

        # If cache miss or invalid, find and load the appropriate file
        filepath = self.get_filepath(dt)

        if filepath is None:
            raise FileNotFoundError(
                f"No radar file found containing the time {dt} for campaign {self.campaign}, scan_type {self.scan_type}"
            )

        if filepath != self.current_file_path or self.radar_data is None:
            self.radar_data = pyart.io.read(filepath)
            self.current_file_path = filepath

        # Final check to ensure the loaded data actually contains dt
        # This is a safeguard, as get_filepath should have ensured this.
        if self.radar_data:  # Ensure radar_data is loaded
            scan_start_time = pyart.util.datetime_from_radar(self.radar_data)
            if (
                hasattr(self.radar_data, "time")
                and "data" in self.radar_data.time
                and self.radar_data.time["data"].size > 0
            ):
                scan_duration_seconds = (
                    self.radar_data.time["data"][-1] - self.radar_data.time["data"][0]
                )
                scan_end_time = scan_start_time + datetime.timedelta(
                    seconds=scan_duration_seconds
                )
                # if not (scan_start_time <= dt <= scan_end_time): #no longer the case, because we're reverting to exactly matching the filename (rounding the start time)
                #     # This case should ideally not be reached if get_filepath is correct.
                #     raise FileNotFoundError(
                #         f"Loaded file {filepath} (scan: {scan_start_time} to {scan_end_time}) "
                #         f"does not actually contain requested dt {dt}. This indicates an issue in logic."
                #     )
            else:  # Should not happen if pyart.io.read was successful and file is valid
                raise ValueError(
                    f"Loaded radar file {filepath} has no valid time data."
                )
        else:  # Should not happen if filepath was found
            raise ValueError("Radar data is None after attempting to load file.")

        return self.radar_data

    def get_filepath(self, dt: datetime.datetime) -> str | None:
        """
        Finds the filepath of a radar scan that contains the given datetime 'dt'.
        It searches all files for the day of 'dt', then checks their actual
        scan start and end times by reading the file metadata.
        Returns the path to the first suitable file found (sorted by filename).
        """
        all_day_files = self.locator.search(
            "ARGUS",
            self.scan_type,
            campaign=self.campaign,
            year=dt.year,
            mon=dt.month,
            day=dt.day,
            hour="**",  # Search all hours for the day
            min="**",
            second="**",
        )

        if not all_day_files:
            return None

        # Sort files by name, which should correspond to start times.
        # This helps in picking the earliest scan if multiple cover 'dt'.
        for f_path in sorted(all_day_files):
            try:
                # Optional: Heuristic pre-check based on filename to speed up.
                # Assumes filename contains start time and a max scan duration.
                filename_base = os.path.basename(f_path)
                match = re.search(r"(\d{4}\d{2}\d{2}-\d{2}\d{2}\d{2})", filename_base)
                if match:
                    file_start_dt_from_name = datetime.datetime.strptime(
                        match.group(1), "%Y%m%d-%H%M%S"
                    )
                    if file_start_dt_from_name.replace(microsecond=0) == dt.replace(
                        microsecond=0
                    ):
                        return f_path

                #     # If file starts significantly after dt, skip.
                #     if file_start_dt_from_name > dt + datetime.timedelta(minutes=5): # Grace for dt slightly before actual scan
                #         continue
                #     # If file (estimated) ends significantly before dt, skip.
                #     max_estimated_scan_duration = datetime.timedelta(minutes=15) # Generous estimate
                #     if file_start_dt_from_name + max_estimated_scan_duration < dt:
                #         continue

                # # Read the radar file to get its precise time coverage
                # radar_obj = pyart.io.read(f_path)

                # actual_scan_start_time = pyart.util.datetime_from_radar(radar_obj)

                # if not (hasattr(radar_obj, 'time') and 'data' in radar_obj.time and radar_obj.time['data'].size > 0):
                #     continue # Skip if no valid time data

                # scan_duration_seconds = radar_obj.time['data'][-1] - radar_obj.time['data'][0]
                # actual_scan_end_time = actual_scan_start_time + datetime.timedelta(seconds=scan_duration_seconds)

                # if actual_scan_start_time <= dt <= actual_scan_end_time:
                #     return f_path # Found a file whose scan duration covers dt

            except FileNotFoundError:
                # Log or handle if a file listed by locator is not found
                continue
            except Exception:
                # Log or handle errors during file reading/processing (e.g., corrupted file)
                continue

        return None  # No suitable file found

    def get_next_time(self, dt, max_gap_hrs=48):
        """
        Get the next available time for the radar data.
        This method might need review if its definition of "next time"
        should also consider scans crossing hour boundaries more explicitly,
        but for now, it finds the next *file start time*.
        """
        # Search all files for the day
        files = self.locator.search(
            "ARGUS",
            self.scan_type,
            campaign=self.campaign,
            year=dt.year,
            mon=dt.month,
            day=dt.day,
            hour="**",  # Search all hours
            min="**",
            second="**",
        )
        files = sorted(files)

        dates = []
        for f_path in files:
            try:
                # Parse start time from filename
                filename_base = os.path.basename(f_path)
                # Assuming YYYYMMDD-HHMMSS is in a part like 'YYYYMMDD-HHMMSS'
                match = re.search(r"(\d{4}\d{2}\d{2}-\d{2}\d{2}\d{2})", filename_base)
                if match:
                    file_start_dt = datetime.datetime.strptime(
                        match.group(1), "%Y%m%d-%H%M%S"
                    )
                    dates.append(file_start_dt)
            except ValueError:
                # Skip if filename format is unexpected for date parsing
                pass

        next_available = None
        for i, date_obj in enumerate(dates):  # Changed 'date' to 'date_obj'
            if date_obj > dt:  # We want the start time of the *next* file
                # check if the gap is less than max_gap_hrs
                if (date_obj - dt).total_seconds() / 3600 > max_gap_hrs:
                    return None  # Gap is too large
                next_available = date_obj
                return next_available  # Return start time of the next file

        # If no file found today after dt, try the next day recursively.
        start_of_next_day = (dt + datetime.timedelta(days=1)).replace(
            hour=0, minute=0, second=0, microsecond=0
        )
        time_to_next_day_hrs = (start_of_next_day - dt).total_seconds() / 3600

        remaining_max_gap_hrs = max_gap_hrs - time_to_next_day_hrs
        if remaining_max_gap_hrs < 0:
            return None

        return self.get_next_time(start_of_next_day, remaining_max_gap_hrs)
