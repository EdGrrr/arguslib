from arguslib.misc.met import download_era5_winds
from arguslib.misc.times import convert_to_london_naive
import numpy as np
import netCDF4
import tqdm
import datetime
import pandas as pd

from ..misc import geo
import xarray as xr

# Need some system to maintain aircraft position inforamtion
# 4 * 60 * 24 = ~6000 slots per day per aircraft, would need to store lon, lat, alt, speed, direction along with some aircarft metadata (where required).

# This global handler can cause issues if AircraftInterface is instantiated more than once,
# as the underlying file handles in the ERA5WindData objects can become stale.
# The fix is to create a local instance of this handler inside assign_era5_winds.
try:
    from csat2.ECMWF import ERA5WindData

    ERA5_LEVELS = ["150hPa", "175hPa", "200hPa", "225hPa", "250hPa", "300hPa", "350hPa"]
    # This global handler is left for any other part of the code that might use it,
    # but the problematic assign_era5_winds will use its own local copy.
    ERA5_DATA_HANDLER = {lvl: ERA5WindData(lvl) for lvl in ERA5_LEVELS}
except ImportError:
    ERA5_DATA_HANDLER = None
    ERA5_LEVELS = {}


def jsonfloat(value):
    output = float(value)
    if np.isfinite(output):
        return output
    else:
        return -9999999


class AircraftPos:
    def __init__(self, time_resolution=15, variables=["lon", "lat", "alt", "geom"]):
        if (time_resolution * (60 // time_resolution)) != 60:
            raise ValueError("Time resolution must be factor of 60")
        self.variables = variables
        self.positions = np.full(
            (60 * 60 * 24 // time_resolution, len(self.variables)), np.nan
        )  # index 0  is at midnight. index 1 is at T, index 2 is 2T
        self.time_resolution = time_resolution  # T in seconds
        self.last_update = None
        self.times = np.arange(0, 60 * 60 * 24, self.time_resolution)

    def add_position(self, dtime, aircraft_data):
        # Add aircraft data to the position array
        add_index = (
            dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        ) // self.time_resolution
        self.positions[add_index] = [
            aircraft_data.get(name, np.nan) for name in self.variables
        ]
        self.last_update = dtime

    def valid_points(self):
        return np.sum(np.isfinite(self.positions[:, 0]))

    def get_current(self, dtime, vname=None):
        # Return the aircraft position (and data) for a given time index)
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution
        if vname is None:
            # Return all variables
            tempdata = self.positions[index]
            return {
                vn: jsonfloat(tempdata[self.variables.index(vn)])
                for vn in self.variables
            }

        # Account for a string/single variable name
        if type(vname) == type("str"):
            vname = [vname]
        return {
            vn: jsonfloat(self.positions[index, self.variables.index(vn)])
            for vn in vname
        }

    def get_trail(
        self,
        dtime,
        tlen=2 * 60 * 60,
        spread_velocity=-1,
        wind_filter=-1,
        winds="era5",
        include_time=False,
        adjust_mps=(0, 0),
    ):
        # The advected flight locations (based on the aircraft
        # measured windspeed for some time after dtime). If
        # 'spread_velocity' is set positive (in m/s), this also
        # returns lon/lat positions for a given edge of the plume
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = (
            daysec // self.time_resolution + 1
        )  # Index immediately after requested dt
        length = (tlen // self.time_resolution) + 2
        startind = int(max(0, index - length))

        times = (daysec - self.times)[
            startind:index
        ]  # Time since the aircraft passed this point

        # Slice data
        lon = self.positions[startind:index, self.variables.index("lon")]
        lat = self.positions[startind:index, self.variables.index("lat")]

        if winds == "aircraft":
            ws = self.positions[startind:index, self.variables.index("ws")]
            wd = self.positions[startind:index, self.variables.index("wd")]
            wind_u, wind_v = (
                -0.51444 * ws * np.sin(np.deg2rad(wd)),
                -0.51444 * ws * np.cos(np.deg2rad(wd)),
            )
        elif winds == "era5":
            try:
                u_idx = self.variables.index("uwind")
                v_idx = self.variables.index("vwind")
                wind_u = self.positions[startind:index, u_idx]
                wind_v = self.positions[startind:index, v_idx]
                if np.all(np.isnan(wind_u)):  # Check if winds are missing
                    raise ValueError("Pre-calculated ERA5 winds are all NaN.")
            except (ValueError, IndexError):
                # This fallback should ideally not be used if assign_era5_winds has been run
                if include_time:
                    return np.full((3, len(lon)), np.nan)
                return np.full((2, len(lon)), np.nan)

        # Apply the adjustment before filtering or any other processing
        if adjust_mps != (0, 0):
            wind_u = wind_u + adjust_mps[0]
            wind_v = wind_v + adjust_mps[1]

        # Use a NaN-aware rolling average to filter the winds
        if wind_filter > 0 and len(wind_u) > 0:
            wind_u = (
                pd.Series(wind_u)
                .rolling(window=int(wind_filter), min_periods=1, center=True)
                .median()
                .to_numpy()
            )
            wind_v = (
                pd.Series(wind_v)
                .rolling(window=int(wind_filter), min_periods=1, center=True)
                .median()
                .to_numpy()
            )

        track_offset_km = np.array([wind_u, wind_v]) * times / 1000

        if not include_time:
            track_pos = np.array(geo.xy_offset_to_ll(lon, lat, *track_offset_km))
        else:
            track_pos = np.concatenate(
                [geo.xy_offset_to_ll(lon, lat, *track_offset_km), np.atleast_2d(times)]
            )

        if spread_velocity > 0:
            # Calculate the cross-track spreading of the trail
            track = self.positions[startind:index, self.variables.index("track")]
            spread_direction = track + 90
            spread_u, spread_v = spread_velocity * (
                np.array(
                    [
                        np.sin(np.deg2rad(spread_direction)),
                        np.cos(np.deg2rad(spread_direction)),
                    ]
                )
            )
            track_leftoffset_km = (
                track_offset_km + np.array([spread_u, spread_v]) * times / 1000
            )
            track_leftpos = np.array(
                geo.xy_offset_to_ll(lon, lat, *track_leftoffset_km)
            )
            track_rightoffset_km = (
                track_offset_km - np.array([spread_u, spread_v]) * times / 1000
            )
            track_rightpos = np.array(
                geo.xy_offset_to_ll(lon, lat, *track_rightoffset_km)
            )
            return track_pos, track_leftpos, track_rightpos

        return track_pos

    def interpolate_position(self, dtime, alt_var="alt_geom"):
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution
        # index before the current pos
        # different to get_trail's index, which indexes up to (and includes) this.

        daysec_us = daysec + dtime.microsecond * 1e-6

        # Collect the values either side of the aircraft
        time = (daysec_us - self.times)[
            index : index + 2
        ]  # Time since the aircraft passed this point
        lon = self.positions[index : index + 2, self.variables.index("lon")]
        lat = self.positions[index : index + 2, self.variables.index("lat")]
        alt = self.positions[index : index + 2, self.variables.index(alt_var)]

        pos = np.array([lon, lat, alt])
        # Handle NaNs that might arise from slicing at the edge
        if np.any(np.isnan(pos)) or np.any(np.isnan(time)):
            return np.array([np.nan, np.nan, np.nan])

        dpos_dtime = (pos[:, 1] - pos[:, 0]) / (time[1] - time[0])
        pos = pos[:, 0] + (0 - time[0]) * dpos_dtime

        return pos

    def get_heading(self, dtime, advection_winds="none", **kwargs):
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution

        daysec_us = daysec + dtime.microsecond * 1e-6

        if advection_winds == "none":
            lon = self.positions[index : index + 2, self.variables.index("lon")]
            lat = self.positions[index : index + 2, self.variables.index("lat")]
        else:
            trail = self.get_trail(
                dtime,
                tlen=60,
                winds=advection_winds,
                include_time=True,
                **kwargs,
            )
            last_time = trail[2, -1]
            if last_time != 0:
                current_pos = self.interpolate_position(dtime)
                lon = np.array([trail[0, -1], current_pos[0]])
                lat = np.array([trail[1, -1], current_pos[1]])
            else:
                lon = trail[0, -2:]
                lat = trail[1, -2:]

        if np.any(np.isnan(lon)) or np.any(np.isnan(lat)):
            return np.nan

        bearing = geo.calculate_bearing(lat[0], lon[0], lat[1], lon[1])
        return bearing

    def get_track(self, dtime, tlen=2 * 60 * 60, include_time=False):
        # Provide the historical locations for the aircraft for some
        # period 'tlen' seconds behind the aircraft for a given time
        # 'dtime'
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution + 1  # index immediately after current
        length = (tlen // self.time_resolution) + 2  # For before+after slots
        startind = int(max(0, index - length))

        times = (daysec - self.times)[startind:index]

        lon = self.positions[startind:index, self.variables.index("lon")]
        lat = self.positions[startind:index, self.variables.index("lat")]
        alt = self.positions[startind:index, self.variables.index("alt_geom")]
        return np.array([lon, lat, alt] + ([times] if include_time else []))

    def get_data(self, dtime, vname, tlen=2 * 60 * 60):
        #'Data' here is constant (it does not vary with
        # location/advection). This could include the initial
        # temperature, altitude or aircraft speed. It can vary with
        # aircraft position/time (unlike aircraft type, for example)
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution + 1  # index immediately after current
        length = (tlen // self.time_resolution) + 2
        startind = int(max(0, index - length))

        times = (daysec - self.times)[startind:index]

        if type(vname) == type("str"):
            vname = [vname]
        return {
            vn: self.positions[startind:index, self.variables.index(vn)] for vn in vname
        }


class Aircraft:
    def __init__(
        self,
        icao24,
        atype=None,
        time_resolution=15,
        variables=["lon", "lat", "alt_geom"],
    ):
        self.icao24 = icao24
        self.atype = atype
        self.pos = AircraftPos(time_resolution=time_resolution, variables=variables)

    def add_position(self, dtime, aircraft_data):
        self.pos.add_position(dtime, aircraft_data)

    def valid_points(self):
        return self.pos.valid_points()

    def __repr__(self):
        return f"{self.icao24} - {self.atype}: {self.valid_points()} points"

    def get_current(self, dtime, vname=None):
        # We don't need to add icao24 here as it is used as an
        # identifier elsewhere
        acdata = self.pos.get_current(dtime, vname)
        acdata["atype"] = self.atype
        return acdata

    def get_track(self, dtime, tlen=2 * 60 * 60, include_time=False):
        return self.pos.get_track(dtime, tlen, include_time=include_time)

    def get_trail(
        self,
        dtime,
        tlen=2 * 60 * 60,
        spread_velocity=-1,
        winds="era5",
        wind_filter=-1,
        include_time=False,
        adjust_mps=(0, 0),
    ):
        return self.pos.get_trail(
            dtime,
            tlen,
            spread_velocity=spread_velocity,
            wind_filter=wind_filter,
            include_time=include_time,
            winds=winds,
            adjust_mps=adjust_mps,
        )

    def get_data(self, dtime, vname, tlen=2 * 60 * 60):
        return self.pos.get_data(dtime, vname, tlen)


class Fleet:
    def __init__(self, time_resolution=15, variables=["lon", "lat", "alt_geom"]):
        self.aircraft = {}
        self.variables = variables
        self.last_update = None
        self.time_resolution = time_resolution

        self.loaded_file = None

    def add_data(self, dtime, acdata):
        if self.last_update and (dtime.hour != self.last_update.hour):
            self.write_output(self.last_update.strftime("%Y%m%d_ADS-B"))
        if self.last_update and (dtime.day != self.last_update.day):
            # Have passed midnight, update the internal arrays!
            self.update_internal()

        for icao24 in acdata.keys():
            try:
                self.aircraft[icao24].add_position(dtime, acdata[icao24][1])
            except KeyError:
                # Add new aircraft
                self.aircraft[icao24] = Aircraft(
                    icao24,
                    acdata[icao24][0],
                    time_resolution=self.time_resolution,
                    variables=self.variables,
                )
                self.aircraft[icao24].add_position(dtime, acdata[icao24][1])
        self.last_update = dtime

    def has_notnull_data(self, var):
        if var not in self.variables:
            return False

        var_pos = self.variables.index(var)
        for ac in self.aircraft.keys():
            if np.any(np.isfinite(self.aircraft[ac].pos.positions[:, var_pos])):
                return True

        return False

    def update_internal(self):
        self.write_output(self.last_update.strftime("%Y%m%d_ADS-B"))
        self.aircraft = {}

    def __repr__(self):
        return (
            f"Fleet: {len(self.aircraft.keys())}"
            + f" aircraft, {np.sum([self.aircraft[ac].valid_points() for ac in self.aircraft.keys()])} points"
        )

    def write_output(self, filename):
        # import pickle
        # with open(filename, 'wb') as f:
        #     pickle.dump(self.aircraft, f)
        aircraft_list = self.aircraft.keys()

        datalength = self.aircraft[list(aircraft_list)[0]].pos.positions.shape[0]
        datastart = 0

        with netCDF4.Dataset(filename + ".nc", "w", format="NETCDF4") as ncdf:
            ncdf.createDimension("time", datalength)
            ncdf.createDimension("aircraft", len(aircraft_list))
            dnames = ["aircraft", "time"]

            for vind, vname in enumerate(self.variables):
                Var = ncdf.createVariable(vname, "f", dnames, zlib=True)
                output = np.array(
                    [
                        self.aircraft[acft].pos.positions[datastart:, vind]
                        for acft in aircraft_list
                    ]
                )
                Var[:] = output.astype("float")

        with open(filename + ".txt", "w") as metafile:
            for acft in aircraft_list:
                metafile.write(f"{acft} {self.aircraft[acft].atype}\n")

    def load_output(self, filename):
        # import pickle
        # with open(filename, 'rb') as f:
        #     self.aircraft = pickle.load(f)

        if self.loaded_file == filename:
            return

        print(f"Loading ADS-B data from: {filename}")

        acft_list = []
        acft_types = {}
        with open(filename + ".txt", "r") as metafile:
            for line in metafile.readlines():
                splitted = line.strip().split(" ")
                acft = splitted[0]
                atype = splitted[1]
                acft_list.append(acft)
                acft_types[acft] = atype

        # determine offset based on the date
        # adsb data is reported in local (UK) time. - in summer, UTC+1.
        # Determine the offset move that many indices EARLIER
        file_dtime = datetime.datetime.strptime(
            (filename.split("/")[-1]).split("_")[0], "%Y%m%d"
        )
        offset = int(
            (
                convert_to_london_naive(file_dtime.replace(hour=12))
                - file_dtime.replace(hour=12)
            ).total_seconds()
            / self.time_resolution
        )

        self.aircraft: dict[str, Aircraft] = {}
        with netCDF4.Dataset(filename + ".nc") as ncdf:
            var_data = []
            for vind, vname in enumerate(self.variables):
                try:
                    var_data.append(ncdf.variables[vname])  # axes: 0: aircraft, 1: time
                except (
                    KeyError
                ):  # Likely a field that has been added programatically to be obtained from another source (i.e. uwind, vwind). Set these to be NaN
                    var_data.append(None)
            first_actual_vardata = 0
            while var_data[first_actual_vardata] is None:
                first_actual_vardata += 1
            var_data = [
                (
                    v
                    if v is not None
                    else np.full_like(var_data[first_actual_vardata], np.nan)
                )
                for v in var_data
            ]
            var_data = np.array(var_data)  # axes: 0: variable, 1: aircraft, 2: time
            if offset != 0:
                var_data = np.concatenate(
                    [
                        var_data[:, :, offset:],
                        np.full((*var_data.shape[:2], offset), np.nan),
                    ],
                    axis=-1,
                )

            for aind, acft_name in tqdm.tqdm(
                enumerate(acft_list), desc="Loading ADS-B data", unit="acft"
            ):
                self.aircraft[acft_name] = Aircraft(
                    acft_name,
                    acft_types[acft_name],
                    self.time_resolution,
                    self.variables,
                )
                self.aircraft[acft_name].pos.positions = var_data[:, aind].T

        self.loaded_file = filename

    def list_current(self):
        return self.aircraft.keys()

    def get_current(self, dtime, vname=None):
        acdata = {}
        for ac in self.aircraft.keys():
            tempdata = self.aircraft[ac].get_current(dtime, vname)
            if np.isfinite(tempdata["lon"] + tempdata["lat"]):
                acdata[ac] = tempdata
        return acdata

    def get_tracks(self, dtime, tlen=2 * 60 * 60, include_time=False):
        """Returns unadvected track position (lon, lat, alt, and potentially time [in seconds befor dtime]) for now and every previous 15 sec until tlen (in min)

        Currently uses aircraft wind - I don't think this is accurate when the aircraft is climbing of descending
        """
        tracks = {}
        for ac in self.aircraft.keys():
            tracks[ac] = self.aircraft[ac].get_track(
                dtime, tlen, include_time=include_time
            )
        return tracks

    def get_trails(
        self,
        dtime,
        tlen=2 * 60 * 60,
        spread_velocity=-1,
        wind_filter=-1,
        include_time=False,
        winds="era5",
        adjust_mps=(0, 0),
    ):
        """Returns trail position (lon, lat, and potentially time [in seconds befor dtime]) for now and every previous 15 sec until tlen (in min)

        Currently uses aircraft wind - I don't think this is accurate when the aircraft is climbing of descending
        """
        if winds == "none":
            # do not advect at all
            return self.get_tracks(dtime, tlen, include_time=include_time)
        if (winds == "era5") and (
            ("uwind" not in self.variables) or (not self.has_notnull_data("uwind"))
        ):
            print(
                f"Warning (Fleet.get_trails): No ERA5 wind data available for {dtime}. Using aircraft ADS-B winds (quite inaccurate)."
            )
            return self.get_trails(
                dtime,
                tlen,
                spread_velocity,
                wind_filter,
                include_time,
                winds="aircraft",
                adjust_mps=adjust_mps,
            )

        trails = {}
        for ac in self.aircraft.keys():
            trails[ac] = self.aircraft[ac].get_trail(
                dtime,
                tlen,
                spread_velocity=spread_velocity,
                wind_filter=wind_filter,
                include_time=include_time,
                winds=winds,
                adjust_mps=adjust_mps,
            )
        return trails

    def get_data(self, dtime, vname, tlen=2 * 60 * 60):
        metadata = {}
        for ac in self.aircraft.keys():
            metadata[ac] = self.aircraft[ac].get_data(dtime, vname, tlen)
        return metadata

    def assign_era5_winds(self, _download_attempted_this_call=False):
        """
        Pre-calculates and assigns ERA5 wind data (uwind, vwind) to all aircraft positions.

        This optimized method uses vectorized NumPy operations to accelerate performance. It first
        stacks all flight data into a single large array, identifies all valid data points at once,
        and then performs a single, efficient 4D interpolation for time-grouped chunks of points.
        The final assignment of wind data is also vectorized, avoiding slow Python loops.

        This is intended to be run once after loading data to prevent slow on-the-fly wind
        lookups during repeated calls to get_trail.

        Raises:
            RuntimeError: If ERA5_DATA_HANDLER is not available, if data has not been loaded,
                          or if a download was attempted but files are still missing.
            ValueError: If the date cannot be parsed from the loaded filename.
        """
        if ERA5_DATA_HANDLER is None:
            raise RuntimeError(
                "ERA5WindData handler is not available. Please check your csat2 installation."
            )

        if self.loaded_file is None:
            raise RuntimeError(
                "Fleet data must be loaded using 'load_output' before assigning wind data."
            )

        # --- FIX: Create a fresh, local instance of the ERA5 data handler ---
        # This avoids using the global, shared handler which causes file handle errors
        # on subsequent calls by ensuring each call to this method gets its own objects.
        try:
            from csat2.ECMWF import ERA5WindData

            local_era5_data_handler = {lvl: ERA5WindData(lvl) for lvl in ERA5_LEVELS}
        except ImportError:
            raise RuntimeError("csat2.ECMWF.ERA5WindData could not be imported.")

        # --- 1. Add uwind and vwind to variables ---
        new_vars = [var for var in ["uwind", "vwind"] if var not in self.variables]
        if new_vars:
            # print(f"Adding {new_vars} to fleet variables.")
            original_var_count = len(self.variables)
            self.variables.extend(new_vars)
            for acft in self.aircraft.values():
                old_positions = acft.pos.positions
                new_shape = (old_positions.shape[0], len(self.variables))
                new_positions = np.full(new_shape, np.nan, dtype=np.float32)
                new_positions[:, :original_var_count] = old_positions
                acft.pos.positions = new_positions
                acft.pos.variables = self.variables

        # --- 2. Parse date from filename ---
        import os
        import datetime

        try:
            fname = os.path.basename(self.loaded_file)
            date_str = fname.split("_")[0]
            base_date = datetime.datetime.strptime(date_str, "%Y%m%d")
        except (ValueError, IndexError):
            raise ValueError(
                f"Could not parse date from filename: {self.loaded_file}. Expected format is 'YYYYMMDD_...'."
            )

        # check if the era5 data is there using the new local handler
        needs_download = False
        for level, handler in local_era5_data_handler.items():
            try:
                handler.get_data_time(base_date)
            except IndexError:
                print(
                    f"ERA5 data for {level} on {base_date.date()} appears to be missing (locator.search failed)."
                )
                needs_download = True
                break
            except Exception as e:
                print(
                    f"Warning: Error checking/loading ERA5 data for {level} on {base_date.date()}: {e}. "
                )
                # Continue, as re-downloading might not fix non-IndexError issues.

        if needs_download:
            if _download_attempted_this_call:
                msg = (
                    f"ERA5 download was already attempted for {base_date.date()} but files are still "
                    f"not found. Aborting wind assignment to prevent infinite loop."
                )
                print(msg)
                raise RuntimeError(msg)

            # print(f"Attempting to download ERA5 wind data for {base_date.date()}...")
            raise ValueError("There is no ERA5 data downloaded for this date.")
            # try:
            #     download_era5_winds(base_date)
            #     print("Download attempt finished. Re-attempting wind assignment.")
            #     return self.assign_era5_winds(_download_attempted_this_call=True)
            # except Exception as e:
            #     print(f"An unexpected error occurred during ERA5 download for {base_date.date()}: {e}. ERA5 winds will be unavailable.")
            #     return

        # --- 3. Vectorized Data Preparation ---
        lon_idx, lat_idx, alt_idx = (
            self.variables.index("lon"),
            self.variables.index("lat"),
            self.variables.index("alt_geom"),
        )
        u_idx, v_idx = self.variables.index("uwind"), self.variables.index("vwind")

        aircraft_list = list(self.aircraft.values())
        if not aircraft_list:
            print("No aircraft in fleet. Aborting wind assignment.")
            return

        all_positions = np.stack([ac.pos.positions for ac in aircraft_list]).astype(
            np.float32
        )
        valid_mask = (
            np.isfinite(all_positions[:, :, lon_idx])
            & np.isfinite(all_positions[:, :, lat_idx])
            & np.isfinite(all_positions[:, :, alt_idx])
        )
        ac_indices, t_indices = np.where(valid_mask)

        if ac_indices.size == 0:
            print("No valid flight data points found. Aborting wind assignment.")
            return

        valid_lons = all_positions[ac_indices, t_indices, lon_idx]
        valid_lats = all_positions[ac_indices, t_indices, lat_idx]
        valid_alts = all_positions[ac_indices, t_indices, alt_idx]

        base_datetime64 = np.datetime64(base_date)
        time_deltas = t_indices * np.timedelta64(self.time_resolution, "s")
        valid_datetimes_np = base_datetime64 + time_deltas

        # --- 4. Group points by ERA5 time interval ---
        print("Grouping flight data by time interval...")
        steps_per_interval = (3 * 3600) // self.time_resolution
        interval_bins = t_indices // steps_per_interval
        num_intervals = (24 * 3600) // (steps_per_interval * self.time_resolution)

        # --- 5. Process each group with a single 4D interpolation ---
        print("Assigning ERA5 wind data to fleet positions...")
        for i in tqdm.tqdm(range(num_intervals), desc="Processing 4D time chunks"):
            chunk_mask = interval_bins == i
            if not np.any(chunk_mask):
                continue

            start_hour = i * 3
            start_t = base_date + datetime.timedelta(hours=start_hour)
            end_t = base_date + datetime.timedelta(hours=start_hour + 3)

            chunk_ac_indices = ac_indices[chunk_mask]
            chunk_t_indices = t_indices[chunk_mask]

            try:
                # This nested function will now use the 'local_era5_data_handler' from the parent scope
                def load_3d_wind_field(dtime):
                    wind_u_3d, wind_v_3d = [], []
                    for level in ERA5_LEVELS:
                        wind_u_da, wind_v_da = local_era5_data_handler[
                            level
                        ].get_data_time(dtime)
                        wind_u_da = geo.xr_add_cyclic_points(wind_u_da)
                        wind_v_da = geo.xr_add_cyclic_points(wind_v_da)
                        wind_u_3d.append(wind_u_da)
                        wind_v_3d.append(wind_v_da)
                    plevels = [float(lvl.removesuffix("hPa")) for lvl in ERA5_LEVELS]
                    u_field = xr.concat(wind_u_3d, "plevel").assign_coords(
                        {"plevel": plevels}
                    )
                    v_field = xr.concat(wind_v_3d, "plevel").assign_coords(
                        {"plevel": plevels}
                    )
                    return u_field, v_field

                u_field_start, v_field_start = load_3d_wind_field(start_t)
                u_field_end, v_field_end = load_3d_wind_field(end_t)

                wind_u_4d = xr.concat(
                    [u_field_start, u_field_end], dim="time"
                ).assign_coords({"time": [start_t, end_t]})
                wind_v_4d = xr.concat(
                    [v_field_start, v_field_end], dim="time"
                ).assign_coords({"time": [start_t, end_t]})

                ds_track = xr.Dataset(
                    {
                        "lon": ("points", valid_lons[chunk_mask] % 360),
                        "lat": ("points", valid_lats[chunk_mask]),
                        "alt": (
                            "points",
                            np.clip(
                                geo.ft_to_hPa(valid_alts[chunk_mask]), 150.0, 350.0
                            ),
                        ),
                        "time": ("points", valid_datetimes_np[chunk_mask]),
                    }
                )

                u_vals = wind_u_4d.interp(
                    lon=ds_track.lon,
                    lat=ds_track.lat,
                    plevel=ds_track.alt,
                    time=ds_track.time,
                ).values
                v_vals = wind_v_4d.interp(
                    lon=ds_track.lon,
                    lat=ds_track.lat,
                    plevel=ds_track.alt,
                    time=ds_track.time,
                ).values

                # --- 6. Vectorized Final Assignment ---
                all_positions[chunk_ac_indices, chunk_t_indices, u_idx] = u_vals
                all_positions[chunk_ac_indices, chunk_t_indices, v_idx] = v_vals

            except Exception as e:
                print(
                    f"\nWarning: Could not process wind chunk for {start_t} to {end_t}: {e}"
                )

        # --- 7. Unstack data back into individual aircraft objects ---
        for i, acft in enumerate(aircraft_list):
            acft.pos.positions = all_positions[i]

        print("ERA5 wind assignment process complete.")
