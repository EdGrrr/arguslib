from arguslib.misc.met import download_era5_winds
import numpy as np
import netCDF4
import tqdm

from ..misc import geo
import xarray as xr

# Need some system to maintain aircraft position inforamtion
# 4 * 60 * 24 = ~6000 slots per day per aircraft, would need to store lon, lat, alt, speed, direction along with some aircarft metadata (where required).

try:
    from csat2.ECMWF import ERA5WindData
    ERA5_LEVELS = [
        "150hPa",
        "175hPa",
        "200hPa",
        "225hPa",
        "250hPa",
        "300hPa",
        "350hPa"
    ]
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
        ) # index 0  is at midnight. index 1 is at T, index 2 is 2T
        self.time_resolution = time_resolution # T in seconds
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
        winds='era5',
        include_time=False,
    ):
        # The advected flight locations (based on the aircraft
        # measured windspeed for some time after dtime). If
        # 'spread_velocity' is set positive (in m/s), this also
        # returns lon/lat positions for a given edge of the plume
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution + 1 # Index immediately after requested dt (slices go up to this value, not including it)
        length = (tlen / self.time_resolution) + 2
        startind = int(max(0, index - length))

        times = (daysec - self.times)[
            startind:index
        ]  # Time since the aircraft passed this point

        gs = self.positions[startind:index, self.variables.index("gs")]
        track = self.positions[startind:index, self.variables.index("track")]
        lon = self.positions[startind:index, self.variables.index("lon")]
        lat = self.positions[startind:index, self.variables.index("lat")]
        
        if winds == 'aircraft':
            ws = self.positions[startind:index, self.variables.index("ws")]
            wd = self.positions[startind:index, self.variables.index("wd")]
            wind_u, wind_v = ( # negative sign because these are "metorological wind directions" (e.g. coming from this direction)
                -0.51444 * ws * np.sin(np.deg2rad(wd)), # knots to m/s
                -0.51444 * ws * np.cos(np.deg2rad(wd)),
            )
        elif winds == 'era5':
            alt_geom = self.positions[startind:index, self.variables.index("alt_geom")]
            try:
                # --- NEW: Fast path using pre-calculated data ---
                u_idx = self.variables.index('uwind')
                v_idx = self.variables.index('vwind')
                wind_u = self.positions[startind:index, u_idx]
                wind_v = self.positions[startind:index, v_idx]

                # If data is all NaN, it means either there are no valid positions or pre-calculation wasn't run or failed.
                has_valid_pos = np.any(np.isfinite(lon) & np.isfinite(lat) & np.isfinite(alt_geom))
                if has_valid_pos and (np.all(np.isnan(wind_u)) or np.all(np.isnan(wind_v))):
                     raise ValueError("Pre-calculated ERA5 winds are all NaN.")

            except (ValueError, IndexError) as e:
                raise ValueError("ERA5 data not on fleet! Using aircraft wind advection, or loading ERA5.")
                # --- OLD: Fallback to slow on-the-fly calculation ---
                # FIXME: this is really slow. I have to do so much indexing every time advect to a different timestep.
                wind_u, wind_v = ERA5_DATA_HANDLER.get_data_time(dtime)
                ds_track = xr.Dataset(
                    {'lon':('points', lon % 360), 'lat':('points', lat % 360)}
                )
                wind_u = wind_u.sel(lon=ds_track.lon, lat=ds_track.lat, method='nearest').values
                wind_v = wind_v.sel(lon=ds_track.lon, lat=ds_track.lat, method='nearest').values
            
            
            
        if wind_filter > 0:

            def wind_conv_filter(wval, wind_filter):
                filtered = np.convolve(
                    np.ones(int(wind_filter)), wval, mode="same"
                ) / np.convolve(
                    np.ones(int(wind_filter)), np.ones(len(wval)), mode="same"
                )
                if wind_filter > len(wval):
                    difflen = int(wind_filter) - len(wval)
                    filtered = filtered[difflen // 2 : (difflen // 2 + len(wval))]
                return filtered

            wind_u = wind_conv_filter(wind_u, wind_filter)
            wind_v = wind_conv_filter(wind_v, wind_filter)

        track_offset_km = np.array([wind_u, wind_v]) * times / 1000


        if not include_time:
            track_pos = np.array(geo.xy_offset_to_ll(lon, lat, *track_offset_km))
        else:
            track_pos = np.concatenate(
                [geo.xy_offset_to_ll(lon, lat, *track_offset_km), np.atleast_2d(times)]
            )

        if spread_velocity > 0:
            # Calculate the cross-track spreading of the trail
            spread_direction = (
                track + 90
            )  # Simple for now - mostly governed by aircraft speed/direction
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
            return track_pos, track_leftpos, track_rightp

        return track_pos
    
    def interpolate_position(self, dtime):
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution 
        # index before the current pos
        # different to get_trail's index, which indexes up to (and includes) this.
        
        daysec_us = daysec + dtime.microsecond*1e-6

        # Collect the values either side of the aircraft
        time = (daysec_us - self.times)[
            index:index+2
        ]  # Time since the aircraft passed this point
        lon = self.positions[index:index+2, self.variables.index("lon")]
        lat = self.positions[index:index+2, self.variables.index("lat")]
        alt = self.positions[index:index+2, self.variables.index("alt_geom")]
        
        pos = np.array([lon, lat, alt])
        dpos_dtime = (pos[:, 1] - pos[:, 0]) / (time[1] - time[0])
        pos = pos[:, 0] + (0 - time[0]) * dpos_dtime
        
        return pos
    
    def get_heading(self, dtime):
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution
        
        daysec_us = daysec + dtime.microsecond*1e-6

        time = (daysec_us - self.times)[
            index:index+2 # time before and after
        ]  # Time since the aircraft passed this point
        lon = self.positions[index:index+2, self.variables.index("lon")]
        lat = self.positions[index:index+2, self.variables.index("lat")]
        
        bearing = geo.calculate_bearing(lat[0], lon[0], lat[1], lon[1])
        return bearing


    def get_track(self, dtime, tlen=2 * 60 * 60, include_time=False):
        # Provide the historical locations for the aircraft for some
        # period 'tlen' seconds behind the aircraft for a given time
        # 'dtime'
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution + 1 # index immediately after current
        length = (tlen / self.time_resolution) + 2  # For before+after slots
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
        index = daysec // self.time_resolution + 1 # index immediately after current
        length = (tlen / self.time_resolution) + 2  # For before+after slots
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
        winds='era5',
        wind_filter=-1,
        include_time=False,
    ):
        return self.pos.get_trail(
            dtime,
            tlen,
            spread_velocity=spread_velocity,
            wind_filter=wind_filter,
            include_time=include_time,
            winds=winds
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

        acft_list = []
        acft_types = {}
        with open(filename + ".txt", "r") as metafile:
            for line in metafile.readlines():
                acft, atype = line.strip().split(" ")
                acft_list.append(acft)
                acft_types[acft] = atype

        self.aircraft: dict[str, Aircraft] = {}
        with netCDF4.Dataset(filename + ".nc") as ncdf:
            var_data = []
            for vind, vname in enumerate(self.variables):
                var_data.append(ncdf.variables[vname])
            var_data = np.array(var_data)

            for aind, acft_name in tqdm.tqdm(
                enumerate(acft_list), desc="Loading ADS-B data", unit="acft"
            ):
                self.aircraft[acft_name] = Aircraft(
                    acft_name,
                    acft_types[acft_name],
                    self.time_resolution,
                    self.variables,
                )
                for vind, vname in enumerate(self.variables):
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
            tracks[ac] = self.aircraft[ac].get_track(dtime, tlen, include_time=include_time)
        return tracks

    def get_trails(
        self,
        dtime,
        tlen=2 * 60 * 60,
        spread_velocity=-1,
        wind_filter=-1,
        include_time=False,
        winds='era5',
    ):
        """Returns trail position (lon, lat, and potentially time [in seconds befor dtime]) for now and every previous 15 sec until tlen (in min)

        Currently uses aircraft wind - I don't think this is accurate when the aircraft is climbing of descending
        """
        trails = {}
        for ac in self.aircraft.keys():
            trails[ac] = self.aircraft[ac].get_trail(
                dtime,
                tlen,
                spread_velocity=spread_velocity,
                wind_filter=wind_filter,
                include_time=include_time,
                winds=winds,
            )
        return trails

    def get_data(self, dtime, vname, tlen=2 * 60 * 60):
        metadata = {}
        for ac in self.aircraft.keys():
            metadata[ac] = self.aircraft[ac].get_data(dtime, vname, tlen)
        return metadata
    
    def assign_era5_winds(self):
        """
        Pre-calculates and assigns ERA5 wind data (uwind, vwind) to all aircraft positions,
        using 4D interpolation for time, pressure level, latitude, and longitude.

        This method chunks flight data into the intervals between available ERA5 data points
        (e.g., 09:00-12:00). For each chunk, it loads the wind fields for the start and
        end times and performs a single 4D interpolation for all points within it.

        This is intended to be run once after loading data to prevent slow
        on-the-fly wind lookups during repeated calls to get_trail.

        Raises:
            RuntimeError: If ERA5_DATA_HANDLER is not available or if data
                        has not been loaded via load_output.
            ValueError: If the date cannot be parsed from the loaded filename.
        """
        if ERA5_DATA_HANDLER is None:
            raise RuntimeError("ERA5WindData handler is not available. Please check your csat2 installation.")

        if self.loaded_file is None:
            raise RuntimeError("Fleet data must be loaded using 'load_output' before assigning wind data.")

        # --- 1. Add uwind and vwind to variables if they don't exist ---
        new_vars = []
        if 'uwind' not in self.variables:
            new_vars.append('uwind')
        if 'vwind' not in self.variables:
            new_vars.append('vwind')

        if new_vars:
            print(f"Adding {new_vars} to fleet variables.")
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
            date_str = fname.split('_')[0]
            base_date = datetime.datetime.strptime(date_str, '%Y%m%d')
        except (ValueError, IndexError):
            raise ValueError(f"Could not parse date from filename: {self.loaded_file}. Expected format is 'YYYYMMDD_...'.")

        # check if the era5 data is there
        for level, handler in ERA5_DATA_HANDLER.items():
            try: handler.get_data_time(base_date)
            except:
                download_era5_winds(base_date)
                return self.assign_era5_winds()
            
                

        # --- 3. Get variable indices ---
        lon_idx = self.variables.index("lon")
        lat_idx = self.variables.index("lat")
        alt_idx = self.variables.index("alt_geom")
        u_idx = self.variables.index("uwind")
        v_idx = self.variables.index("vwind")

        aircraft_list = list(self.aircraft.values())
        if not aircraft_list:
            print("No aircraft in fleet. Aborting wind assignment.")
            return
            
        num_timesteps = aircraft_list[0].pos.positions.shape[0]

        # --- 4. Group points by ERA5 time interval, then process in 4D chunks. ---
        print("Collecting and grouping all flight positions by time interval...")
        avail_era5_times = [base_date + datetime.timedelta(hours=hrs) for hrs in range(0, 24, 3)]
        
        # Create the time intervals (e.g., (00:00, 03:00), (03:00, 06:00), ...)
        time_intervals = list(zip(avail_era5_times, avail_era5_times[1:]))
        
        # Prepare the data structure for grouping
        points_by_interval = {interval: {'lons': [], 'lats': [], 'alts': [], 'times': [], 'ac_indices': [], 't_indices': []} for interval in time_intervals}

        for ac_idx, acft in enumerate(tqdm.tqdm(aircraft_list, desc="Grouping flight data")):
            for t_idx in range(num_timesteps):
                lon = acft.pos.positions[t_idx, lon_idx]
                lat = acft.pos.positions[t_idx, lat_idx]
                alt = acft.pos.positions[t_idx, alt_idx]

                if np.isfinite(lon) and np.isfinite(lat) and np.isfinite(alt):
                    current_dtime = base_date + datetime.timedelta(seconds=(t_idx * self.time_resolution))
                    
                    # Find which time interval this point falls into
                    for start_t, end_t in time_intervals:
                        if start_t <= current_dtime < end_t:
                            group = points_by_interval[(start_t, end_t)]
                            group['lons'].append(lon)
                            group['lats'].append(lat)
                            group['alts'].append(geo.ft_to_hPa(alt))
                            group['times'].append(current_dtime)
                            group['ac_indices'].append(ac_idx)
                            group['t_indices'].append(t_idx)
                            break

        # 4b. Process each group with a single large 4D interpolation call.
        print("Assigning ERA5 wind data to fleet positions...")
        for (start_t, end_t), group in tqdm.tqdm(points_by_interval.items(), desc="Processing 4D time chunks"):
            if not group['lons']:
                continue

            try:
                # --- Load data for BOTH start and end of the time interval ---
                def load_3d_wind_field(dtime):
                    wind_u_3d, wind_v_3d = [], []
                    for level in ERA5_LEVELS:
                        wind_u_da, wind_v_da = ERA5_DATA_HANDLER[level].get_data_time(dtime)
                        wind_u_da = geo.xr_add_cyclic_points(wind_u_da)
                        wind_v_da = geo.xr_add_cyclic_points(wind_v_da)
                        wind_u_3d.append(wind_u_da)
                        wind_v_3d.append(wind_v_da)
                    
                    u_field = xr.concat(wind_u_3d, "plevel").assign_coords({"plevel":[float(lvl.removesuffix("hPa")) for lvl in ERA5_LEVELS]})
                    v_field = xr.concat(wind_v_3d, "plevel").assign_coords({"plevel":[float(lvl.removesuffix("hPa")) for lvl in ERA5_LEVELS]})
                    return u_field, v_field

                u_field_start, v_field_start = load_3d_wind_field(start_t)
                u_field_end, v_field_end = load_3d_wind_field(end_t)

                # --- Combine 3D fields into a 4D field for interpolation ---
                wind_u_4d = xr.concat([u_field_start, u_field_end], dim='time').assign_coords({'time': [start_t, end_t]})
                wind_v_4d = xr.concat([v_field_start, v_field_end], dim='time').assign_coords({'time': [start_t, end_t]})

                # Create a single large Dataset of points for interpolation, including their exact times
                ds_track = xr.Dataset({
                    'lon': ('points', np.array(group['lons']) % 360),
                    'lat': ('points', np.array(group['lats'])),
                    'alt': ('points', np.clip(group['alts'], 150., 350.)),
                    'time': ('points', np.array(group['times'])) # The exact timestamp of each point
                })
                
                # --- Perform one big, efficient 4D interpolation ---
                u_vals = wind_u_4d.interp(lon=ds_track.lon, lat=ds_track.lat, plevel=ds_track.alt, time=ds_track.time).values
                v_vals = wind_v_4d.interp(lon=ds_track.lon, lat=ds_track.lat, plevel=ds_track.alt, time=ds_track.time).values

                # Assign all results back to the main data structure
                for i in range(len(group['lons'])):
                    ac_idx = group['ac_indices'][i]
                    t_idx = group['t_indices'][i]
                    aircraft_list[ac_idx].pos.positions[t_idx, u_idx] = u_vals[i]
                    aircraft_list[ac_idx].pos.positions[t_idx, v_idx] = v_vals[i]

            except Exception as e:
                print(f"Warning: Could not process wind chunk for {start_t} to {end_t}: {e}")
