import numpy as np
import netCDF4
import tqdm

from ..misc import geo

# Need some system to maintain aircraft position inforamtion
# 4 * 60 * 24 = ~6000 slots per day per aircraft, would need to store lon, lat, alt, speed, direction along with some aircarft metadata (where required).


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
        )
        self.time_resolution = time_resolution
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

    def get_trail(self, dtime, tlen=2 * 60 * 60, spread_velocity=-1, wind_filter=-1):
        # The advected flight locations (based on the aircraft
        # measured windspeed for some time after dtime). If
        # 'spread_velocity' is set positive (in m/s), this also
        # returns lon/lat positions for a given edge of the plume
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution
        length = (tlen / self.time_resolution) + 2  # For before+after slots
        startind = int(max(0, index - length))

        times = (daysec - self.times)[startind:index]

        gs = self.positions[startind:index, self.variables.index("gs")]
        track = self.positions[startind:index, self.variables.index("track")]
        ws = self.positions[startind:index, self.variables.index("ws")]
        wd = self.positions[startind:index, self.variables.index("wd")]
        wind_u, wind_v = (
            0.51444 * ws * np.sin(np.deg2rad(wd)),
            0.51444 * ws * np.cos(np.deg2rad(wd)),
        )
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

        acft_u, acft_v = (
            0.51444 * gs * np.sin(np.deg2rad(track)),
            0.51444 * gs * np.cos(np.deg2rad(track)),
        )

        track_offset_km = np.array([wind_u, wind_v]) * times / 1000

        lon = self.positions[startind:index, self.variables.index("lon")]
        lat = self.positions[startind:index, self.variables.index("lat")]

        track_pos = np.array(geo.xy_offset_to_ll(lon, lat, *track_offset_km))

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

    def get_track(self, dtime, tlen=2 * 60 * 60):
        # Provide the historical locations for the aircraft for some
        # period 'tlen' seconds behind the aircraft for a given time
        # 'dtime'
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution
        length = (tlen / self.time_resolution) + 2  # For before+after slots
        startind = int(max(0, index - length))

        times = (daysec - self.times)[startind:index]

        lon = self.positions[startind:index, self.variables.index("lon")]
        lat = self.positions[startind:index, self.variables.index("lat")]
        alt = self.positions[startind:index, self.variables.index("alt_geom")]

        return np.array([lon, lat, alt])

    def get_data(self, dtime, vname, tlen=2 * 60 * 60):
        #'Data' here is constant (it does not vary with
        # location/advection). This could include the initial
        # temperature, altitude or aircraft speed. It can vary with
        # aircraft position/time (unlike aircraft type, for example)
        daysec = dtime.hour * 3600 + dtime.minute * 60 + dtime.second
        offset = daysec % self.time_resolution
        index = daysec // self.time_resolution
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

    def get_track(self, dtime, tlen=2 * 60 * 60):
        return self.pos.get_track(dtime, tlen)

    def get_trail(self, dtime, tlen=2 * 60 * 60, spread_velocity=-1, wind_filter=-1):
        return self.pos.get_trail(
            dtime, tlen, spread_velocity=spread_velocity, wind_filter=wind_filter
        )

    def get_data(self, dtime, vname, tlen=2 * 60 * 60):
        return self.pos.get_data(dtime, vname, tlen)


class Fleet:
    def __init__(self, time_resolution=15, variables=["lon", "lat", "alt_geom"]):
        self.aircraft = {}
        self.variables = variables
        self.last_update = None
        self.time_resolution = time_resolution

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

        self.aircraft = {}
        with netCDF4.Dataset(filename + ".nc") as ncdf:
            var_data = []
            for vind, vname in enumerate(self.variables):
                var_data.append(ncdf.variables[vname])
            var_data = np.array(var_data)

            for aind, acft_name in tqdm.tqdm(enumerate(acft_list)):
                self.aircraft[acft_name] = Aircraft(
                    acft_name,
                    acft_types[acft_name],
                    self.time_resolution,
                    self.variables,
                )
                for vind, vname in enumerate(self.variables):
                    self.aircraft[acft_name].pos.positions = var_data[:, aind].T

    def list_current(self):
        return self.aircraft.keys()

    def get_current(self, dtime, vname=None):
        acdata = {}
        for ac in self.aircraft.keys():
            tempdata = self.aircraft[ac].get_current(dtime, vname)
            if np.isfinite(tempdata["lon"] + tempdata["lat"]):
                acdata[ac] = tempdata
        return acdata

    def get_tracks(self, dtime, tlen=2 * 60 * 60):
        tracks = {}
        for ac in self.aircraft.keys():
            tracks[ac] = self.aircraft[ac].get_track(dtime, tlen)
        return tracks

    def get_trails(self, dtime, tlen=2 * 60 * 60, spread_velocity=-1, wind_filter=-1):
        """Returns trail position for now and every previous 15 sec until tlen (in min)

        Currently uses aircraft wind - I don't think this is accurate when the aircraft is climbing of descending
        """
        trails = {}
        for ac in self.aircraft.keys():
            trails[ac] = self.aircraft[ac].get_trail(
                dtime, tlen, spread_velocity=spread_velocity, wind_filter=wind_filter
            )
        return trails

    def get_data(self, dtime, vname, tlen=2 * 60 * 60):
        metadata = {}
        for ac in self.aircraft.keys():
            metadata[ac] = self.aircraft[ac].get_data(dtime, vname, tlen)
        return metadata
