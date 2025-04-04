from matplotlib.axes import Axes

from ..misc.geo import haversine, bearing
import numpy as np
from pathlib import Path

default_calibration_file = (
    Path(__file__).parent.parent / "instruments/cam1_calibration.yml"
)


def ead_to_xyz(target_elevation, target_azimuth, target_distance):
    target_z = target_distance * np.sin(np.deg2rad(target_elevation))
    target_hor_range = target_distance * np.cos(np.deg2rad(target_elevation))
    target_x = target_hor_range * np.sin(np.deg2rad(target_azimuth))
    target_y = target_hor_range * np.cos(np.deg2rad(target_azimuth))
    return np.array([target_x, target_y, target_z])


def xyz_to_ead(target_x, target_y, target_z):
    target_distance = np.sqrt(target_x**2 + target_y**2 + target_z**2)
    target_elevation = np.rad2deg(
        np.arctan2(target_z, np.sqrt(target_x**2 + target_y**2))
    )
    target_azimuth = np.rad2deg(np.arctan2(target_x, target_y))
    return np.array([target_elevation, target_azimuth, target_distance])


def rotation_matrix_i_to_g(elevation, azimuth, roll):
    # Convert angles to radians
    azimuth = np.deg2rad(azimuth)
    elevation = np.deg2rad(90 - elevation)
    roll = np.deg2rad(roll)

    # Azimuth rotation (Rₐ) - global z-axis
    R_a = np.array(
        [
            [np.cos(azimuth), -np.sin(azimuth), 0],
            [np.sin(azimuth), np.cos(azimuth), 0],
            [0, 0, 1],
        ]
    )

    # Elevation rotation (Rₑ) - global y-axis
    R_e = np.array(
        [
            [np.cos(elevation), 0, np.sin(elevation)],
            [0, 1, 0],
            [-np.sin(elevation), 0, np.cos(elevation)],
        ]
    )

    # Roll rotation (Rᵣ) - local z-axis
    R_r = np.array(
        [[np.cos(roll), -np.sin(roll), 0], [np.sin(roll), np.cos(roll), 0], [0, 0, 1]]
    )

    # Combined rotation matrix
    R = R_a @ R_e @ R_r
    return R


class Position:
    # Position as
    # - lon, lat, alt (assuming Earth is locally flat)
    def __init__(self, lon, lat, alt):
        self.lon = lon
        self.lat = lat
        self.alt = alt

    # Elevation, Azimuth, distance
    #  - Elevation from local horizontal
    #  - Azimuth relative to local north (compass bearing)
    #  - Distance in km

    def target_ead(self, target_position):
        # Approximate - can use pyproj for exact (or with large distances)
        distance = haversine(
            self.lon, self.lat, target_position.lon, target_position.lat
        )
        alt_diff = target_position.alt - self.alt
        target_distance = np.sqrt(distance**2 + alt_diff**2)
        target_elevation = np.rad2deg(np.arctan2(alt_diff, distance))
        target_azimuth = bearing(
            self.lon, self.lat, target_position.lon, target_position.lat
        )
        return np.array([target_elevation, target_azimuth, target_distance])

    # XYZ
    #  - X East-West distace in km
    #  - Y North-south distance in km
    #  - Z altitude in km

    def target_xyz(self, target_position):
        # Assume Earth is locally flat, but also spherical
        distance = haversine(
            self.lon, self.lat, target_position.lon, target_position.lat
        )
        target_z = target_position.alt - self.alt
        target_bearing = bearing(
            self.lon, self.lat, target_position.lon, target_position.lat
        )
        target_x = distance * np.sin(np.deg2rad(target_bearing))
        target_y = distance * np.cos(np.deg2rad(target_bearing))
        return np.array([target_x, target_y, target_z])

    def ead_to_lla(self, target_elevation, target_azimuth, target_distance):
        return self.xyz_to_lla(
            *ead_to_xyz(target_elevation, target_azimuth, target_distance)
        )

    def xyz_to_lla(self, target_x, target_y, target_z):
        lon = self.lon + target_x / (111.111 * np.cos(self.lat))
        lat = self.lat + target_y / (111.111)
        return Position(lon, lat, target_z + self.alt)

    def __repr__(self):
        return f"{self.lon:.3f} {self.lat:.3f} {self.alt:.3f}"


class PlottableInstrument:
    def __init__(self, **attrs):
        self.attrs = attrs

    def show(self, dt, ax=None, **kwargs) -> Axes:  # or consistently a 2D list of axes
        raise NotImplementedError("Visualisation not implemented for this instrument")

    def annotate_positions(
        self, positions: list[Position], dt, ax, **kwargs
    ) -> Axes:  # or consistently a 2D list of axes
        raise NotImplementedError(
            "Annotating positions not implemented for this instrument"
        )


class Instrument(PlottableInstrument):
    def __init__(self, position: Position, rotation: list, **attrs):
        """Physical instruments with coordinate transforms and affiliated data loaders."""
        self.position = position
        if isinstance(rotation, int):
            rotation = [90, 0, rotation]
        self.rotation = rotation

        self.data_loader = None
        super().__init__(**attrs)

    # lla and xyz are functions of the instrument position, rather than the instrument itself
    # ead on the other hand is an instrument property, as it depends on the instrument rotation.

    def iead_to_gead(self, elevation, azimuth, dist):
        # convert to instrument-relative xyz coordinates "view coordinates"
        ixyz = ead_to_xyz(elevation, azimuth, dist)

        # apply the rotations in this order:
        R = rotation_matrix_i_to_g(self.rotation[0], self.rotation[1], self.rotation[2])

        # rotate around the z axis (roll)
        gxyz = R @ ixyz

        # convert to global ead
        global_elevation, global_azimuth, dist = xyz_to_ead(*gxyz)

        return global_elevation, global_azimuth, dist

    def gead_to_iead(self, elevation, azimuth, dist):
        gxyz = ead_to_xyz(elevation, azimuth, dist)
        # apply the inverse rotation
        R = rotation_matrix_i_to_g(
            self.rotation[0], self.rotation[1], self.rotation[2]
        ).T
        ixyz = R @ gxyz
        # convert to instrument-relative ead coordinates
        instrument_elevation, instrument_azimuth, dist = xyz_to_ead(*ixyz)

        return instrument_elevation, instrument_azimuth, dist

    def target_iead(self, target_position: Position):
        ead = self.position.target_ead(target_position)
        return self.gead_to_iead(*ead)

    def iead_to_lla(self, elevation, azimuth, dist):
        return self.position.ead_to_lla(*self.iead_to_gead(elevation, azimuth, dist))

    def iead_to_xyzworld(self, elevation, azimuth, dist):
        return ead_to_xyz(*self.iead_to_gead(elevation, azimuth, dist))

    def xyzworld_to_iead(self, target_x, target_y, target_z):
        gead = self.position.xyz_to_ead(target_x, target_y, target_z)
        return self.gead_to_iead(*gead)

    def initialise_data_loader(self):
        raise NotImplementedError(
            "Data loader initialisation not implemented for this instrument"
        )

    def get_data_time(self, dt):
        if self.data_loader is None:
            self.initialise_data_loader()
        return self.data_loader.get_data_time(dt)

    def show(self, dt, ax=None, **kwargs):
        if self.data_loader is None:
            self.initialise_data_loader()
        return self._show(dt, ax=ax, **kwargs)

    def _show(self, dt, ax=None, **kwargs):
        raise NotImplementedError("Visualisation not implemented for this instrument")

    def annotate_positions(self, positions: list[Position], ax, **kwargs):
        raise NotImplementedError(
            "Annotating positions not implemented for this instrument"
        )
