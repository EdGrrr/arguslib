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

    # Elevation rotation (Rₑ) - global x-axis
    R_e = np.array(
        [
            [1, 0, 0],
            [0, np.cos(elevation), -np.sin(elevation)],
            [0, np.sin(elevation), np.cos(elevation)],
        ]
    )

    # Roll rotation (Rᵣ) - local z-axis
    R_r = np.array(
        [[np.cos(roll), np.sin(roll), 0], [-np.sin(roll), np.cos(roll), 0], [0, 0, 1]]
    )

    # Combined rotation matrix
    R = R_a @ R_e @ R_r
    return R


class Position:
    """Represents a geographical point in longitude, latitude, and altitude.

    This class provides the core functionality for geolocation within ``arguslib``.
    It assumes a locally flat Earth for calculations, which is suitable for the
    typical ranges of ground-based instruments.

    All distance and altitude units are in kilometers.

    Attributes:
        lon (float): Longitude in degrees.
        lat (float): Latitude in degrees.
        alt (float): Altitude in kilometers.
    """

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
        """Calculates the elevation, azimuth, and distance to a target position.

        Args:
            target_position (Position): The target position.

        Returns:
            np.ndarray: An array containing [elevation, azimuth, distance].
                        Elevation and azimuth are in degrees, distance is in km.
        """
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
        """Calculates the relative X, Y, Z coordinates to a target position.

        Args:
            target_position (Position): The target position.

        Returns:
            np.ndarray: An array containing [X, Y, Z] distances in km, where
                        X is East-West, Y is North-South, and Z is vertical.
        """
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
        """Converts elevation, azimuth, and distance to a new Position.

        Args:
            target_elevation (float): Elevation from the horizontal in degrees.
            target_azimuth (float): Azimuth from North in degrees.
            target_distance (float): Distance in km.

        Returns:
            Position: The new geographical position.
        """
        return self.xyz_to_lla(
            *ead_to_xyz(target_elevation, target_azimuth, target_distance)
        )

    def xyz_to_lla(self, target_x, target_y, target_z):
        # Assume Earth is locally flat, but also spherical
        lon = self.lon + target_x / (111.111 * np.cos(self.lat))
        lat = self.lat + target_y / (111.111)
        return Position(lon, lat, target_z + self.alt)

    def __repr__(self):
        return f"{self.lon:.3f} {self.lat:.3f} {self.alt:.3f}"


class PlottableInstrument:
    """An abstract base class for any object that can be visualized.

    This class defines the core API for visualization within ``arguslib``. Any
    instrument or interface that can be displayed on a Matplotlib plot or
    rendered into an image should implement these methods.

    Attributes:
        attrs (dict): A dictionary of attributes describing the instrument.
    """

    def __init__(self, **attrs):
        self.attrs = attrs

    def show(self, dt, ax=None, **kwargs) -> Axes:  # or consistently a 2D list of axes
        """Renders the instrument's primary visualization for a given time.

        Args:
            dt (datetime.datetime): The timestamp for the visualization.
            ax (matplotlib.axes.Axes, optional): The axis to plot on. If None,
                a new figure and axis are typically created. For some instruments
                that render directly to images (like `DirectCamera`), this may
                be ignored or required to be None.
            **kwargs: Additional keyword arguments specific to the instrument's
                plotting implementation.

        Returns:
            matplotlib.axes.Axes: The axis (or array of axes) on which the
            visualization was drawn. May be None for non-Matplotlib backends.
        """
        raise NotImplementedError("Visualisation not implemented for this instrument")

    def annotate_positions(
        self, positions: list[Position], dt, ax, **kwargs
    ) -> Axes:  # or consistently a 2D list of axes
        """Annotates one or more geographical positions on the instrument's plot.

        Args:
            positions (list[Position]): A list of `Position` objects to annotate.
            dt (datetime.datetime): The datetime for which the view is valid.
            ax (matplotlib.axes.Axes): The axis (or axes) to plot on.
            **kwargs: Keyword arguments passed to the underlying plotting function
                (e.g., `ax.plot` or `ax.scatter`).

        Returns:
            matplotlib.axes.Axes: The updated axis (or axes).
        """
        raise NotImplementedError(
            "Annotating positions not implemented for this instrument"
        )


class Instrument(PlottableInstrument):
    """Represents a physical instrument with a specific location and orientation.

    This class extends `PlottableInstrument` by adding coordinate transformation
    capabilities. It holds a `Position` and a `rotation` vector, allowing it
    to convert between its own local coordinate system (instrument-relative
    elevation/azimuth) and the global coordinate system (world-relative
    elevation/azimuth).

    It also introduces the concept of a `data_loader`, which is responsible
    for fetching the instrument's data (e.g., images, scans) for a given time.

    Attributes:
        position (Position): The geographical location of the instrument.
        rotation (list): A list of three angles [elevation, azimuth, roll]
            defining the instrument's orientation.
        data_loader: An object responsible for loading data for this instrument.
            It is typically initialized on the first call to `get_data_time` or `show`.
    """

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
        R = rotation_matrix_i_to_g(self.rotation[0], self.rotation[1], self.rotation[2])
        inv_R = np.linalg.inv(R)
        ixyz = inv_R @ gxyz
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

    def get_data_time(self, dt, **kwargs):
        if self.data_loader is None:
            self.initialise_data_loader()
        return self.data_loader.get_data_time(dt, **kwargs)

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
