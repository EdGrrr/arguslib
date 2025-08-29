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
    # Clamp the argument to arcsin to avoid NaN from floating point errors
    elevation_ratio = np.clip(target_z / target_distance, -1.0, 1.0)
    target_elevation = np.rad2deg(np.arcsin(elevation_ratio))
    target_azimuth = np.rad2deg(np.arctan2(target_x, target_y))
    return np.array([target_elevation, target_azimuth, target_distance])



def rotation_matrix_i_to_g(elevation, azimuth, roll):
    """
    Creates a rotation matrix to transform coordinates from the instrument's
    local frame (i) to the global frame (g).

    The instrument frame is defined as:
    - Z-axis: Along the optical axis (pointing direction).
    - Y-axis: "Up" direction of the image sensor.
    - X-axis: "Right" direction of the image sensor.

    The global frame is an East-North-Up (ENU) system.

    Args:
        elevation (float): The instrument's pointing elevation in degrees from the horizon.
        azimuth (float): The instrument's pointing azimuth in degrees from North.
        roll (float): The instrument's roll in degrees around its optical axis.

    Returns:
        np.ndarray: The 3x3 rotation matrix R_i_to_g.
    """
    el_rad = np.deg2rad(elevation)
    az_rad = np.deg2rad(azimuth)
    roll_rad = np.deg2rad(roll)

    # General case for non-vertical cameras
    if np.abs(el_rad - np.pi / 2) > 1e-6:
        # Define the instrument's optical axis (the local z-axis)
        i_z = np.array([
            np.cos(el_rad) * np.sin(az_rad),
            np.cos(el_rad) * np.cos(az_rad),
            np.sin(el_rad)
        ])
        
        # Define the instrument's "right" vector (local x-axis) before roll
        g_up = np.array([0, 0, 1])
        i_x_preroll = np.cross(i_z, g_up)
        # Handle nadir-pointing case
        if np.linalg.norm(i_x_preroll) < 1e-6:
            i_x_preroll = np.array([1, 0, 0])
        i_x_preroll /= np.linalg.norm(i_x_preroll)

        # Define the instrument's "up" vector (local y-axis) before roll
        i_y_preroll = np.cross(i_z, i_x_preroll)
    
    # Special case for a vertically-pointing (zenith) camera
    else:
        # The optical axis is straight up.
        i_z = np.array([0, 0, 1])
        
        # The pre-roll orientation is fixed to North-up, East-right, IGNORING azimuth.
        i_x_preroll = np.array([1, 0, 0])  # East
        i_y_preroll = np.array([0, 1, 0])  # North
    
    # Apply roll to the pre-roll basis vectors to get the final orientation
    cos_r = np.cos(roll_rad)
    sin_r = np.sin(roll_rad)

    # **CORRECTED SIGNS:** This implements a clockwise rotation to match the original system
    i_x_final = i_x_preroll * cos_r - i_y_preroll * sin_r
    i_y_final = i_x_preroll * sin_r + i_y_preroll * cos_r

    # Construct the rotation matrix from the final basis vectors
    R = np.column_stack([i_x_final, i_y_final, i_z])
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
        """
        Calculates the elevation, azimuth, and distance to a target position.
        Vectorized to handle a list of target positions.
        """
        # Check if input is a list or a single object
        is_list = isinstance(target_position, (list, tuple, np.ndarray)) and len(target_position) > 0
        
        targets = target_position if is_list else [target_position]

        # Extract target coordinates into numpy arrays
        target_lons = np.array([p.lon for p in targets])
        target_lats = np.array([p.lat for p in targets])
        target_alts = np.array([p.alt for p in targets])

        # Perform vectorized calculations
        distance = haversine(self.lon, self.lat, target_lons, target_lats)
        alt_diff = target_alts - self.alt
        
        # Avoid division by zero for elevation calculation when distance is zero
        distance_safe = np.where(distance == 0, 1e-9, distance)
        
        target_distance = np.sqrt(distance**2 + alt_diff**2)
        target_elevation = np.rad2deg(np.arctan2(alt_diff, distance_safe))
        target_azimuth = bearing(self.lon, self.lat, target_lons, target_lats)
        
        # Stack results into an (N, 3) array
        result = np.vstack([target_elevation, target_azimuth, target_distance]).T
        
        # Return a single array or the full array to match input type
        return result[0] if not is_list else result

    # XYZ
    #  - X East-West distace in km
    #  - Y North-south distance in km
    #  - Z altitude in km

    def target_xyz(self, target_position):
        """
        Calculates the relative X, Y, Z coordinates to a target position.
        Vectorized to handle a list of target positions.
        """
        # Check if input is a list or a single object
        is_list = isinstance(target_position, (list, tuple, np.ndarray)) and len(target_position) > 0
        
        targets = target_position if is_list else [target_position]
    
        # Extract target coordinates into numpy arrays
        target_lons = np.array([p.lon for p in targets])
        target_lats = np.array([p.lat for p in targets])
        target_alts = np.array([p.alt for p in targets])

        # Perform vectorized calculations
        distance = haversine(self.lon, self.lat, target_lons, target_lats)
        target_z = target_alts - self.alt
        target_bearing = bearing(self.lon, self.lat, target_lons, target_lats)
        
        target_x = distance * np.sin(np.deg2rad(target_bearing))
        target_y = distance * np.cos(np.deg2rad(target_bearing))

        # Stack results into an (N, 3) array
        result = np.vstack([target_x, target_y, target_z]).T
        
        # Return a single array or the full array to match input type
        return result[0] if not is_list else result

    def ead_to_lla(self, target_elevation, target_azimuth, target_distance):
        """
        Converts elevation, azimuth, and distance to a new Position
        using a spherical Earth model. Handles elevations > 90 degrees.
        """
        from ..misc.geo import destination_point

        # Normalize elevation and azimuth for cases where elevation > 90 degrees.
        # This maps "over-the-top" coordinates to their forward-facing equivalent.
        if target_elevation > 90.0:
            target_azimuth = (target_azimuth + 180.0) % 360.0
            target_elevation = 180.0 - target_elevation

        # Calculate the horizontal distance traveled over the Earth's surface
        horiz_dist_km = target_distance * np.cos(np.deg2rad(target_elevation))
        
        # Calculate the change in altitude
        alt_diff_km = target_distance * np.sin(np.deg2rad(target_elevation))
        
        # Calculate the new lat/lon using the accurate spherical model
        new_lon, new_lat = destination_point(self.lon, self.lat, target_azimuth, horiz_dist_km)
        
        # Calculate the new altitude
        new_alt = self.alt + alt_diff_km
        
        return Position(new_lon.item(), new_lat.item(), new_alt.item())


    def xyz_to_lla(self, target_x, target_y, target_z):
        # Assume Earth is locally flat, but also spherical
        lon = self.lon + target_x / (111.111 * np.cos(np.deg2rad(self.lat)))
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
        elevation_from_image = 90.0 - elevation

        # convert to instrument-relative xyz coordinates "view coordinates"
        ixyz = ead_to_xyz(elevation_from_image, azimuth, dist)

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
        
        # Convert elevation from 'angle from instrument horizon' 
        # to 'angle from optical axis' (90 - elev)
        instrument_elevation = 90 - instrument_elevation

        return instrument_elevation, instrument_azimuth, dist

    def target_iead(self, target_position: Position):
        """
        Calculates instrument-relative EAD from a target position.
        Vectorized to handle a list of target positions.
        """
        # This call is now vectorized
        eads = self.position.target_ead(target_position)

        # Handle both single (1D array) and multiple (2D array) positions
        if eads.ndim == 1:
            return np.array(self.gead_to_iead(*eads))
        else:
            elevations, azimuths, dists = eads[:, 0], eads[:, 1], eads[:, 2]
            # gead_to_iead is already vector-ready due to numpy operations
            iead_tuple = self.gead_to_iead(elevations, azimuths, dists)
            # Return as a stacked (N, 3) numpy array
            return np.vstack(iead_tuple).T

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