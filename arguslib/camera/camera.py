"""
Defines the Camera instrument, providing methods for geolocation and visualization of all-sky or perspective camera imagery.
"""

# %%
import datetime
from .calibration import PerspectiveProjection, Projection, unit
from ..instruments.instruments import (
    Instrument,
    Position,
    ead_to_xyz,
    xyz_to_ead,
    rotation_matrix_i_to_g,
)
from ..misc.geo import haversine, destination_point
from ..misc.plotting import (
    get_pixel_transform,
    make_camera_axes,
    plot_range_rings,
)
from .loaders import get_data_loader_class
from .locator import CameraData


import numpy as np
from tqdm import trange

from pathlib import Path
from typing_extensions import override
import importlib.resources as il_res
from ..config import resolve_config_resource

DEFAULT_CALIBRATION_NAME = "cobalt_3-7_calibration.yml"


def _resolve_calibration_resource(spec) -> object:
    """
    Resolve a calibration spec (None, str, or Path) using config.py helpers.
    """
    if spec is None:
        return resolve_config_resource(DEFAULT_CALIBRATION_NAME)
    return resolve_config_resource(spec)


class Camera(Instrument):
    """Represents a camera instrument, handling intrinsic and extrinsic calibration.

    This class provides the core functionality for a single camera, including
    loading configuration, converting between world coordinates and pixel
    coordinates, and rendering the camera's view on a plot. It supports both
    'allsky' fisheye cameras and standard 'perspective' cameras.

    Attributes:
        intrinsic (Projection): The intrinsic calibration model.
        scale_factor (float): A scaling factor applied to image dimensions.
        camera_type (str): The type of camera, e.g., 'allsky' or 'perspective'.
        image_size_px (np.ndarray): The dimensions of the camera image in pixels.
    """

    def __init__(
        self,
        intrinsic_calibration: Projection,
        image_size_px: list[int],
        camera_type: str,
        *args,
        scale_factor=1,
        time_offset_s=0.0,
        data_loader_class=None,
        **kwargs,
    ):
        self.intrinsic = intrinsic_calibration
        self.scale_factor = scale_factor
        self.camera_type = camera_type
        self.time_offset_s = time_offset_s

        self.image_size_px = image_size_px
        self.image_size_px = np.array(self.image_size_px) * scale_factor

        if data_loader_class is None:
            self._data_loader_class = CameraData
        else:
            self._data_loader_class = data_loader_class
        super().__init__(*args, **kwargs)

    @classmethod
    def from_filename(cls, filename, *args, **kwargs):
        """
        Accepts either a filesystem Path/str or an importlib.resources Traversable.
        """
        # Check for Traversable-like object (has open/joinpath) from importlib.resources
        if hasattr(filename, "open") and hasattr(filename, "joinpath"):
            # Use as_file to materialize on filesystem if needed
            with il_res.as_file(filename) as tmp_path:
                return cls(Projection.fromfile(tmp_path), *args, **kwargs)
        # Otherwise treat as filesystem path
        return cls(Projection.fromfile(Path(filename)), *args, **kwargs)

    @classmethod
    def from_config(
        cls, campaign, camstr, **kwargs
    ):  # TODO: make from_config a method of Instrument...?
        from arguslib.config import load_config

        cameras = load_config("cameras.yml")

        camera_config = cameras[campaign][camstr]

        loader_name = camera_config.get("data_loader", campaign)
        LoaderClass = get_data_loader_class(loader_name)

        kwargs["data_loader_class"] = LoaderClass

        if camera_config["camera_type"] == "allsky":
            # Resolve calibration file across config locations and packaged defaults
            calib_resource = _resolve_calibration_resource(
                camera_config.get("calibration_file")
            )
            # Propagate image_size_px if provided (expected order: [width, height])
            img_size = camera_config.get("image_size_px", None)

            time_offset_s = camera_config.get("time_offset_s", 0.0)

            kwargs = (
                {
                    "filename": calib_resource,
                    "position": Position(*camera_config["position"]),
                    "rotation": camera_config["rotation"],
                    "time_offset_s": time_offset_s,
                    "camera_type": camera_config["camera_type"],
                }
                | ({"image_size_px": img_size} if img_size is not None else {})
                | kwargs
            )

            kwargs |= {"campaign": campaign, "camstr": camstr}
            return cls.from_filename(**kwargs)

        else:  # Then this is a perspective camera
            # Propagate image_size_px if provided (expected order: [width, height])
            img_size = camera_config.get("image_size_px", None)

            kwargs |= {"campaign": campaign, "camstr": camstr}
            return cls(
                PerspectiveProjection(
                    camera_config["focal_lengths"],
                    camera_config["principal_point"],
                    camera_config.get("distortion_coeffs", None),
                ),
                position=Position(*camera_config["position"]),
                rotation=camera_config["rotation"],
                camera_type=camera_config["camera_type"],
                **({"image_size_px": img_size} if img_size is not None else {}),
                **kwargs,
            )

    # view here is in the camera space - we need to get the xyz in camera projection, not the world projection
    def target_pix(self, target_position: Position):
        """
        Calculates pixel coordinates for a target position.
        Vectorized to handle a list of positions.
        """
        # This call is now vectorized
        ieads = self.target_iead(target_position)

        # Handle single vs. multiple positions
        if ieads.ndim == 1:
            return self.iead_to_pix(*ieads)
        else:
            elevations, azimuths, dists = ieads[:, 0], ieads[:, 1], ieads[:, 2]
            # iead_to_pix is vector-ready
            return self.iead_to_pix(elevations, azimuths, dists)

    def pix_to_iead(self, pix_x, pix_y, distance=None, altitude=None):
        xyz = self.intrinsic.image_to_view(
            np.array([pix_x * self.scale_factor, pix_y * self.scale_factor]).T,
        )

        uv = unit(xyz)
        if distance is None and altitude is None:
            raise Exception("Must specify either distance or altitude")

        # This inner function calculates the EAD and performs the critical elevation conversion
        def calculate_and_convert_ead(dist):
            final_xyz = uv * dist
            # xyz_to_ead returns [elevation_from_plane, azimuth, distance]
            ead = xyz_to_ead(*final_xyz)
            # Elevation: 90° - angle_from_xy_plane = angle_from_z_axis
            ead[0] = 90.0 - ead[0]
            # Azimuth: shift so 0° is at image-top (instead of +X/right)
            # ead[1] = (ead[1] + 90.0) % 360.0
            return ead

        if distance is not None and altitude is None:
            return calculate_and_convert_ead(distance)
        elif distance is None and altitude is not None:
            R_i_to_g = rotation_matrix_i_to_g(*self.rotation)
            uv_global = R_i_to_g @ uv

            if uv_global[2] <= 1e-6:
                return np.array([np.nan, np.nan, np.nan])

            distance = (altitude - self.position.alt) / uv_global[2]
            return calculate_and_convert_ead(distance)
        else:
            raise ValueError("Cannot specify both distance and altitude.")

    def get_bearing_to_image_top(self):
        # Use instrument -Y (image-top) rotated to global to determine bearing
        R_i_to_g = rotation_matrix_i_to_g(*self.rotation)
        dir_i_top = np.array([0.0, 1.0, 0.0])  # image-top in instrument frame
        dir_g = unit(R_i_to_g @ dir_i_top).ravel()

        # Bearing from North, clockwise: atan2(East, North); +90° per your convention
        top_bearing_deg = (np.degrees(np.arctan2(dir_g[0], dir_g[1])) + 360.0) % 360.0
        # want to remove the +90.
        return top_bearing_deg

    def iead_to_pix(self, elevation, azimuth, dist=10):
        # Convert elevation: angle_from_xy_plane = 90 degrees - angle_from_z_axis
        elevation_for_xyz = 90.0 - elevation
        # Undo the 90° top-zero shift
        azimuth_for_xyz = (np.asarray(azimuth)) % 360.0

        view_vector = ead_to_xyz(
            elevation_for_xyz, azimuth_for_xyz, dist
        ).T  # ead_to_xyz seems to transpose them...

        # Transpose view_vector to (N, 3) before passing to the projection function
        return self.intrinsic.view_to_image(view_vector) / self.scale_factor

    def radar_beam(self, target_elevation, target_azimuth, radar):
        # TODO: this should maybe belong somewhere else. is it still used?
        positions = radar.beam(target_elevation, target_azimuth, [1, 10])
        pix = (
            np.array([self.target_pix(k) for k in positions.reshape(-1)])
            .reshape(2, -1, 2)
            .squeeze()
        )  # -1 and squeeze allows for either 5 or 1 position depending on radar beamwidth
        return pix

    @override
    def initialise_data_loader(self):
        from .locator import CameraData

        LoaderClass = self._data_loader_class

        self.data_loader = LoaderClass(self.attrs["campaign"], self.attrs["camstr"])

    @override
    def _show(
        self,
        dt,
        ax=None,
        fail_if_no_data=True,
        imshow_kw={},
        brightness_adjust=1.0,
        allow_timestamp_updates=True,
        **kwargs,
    ):
        """Renders the camera image for a given datetime on a Matplotlib axis.

        This method fetches the camera frame closest to the specified datetime,
        applies transformations for correct orientation and projection, and
        displays it. It can create a new plot or draw on an existing one.
        For all-sky cameras, this typically produces a polar plot.

        Args:
            dt (datetime.datetime): The timestamp for which to show the image.
            ax (matplotlib.axes.Axes, optional): The axis to plot on. If None,
                a new figure and axis are created. Defaults to None.
            fail_if_no_data (bool, optional): If True, raises a FileNotFoundError
                if no image data is available. If False, returns the axis without
                plotting an image. Defaults to True.
            imshow_kw (dict, optional): Keyword arguments passed to `ax.imshow`.
                Defaults to {}.
            brightness_adjust (float, optional): A factor to adjust image brightness.
                Values > 1 increase brightness. Defaults to 1.0.
            **kwargs:
                theta_behaviour (str): Controls the polar axis orientation.
                    Can be 'bearing' (North at top), 'pixels' (image y-axis up),
                    or 'unflipped_ordinal_aligned'.
                lr_flip (bool): Whether to flip the image left-to-right.

        Returns:
            matplotlib.axes.Axes: The axis on which the image was plotted.
        """
        defaults = {"theta_behaviour": "bearing", "lr_flip": True}

        if "theta_behaviour" in kwargs and "lr_flip" not in kwargs:
            # if theta_behaviour is set, assume lr_flip should be false, unless it's bearing
            if kwargs["theta_behaviour"] != "bearing":
                defaults["lr_flip"] = False

        if "lr_flip" in kwargs and "theta_behaviour" not in kwargs:
            # if lr_flip is set, but we are inferring theta_behaviour, assume it's bearing if flipped, ordinal aligned if not flipped
            if kwargs["lr_flip"]:
                defaults["theta_behaviour"] = "bearing"
            else:
                defaults["theta_behaviour"] = "unflipped_ordinal_aligned"

        kwargs = defaults | kwargs

        lr_flip = kwargs.pop("lr_flip")
        theta_behaviour = kwargs.pop("theta_behaviour")

        if ax is None:
            # If no axis is provided, make_camera_axes will create one,
            # potentially polar, and set theta_offset/direction.
            ax = make_camera_axes(
                self, theta_behaviour=theta_behaviour, dt=dt, **kwargs
            )
        else:
            # If an axis is provided, check if it's polar and apply theta settings.
            if hasattr(ax, "set_theta_zero_location"):  # It's a polar axis
                if theta_behaviour == "pixels":
                    ax.set_theta_offset(np.deg2rad(self.rotation[-1]))
                    # total_rotation_deg = self.rotation[1] + self.rotation[2]
                    # ax.set_theta_offset(np.deg2rad(total_rotation_deg))
                    # Ensure default direction if not specified or handled by lr_flip logic later
                    if not hasattr(
                        ax, "_theta_direction_set_by_camera_show"
                    ):  # Avoid double-setting
                        ax.set_theta_direction(1)  # Example default, adjust if needed
                        ax._theta_direction_set_by_camera_show = True
                elif theta_behaviour == "bearing":
                    ax.set_theta_offset(np.pi / 2)
                    ax.set_theta_direction(-1)
                elif theta_behaviour == "unflipped_ordinal_aligned":
                    ax.set_theta_offset(np.pi / 2)
                    # ax.set_theta_direction(1) # Example default
        is_polar = hasattr(ax, "set_theta_zero_location")

        # if polar axes, assume it's a camera axes with
        transform = get_pixel_transform(self, ax, lr_flip=lr_flip)
        # if is_polar:
        #     transform = get_pixel_transform(self, ax, lr_flip=lr_flip)
        # else:
        #     transform = ax.transData

        ax.transData = transform

        if is_polar:
            ax.set_rticks([])

        ax.grid(False)
        if is_polar:
            ax.set_rticks([])
        ax.grid(False)

        plot_range_rings(self, self, dt, ax=ax)

        try:
            img, timestamp = self.get_data_time(dt, return_timestamp=True)
            if allow_timestamp_updates:
                ax.get_figure().timestamp = timestamp + datetime.timedelta(
                    seconds=self.time_offset_s
                )
        except FileNotFoundError as e:
            if fail_if_no_data:
                raise e
            else:
                return ax

        img = np.clip(np.uint16(img) * brightness_adjust, 0, 255).astype(np.uint8)

        # Draw the image in pixel coordinates with correct extent and aspect
        # Note: origin='upper' because pixel (0,0) is top-left
        ax.imshow(
            img[:, :, ::-1],
            origin="upper",
            aspect="equal",
            **imshow_kw,
        )

        return ax

    @override
    def annotate_positions(
        self, positions, dt, ax, *args, plotting_method=None, max_range_km=90, **kwargs
    ):
        """
        Annotates one or more geographical positions on the camera image.
        This version is vectorized for performance.
        """
        if not isinstance(positions, (list, tuple, np.ndarray)) or len(positions) == 0:
            return ax

        # --- Vectorized Geolocation ---
        # NOTE: This assumes `positions` can be iterated over for list comprehension.
        # This is safe for lists, tuples, and numpy arrays.
        target_lons = np.array([p.lon for p in positions])
        target_lats = np.array([p.lat for p in positions])

        dists = haversine(
            self.position.lon, self.position.lat, target_lons, target_lats
        )
        ieads = self.target_iead(positions)

        if ieads.ndim == 1:
            ieads = ieads.reshape(1, -1)
            dists = np.array([dists]) if np.isscalar(dists) else dists

        # --- Vectorized Pixel Conversion ---
        pl_track = self.target_pix(positions)

        # --- Filtering and Plotting ---
        behind_camera = ieads[:, 0] > 90
        in_range = dists < max_range_km
        plot_filter = in_range & ~behind_camera

        if not np.any(plot_filter):
            return ax

        filtered_track = pl_track[plot_filter]

        if hasattr(ax, "set_xlim"):
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

        if plotting_method is None:
            ax.plot(filtered_track[:, 0], filtered_track[:, 1], *args, **kwargs)
        else:
            getattr(ax, plotting_method)(
                filtered_track[:, 0], filtered_track[:, 1], *args, **kwargs
            )

        if hasattr(ax, "set_xlim"):
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        return ax

    def distance_calibration_img(self, height_km=10):
        X, Y = np.meshgrid(
            np.arange(self.image_size_px[0]), np.arange(self.image_size_px[1])
        )
        grid = np.full_like(X, np.nan)

        for xx in trange(self.image_size_px[0]):
            for yy in range(self.image_size_px[1]):
                elev, _, dist = self.pix_to_iead(
                    X[xx, yy].item(), Y[xx, yy].item(), altitude=height_km
                )
                dist_in_plane = dist * np.tan(np.deg2rad(elev))
                grid[xx, yy] = dist_in_plane
        return grid


# %%
