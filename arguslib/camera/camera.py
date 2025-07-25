"""
Defines the Camera instrument, providing methods for geolocation and visualization of all-sky or perspective camera imagery.
"""
# %%
import datetime
from .calibration import PerspectiveProjection, Projection, unit
from ..instruments.instruments import (
    Instrument,
    Position,
    default_calibration_file,
    ead_to_xyz,
    xyz_to_ead,
    rotation_matrix_i_to_g
)
from ..misc.geo import haversine
from ..misc.plotting import (
    get_pixel_transform,
    make_camera_axes,
    plot_range_rings,
)


import numpy as np
from tqdm import trange

from pathlib import Path
from typing import override


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
        *args,
        scale_factor=1,
        camera_type="allsky",
        **kwargs,
    ):
        self.intrinsic = intrinsic_calibration
        self.scale_factor = scale_factor
        self.camera_type = camera_type

        self.image_size_px = (
            [4608, 2592] if self.camera_type == "perspective" else [3040, 3040]
        )
        self.image_size_px = np.array(self.image_size_px) * scale_factor
        super().__init__(*args, **kwargs)

    @classmethod
    def from_filename(cls, filename, *args, **kwargs):
        return cls(Projection.fromfile(filename), *args, **kwargs)

    @classmethod
    def from_config(
        cls, campaign, camstr, **kwargs
    ):  # TODO: make from_config a method of Instrument...?
        import yaml

        # look for a config file - try ~/.config/arguslib/cameras.yml, then ~/.arguslib/cameras.yml, then /etc/arguslib/cameras.yml
        # read from all that exists
        config_paths = [
            Path("~/.config/arguslib/cameras.yml").expanduser(),
            Path("~/.arguslib/cameras.yml").expanduser(),
            Path("/etc/arguslib/cameras.yml"),
        ]
        configs = []
        for config_file in config_paths:
            if not config_file.exists():
                continue
            with open(config_file, "r") as f:
                configs.append(yaml.safe_load(f))

        if not configs:
            raise FileNotFoundError("No camera configuration file found")

        cameras = {}
        for config in configs[::-1]:
            cameras.update(config)

        camera_config = cameras[campaign][camstr]

        if "calibration_file" in camera_config:  # Then this is an allsky camera
            if camera_config["calibration_file"] is None:
                camera_config["calibration_file"] = default_calibration_file
            kwargs = {
                "filename": Path(camera_config["calibration_file"])
                .expanduser()
                .absolute(),
                "position": Position(*camera_config["position"]),
                "rotation": camera_config["rotation"],
            } | kwargs  # will ignore config if kwargs contains any of the keys in camera_config

            kwargs |= {"campaign": campaign, "camstr": camstr}
            return cls.from_filename(**kwargs)

        else:  # Then this is a perspective camera
            kwargs |= {"campaign": campaign, "camstr": camstr}
            return cls(
                PerspectiveProjection(
                    camera_config["focal_lengths"],
                    camera_config["principal_point"],
                    camera_config.get("distortion_coeffs", None),
                ),
                position=Position(*camera_config["position"]),
                rotation=camera_config["rotation"],
                camera_type="perspective",
                **kwargs,
            )

    # view here is in the camera space - we need to get the xyz in camera projection, not the world projection
    def target_pix(self, target_position: Position):
        return self.iead_to_pix(*self.target_iead(target_position))


    def pix_to_iead(self, pix_x, pix_y, distance=None, altitude=None):
        xyz = self.intrinsic.image_to_view(
            [pix_x * self.scale_factor, pix_y * self.scale_factor]
        )
        
        uv = unit(xyz)
        if distance is None and altitude is None:
            raise Exception("Must specify either distance or altitude")
        
        # This inner function calculates the EAD and performs the critical elevation conversion
        def calculate_and_convert_ead(dist):
            final_xyz = uv * dist
            # xyz_to_ead returns [elevation_from_plane, azimuth, distance]
            ead = xyz_to_ead(*final_xyz)
            # Convert elevation: 90 degrees - angle_from_xy_plane = angle_from_z_axis
            ead[0] = 90.0 - ead[0]
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

    def iead_to_pix(self, elevation, azimuth, dist=10):
        # Convert elevation: angle_from_xy_plane = 90 degrees - angle_from_z_axis
        elevation_for_xyz = 90.0 - elevation
        
        view_vector = ead_to_xyz(elevation_for_xyz, azimuth, dist)
        
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
        from ..camera import CameraData

        self.data_loader = CameraData(self.attrs["campaign"], self.attrs["camstr"])

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
            ax = make_camera_axes(self, theta_behaviour=theta_behaviour, dt=dt, **kwargs)
        else:
            # If an axis is provided, check if it's polar and apply theta settings.
            if hasattr(ax, "set_theta_zero_location"): # It's a polar axis
                if theta_behaviour == "pixels":
                    ax.set_theta_offset(np.deg2rad(self.rotation[-1]))
                    # total_rotation_deg = self.rotation[1] + self.rotation[2]
                    # ax.set_theta_offset(np.deg2rad(total_rotation_deg))
                    # Ensure default direction if not specified or handled by lr_flip logic later
                    if not hasattr(ax, '_theta_direction_set_by_camera_show'): # Avoid double-setting
                        ax.set_theta_direction(1) # Example default, adjust if needed
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

        
        plot_range_rings(self, dt, ax=ax)

        try:
            img, timestamp = self.get_data_time(dt, return_timestamp=True)
            if allow_timestamp_updates:
                ax.get_figure().timestamp = timestamp
        except FileNotFoundError as e:
            if fail_if_no_data:
                raise e
            else:
                return ax

        # new_vmin = np.round(brightness_adjust * 255)
        img = np.clip(np.uint16(img) * brightness_adjust, 0, 255).astype(np.uint8)
        # img[img <= new_vmin] = new_vmin
        # img = np.round(255 * (img.astype(float) - new_vmin) / (255 - new_vmin)).astype(
        #     int
        # )

        ax.imshow(img[:, :, ::-1], origin="upper", **imshow_kw)

        return ax

    @override
    def annotate_positions(
        self, positions, dt, ax, *args, plotting_method=None, max_range_km=90, **kwargs
    ):
        """Annotates one or more geographical positions on the camera image.

        Converts a list of `Position` objects (lat, lon, alt) into pixel
        coordinates for the given camera view and plots them on the specified
        axis.

        Args:
            positions (list[Position]): A list of `Position` objects to annotate.
            dt (datetime.datetime): The datetime for which the camera view is valid.
                This can affect calibration.
            ax (matplotlib.axes.Axes): The axis to plot the annotations on.
            *args: Positional arguments passed to the plotting function (e.g., `ax.plot`).
            plotting_method (str, optional): The name of the Matplotlib plotting
                method to use (e.g., 'scatter', 'plot'). If None, `ax.plot` is used.
                Defaults to None.
            max_range_km (float, optional): The maximum horizontal distance (in km)
                from the camera at which to plot positions. Defaults to 90.
            **kwargs: Keyword arguments passed to the plotting function.
        """
        # TODO: this should take dt into account, mostly because the calibration may change for the same camera at different times...
        lats = [p.lat for p in positions]
        lons = [p.lon for p in positions]

        dists = np.array(
            [
                haversine(self.position.lon, self.position.lat, lon, lat)
                for lon, lat in zip(lons, lats)
            ]
        )

        behind_camera = np.array([self.target_iead(p)[0] > 90 for p in positions])

        if (dists[~np.isnan(dists)] > max_range_km).all():
            return

        pl_track = np.array([self.target_pix(p) for p in positions])

        # if nonpolar axes, need to preserve the limits
        if hasattr(ax, "set_xlim"):
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

        if plotting_method is None:
            ax.plot(
                pl_track.T[0][(dists < max_range_km) & ~behind_camera],
                pl_track.T[1][(dists < max_range_km) & ~behind_camera],
                *args,
                **kwargs,
            )
        else:
            getattr(ax, plotting_method)(
                pl_track.T[0][(dists < max_range_km) & ~behind_camera],
                pl_track.T[1][(dists < max_range_km) & ~behind_camera],
                *args,
                **kwargs,
            )

        if hasattr(ax, "set_xlim"):
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            
        return ax
    
    def distance_calibration_img(self, height_km=10):
        X, Y = np.meshgrid(np.arange(self.image_size_px[0]), np.arange(self.image_size_px[1]))
        grid = np.full_like(X, np.nan)
        
        for xx in trange(self.image_size_px[0]):
            for yy in range(self.image_size_px[1]):
                elev, _, dist = self.pix_to_iead(X[xx, yy].item(), Y[xx, yy].item(), altitude=height_km)
                dist_in_plane = dist * np.tan(np.deg2rad(elev))
                grid[xx, yy] = dist_in_plane
        return grid
# %%
