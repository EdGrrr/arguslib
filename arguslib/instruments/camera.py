from arguslib.instruments.calibration import Projection, unit
from arguslib.instruments.instruments import (
    Instrument,
    Position,
    default_calibration_file,
)
from arguslib.misc.geo import haversine
from arguslib.misc.plotting import (
    get_pixel_transform,
    make_camera_axes,
    plot_range_rings,
)


import numpy as np
import yaml


from pathlib import Path
from typing import override


class Camera(Instrument):
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
        if camera_config["calibration_file"] is None:
            camera_config["calibration_file"] = default_calibration_file

        kwargs = {
            "filename": Path(camera_config["calibration_file"]).expanduser().absolute(),
            "position": Position(*camera_config["position"]),
            "rotation": camera_config["rotation"],
        } | kwargs
        # will ignore config if kwargs contains any of the keys in camera_config

        kwargs |= {"campaign": campaign, "camstr": camstr}
        return cls.from_filename(**kwargs)

    # view here is in the camera space - we need to get the xyz in camera projection, not the world projection
    def target_pix(self, target_position: Position):
        return self.iead_to_pix(*self.target_iead(target_position))

    def pix_to_iead(self, pix_x, pix_y, distance):
        xyz = self.intrinsic.image_to_view(
            [pix_x * self.scale_factor, pix_y * self.scale_factor]
        )
        # xyz for projection, convert to ead for projection
        return self.position.xyz_to_ead(distance * unit(xyz))

    def iead_to_pix(self, elevation, azimuth, dist=10):
        return (
            self.intrinsic.view_to_image(
                self.position.ead_to_xyz(elevation, azimuth, dist)
            )
            / self.scale_factor
        )

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
        from ..video import CameraData

        self.data_loader = CameraData(self.attrs["campaign"], self.attrs["camstr"])

    @override
    def _show(self, dt, ax=None, **kwargs):
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
            ax = make_camera_axes(self, theta_behaviour=theta_behaviour, **kwargs)

        is_polar = hasattr(ax, "set_theta_zero_location")

        # if polar axes, assume it's a camera axes with
        if is_polar:
            transform = get_pixel_transform(self, ax, lr_flip=lr_flip)
        else:
            transform = ax.transData

        ax.transData = transform

        img = self.get_data_time(dt)
        ax.imshow(img[:, :, ::-1], origin="upper")
        plot_range_rings(self, ax=ax)

        if is_polar:
            ax.set_rticks([])
        return ax

    @override
    def annotate_positions(
        self, positions, dt, ax, *args, plotting_method=None, **kwargs
    ):
        # TODO: this should take dt into account, mostly because the calibration may change for the same camera at different times...
        lats = [p.lat for p in positions]
        lons = [p.lon for p in positions]

        dists = np.array(
            [
                haversine(self.position.lon, self.position.lat, lon, lat)
                for lon, lat in zip(lons, lats)
            ]
        )

        if (dists[~np.isnan(dists)] > 90).all():
            return

        pl_track = np.array([self.target_pix(p) for p in positions])
        if plotting_method is None:
            ax.plot(
                pl_track.T[0][dists < 90], pl_track.T[1][dists < 90], *args, **kwargs
            )
        else:
            getattr(ax, plotting_method)(
                pl_track.T[0][dists < 90],
                pl_track.T[1][dists < 90],
                *args,
                **kwargs,
            )
        return ax
