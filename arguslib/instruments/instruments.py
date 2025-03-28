from typing import override

from arguslib.arguslib.misc.plotting import (
    get_pixel_transform,
    make_camera_axes,
    plot_range_rings,
    plot_rhi_beam,
)
from .calibration import Projection, unit
from ..misc.geo import haversine, bearing
import numpy as np
from pathlib import Path

default_calibration_file = (
    Path(__file__).parent.parent / "instruments/cam1_calibration.yml"
)


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

    def ead_to_xyz(self, target_elevation, target_azimuth, target_distance):
        target_z = target_distance * np.sin(np.deg2rad(target_elevation))
        target_hor_range = target_distance * np.cos(np.deg2rad(target_elevation))
        target_x = target_hor_range * np.sin(np.deg2rad(target_azimuth))
        target_y = target_hor_range * np.cos(np.deg2rad(target_azimuth))
        return np.array([target_x, target_y, target_z])

    def ead_to_lla(self, target_elevation, target_azimuth, target_distance):
        return self.xyz_to_lla(
            *self.ead_to_xyz(target_elevation, target_azimuth, target_distance)
        )

    def xyz_to_ead(self, target_x, target_y, target_z):
        target_distance = np.sqrt(target_x**2 + target_y**2 + target_z**2)
        target_elevation = np.rad2deg(
            np.arctan2(target_z, np.sqrt(target_x**2 + target_y**2))
        )
        target_azimuth = np.rad2deg(np.arctan2(target_x, target_y))
        return np.array([target_elevation, target_azimuth, target_distance])

    def xyz_to_lla(self, target_x, target_y, target_z):
        lon = self.lon + target_x / (111.111 * np.cos(self.lat))
        lat = self.lat + target_y / (111.111)
        return Position(lon, lat, target_z + self.alt)

    def __repr__(self):
        return f"{self.lon:.3f} {self.lat:.3f} {self.alt:.3f}"


class Instrument:
    def __init__(self, position: Position, rotation: list, **attrs):
        self.position = position
        self.rotation = rotation

        self.attrs = attrs
        self.data_loader = None

    # lla and xyz are functions of the instrument position, rather than the instrument itself
    # ead on the other hand is an instrument property, as it depends on the instrument rotation.

    def iead_to_gead(self, elevation, azimuth, dist):
        """Instrument El, az, dist to global values, based on an
        instrument facing directly up/north.  For the allsky cameras,
        this is probably more straightforwrd. A little trickier for
        perspective cameras.

        """
        return elevation, azimuth + self.rotation, dist

    def gead_to_iead(self, elevation, azimuth, dist):
        return elevation, azimuth - self.rotation, dist

    def target_iead(self, target_position: Position):
        ead = self.position.target_ead(target_position)
        return self.gead_to_iead(*ead)

    def iead_to_lla(self, elevation, azimuth, dist):
        return self.position.ead_to_lla(*self.iead_to_gead(elevation, azimuth, dist))

    def iead_to_xyzworld(self, elevation, azimuth, dist):
        return self.position.ead_to_xyz(*self.iead_to_gead(elevation, azimuth, dist))

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


class Radar(Instrument):
    def __init__(self, beamwidth, *args, **kwargs):
        self.beamwidth = beamwidth
        super().__init__(*args, **kwargs)

    def beam(self, radar_elevation, radar_azimuth, radar_distances):
        if self.beamwidth:
            return np.array(
                [
                    [
                        self.iead_to_lla(radar_elevation, radar_azimuth, radar_dist),
                        self.iead_to_lla(
                            radar_elevation + self.beamwidth / 2,
                            radar_azimuth,
                            radar_dist,
                        ),
                        self.iead_to_lla(
                            radar_elevation,
                            radar_azimuth + self.beamwidth / 2,
                            radar_dist,
                        ),
                        self.iead_to_lla(
                            radar_elevation - self.beamwidth / 2,
                            radar_azimuth,
                            radar_dist,
                        ),
                        self.iead_to_lla(
                            radar_elevation,
                            radar_azimuth - self.beamwidth / 2,
                            radar_dist,
                        ),
                    ]
                    for radar_dist in radar_distances
                ]
            )
        else:
            return np.array(
                [
                    [self.iead_to_lla(radar_elevation, radar_azimuth, radar_dist)]
                    for radar_dist in radar_distances
                ]
            )

    @classmethod
    def from_config(
        cls, campaign, **kwargs
    ):  # TODO: make from_config a method of Instrument...?
        import yaml

        # look for a config file - try ~/.config/arguslib/radars.yml, then ~/.arguslib/radars.yml, then /etc/arguslib/radars.yml
        # read from all that exists
        config_paths = [
            Path("~/.config/arguslib/radars.yml").expanduser(),
            Path("~/.arguslib/radars.yml").expanduser(),
            Path("/etc/arguslib/radars.yml"),
        ]
        configs = []
        for config_file in config_paths:
            if not config_file.exists():
                continue
            with open(config_file, "r") as f:
                configs.append(yaml.safe_load(f))

        if not configs:
            raise FileNotFoundError("No radar configuration file found")

        radars = {}
        for config in configs[::-1]:
            radars.update(config)

        radar_config = radars[campaign]

        kwargs = {
            "beamwidth": radar_config["beamwidth"],
            "position": Position(*radar_config["position"]),
            "rotation": radar_config["rotation"],
        } | kwargs
        # will ignore config if kwargs contains any of the keys in camera_config

        kwargs |= {"campaign": campaign}

        return cls(**kwargs)

    @override
    def initialise_data_loader(self):
        from ..radar import RadarData

        self.data_loader = RadarData(self.attrs["campaign"], "rhi")

    @override
    def _show(self, dt, var, ax=None, **kwargs):
        import pyart
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        radar = self.data_loader.get_pyart_radar(dt)
        display = pyart.graph.RadarDisplay(radar)

        kwargs = {
            "vmin": -60,
            "vmax": 40,
        } | kwargs
        display.plot(var, ax=ax, **kwargs)
        ax.set_aspect("equal")

        elevs = radar.elevation["data"]
        plot_rhi_beam(ax, elevs[0])
        plot_rhi_beam(ax, elevs[-1])

        return ax
