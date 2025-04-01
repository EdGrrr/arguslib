from arguslib.instruments.instruments import Instrument, Position
from arguslib.misc.plotting import plot_beam


import numpy as np
import pyart
import yaml


from pathlib import Path
from typing import override


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

        # elevs = radar.elevation["data"]

        elev_azi_start = radar.elevation["data"][0], radar.azimuth["data"][0]
        elev_azi_end = radar.elevation["data"][-1], radar.azimuth["data"][-1]
        plot_beam(self, self, elev_azi_start, dt=dt, ax=ax, **kwargs)
        plot_beam(self, self, elev_azi_end, dt=dt, ax=ax, **kwargs)
        # plot_rhi_beam(ax, elevs[0])
        # plot_rhi_beam(ax, elevs[-1])

        return ax

    def annotate_positions(
        self, positions, dt, ax, *args, plotting_method=None, **kwargs
    ):
        # project the positions to the xy plane...
        xlims = ax.get_xlim()

        furthest_distance = np.max(np.abs(xlims))

        if dt is None:
            raise ValueError(
                "dt must be provided for radar positions to get the azimuth"
            )

        eads = [self.target_iead(p) for p in positions]

        elevs = np.array([ead[0] for ead in eads])
        azimuths = np.array([ead[1] for ead in eads])
        dists = np.array([ead[2] for ead in eads])

        azimuth = self.data_loader.get_pyart_radar(dt).azimuth["data"][0]

        theta_seps = azimuths - azimuth
        offsets = dists * np.sin(
            np.deg2rad(theta_seps)
        )  # dist btwn object and radar plane
        # filtered to be less than 5km

        ys = dists * np.sin(np.deg2rad(elevs))
        xs = dists * np.cos(np.deg2rad(elevs)) * np.cos(np.deg2rad(theta_seps))

        plot_filter = (np.abs(offsets) < 5) & (xs < furthest_distance)
        if not plot_filter.any():
            return
        if plotting_method is None:
            ax.plot(xs[plot_filter], ys[plot_filter], *args, **kwargs)
        else:
            getattr(ax, plotting_method)(
                xs[plot_filter], ys[plot_filter], *args, **kwargs
            )

        ax.set_xlim(xlims)
