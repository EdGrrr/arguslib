from arguslib.instruments.instruments import Instrument, Position
from arguslib.misc.plotting import TimestampedFigure, plot_beam
from arguslib.misc.interpolation import interpolate_to_intersection

import datetime
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
                        self.position.ead_to_lla(radar_elevation, radar_azimuth, radar_dist),
                        self.position.ead_to_lla(
                            radar_elevation + self.beamwidth / 2,
                            radar_azimuth,
                            radar_dist,
                        ),
                        self.position.ead_to_lla(
                            radar_elevation,
                            radar_azimuth + self.beamwidth / 2,
                            radar_dist,
                        ),
                        self.position.ead_to_lla(
                            radar_elevation - self.beamwidth / 2,
                            radar_azimuth,
                            radar_dist,
                        ),
                        self.position.ead_to_lla(
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
                    [self.ead_to_lla(radar_elevation, radar_azimuth, radar_dist)]
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
        from . import RadarData

        self.data_loader = RadarData(self.attrs["campaign"], "rhi")

    @override
    def _show(self, dt, var, ax=None, kwargs_beam={}, **kwargs):
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(FigureClass=TimestampedFigure, timestamp=dt)
        radar = self.data_loader.get_pyart_radar(dt)
        display = pyart.graph.RadarDisplay(radar)

        kwargs = {
            "vmin": -60,
            "vmax": 40,
            "reverse_xaxis": False,  # explicitly False to avoid flipping when all distances are negative.
        } | kwargs
        display.plot(var, ax=ax, **kwargs)
        ax.set_aspect("equal")

        # elevs = radar.elevation["data"]

        elev_azi_start = radar.elevation["data"][0], radar.azimuth["data"][0]
        elev_azi_end = radar.elevation["data"][-1], radar.azimuth["data"][-1]
        plot_beam(
            self, self, elev_azi_start, dt=dt, ax=ax, color="darkgreen", **kwargs_beam
        )
        plot_beam(
            self, self, elev_azi_end, dt=dt, ax=ax, color="limegreen", **kwargs_beam
        )

        return ax
    
    

    def annotate_positions(
        self, positions, dt, ax, *args, plotting_method=None, label=None, **kwargs
    ):
        interpolation_data = kwargs.pop('interpolation_data', None)
        # project the positions to the xy plane...
        xlims = ax.get_xlim()

        if dt is None:
            raise ValueError(
                "dt must be provided for radar positions to get the azimuth"
            )

        eads = [self.position.target_ead(p) for p in positions]

        elevs = np.array([ead[0] for ead in eads])
        azimuths = np.array([ead[1] for ead in eads])
        dists = np.array([ead[2] for ead in eads])

        azimuth = self.data_loader.get_pyart_radar(dt).azimuth["data"][0]

        theta_seps = azimuths % 360 - azimuth % 360
        offsets = dists * np.sin(
            np.deg2rad(theta_seps)
        )  # dist btwn object and radar plane
        # filtered to be less than 5km

        ys = dists * np.sin(np.deg2rad(elevs))
        xs = dists * np.cos(np.deg2rad(elevs)) * np.cos(np.deg2rad(theta_seps))

        # xaxis runs west to east, so if azimuth % 360 > 180, then the x axis is flipped
        if azimuth > 90 and azimuth < 270:
            xs = -xs
        plot_filter = (np.abs(offsets) < 2.5) & (xs > xlims[0]) & (xs < xlims[1])
        if not plot_filter.any():
            return
        if plotting_method is None:
            ax.plot(xs[plot_filter], ys[plot_filter], *args,label=label, **kwargs)
        elif plotting_method == "sized_plot":
            # vary the size based on the offset
            sizes = 1 / np.clip(np.abs(offsets[plot_filter]), 0, 5) * 10
            ax.plot(
                xs[plot_filter],
                ys[plot_filter],
                *args,
                label=label,
                **kwargs,
            )
            ax.scatter(
                xs[plot_filter],
                ys[plot_filter],
                s=sizes,
                *args,
                **kwargs,
            )
        else:
            getattr(ax, plotting_method)(
                xs[plot_filter], ys[plot_filter], *args, **kwargs
            )

        ax.set_xlim(xlims)
        
        
    def annotate_intersections(self, positions, times, dt, ax, **kwargs):
        """Calculates and annotates where a given path intersects the radar scan."""

        # 1. This logic is the same as the start of annotate_positions
        if dt is None:
            raise ValueError("dt must be provided for radar positions")

        eads = [self.position.target_ead(p) for p in positions]
        elevs = np.array([ead[0] for ead in eads])
        azimuths = np.array([ead[1] for ead in eads])
        dists = np.array([ead[2] for ead in eads])
        azimuth = self.data_loader.get_pyart_radar(dt).azimuth["data"][0]
        theta_seps = azimuths % 360 - azimuth % 360
        offsets = dists * np.sin(np.deg2rad(theta_seps))
        ys = dists * np.sin(np.deg2rad(elevs))
        xs = dists * np.cos(np.deg2rad(elevs)) * np.cos(np.deg2rad(theta_seps))
        if azimuth > 90 and azimuth < 270:
            xs = -xs

        # Use the interpolation helper to find the intersection points
        intersections = interpolate_to_intersection(
            offsets=offsets,
            coords_to_interpolate={'x': xs, 'y': ys},
            data_to_interpolate={'time': times}
        )

        # Get the current x-axis limits from the plot before looping
        xlims = ax.get_xlim()

        # Annotate the results on the plot
        for point in intersections:
            # *** ADDED CHECK HERE ***
            # Only annotate if the point is within the visible x-axis of the plot.
            if xlims[0] <= point['x'] <= xlims[1]:
                # This logic is now inside the 'if' block.
                # midnight = dt.replace(hour=0, minute=0, second=0, microsecond=0)
                # dt_seconds_from_midnight = (dt - midnight).total_seconds()
                seconds_at_intersect = point['time']
                # time_delta_seconds = dt_seconds_from_midnight - seconds_at_intersect
                label = f"{seconds_at_intersect/60:.0f}min"
                
                ax.scatter(point['x'], point['y'], **kwargs)
                ax.text(point['x'], point['y'], label, fontsize=8, color='white',
                    bbox=dict(facecolor=kwargs.get('color', 'magenta'), alpha=0.6, pad=1), ha='left', va='bottom')
        
        ax.set_xlim(xlims)

