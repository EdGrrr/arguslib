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
        # **THE FIX**: Check for list, tuple, or numpy array, and check for length.
        if not isinstance(positions, (list, tuple, np.ndarray)) or len(positions) == 0:
            return

        interpolation_data = kwargs.pop('interpolation_data', None)
        # project the positions to the xy plane...
        xlims = ax.get_xlim()

        if dt is None:
            raise ValueError(
                "dt must be provided for radar positions to get the azimuth"
            )

        # Vectorized call to get elevation, azimuth, and distance for all points
        eads = self.position.target_ead(positions)
        if eads.ndim == 1:
            eads = eads.reshape(1, -1)
        
        elevs, azimuths, dists = eads[:, 0], eads[:, 1], eads[:, 2]

        azimuth = self.data_loader.get_pyart_radar(dt).azimuth["data"][0]

        # Calculate angular separation and perpendicular offset from the scan plane
        theta_seps = azimuths % 360 - azimuth % 360
        offsets = dists * np.sin(np.deg2rad(theta_seps))

        # Calculate height (y) and horizontal range (x) for the RHI plot
        ys = dists * np.sin(np.deg2rad(elevs))
        xs = dists * np.cos(np.deg2rad(elevs)) * np.cos(np.deg2rad(theta_seps))

        # xaxis runs west to east, so if azimuth % 360 > 180, then the x axis is flipped
        if azimuth > 90 and azimuth < 270:
            xs = -xs
        plot_filter = (np.abs(offsets) < 2.5) & (xs > xlims[0]) & (xs < xlims[1])
        if not plot_filter.any():
            return
            
        filtered_xs = xs[plot_filter]
        filtered_ys = ys[plot_filter]
        
        if plotting_method is None:
            ax.plot(filtered_xs, filtered_ys, *args,label=label, **kwargs)
        elif plotting_method == "sized_plot":
            # vary the size based on the offset
            filtered_offsets = offsets[plot_filter]
            sizes = 1 / np.clip(np.abs(filtered_offsets), 0.1, 5) * 10
            ax.plot(
                filtered_xs,
                filtered_ys,
                *args,
                label=label,
                **kwargs,
            )
            ax.scatter(
                filtered_xs,
                filtered_ys,
                s=sizes,
                *args,
                **kwargs,
            )
        else:
            getattr(ax, plotting_method)(
                filtered_xs, filtered_ys, *args, **kwargs
            )

        ax.set_xlim(xlims)
        
        
    def annotate_intersections(self, positions, ages, dt, ax, **kwargs):
        """Calculates and annotates where a given path intersects the radar scan."""
        
        (dt_start, dt_end) = kwargs.pop('time_bounds', (None, None))

        if dt is None:
            raise ValueError("dt must be provided for radar positions")

        # A single vectorized call now handles all positions efficiently.
        eads = self.position.target_ead(positions)
        
        # Use direct numpy slicing for clarity and performance.
        elevs = eads[:, 0]
        azimuths = eads[:, 1]
        dists = eads[:, 2]

        pyart_radar = self.data_loader.get_pyart_radar(dt)
        
        
        azimuth = pyart_radar.azimuth["data"][0]
        theta_seps = azimuths % 360 - azimuth % 360
        offsets = dists * np.sin(np.deg2rad(theta_seps))
        ys = dists * np.sin(np.deg2rad(elevs))
        xs = dists * np.cos(np.deg2rad(elevs)) * np.cos(np.deg2rad(theta_seps))
        if azimuth > 90 and azimuth < 270:
            xs = -xs
        
        if dt_start is not None:
            # get the elev at the two dtimes
            radar_elevs = pyart_radar.elevation["data"]
            radar_times = pyart_radar.time["data"]
            try:
                start_elev = np.interp(dt_start.hour *60*60 +dt_start.minute*60 + dt_start.second + dt_start.microsecond*1e-6, radar_times, radar_elevs)
            except:
                start_elev = radar_elevs[np.argmin(radar_times)]
            try:
                end_elev = np.interp(dt_end.hour *60*60 +dt_end.minute*60 + dt_end.second + dt_end.microsecond*1e-6, radar_times, radar_elevs)
            except:
                end_elev = radar_elevs[np.argmax(radar_times)]
            
            min_elev = min(start_elev, end_elev)
            max_elev = max(start_elev, end_elev)
        else:
            min_elev, max_elev = None, None

        # Use the interpolation helper to find the intersection points
        intersections = interpolate_to_intersection(
            offsets=offsets,
            coords_to_interpolate={'x': xs, 'y': ys},
            data_to_interpolate={'age': ages}
        )

        # Get the current x-axis limits from the plot before looping
        xlims = ax.get_xlim()

        reverse_elevation = azimuth > 90 and azimuth < 270
        elevs = np.array([np.rad2deg(np.arctan2(point['y'], point['x'])) for point in intersections])
        elevs = 180-elevs if reverse_elevation else elevs

        if min_elev is not None:
            intersections_in_elev_range = [
                point for point, elev in zip(intersections, elevs) if min_elev <= elev <= max_elev
            ]
        else:
            intersections_in_elev_range = intersections

        valid_intersections = [
            point for point in intersections_in_elev_range if xlims[0] <= point['x'] <= xlims[1]
        ]

        
        intersect_positions = []
        if valid_intersections:
            # Extract data for plotting (vectorized)
            plot_xs = np.array([point['x'] for point in valid_intersections])
            plot_ys = np.array([point['y'] for point in valid_intersections])
            ages_at_intersect = np.array([point['age'] for point in valid_intersections])

            # Calculate labels
            labels = [f"{age/60:.0f}min" for age in ages_at_intersect]

            # Plot intersections (single scatter call)
            ax.scatter(plot_xs, plot_ys, **kwargs)

            # Annotate intersections
            for i in range(len(plot_xs)):
                ax.text(
                    plot_xs[i],
                    plot_ys[i],
                    labels[i],
                    fontsize=8,
                    color='white',
                    bbox=dict(facecolor=kwargs.get('color', 'magenta'), alpha=0.6, pad=1),
                    ha='left',
                    va='bottom'
                )
            
            elevs = np.rad2deg(np.arctan2(plot_ys, plot_xs))
            elevs = 180 - elevs if reverse_elevation else elevs
            
            intersect_positions.append(self.position.ead_to_lla(
                elevs,
                np.ones_like(elevs)*azimuth,
                np.sqrt(plot_xs**2 + plot_ys**2)
            ))
            
            
        
        ax.set_xlim(xlims)
        return intersect_positions