import matplotlib.pyplot as plt
from pyart.util import datetime_from_radar
import datetime

from arguslib.instruments.instruments import PlottableInstrument

from ..instruments.radar import Radar

from ..instruments.camera import Camera
from ..misc.plotting import TimestampedFigure, plot_beam


class RadarInterface(PlottableInstrument):
    def __init__(self, radar: Radar, camera: PlottableInstrument):
        self.radar = radar
        self.camera = camera

        if self.radar.data_loader is None:
            radar.initialise_data_loader()

        attrs = {
            "camera": self.camera.attrs,
            "radar": self.radar.attrs,
        }
        super().__init__(**attrs)

    @classmethod
    def from_campaign(cls, campaign, camstr):
        return cls(
            Radar.from_config(campaign),
            Camera.from_config(campaign, camstr),
        )

    def show_camera(self, dt, show_legend=False, ax=None, kwargs_beam={}, **kwargs):
        radar = self.radar.data_loader.get_pyart_radar(dt)
        dt_radar = datetime.datetime.fromisoformat(
            datetime_from_radar(radar).isoformat()
        )
        if dt_radar.replace(microsecond=0) != dt.replace(microsecond=0):
            raise ValueError(
                f"dt ({dt}) does not match any radar data. (Try {dt_radar.replace(microsecond=0)})"
            )

        ax = self.camera.show(dt_radar, replace_ax=ax, **kwargs)

        elev_azi_start = radar.elevation["data"][0], radar.azimuth["data"][0]
        elev_azi_end = radar.elevation["data"][-1], radar.azimuth["data"][-1]

        kwargs_beam = {"lw": 0.7} | kwargs_beam
        plot_beam(
            self.camera,
            self.radar,
            elev_azi_start,
            dt=dt,
            ax=ax,
            color="darkgreen",
            **kwargs_beam,
        )
        plot_beam(
            self.camera,
            self.radar,
            elev_azi_end,
            dt=dt,
            ax=ax,
            color="limegreen",
            **kwargs_beam,
        )

        # get the cross-scan bounding box at 10km.
        pos_start_10km = self.radar.position.ead_to_lla(*elev_azi_start, 10)
        # displace by pm5km
        orthogonal_direction = (elev_azi_start[1] + 90) % 360
        position_corners = [
            pos_start_10km.ead_to_lla(0, orthogonal_direction, 2.5),
            pos_start_10km.ead_to_lla(0, orthogonal_direction, -2.5),
        ]

        pos_end_10km = self.radar.position.ead_to_lla(*elev_azi_end, 10)
        position_corners += [
            pos_end_10km.ead_to_lla(0, orthogonal_direction, -2.5),
            pos_end_10km.ead_to_lla(0, orthogonal_direction, 2.5),
        ]
        position_corners += [position_corners[0]]
        self.camera.annotate_positions(
            position_corners, dt, ax=ax, color="yellow", linewidth=1
        )

        if show_legend:
            try:
                ax.legend()
            except AttributeError:
                # probably a multicam, so just place the legend on the top left of the figure
                # but only include the first two items in the legend
                fig = plt.gcf()
                handles = fig.axes[0].get_legend_handles_labels()[0]
                labels = fig.axes[0].get_legend_handles_labels()[1]
                handles = handles[:2]
                labels = labels[:2]
                fig.legend(handles=handles, labels=labels, loc="upper left", fontsize=8)

        return ax

    def show(self, dt, ax=None, var="DBZ", kwargs_camera={}, kwargs_beam={}, **kwargs):
        if ax is not None:
            raise ValueError("We need to start with a clean figure")
        fig, (ax_cam, ax_radar) = plt.subplots(
            1,
            2,
            figsize=(10, 4.2),
            dpi=300,
            width_ratios=[0.8, 1.2],
            FigureClass=TimestampedFigure,
            timestamp=dt,
            # constrained_layout=True,
        )
        ax_cam = self.show_camera(
            dt, ax=ax_cam, kwargs_beam=kwargs_beam, **kwargs_camera
        )

        ax_radar = self.radar.show(
            dt, ax=ax_radar, var=var, kwargs_beam=kwargs_beam, **kwargs
        )

        return ax_cam, ax_radar

    def annotate_positions(
        self, positions, dt, ax, cam_kwargs={}, radar_kwargs={}, **kwargs
    ):
        ax_cam, ax_radar = ax

        ax_cam = self.camera.annotate_positions(
            positions, dt, ax=ax_cam, **(cam_kwargs | kwargs)
        )
        ax_radar = self.radar.annotate_positions(
            positions, dt, ax=ax_radar, **(radar_kwargs | kwargs)
        )

        return ax_cam, ax_radar
