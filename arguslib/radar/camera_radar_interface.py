import matplotlib.pyplot as plt
from pyart.util import datetime_from_radar
import datetime
import numpy as np

from ..instruments import Camera, Radar
from ..video.locator import CameraData
from .locator import RadarData
from ..misc.plotting import plot_beam


class CameraRadarInterface:
    def __init__(self, radar, camera, camera_data, radar_data):
        self.radar = radar
        self.camera = camera
        self.camera_data = camera_data
        self.radar_data = radar_data

    @classmethod
    def from_campaign(cls, campaign, camstr):
        return cls(
            Radar.from_config(campaign),
            Camera.from_config(campaign, camstr),
            CameraData(campaign, camstr),
            RadarData(campaign, "rhi"),
        )

    def show_camera(self, dt, ax=None, **kwargs):
        radar = self.radar_data.get_pyart_radar(dt)
        dt_radar = datetime.datetime.fromisoformat(
            datetime_from_radar(radar).isoformat()
        )

        if ax is None:
            _, ax = plt.subplots()

        img = self.camera_data.get_data_time(dt_radar)
        ax.imshow(img[:, :, ::-1], origin="lower")

        elev_azi_start = radar.elevation["data"][0], radar.azimuth["data"][0]
        elev_azi_end = radar.elevation["data"][-1], radar.azimuth["data"][-1]

        kwargs = {"c": "limegreen", "lw": 0.7} | kwargs
        plot_beam(self.camera, self.radar, elev_azi_start, ax=ax, **kwargs)
        plot_beam(self.camera, self.radar, elev_azi_end, ax=ax, **kwargs)

        ax.legend()
        return ax

    def show(self, dt):
        fig, (ax_cam, ax_radar) = plt.subplots(
            1, 2, figsize=(10, 4), dpi=300, width_ratios=[1, 2]
        )
        self.show_camera(dt, ax=ax_cam)

        self.radar_data.plot(dt, "DBZ", ax=ax_radar, vmin=-60, vmax=40)

        # plot the sweep extremes
        xlims = ax_radar.get_xlim()
        ylims = ax_radar.get_ylim()

        radar = self.radar_data.get_pyart_radar(dt)
        elev_azi_start = radar.elevation["data"][0], radar.azimuth["data"][0]
        elev_azi_end = radar.elevation["data"][-1], radar.azimuth["data"][-1]
        # from 0,0, to xlims[1], xlims[1]*np.tan(np.deg2rad(90-elev_azi[0])) or ylims[1]*np.tan(np.deg2rad(elev_azi[0])), ylims[1]
        # plot start
        endpoint_start_x = xlims[1], xlims[1] * np.tan(np.deg2rad(elev_azi_start[0]))
        endpoint_start_y = (
            ylims[1] * np.tan(np.deg2rad(90 - elev_azi_start[0])),
            ylims[1],
        )
        endpoint = (
            endpoint_start_y
            if endpoint_start_x[1] > endpoint_start_y[0]
            else endpoint_start_x
        )
        ax_radar.plot([0, endpoint[0]], [0, endpoint[1]], c="limegreen", lw=0.7)

        # plot end
        endpoint_end_x = xlims[1], xlims[1] * np.tan(np.deg2rad(elev_azi_end[0]))
        endpoint_end_y = ylims[1] * np.tan(np.deg2rad(90 - elev_azi_end[0])), ylims[1]
        endpoint = (
            endpoint_end_y if endpoint_end_x[1] > endpoint_end_y[0] else endpoint_end_x
        )
        ax_radar.plot([0, endpoint[0]], [0, endpoint[1]], c="limegreen", lw=0.7)

        ax_radar.set_xlim(xlims)
        ax_radar.set_ylim(ylims)

        return fig, (ax_cam, ax_radar)
