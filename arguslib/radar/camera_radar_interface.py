import matplotlib.pyplot as plt
from pyart.util import datetime_from_radar
import datetime
import numpy as np

from ..instruments import Camera, Radar
from ..video.locator import CameraData
from .locator import RadarData
from ..misc.plotting import plot_beam


class CameraRadarInterface:
    def __init__(self, radar, camera):
        self.radar = radar
        self.camera = camera

        if self.radar.data_loader is None:
            radar.initialise_data_loader()

    @classmethod
    def from_campaign(cls, campaign, camstr):
        return cls(
            Radar.from_config(campaign),
            Camera.from_config(campaign, camstr),
            # CameraData(campaign, camstr),
            # RadarData(campaign, "rhi"),
        )

    def show_camera(self, dt, ax=None, **kwargs):
        radar = self.radar.data_loader.get_pyart_radar(dt)
        dt_radar = datetime.datetime.fromisoformat(
            datetime_from_radar(radar).isoformat()
        )

        ax = self.camera.show(dt_radar, replace_ax=ax, **kwargs)

        elev_azi_start = radar.elevation["data"][0], radar.azimuth["data"][0]
        elev_azi_end = radar.elevation["data"][-1], radar.azimuth["data"][-1]

        kwargs = {"c": "limegreen", "lw": 0.7} | kwargs
        plot_beam(self.camera, self.radar, elev_azi_start, ax=ax, **kwargs)
        plot_beam(self.camera, self.radar, elev_azi_end, ax=ax, **kwargs)

        ax.legend()
        return ax

    def show(self, dt, var="DBZ", **kwargs):
        fig, (ax_cam, ax_radar) = plt.subplots(
            1, 2, figsize=(10, 4), dpi=300, width_ratios=[1, 2]
        )
        self.show_camera(dt, ax=ax_cam)

        self.radar.show(dt, ax=ax_radar, var=var, **kwargs)

        return fig, (ax_cam, ax_radar)
