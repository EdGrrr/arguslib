import matplotlib.pyplot as plt
import numpy as np

from arguslib.arguslib.misc.plotting import plot_range_rings

from ..instruments import Camera, Position
from .fleet import Fleet
from ..video.locator import CameraData
from ..misc import geo


class CameraAircraftInterface:
    def __init__(self, camera: Camera, fleet: Fleet, camera_data: CameraData):
        self.camera = camera
        self.fleet = fleet
        self.camera_data = camera_data

    @classmethod
    def from_campaign(cls, campaign, camstr):
        return cls(
            Camera.from_config(campaign, camstr),
            Fleet(
                variables=[
                    "lon",
                    "lat",
                    "alt_baro",
                    "alt_geom",
                    "geom_rate",
                    "tas",
                    "gs",
                    "ws",
                    "track",
                    "true_heading",
                    "wd",
                    "oat",
                ]
            ),
            CameraData(campaign, camstr),
        )

    def get_trails(self, time, **kwargs):
        kwargs = {"wind_filter": 10, "tlen": 3600} | kwargs
        return self.fleet.get_trails(time, **kwargs)

    def get_image(self, time):
        return self.camera_data.get_data_time(time)

    def show(self, time, ax=None, tlen=3600, color_icao=False):
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        img = self.get_image(time)
        ax.imshow(img[:, :, ::-1], origin="lower")
        plot_range_rings(self.camera, ax=ax)
        self.plot_trails(time, ax=ax, tlen=tlen, color_icao=color_icao)
        return ax

    def plot_trails(self, time, ax, color_icao=False, **kwargs):
        kwargs = {"wind_filter": 10, "tlen": 3600} | kwargs
        trail_latlons = self.get_trails(time, **kwargs)
        trail_alts_geom = self.fleet.get_data(time, "alt_geom", tlen=kwargs["tlen"])

        current_data = self.fleet.get_current(time, ["lon", "lat", "alt_geom"])

        for acft in trail_latlons.keys():
            if (
                np.isnan(trail_latlons[acft])
                | (trail_alts_geom[acft]["alt_geom"] < 26000)
            ).all():
                continue

            lons = trail_latlons[acft][0]
            lats = trail_latlons[acft][1]
            alts_m = trail_alts_geom[acft]["alt_geom"] / (3.33 * 1000)

            if current_data[acft]["lat"] != -9999999:
                lons = np.append(lons, current_data[acft]["lon"])
                lats = np.append(lats, current_data[acft]["lat"])
                alts_m = np.append(
                    alts_m, current_data[acft]["alt_geom"] / (3.33 * 1000)
                )

            # TODO: I want a function on Camera that will take a lon, lat, alt_m (or iterables of them) and return the pixel coordinates
            dists = geo.haversine(
                self.camera.position.lon, self.camera.position.lat, lons, lats
            )

            if (dists[~np.isnan(dists)] > 90).all():
                continue

            pl_track = np.array(
                [
                    self.camera.target_pix(Position(lon, lat, alt_m))
                    for lat, lon, alt_m in zip(lats, lons, alts_m)
                ]
            )
            c = "r" if not color_icao else f"#{acft}"
            ax.plot(pl_track.T[0][dists < 90], pl_track.T[1][dists < 90], c=c, lw=1)
            if dists[-1] < 90:
                ax.plot(pl_track.T[0][-1], pl_track.T[1][-1], "ro", markersize=2)


# %%
# %%
