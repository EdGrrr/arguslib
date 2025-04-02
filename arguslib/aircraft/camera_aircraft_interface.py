import numpy as np

from arguslib.instruments.instruments import PlottableInstrument

from ..instruments.camera import Camera

from ..instruments import Position
from .fleet import Fleet


class CameraAircraftInterface(PlottableInstrument):
    def __init__(self, camera: PlottableInstrument, fleet: Fleet = None):
        self.camera = camera
        self.fleet = fleet  # TODO: loading this data should be easier/automatic. Maybe an AircraftInterface base class to house functionality for this and the radar.
        if self.fleet is None:
            self.fleet = Fleet(
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
            )

        attrs = {"camera": self.camera.attrs}
        super().__init__(**attrs)

    @classmethod
    def from_campaign(cls, campaign, camstr):
        return cls(
            Camera.from_config(campaign, camstr),
        )

    def get_trails(self, time, **kwargs):
        kwargs = {"wind_filter": 10, "tlen": 3600} | kwargs
        return self.fleet.get_trails(time, **kwargs)

    def show(self, dt, ax=None, tlen=3600, color_icao=False, trail_kwargs={}, **kwargs):
        ax = self.camera.show(dt, ax=ax, **kwargs)

        self.plot_trails(dt, ax=ax, tlen=tlen, color_icao=color_icao, **trail_kwargs)
        return ax

    def annotate_positions(self, positions, dt, ax, *args, **kwargs):
        return self.camera.annotate_positions(positions, dt, ax, *args, **kwargs)

    def plot_trails(
        self,
        dt,
        ax,
        color_icao=True,
        plot_kwargs={},
        plot_trails_kwargs={},
        plot_plane_kwargs={},
        **kwargs,
    ):
        kwargs = {"wind_filter": 10, "tlen": 3600} | kwargs
        trail_latlons = self.get_trails(dt, **kwargs)
        trail_alts_geom = self.fleet.get_data(dt, "alt_geom", tlen=kwargs["tlen"])

        current_data = self.fleet.get_current(dt, ["lon", "lat", "alt_geom"])

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

            positions = [
                Position(lon, lat, alt_m) for lon, lat, alt_m in zip(lons, lats, alts_m)
            ]
            self.camera.annotate_positions(
                positions,
                dt,
                ax,
                color="r" if not color_icao else f"#{acft}",
                lw=1,
                **(plot_kwargs | plot_trails_kwargs),
            )
            self.camera.annotate_positions(
                positions[-1:],
                dt,
                ax,
                color="r",
                marker="o",
                markersize=2,
                **(plot_kwargs | plot_plane_kwargs),
            )


# %%
# %%
