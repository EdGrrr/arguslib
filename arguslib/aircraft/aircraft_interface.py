from arguslib.misc.geo import ft_to_km
import numpy as np

from arguslib.instruments.instruments import PlottableInstrument

from ..instruments.camera import Camera

from ..instruments import Position
from .fleet import Fleet


class AircraftInterface(PlottableInstrument):
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
        label_acft=False,
        icao_include: list = None,
        plot_kwargs={},
        plot_trails_kwargs={},
        plot_plane_kwargs={},
        advection_winds="era5",
        **kwargs,
    ):
        kwargs = {"wind_filter": 10, "tlen": 3600} | kwargs
        
        if ax is None:
            # ax is None - which is indicative of a DirectCamera - i.e. matplotlib avoidant
            timestamp = self.camera.data_loader.current_image_time
        else:
            try:
                timestamp = ax.get_figure().timestamp
            except AttributeError:
                timestamp = ax[-1].get_figure().timestamp
            
        kwargs['winds'] = advection_winds
        trail_latlons = self.get_trails(timestamp, **kwargs)
        trail_alts_geom = self.fleet.get_data(timestamp, "alt_geom", tlen=kwargs["tlen"])

        if icao_include is not None:
            trail_latlons = {icao: trail_latlons[icao] for icao in icao_include}

        for acft in trail_latlons.keys():
            if (
                np.isnan(trail_latlons[acft])
                | (trail_alts_geom[acft]["alt_geom"] < 26000)
            ).all():
                continue

            lons = trail_latlons[acft][0]
            lats = trail_latlons[acft][1]
            alts_km = ft_to_km(trail_alts_geom[acft]["alt_geom"])

            current_pos = self.fleet.aircraft[acft].pos.interpolate_position(timestamp)
            lons = np.append(lons, current_pos[0])
            lats = np.append(lats, current_pos[1])
            alts_km = np.append(
                alts_km, ft_to_km(current_pos[2])
            )

            positions = [
                Position(lon, lat, alt_m) for lon, lat, alt_m in zip(lons, lats, alts_km)
            ]
            self.camera.annotate_positions(
                positions,
                timestamp,
                ax,
                color="r" if not color_icao else f"#{acft}",
                lw=1,
                label=f"{acft}" if label_acft else None,
                **(plot_kwargs | plot_trails_kwargs),
            )
            self.camera.annotate_positions(
                positions[-1:],
                timestamp,
                ax,
                color="r",
                marker="o",
                markersize=2,
                **(plot_kwargs | plot_plane_kwargs),
            )


# %%
# %%
