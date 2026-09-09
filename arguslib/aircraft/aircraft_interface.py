from matplotlib.lines import Line2D
import numpy as np
import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Union

from arguslib.config import load_config
from tqdm import tqdm


from arguslib.misc.plotting import get_fig_from_ax_or_axs, get_timestamp_from_ax

from ..misc.geo import ft_to_km, haversine
from ..protocols import DirectRenderable, ProvidesRadarScanTime, SupportsIntersections

from ..instruments.instruments import PlottableInstrument
from ..instruments import Position
from .fleet import Fleet

if TYPE_CHECKING:
    # This import is only for static type checkers, preventing runtime circular imports.
    from arguslib.radar.radar_interface import RadarInterface


def adjust_trail_positions(positions: list[Position], adjust_km):
    if adjust_km[0] == 0 and adjust_km[1] == 0:
        return positions
    positions = [p.xyz_to_lla(adjust_km[0], adjust_km[1], 0) for p in positions]
    return positions


class AircraftInterface(PlottableInstrument):
    """An interface for visualizing aircraft flight tracks on a plottable instrument.

    This class acts as a wrapper around another `PlottableInstrument` (like a
    `Camera` or `CameraArray`) and a `Fleet` object containing flight data.
    Its primary role is to orchestrate the plotting of the underlying
    instrument's view and then overlaying the aircraft trails and positions
    on top of it.

    Because it inherits from `PlottableInstrument`, it can be used interchangeably
    wherever a plottable object is expected, allowing for powerful composition
    (e.g., wrapping an `AircraftInterface` inside a `RadarInterface`).

    Attributes:
        camera (PlottableInstrument): The underlying instrument to draw on.
            Despite the name, this can be any `PlottableInstrument`.
        fleet (Fleet): The object managing the aircraft data.
    """

    def __init__(
        self, camera: PlottableInstrument, fleet: Fleet = None, **fleet_kwargs
    ):
        self.camera = camera
        self.fleet = fleet  # TODO: loading this data should be easier/automatic. Maybe an AircraftInterface base class to house functionality for this and the radar.
        if self.fleet is None:
            default_fleet_kwargs = {
                "variables": [
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
            }
            self.fleet = Fleet(**(default_fleet_kwargs | fleet_kwargs))

        attrs = {"camera": self.camera.attrs}
        super().__init__(**attrs)

    @classmethod
    def from_campaign(cls, campaign, camstr):
        from ..camera.camera import Camera

        return cls(
            Camera.from_config(campaign, camstr),
        )

    def load_flight_data(
        self,
        date_or_dt: Union[datetime.date, datetime.datetime],
        adsb_data_dir: Union[str, Path] = None,
        force_reload: bool = False,
        **fleet_load_kwargs,
    ):
        """
        Loads ADS-B flight data for the specified date from the given directory
        and assigns ERA5 wind data to the fleet.

        Args:
            date_or_dt: The date (or datetime object) for which to load data.
            adsb_data_dir: The directory containing the ADS-B data files
                             (e.g., YYYYMMDD_ADS-B.nc and YYYYMMDD_ADS-B.txt).

        Raises:
            TypeError: If date_or_dt is not a datetime.date or datetime.datetime object.
            FileNotFoundError: If the ADS-B data directory or necessary files are not found.
        """

        if adsb_data_dir is None:
            try:
                adsb_data_dir = load_config("adsb_path.txt")
            except FileNotFoundError:
                raise FileNotFoundError(
                    "ADS-B data directory not specified and 'adsb_path.txt' not found."
                )

        if isinstance(date_or_dt, datetime.datetime):
            date_to_load = date_or_dt.date()
        elif isinstance(date_or_dt, datetime.date):
            date_to_load = date_or_dt
        else:
            raise TypeError(
                "date_or_dt must be a datetime.date or datetime.datetime object."
            )

        adsb_dir_path = Path(adsb_data_dir)
        if not adsb_dir_path.is_dir():
            raise FileNotFoundError(f"ADS-B data directory not found: {adsb_dir_path}")

        adsb_file_basename = date_to_load.strftime("%Y%m%d") + "_ADS-B"
        # fleet.load_output expects the base path and appends .nc and .txt itself.
        adsb_file_path_base = adsb_dir_path / adsb_file_basename

        if not (
            adsb_file_path_base.with_suffix(".nc").exists()
            and adsb_file_path_base.with_suffix(".txt").exists()
        ):
            raise FileNotFoundError(
                f"Required ADS-B files not found: {adsb_file_path_base.with_suffix('.nc')} "
                f"or {adsb_file_path_base.with_suffix('.txt')}"
            )

        current_loaded_file = self.fleet.loaded_file
        self.fleet.load_output(
            str(adsb_file_path_base), force_reload=force_reload, **fleet_load_kwargs
        )

        if force_reload or (
            not self.fleet.has_notnull_data("uwind")
            and self.fleet.loaded_file != current_loaded_file
        ):
            print("Attempting to assign ERA5 wind data to fleet...")
            try:
                self.fleet.assign_era5_winds()  # This method has its own error handling and print statements
                print("Flight data loading and ERA5 wind assignment process complete.")
            except ValueError:
                print(
                    "Error occurred during flight data loading or ERA5 wind assignment. Will use aircraft ADS-B wind as a fallback."
                )

    def get_trails(self, time, **kwargs):
        kwargs = {"tlen": 3600, "include_time": True} | kwargs
        return self.fleet.get_trails(time, **kwargs)

    def show(self, dt, ax=None, tlen=3600, color_icao=False, trail_kwargs={}, **kwargs):
        ax = self.camera.show(dt, ax=ax, **kwargs)

        # Use the protocol for type-safe, decoupled access to radar time bounds
        if isinstance(self.camera, ProvidesRadarScanTime):
            self.start_time, self.end_time = self.camera.get_scan_time_bounds(dt)

        self.plot_trails(dt, ax=ax, tlen=tlen, color_icao=color_icao, **trail_kwargs)
        return ax

    def annotate_positions(self, positions, dt, ax, *args, **kwargs):
        return self.camera.annotate_positions(positions, dt, ax, *args, **kwargs)

    def plot_trails(
        self,
        dt,
        ax,
        adjust_km=(0, 0),
        adjust_mps=(0, 0),
        color_icao=True,
        label_acft=False,
        icao_include: list = None,
        plot_kwargs=None,
        plot_plane_kwargs=None,
        advection_winds=None,
        plot_baseline=True,
        plot_intersections=None,
        **kwargs,
    ):
        plot_kwargs = dict(plot_kwargs or {})
        plot_plane_kwargs = dict(plot_plane_kwargs or {})
        kwargs = {
            "tlen": 3600,
            "adjust_mps": adjust_mps,
            "winds": advection_winds,
        } | kwargs

        if plot_intersections is None:
            plot_intersections = isinstance(self.camera, SupportsIntersections)

        if not self.fleet.loaded_file:
            print(
                f"Warning (AircraftInterface.plot_trails): No ADS-B data seems to have been loaded (fleet.loaded_file is None). "
                f"Call 'load_flight_data()' on the AircraftInterface instance. No trails will be plotted for {dt}."
            )
            return
        if not self.fleet.aircraft:
            # This means loaded_file might be set, but the file contained no aircraft or was empty.
            print(
                f"Warning (AircraftInterface.plot_trails): ADS-B data loaded from '{self.fleet.loaded_file}', "
                f"but fleet.aircraft is empty. No trails will be plotted for {dt}."
            )
            return

        if ax is None:
            timestamp = self.camera.data_loader.current_image_time
        else:
            timestamp = dt

        loaded_date = datetime.datetime.strptime(
            self.fleet.loaded_file.split("/")[-1], "%Y%m%d_ADS-B"
        )
        if timestamp.replace(hour=0, minute=0, second=0, microsecond=0) != loaded_date:
            raise ValueError(
                f"Plotting timestamp {timestamp} does not match loaded date {loaded_date}"
            )

        if plot_baseline:
            dict_positions = self.get_trail_positions(
                timestamp, icao_include=icao_include, **kwargs
            )
            for acft, (positions, _ages) in dict_positions.items():
                positions = adjust_trail_positions(positions, adjust_km)
                self.camera.annotate_trail(
                    positions,
                    timestamp,
                    ax,
                    trail_kwargs=plot_kwargs
                    | {
                        "color": f"#{acft}" if color_icao else "red",
                        "label": (
                            acft if (label_acft and not plot_intersections) else None
                        ),
                    },
                    plane_kwargs=plot_kwargs
                    | {"color": "r", "marker": "o", "markersize": 2}
                    | plot_plane_kwargs,
                )

        if plot_intersections:
            self._plot_intersections(
                dt,
                ax,
                adjust_km=adjust_km,
                icao_include=icao_include,
                color_icao=color_icao,
                label_acft=label_acft,
                plot_kwargs=plot_kwargs,
                **kwargs,
            )

        if ax is None:
            return
        get_fig_from_ax_or_axs(ax).canvas.draw()

        if isinstance(ax, tuple):  # the complicated case of likely a radar interface.
            axes = []
            for a in ax:
                if not hasattr(a, "flat"):
                    axes.append(a)
                else:
                    axes.extend(a.flat)
        elif not hasattr(ax, "flat"):
            axes = [ax]
        else:
            axes = ax.flat

        for ax_to_clean in axes:
            boundary_path = ax_to_clean.patch.get_path()
            transform = ax_to_clean.patch.get_transform()
            transformed_boundary = boundary_path.transformed(transform)
            # data_bbox = Bbox.from_bounds(*transformed_boundary.get_extents().bounds)

            artists = ax_to_clean.get_children()
            trail_artists = [
                a
                for a in artists
                if isinstance(a, Line2D)
                and a.get_label()
                != "_nolegend_"  # Exclude behind the scenes lines (like axis lines)
            ]
            for artist in trail_artists:
                artist_bbox = artist.get_window_extent()

                # Check for intersection with the precise boundary path
                if not transformed_boundary.intersects_bbox(artist_bbox):
                    artist.remove()

    def _plot_intersections(
        self,
        dt,
        ax,
        adjust_km=(0, 0),
        icao_include=None,
        color_icao=True,
        label_acft=False,
        plot_kwargs=None,
        chunk_size=10,
        **kwargs,
    ):
        plot_kwargs = plot_kwargs or {}
        start, end = self.start_time.timestamp(), self.end_time.timestamp()
        step = chunk_size / 2
        midpoints = np.arange(start + step, end, step)
        edges = np.arange(start, end + step, step)

        plotted_icaos = set()
        for t, ti, tf in tqdm(
            zip(midpoints, edges[:-2], edges[2:]),
            total=len(midpoints),
            desc="Processing radar intersections in time chunks...",
        ):
            dict_positions = self.get_trail_positions(
                datetime.datetime.fromtimestamp(t), icao_include=icao_include, **kwargs
            )
            for acft, (positions, ages) in dict_positions.items():
                if acft in plotted_icaos:
                    continue
                positions = adjust_trail_positions(positions, adjust_km)
                if self.camera.annotate_intersections(
                    positions,
                    ages,
                    dt,
                    ax,
                    time_bounds=(
                        datetime.datetime.fromtimestamp(ti),
                        datetime.datetime.fromtimestamp(tf),
                    ),
                    **(
                        plot_kwargs
                        | {
                            "color": f"#{acft}" if color_icao else "red",
                            "label": acft if label_acft else None,
                            "marker": "x",
                            "s": 25,
                        }
                    ),
                ):
                    plotted_icaos.add(acft)

    def get_trail_positions(self, timestamp, icao_include=None, **kwargs):
        trail_latlons = self.get_trails(timestamp, **kwargs)
        trail_alts_geom = self.fleet.get_data(
            timestamp,
            "alt_geom",
            tlen=kwargs["tlen"],
        )

        # Here we want to filter the trails (if possible) for othly those that pass within 30km of Chilbolton
        # Using trackerlib

        radar = getattr(
            self.camera, "radar", self.camera
        )  # RadarInterface -> .radar, Radar -> itself
        try:
            cao = radar.position
        except AttributeError:
            cao = None

        if cao is not None:
            dist_limit = 30  # radar range limit in km

            trail_array = self.fleet.get_trails_arr(
                timestamp, kwargs["tlen"], winds=kwargs.get("winds", "era5")
            )
            dists = haversine(
                trail_array[:, :, 0], trail_array[:, :, 1], cao.lon, cao.lat
            )
            valid_inds = np.where((dists < dist_limit).sum(axis=1) > 0)
            valid_ids = np.array(self.fleet.get_ids())[valid_inds].tolist()

            if icao_include is not None:
                icao_include = [i for i in icao_include if i in set(valid_ids)]
            else:
                icao_include = valid_ids

        if icao_include is not None:
            trail_latlons = {icao: trail_latlons[icao] for icao in icao_include}

        acfts = []
        positions_lists = []
        ages_lists = []
        for acft in trail_latlons.keys():
            if (
                np.isnan(trail_latlons[acft])
                | (trail_alts_geom[acft]["alt_geom"] < 26000)
            ).all():
                continue

            lons = trail_latlons[acft][0]
            lats = trail_latlons[acft][1]
            alts_km = ft_to_km(trail_alts_geom[acft]["alt_geom"])
            # Get the times array!
            ages = trail_latlons[acft][2]  # Assuming get_trails returns this

            current_pos = self.fleet.interpolate_position(acft, timestamp)

            # differentce between the current timestamp and the last point in the get_trails...

            current_time = 0.0  # Or get a more precise time if available

            lons = np.append(lons, current_pos[0])
            lats = np.append(lats, current_pos[1])
            alts_km = np.append(alts_km, ft_to_km(current_pos[2]))
            # Append the current time to the ages array
            ages = np.append(ages, current_time)

            positions = [
                Position(lon, lat, alt_m)
                for lon, lat, alt_m in zip(lons, lats, alts_km, strict=True)
            ]
            acfts.append(acft)
            positions_lists.append(positions)
            ages_lists.append(ages)
        return dict(zip(acfts, zip(positions_lists, ages_lists)))

    def to_image_array(self, time=True):
        """
        If the underlying camera is a DirectCamera, this method calls its
        to_image_array() method. This is typically called after self.show()
        has prepared the image (including any trails).
        """
        if isinstance(self.camera, DirectRenderable):
            return self.camera.to_image_array(time=time)
        else:
            raise NotImplementedError(
                "to_image_array is only available when the underlying instrument is DirectRenderable."
            )

    @property
    def image(self):
        """
        If the underlying camera is a DirectCamera, this property accesses its image property.
        """
        if isinstance(self.camera, DirectRenderable):
            return self.camera.image
        else:
            raise NotImplementedError(
                "image property is only available when the underlying instrument is DirectRenderable."
            )

    def annotate_trail(self, positions, dt, ax, **kwargs):
        return self.camera.annotate_trail(positions, dt, ax, **kwargs)


class AutomaticADSBAircraftInterface(AircraftInterface):
    def __init__(self, camera: PlottableInstrument, wind_filter=10, **kwargs):
        self.wind_filter = wind_filter
        super().__init__(camera, **kwargs)

    def show(self, dt, *args, **kwargs):
        self.load_flight_data(dt, wind_filter=self.wind_filter)
        return super().show(dt, *args, **kwargs)


# %%
# %%
