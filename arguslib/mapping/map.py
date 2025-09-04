import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from typing_extensions import override
import cartopy.feature as cfeature

from arguslib.misc.plotting import TimestampedFigure

from ..instruments import PlottableInstrument, Position


class MapInstrument(PlottableInstrument):
    """
    An instrument for creating and displaying maps using Cartopy.

    This instrument acts as a plottable background for other instruments,
    providing a geographical context.
    """

    def __init__(self, projection=None, extent=None, **attrs):
        """
        Initializes the map instrument.

        Args:
            projection: A Cartopy CRS projection object. Defaults to PlateCarree.
            extent (list, optional): The map extent as [lon_min, lon_max, lat_min, lat_max].
            **attrs: Additional attributes for the instrument.
        """
        if projection is None:
            projection = ccrs.PlateCarree()
        self.projection = projection
        self.extent = extent
        super().__init__(**attrs)

    @classmethod
    def from_config(cls, map_name: str, **kwargs):
        """
        Creates a MapInstrument from a configuration file.

        Args:
            map_name (str): The name of the map configuration in 'maps.yml'.
            **kwargs: Additional keyword arguments to override config values.

        Returns:
            MapInstrument: An instance of the map instrument.
        """
        from arguslib.config import load_config
        import cartopy.crs as ccrs

        maps_config = load_config("maps.yml")
        map_details = maps_config.get(map_name, {})

        # Allow kwargs to override config
        config = map_details | kwargs

        # Dynamically get projection from cartopy
        proj_name = config.pop(
            "projection", "PlateCarree"
        )  # Pop to avoid passing it twice
        proj_class = getattr(ccrs, proj_name, ccrs.PlateCarree)
        projection = proj_class()

        # Pop extent to avoid passing it twice, with a default of None
        extent = config.pop("extent", None)

        return cls(projection=projection, extent=extent, **config)

    @override
    def show(self, dt, ax=None, **kwargs):
        """
        Renders the map, creating a new figure if no axes are provided.

        Args:
            dt (datetime.datetime): The timestamp for the visualization.
                (can be ignored for a static map but required by the interface).
            ax (matplotlib.axes.Axes, optional): The axis to plot on. If None,
                a new figure and axis with the specified projection are created.
            **kwargs: Additional keyword arguments for map features (e.g., `resolution`).

        Returns:
            matplotlib.axes.Axes: The axis on which the map was drawn.
        """
        if ax is None:
            fig, ax = plt.subplots(
                subplot_kw={"projection": self.projection},
                FigureClass=TimestampedFigure,
                timestamp=dt,
            )

        # Add common map features first
        ax.coastlines(resolution=kwargs.get("resolution", "10m"))
        ax.gridlines(draw_labels=True, linestyle="--", alpha=0.5)
        ax.add_feature(cfeature.LAND, facecolor="lightgray", alpha=0.5)
        ax.add_feature(cfeature.OCEAN, facecolor="lightblue", alpha=0.5)
        ax.add_feature(cfeature.LAKES, facecolor="lightblue", alpha=0.5)
        # Set the extent after drawing features, as they can sometimes override it.
        if self.extent:
            ax.set_extent(self.extent, crs=ccrs.PlateCarree())

        return ax

    @override
    def annotate_positions(self, positions: list[Position], dt, ax, **kwargs):
        """
        Annotates one or more geographical positions on the map.

        Args:
            positions (list[Position]): A list of `Position` objects to annotate.
            dt (datetime.datetime): The datetime for which the view is valid (can be ignored).
            ax (matplotlib.axes.Axes): The axis to plot on. Must be a Cartopy GeoAxes.
            **kwargs: Keyword arguments passed to `ax.plot` or `ax.scatter`.
                      A special `plotting_method` kwarg can be 'scatter'.

        Returns:
            matplotlib.axes.Axes: The updated axis.
        """
        if not positions:
            return ax

        if not hasattr(ax, "coastlines"):
            raise TypeError(
                "The provided axes `ax` must be a Cartopy GeoAxes for map annotations."
            )

        lons = [p.lon for p in positions]
        lats = [p.lat for p in positions]

        plot_kwargs = {"transform": ccrs.PlateCarree()} | kwargs
        plotting_method = plot_kwargs.pop("plotting_method", "plot")

        if plotting_method == "scatter":
            ax.scatter(lons, lats, **plot_kwargs)
        else:  # default to plot
            ax.plot(lons, lats, **plot_kwargs)

        # Restore the original extent. The CRS of the extent tuple is the native
        # projection of the axes, so we pass that to set_extent.
        ax.set_extent(self.extent, crs=ccrs.PlateCarree())

        return ax
