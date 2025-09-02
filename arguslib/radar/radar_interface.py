"""
Provides an interface for visualizing radar data alongside a plottable instrument like a camera.
"""

import matplotlib.pyplot as plt
from pyart.util import datetime_from_radar
import datetime
import numpy as np

from arguslib.instruments.instruments import PlottableInstrument

from .radar import Radar
from ..camera.camera import Camera
from ..misc.plotting import TimestampedFigure

from .radar_overlay_interface import RadarOverlayInterface


class RadarInterface(PlottableInstrument):
    """Combines a radar and a plottable instrument (e.g., Camera) for synchronized visualization.

    This class facilitates the creation of plots that show a camera's view
    (or another instrument's view) side-by-side with a corresponding radar scan.
    It also handles overlaying radar-derived information, such as the scan
    volume or individual beams, onto the camera's display.

    The core functionality is provided by the `show` method, which generates
    the combined plot.

    Attributes:
        radar (Radar): The radar instrument instance.
        camera (PlottableInstrument): The camera or other instrument to plot alongside the radar.
    """

    def __init__(self, radar: Radar, camera: PlottableInstrument):
        self.radar = radar
        self.camera = camera

        self._overlay_interface = RadarOverlayInterface(radar, self.camera)

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
        # Uses the overlay interface to show the target, which then delegates
        return self._overlay_interface.show(
            dt, ax=ax, allow_timestamp_updates=False, **kwargs
        )

    def show(
        self,
        dt,
        ax=None,
        var="DBZ",
        kwargs_camera=None,
        kwargs_radar_scan=None,
        kwargs_radar_beams=None,
        annotate_beams=True,
        beam_type="start_end",
        ranges_km_for_beams=None,
        annotate_scan_box=True,
        kwargs_scan_box=None,
        show_legend=False,
        **kwargs,
    ):
        """Displays the camera view and radar scan side-by-side for a specific time.

        This is the primary method for this class. It creates a figure with two
        subplots: one for the camera/target instrument's view and one for the
        radar scan (e.g., an RHI or PPI plot). It can also orchestrate the
        annotation of radar beams and scan boundaries on the camera view.

        Args:
            dt (datetime.datetime): The datetime for the visualization. This time
                must correspond to an available radar scan.
            ax (tuple[Axes, Axes], optional): A tuple of two Matplotlib axes
                `(ax_camera, ax_radar)` to plot on. If None, new axes are created.
                Defaults to None.
            var (str, optional): The radar variable to plot (e.g., 'DBZ').
                Defaults to "DBZ".
            kwargs_camera (dict, optional): Keyword arguments passed to the
                `camera.show()` method. Defaults to None.
            kwargs_radar_scan (dict, optional): Keyword arguments passed to the
                `radar.show()` method. Defaults to None.
            kwargs_radar_beams (dict, optional): Keyword arguments for plotting
                radar beams on the camera, passed to `annotate_radar_beams`.
                Defaults to None.
            annotate_beams (bool, optional): If True, overlays radar beams on the
                camera view. Defaults to True.
            beam_type (str, optional): Type of beams to show ('start_end', 'active').
                Defaults to 'start_end'.
            ranges_km_for_beams (list, optional): Distances along the beam to plot.
                Defaults to None.
            annotate_scan_box (bool, optional): If True, overlays the scan extent
                on the camera view. Defaults to True.
            kwargs_scan_box (dict, optional): Keyword arguments for plotting the
                scan box. Defaults to None.
            show_legend (bool, optional): If True, attempts to display a legend
                on the camera plot. Defaults to False.
            **kwargs: Additional keyword arguments passed to `radar.show()`.
        """
        radar = self.radar.data_loader.get_pyart_radar(dt)
        dt_radar = datetime.datetime.fromisoformat(
            datetime_from_radar(radar).isoformat()
        )
        if dt_radar.replace(microsecond=0) != dt.replace(microsecond=0):
            raise ValueError(
                f"Requested dt ({dt}) does not match radar scan time ({dt_radar.replace(microsecond=0)}).\nEnsure dt corresponds to an actual radar scan time."
            )
        # Use dt_radar for all plotting operations to ensure consistency
        current_dt = dt_radar

        _kwargs_camera = kwargs_camera or {}
        _kwargs_radar_scan = kwargs_radar_scan or {}
        _kwargs_radar_beams = kwargs_radar_beams or {}
        _kwargs_scan_box = kwargs_scan_box or {}

        if ax is not None:
            if not isinstance(ax, (tuple, list)) or len(ax) < 2:
                raise ValueError(
                    "If 'ax' is provided, it must be a tuple/list of two axes (ax_target, ax_radar)."
                )
            ax_cam, ax_radar_plot = ax
            fig = ax_cam.figure
            if isinstance(fig, TimestampedFigure):
                fig.timestamp = current_dt
        else:
            fig = plt.figure(
                FigureClass=TimestampedFigure,
                timestamp=current_dt,
                figsize=(10, 4.9),
                dpi=300,
                constrained_layout=True,
            )
            gs = fig.add_gridspec(1, 2)

            camera_subplot_kwargs = {}
            # Check if the camera is 'allsky' to set polar projection by default
            if (
                hasattr(self.camera, "camera_type")
                and self.camera.camera_type == "allsky"
            ):
                camera_subplot_kwargs["projection"] = "polar"

            ax_cam = fig.add_subplot(gs[0], **camera_subplot_kwargs)
            ax_radar_plot = fig.add_subplot(gs[1])  # Radar plot is typically Cartesian

            # Note: Detailed polar axis setup (theta_offset, theta_direction)
            # will be handled by Camera.show() based on its kwargs (_kwargs_camera),
            # now that it will correctly configure a provided polar axis.

        # Pass the (potentially polar) ax_cam to the camera's show method.
        # _kwargs_camera can include theta_behaviour, lr_flip etc.
        ax_cam = self.show_camera(current_dt, ax=ax_cam, **_kwargs_camera)

        ax_radar = self.radar.show(
            current_dt, ax=ax_radar_plot, var=var, **(_kwargs_radar_scan | kwargs)
        )

        if annotate_beams:
            self._overlay_interface.annotate_radar_beams(
                current_dt,
                ax=ax_cam,
                ranges_km=ranges_km_for_beams,
                beam_type=beam_type,
                **_kwargs_radar_beams,
            )

        if annotate_scan_box:
            self._overlay_interface.annotate_scan_box(
                current_dt,
                ax=ax_cam,
                **(_kwargs_scan_box),  # Pass specific scan box kwargs
            )

        if show_legend and hasattr(ax_cam, "legend"):
            # Attempt to add a legend to the target instrument's plot
            # This might need adjustment if the target is a CameraArray (multiple axes)
            # or DirectCamera (no axes).
            try:
                handles, labels = [], []
                # Consolidate legend items from target instrument if possible
                if hasattr(ax_cam, "get_legend_handles_labels"):
                    h, l = ax_cam.get_legend_handles_labels()
                    handles.extend(h)
                    labels.extend(l)

                # If it's a CameraArray, it returns a list of axes.
                # We might want to create a figure-level legend.
                if (
                    isinstance(self.camera, PlottableInstrument)
                    and hasattr(self.camera, "cameras")
                    and isinstance(ax_cam, np.ndarray)
                ):  # Heuristic for CameraArray
                    # For CameraArray, collect unique legend items from all subplots
                    all_handles, all_labels = [], []
                    for sub_ax in ax_cam.ravel():
                        if hasattr(sub_ax, "get_legend_handles_labels"):
                            h, l = sub_ax.get_legend_handles_labels()
                            for handle, label in zip(h, l):
                                if label not in all_labels:  # Add unique items
                                    all_labels.append(label)
                                    all_handles.append(handle)
                    if all_handles:
                        fig.legend(
                            all_handles, all_labels, loc="upper left", fontsize=8
                        )
                elif handles:  # For single axis target
                    ax_cam.legend(handles, labels)

            except Exception as e:
                print(f"Could not generate legend for target instrument: {e}")

        return ax_cam, ax_radar

    def annotate_positions(
        self, positions, dt, ax, cam_kwargs={}, radar_kwargs={}, **kwargs
    ):
        """Annotates geographical positions on both the camera and radar plots.

        Args:
            positions (list[Position]): A list of `Position` objects to annotate.
            dt (datetime.datetime): The datetime for the annotation.
            ax (tuple[Axes, Axes]): A tuple of the two axes `(ax_camera, ax_radar)`
                to plot on.
            cam_kwargs (dict, optional): Keyword arguments passed specifically to
                `camera.annotate_positions`. Defaults to {}.
            radar_kwargs (dict, optional): Keyword arguments passed specifically to
                `radar.annotate_positions`. Defaults to {}.
            **kwargs: Keyword arguments passed to both annotation methods.

        Returns:
            tuple[Axes, Axes]: The updated camera and radar axes.
        """
        ax_cam, ax_radar = ax

        ax_cam = self.camera.annotate_positions(
            positions, dt, ax=ax_cam, **(cam_kwargs | kwargs)
        )
        ax_radar = self.radar.annotate_positions(
            positions, dt, ax=ax_radar, **(radar_kwargs | kwargs)
        )

        return ax_cam, ax_radar

    def annotate_intersections(self, positions, ages, dt, ax, **kwargs):
        intersect_positions = self.radar.annotate_intersections(
            positions, ages, dt, ax[1], **kwargs
        )
        kwargs.pop("time_bounds", None)
        kwargs["plotting_method"] = "scatter"
        self.camera.annotate_positions(intersect_positions, dt, ax=ax[0], **kwargs)
        return len(intersect_positions) > 0
