import matplotlib.pyplot as plt
from pyart.util import datetime_from_radar
import datetime
import numpy as np

from arguslib.instruments.instruments import PlottableInstrument

from ..instruments.radar import Radar
from ..instruments.camera import Camera
from ..misc.plotting import TimestampedFigure

from .radar_overlay_interface import RadarOverlayInterface


class RadarInterface(PlottableInstrument):
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
        return self._overlay_interface.show(dt, ax=ax, **kwargs)

    def show(self, dt, ax=None, var="DBZ",
             kwargs_camera=None,
             kwargs_radar_scan=None,
             kwargs_radar_beams=None,
             annotate_beams=True,
             beam_type='start_end',
             ranges_km_for_beams=None,
             annotate_scan_box=False,
             kwargs_scan_box=None,
             show_legend=False,
             **kwargs):

        radar = self.radar.data_loader.get_pyart_radar(dt)
        dt_radar = datetime.datetime.fromisoformat(
            datetime_from_radar(radar).isoformat()
        )
        if dt_radar.replace(microsecond=0) != dt.replace(microsecond=0):
            raise ValueError(
        f"Requested dt ({dt}) does not match radar scan time ({dt_radar.replace(microsecond=0)}).\nEnsure dt corresponds to an actual radar scan time.")
        # Use dt_radar for all plotting operations to ensure consistency
        current_dt = dt_radar
        
        _kwargs_camera = kwargs_camera or {}
        _kwargs_radar_scan = kwargs_radar_scan or {}
        _kwargs_radar_beams = kwargs_radar_beams or {}
        _kwargs_scan_box = kwargs_scan_box or {}
        
        
        if ax is not None:
            if not isinstance(ax, (tuple, list)) or len(ax) < 2:
                 raise ValueError("If 'ax' is provided, it must be a tuple/list of two axes (ax_target, ax_radar).")
            ax_cam, ax_radar_plot = ax
            fig = ax_cam.figure
            if isinstance(fig, TimestampedFigure): fig.timestamp = current_dt
        else:
            fig, (ax_cam, ax_radar_plot) = plt.subplots(
            1,
            2,
            figsize=(10, 4.2),
            dpi=300,
            width_ratios=[0.8, 1.2],
            FigureClass=TimestampedFigure,
            timestamp=current_dt,
        )
        ax_cam = self.show_camera(
            current_dt, ax=ax_cam, **_kwargs_camera
        )

        ax_radar = self.radar.show(
            current_dt, ax=ax_radar_plot, var=var, **(_kwargs_radar_scan | kwargs)
        )

        
        if annotate_beams:
            self._overlay_interface.annotate_radar_beams(
                current_dt, ax=ax_cam,
                ranges_km=ranges_km_for_beams,
                beam_type=beam_type,
                **_kwargs_radar_beams
            )

        if annotate_scan_box:
            self._overlay_interface.annotate_scan_box(
                current_dt, ax=ax_cam,
                **(_kwargs_scan_box) # Pass specific scan box kwargs
            )

        if show_legend and hasattr(ax_cam, 'legend'):
            # Attempt to add a legend to the target instrument's plot
            # This might need adjustment if the target is a CameraArray (multiple axes)
            # or DirectCamera (no axes).
            try:
                handles, labels = [], []
                # Consolidate legend items from target instrument if possible
                if hasattr(ax_cam, 'get_legend_handles_labels'):
                    h, l = ax_cam.get_legend_handles_labels()
                    handles.extend(h)
                    labels.extend(l)
                
                # If it's a CameraArray, it returns a list of axes.
                # We might want to create a figure-level legend.
                if isinstance(self.camera, PlottableInstrument) and \
                   hasattr(self.camera, 'cameras') and \
                   isinstance(ax_cam, np.ndarray): # Heuristic for CameraArray
                    # For CameraArray, collect unique legend items from all subplots
                    all_handles, all_labels = [], []
                    for sub_ax in ax_cam.ravel():
                        if hasattr(sub_ax, 'get_legend_handles_labels'):
                            h, l = sub_ax.get_legend_handles_labels()
                            for handle, label in zip(h,l):
                                if label not in all_labels: # Add unique items
                                    all_labels.append(label)
                                    all_handles.append(handle)
                    if all_handles:
                        fig.legend(all_handles, all_labels, loc="upper left", fontsize=8)
                elif handles: # For single axis target
                    ax_cam.legend(handles, labels)

            except Exception as e:
                print(f"Could not generate legend for target instrument: {e}")

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
