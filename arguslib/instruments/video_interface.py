import datetime
from typing import Tuple, Optional, Any

from arguslib.aircraft.aircraft_interface import AircraftInterface
import cv2
import numpy as np

from .direct_camera import DirectCamera
from .instruments import PlottableInstrument, Position
from ..radar.radar_overlay_interface import RadarOverlayInterface # Allow RadarOverlayInterface


class VideoInterface(PlottableInstrument):
    """
    A PlottableInstrument that uses a DirectCamera to produce video frames
    and provides utilities for writing these frames to a video file.
    """

    def __init__(self, instrument: PlottableInstrument):
        """
        Args:
            instrument: A PlottableInstrument that behaves like a DirectCamera.
                        It must implement show(dt, ax=None, ...), to_image_array(),
                        and annotate_positions(..., ax=None, ...).
                        Typically DirectCamera, AircraftInterface (wrapping DirectCamera),
                        or RadarOverlayInterface (wrapping DirectCamera).
        """
        if isinstance(instrument, DirectCamera):
            self.direct_camera_delegate = instrument
            attrs_source = instrument
        elif isinstance(instrument, AircraftInterface):
            if not (hasattr(instrument, 'to_image_array') and callable(getattr(instrument, 'to_image_array'))):
                raise TypeError(
                    "If AircraftInterface is passed to VideoInterface, "
                    "it must support to_image_array (i.e., its underlying camera must be a DirectCamera)."
                )
            self.direct_camera_delegate = instrument
            attrs_source = instrument.camera # Get attrs from the base camera
        elif isinstance(instrument, RadarOverlayInterface):
            if not (hasattr(instrument, 'to_image_array') and callable(getattr(instrument, 'to_image_array'))):
                raise TypeError(
                    "If RadarOverlayInterface is passed to VideoInterface, "
                    "it must support to_image_array (i.e., its target_instrument must be a DirectCamera or wrap one)."
                )
            self.direct_camera_delegate = instrument
            attrs_source = instrument.target_instrument # Get attrs from the base target
        else:
            if not (hasattr(instrument, 'show') and callable(getattr(instrument, 'show')) and \
                    hasattr(instrument, 'to_image_array') and callable(getattr(instrument, 'to_image_array')) and \
                    hasattr(instrument, 'annotate_positions') and callable(getattr(instrument, 'annotate_positions'))):
                raise TypeError(
                    f"Instrument passed to VideoInterface must behave like a DirectCamera. Got {type(instrument)}")
            self.direct_camera_delegate = instrument
            attrs_source = instrument

        super().__init__(**getattr(attrs_source, 'attrs', {}))

        self.video_writer: Optional[cv2.VideoWriter] = None
        self.output_path: Optional[str] = None
        self.fps: Optional[int] = None
        self.resolution: Optional[Tuple[int, int]] = None  # (width, height)

    def show(self, dt: datetime.datetime, ax: Any = None, **kwargs: Any) -> None:
        """
        Prepares the direct camera's internal image for the given datetime.
        Consistent with DirectCamera, this method updates the internal state
        and returns None as no Matplotlib axes are involved.

        Args:
            dt: The datetime for which to prepare the frame.
            ax: Should be None. Included for PlottableInstrument compatibility.
            **kwargs: Additional arguments passed to the direct_camera.show() method
                      (e.g., brightness_adjust).
        """
        if ax is not None:
            raise ValueError("VideoInterface using DirectCamera should not have Matplotlib axes.")
        # This call prepares/updates self.direct_camera.data_loader.image
        self.direct_camera_delegate.show(dt, ax=None, **kwargs)
        return None

    def annotate_positions(self, positions: list[Position], dt: datetime.datetime, ax: Any = None, *args: Any, **kwargs: Any) -> None:
        """
        Delegates annotation to the underlying direct camera.
        'ax' is expected to be None for DirectCamera.
        """
        if ax is not None:
            # This check is also in DirectCamera, but good practice at interface level.
            raise ValueError("VideoInterface using DirectCamera should not have Matplotlib axes.")
        # DirectCamera.annotate_positions will draw on its internal image and return None.
        return self.direct_camera_delegate.annotate_positions(positions, dt, ax, *args, **kwargs)

    def start_video_output(self, output_path: str, fps: int, resolution: Tuple[int, int]) -> None:
        """
        Initializes video writing to a file.

        Args:
            output_path: Path to the output video file.
            fps: Frames per second for the video.
            resolution: Tuple (width, height) for the video frames.
        """
        if self.video_writer is not None:
            self.finish_video_output()  # Close existing video if any

        self.output_path = output_path
        self.fps = fps
        self.resolution = resolution  # (width, height)

        fourcc = cv2.VideoWriter_fourcc(*'avc1')  # Or 'avc1', 'XVID', etc.
        self.video_writer = cv2.VideoWriter(self.output_path, fourcc, self.fps, self.resolution)
        if not self.video_writer.isOpened():
            raise IOError(f"Could not open video writer for path: {self.output_path}")

    def add_frame_to_video(self, dt: datetime.datetime, show_kwargs: Optional[dict] = None, time_overlay: bool = True, image_array: Optional[np.ndarray] = None) -> None:
        """
        Adds a frame for the given datetime 'dt' to the video.
        If image_array is provided (in BGR format), it's used directly.
        Otherwise, a frame is generated using the direct_camera.

        Args:
            dt: Datetime for the frame.
            show_kwargs: Arguments to pass to self.direct_camera.show().
            time_overlay: If True, adds a timestamp overlay to the frame (if generated).
            image_array: Optional pre-rendered BGR numpy array.
        """
        if self.video_writer is None or self.resolution is None:
            raise RuntimeError("Video output not started. Call start_video_output() first.")

        if image_array is None:
            self.show(dt, **(show_kwargs or {}))  # Updates direct_camera's internal image
            # direct_camera.to_image_array returns RGB
            frame_rgb = self.direct_camera_delegate.to_image_array(time=time_overlay)
            frame_bgr = frame_rgb[:, :, ::-1]  # Convert RGB to BGR for OpenCV
        else:
            frame_bgr = image_array  # Assume provided image_array is already BGR

        # Resize if necessary (e.g. if image_array has different dimensions)
        if frame_bgr.shape[1] != self.resolution[0] or frame_bgr.shape[0] != self.resolution[1]:
            frame_bgr = cv2.resize(frame_bgr, self.resolution)

        self.video_writer.write(frame_bgr)

    def finish_video_output(self) -> None:
        """Releases the video writer and resets video properties."""
        if self.video_writer is not None:
            self.video_writer.release()
            print(f"Video saved to {self.output_path}")
        self.video_writer = None
        self.output_path = None
        self.fps = None
        self.resolution = None

    def generate_video(self, output_path: str, start_dt: datetime.datetime,
                       end_dt: datetime.datetime, step_timedelta: datetime.timedelta,
                       fps: int, show_kwargs: Optional[dict] = None, time_overlay: bool = True) -> None:
        """Generates a complete video file by iterating through time."""
        _show_kwargs = show_kwargs or {}
        self.show(start_dt, **_show_kwargs) # Prepare first frame by calling delegate's show
        first_frame_rgb = self.direct_camera_delegate.to_image_array(time=time_overlay) # Get image from delegate
        height, width, _ = first_frame_rgb.shape
        self.start_video_output(output_path, fps, (width, height))

        current_dt = start_dt
        while current_dt <= end_dt:
            self.add_frame_to_video(current_dt, show_kwargs=_show_kwargs, time_overlay=time_overlay)
            print(f"Processed frame for {current_dt.isoformat(timespec='seconds')}", end='\r')
            current_dt += step_timedelta
        print(f"\nFinished processing frames for {output_path}")
        self.finish_video_output()
