import datetime as dtmod
from dataclasses import dataclass
from typing import Tuple, Protocol
from .video import Video  # Core dependency


class FrameSource(Protocol):
    """
    Interface for any class that can provide an image frame for a given datetime.
    """

    def time_bounds(self) -> Tuple[dtmod.datetime, dtmod.datetime]: ...
    def get_data_time(self, dt: dtmod.datetime, return_timestamp: bool = False): ...
    def __len__(self) -> int: ...


@dataclass
class VideoFrameSource:
    """
    A FrameSource that reads from a local video file.
    """

    path: str  # local path only
    timestamp_timezone: str = "UTC"

    def __post_init__(self):
        self._video = Video(self.path, timestamp_timezone=self.timestamp_timezone)

    def time_bounds(self):
        return self._video.time_bounds

    def get_data_time(self, dt, return_timestamp=False):
        return self._video.get_data_time(dt, return_timestamp=return_timestamp)

    def __len__(self):
        return len(self._video)
