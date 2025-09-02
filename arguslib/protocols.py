"""
Defines structural types (Protocols) for use in type hinting across arguslib.

This helps decouple interfaces from concrete implementations, preventing circular
imports and improving static analysis.
"""

import datetime
import numpy as np
from typing import Protocol, runtime_checkable


@runtime_checkable
class ProvidesRadarScanTime(Protocol):
    """A protocol for instruments that can provide radar scan start and end times."""

    def get_scan_time_bounds(
        self, dt: datetime.datetime
    ) -> tuple[datetime.datetime, datetime.datetime]:
        """Returns the UTC start and end time of the radar scan corresponding to dt."""
        ...


@runtime_checkable
class DirectRenderable(Protocol):
    """A protocol for instruments that render directly to an image array.

    This contract is for instruments that do not use Matplotlib for plotting.
    Their `show` and `annotate_positions` methods operate on an internal image
    buffer and do not return an Axes object. The final image can be retrieved
    via `to_image_array()`.
    """

    def show(self, dt: datetime.datetime, ax=None, **kwargs) -> None: ...

    def annotate_positions(self, positions, dt, ax=None, **kwargs) -> None: ...

    def to_image_array(self, time: bool = True) -> np.ndarray: ...
