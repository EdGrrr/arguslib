# In a new file, e.g., arguslib/misc/interpolation.py
import numpy as np
from typing import List, Dict, Any
import datetime as dt


def interpolate_to_intersection(
    offsets: np.ndarray,
    coords_to_interpolate: Dict[str, np.ndarray],
    data_to_interpolate: Dict[str, np.ndarray],
) -> List[Dict[str, Any]]:
    """
    Finds where a path crosses a plane (where offset is zero) and linearly
    interpolates associated coordinate and data values at that point.

    Args:
        offsets: An array of distances from the intersection plane. The function
                 looks for sign changes in this array.
        coords_to_interpolate: A dictionary of coordinate arrays (e.g., {'x': xs, 'y': ys})
                               that define the path's position in some space.
        data_to_interpolate: A dictionary of other data arrays (e.g., {'time': times})
                             associated with the path.

    Returns:
        A list of dictionaries, where each dictionary represents one
        intersection and contains the interpolated values for all provided
        coordinates and data.
    """
    intersections = []

    # Find indices *before* a sign change in the offset
    shift_indices = np.where(np.diff(np.sign(offsets)) != 0)[0]

    for idx in shift_indices:
        # The segment that crosses the plane is between idx and idx + 1
        offset1, offset2 = offsets[idx], offsets[idx + 1]

        # Skip if data is invalid (e.g., NaN)
        if np.isnan(offset1) or np.isnan(offset2) or offset1 == offset2:
            continue

        # Calculate interpolation fraction 't' where the offset is 0
        # t = 0 corresponds to point1, t = 1 corresponds to point2
        t = offset1 / (offset1 - offset2)

        # We only want to interpolate if the crossing is within the segment
        if not (0 <= t <= 1):
            continue

        intersection_point = {}

        # Interpolate coordinates
        for key, values in coords_to_interpolate.items():
            val1, val2 = values[idx], values[idx + 1]
            intersection_point[key] = val1 + t * (val2 - val1)

        # Interpolate other data (e.g., time)
        for key, values in data_to_interpolate.items():
            val1, val2 = values[idx], values[idx + 1]

            # Check for invalid values first. np.isnan doesn't work on datetimes.
            is_val1_nan = isinstance(val1, float) and np.isnan(val1)
            is_val2_nan = isinstance(val2, float) and np.isnan(val2)

            if is_val1_nan or is_val2_nan:
                intersection_point[key] = np.nan
                continue

            # Handle datetimes
            if isinstance(val1, dt.datetime) and isinstance(val2, dt.datetime):
                val1_ts, val2_ts = val1.timestamp(), val2.timestamp()
                interpolated_ts = val1_ts + t * (val2_ts - val1_ts)
                tz = (
                    getattr(val1, "tzinfo", None)
                    or getattr(val2, "tzinfo", None)
                    or dt.timezone.utc
                )
                intersection_point[key] = dt.datetime.fromtimestamp(
                    interpolated_ts, tz=tz
                )
            # Handle numeric types
            elif isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
                intersection_point[key] = val1 + t * (val2 - val1)
            # Handle mixed or other unsupported types by assigning NaN
            else:
                intersection_point[key] = np.nan

        intersections.append(intersection_point)

    return intersections
