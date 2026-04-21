from arguslib.camera.undistorted_camera import UndistortedCamera
from arguslib.camera.camera_array import CameraArray
from arguslib.radar.locator import initialise_locator
from collections import defaultdict
import calendar
import datetime
import pytz
import os
import re

def convert_to_london_naive(dt: datetime.datetime) -> datetime.datetime:
    """
    Converts a datetime object to a naive datetime in the 'Europe/London' timezone.

    Args:
        dt: The input datetime object. If naive, it's assumed to be in UTC.
            If aware, it's converted from its current timezone.

    Returns:
        A naive datetime object representing the equivalent time in London.
    """
    # Define the target timezone
    london_tz = pytz.timezone("Europe/London")

    # 1. If the datetime is naive, assume it's UTC and make it timezone-aware.
    if dt.tzinfo is None or dt.tzinfo.utcoffset(dt) is None:
        dt_utc = pytz.utc.localize(dt)
    # 2. If it's already aware, we'll just use it as is for conversion.
    else:
        dt_utc = dt

    # 3. Convert the UTC datetime to the London timezone.
    #    The result will be a timezone-aware datetime in the London timezone.
    dt_local_aware = dt_utc.astimezone(london_tz)

    # 4. Remove the timezone information to make the datetime naive.
    dt_local_naive = dt_local_aware.replace(tzinfo=None)

    return dt_local_naive


def find_active_year_months(campaign):
    """Scans for all years and months that have data."""
    locator = initialise_locator()
    active_year_months = defaultdict(list)
    # Scan a reasonable range of years
    for year in range(datetime.date.today().year - 10, datetime.date.today().year + 1):
        for month in range(1, 13):
            try:
                # A quick search for any file in the month
                files = locator.search(
                    "ARGUS",
                    "rhi",
                    campaign=campaign,
                    year=year,
                    mon=month,
                    day="**",
                    hour="**",
                )
                if files:
                    active_year_months[year].append(month)
            except Exception:
                continue
    return active_year_months

def get_camera_configs():
    """Defines and caches available camera configurations."""
    # You can expand this dictionary with more camera setups
    return {
        "Single Camera (3-7)": UndistortedCamera.from_config("COBALT", "3-7"),
        # "Camera with Aircraft Tracks (3-7)": AutomaticADSBAircraftInterface(UndistortedCamera.from_config("COBALT", "3-7")),
        "Multicam 2 (5 cameras)": CameraArray(
            [
                UndistortedCamera.from_config("COBALT", "3-7"),
                UndistortedCamera.from_config("COBALT", "5-1"),
                UndistortedCamera.from_config("COBALT", "5-2"),
                UndistortedCamera.from_config("COBALT", "5-3"),
                UndistortedCamera.from_config("COBALT", "5-4"),
            ],
            (3, 3),
        ),
        # Add other CameraArray or single Camera objects here
    }


def find_available_dates_in_month(campaign, year, month):
    """Scans a whole month to find which dates have radar data."""
    locator = initialise_locator()
    available_dates = []
    num_days = calendar.monthrange(year, month)[1]

    for day in range(1, num_days + 1):
        # This is the potentially slow part, but caching helps.
        try:
            files = locator.search(
                "ARGUS",
                "rhi",
                campaign=campaign,
                year=year,
                mon=month,
                day=day,
                hour="**",
            )
            if files:
                available_dates.append(datetime.date(year, month, day))
        except Exception:
            # Ignore errors during search for a single day
            continue
    return available_dates

def find_next_available_date(campaign, start_date):
    """Finds the next date with data, starting from the day after start_date."""
    locator = initialise_locator()
    current_date = start_date + datetime.timedelta(days=1)
    # Limit search to 1 year to prevent infinite loops
    for _ in range(365):
        try:
            files = locator.search(
                "ARGUS",
                "rhi",
                campaign=campaign,
                year=current_date.year,
                mon=current_date.month,
                day=current_date.day,
                hour="**",
            )
            if files:
                return current_date
        except Exception:
            pass  # Ignore errors and continue to the next day
        current_date += datetime.timedelta(days=1)
    return None

def get_available_scans_for_day(
    campaign, scan_type, selected_date, start_time, end_time
):
    """
    Finds all available radar scan start times for a given day within a specified
    UTC time window.
    """
    if not selected_date:
        return []

    locator = initialise_locator()

    try:
        all_day_files = locator.search(
            "ARGUS",
            scan_type,
            campaign=campaign,
            year=selected_date.year,
            mon=selected_date.month,
            day=selected_date.day,
            hour="**",
        )
    except Exception as e:
        return []

    if not all_day_files:
        return []

    # --- Time Filtering (UTC) ---
    utc_start = datetime.datetime.combine(selected_date, start_time)
    utc_end = datetime.datetime.combine(selected_date, end_time)

    scan_times = []
    for f_path in sorted(all_day_files):
        filename_base = os.path.basename(f_path)
        match = re.search(r"(\d{4}\d{2}\d{2}-\d{2}\d{2}\d{2})", filename_base)
        if match:
            try:
                # The file time is already in UTC, and it's "naive"
                scan_dt_utc = datetime.datetime.strptime(match.group(1), "%Y%m%d-%H%M%S")

                # Check if the scan is within our UTC time bounds
                if utc_start <= scan_dt_utc <= utc_end:
                    scan_times.append(scan_dt_utc)
            except ValueError:
                continue

    return scan_times

