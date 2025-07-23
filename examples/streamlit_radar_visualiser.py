import streamlit as st
import datetime as dt
from datetime import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
from pathlib import Path
import random
import calendar
from collections import defaultdict

# --- Important: Make sure arguslib is in your Python Path ---
# You might need to do something like this at the start of your script
# if the `arguslib` directory is in the same folder as this app:
# import sys
# sys.path.append(os.path.abspath('.'))

# --- Importing your Arguslib Classes ---
# Assuming the file structure allows for these imports
try:
    from arguslib.camera.camera_array import CameraArray
    from arguslib.camera.undistorted_camera import UndistortedCamera
    from arguslib.radar.radar import Radar
    from arguslib.radar.radar_interface import RadarInterface
    # We need the locator to find scan times
    from arguslib.radar.locator import RadarData, initialise_locator
except ImportError as e:
    st.error(f"Failed to import arguslib modules: {e}")
    st.info("Please ensure that the 'arguslib' directory is in your Python path or in the same directory as this Streamlit app.")
    st.stop()


# --- App Configuration ---
# Use Streamlit's caching to avoid reloading objects on every interaction.
# This is crucial for performance.
@st.cache_resource
def get_radar(campaign_name):
    """Loads and caches the Radar object."""
    try:
        return Radar.from_config(campaign_name)
    except FileNotFoundError as e:
        st.error(f"Could not find radar configuration for '{campaign_name}'. Please check your config files.")
        st.stop()
    except Exception as e:
        st.error(f"An error occurred while loading the radar: {e}")
        st.stop()


@st.cache_resource
def get_camera_configs():
    """Defines and caches available camera configurations."""
    # You can expand this dictionary with more camera setups
    return {
        "Single Camera (3-7)": UndistortedCamera.from_config("COBALT", "3-7"),
        "Multicam 2 (5 cameras)": CameraArray(
            [
                UndistortedCamera.from_config("COBALT", "3-7"),
                UndistortedCamera.from_config("COBALT", "5-1"),
                UndistortedCamera.from_config("COBALT", "5-2"),
                UndistortedCamera.from_config("COBALT", "5-3"),
                UndistortedCamera.from_config("COBALT", "5-4"),
            ],
            (3, 3)
        ),
        # Add other CameraArray or single Camera objects here
    }

@st.cache_data
def find_active_year_months(campaign):
    """Scans for all years and months that have data."""
    locator = initialise_locator()
    active_year_months = defaultdict(list)
    # Scan a reasonable range of years
    for year in range(dt.date.today().year - 10, dt.date.today().year + 1):
        for month in range(1, 13):
            try:
                # A quick search for any file in the month
                files = locator.search(
                    "ARGUS", "rhi", campaign=campaign,
                    year=year, mon=month, day="**", hour="**"
                )
                if files:
                    active_year_months[year].append(month)
            except Exception:
                continue
    return active_year_months

@st.cache_data
def find_available_dates_in_month(campaign, year, month):
    """Scans a whole month to find which dates have radar data."""
    locator = initialise_locator()
    available_dates = []
    num_days = calendar.monthrange(year, month)[1]

    for day in range(1, num_days + 1):
        # This is the potentially slow part, but caching helps.
        try:
            files = locator.search(
                "ARGUS", "rhi", campaign=campaign,
                year=year, mon=month, day=day, hour="**"
            )
            if files:
                available_dates.append(dt.date(year, month, day))
        except Exception:
            # Ignore errors during search for a single day
            continue
    return available_dates

# This function is not cached with the main resource cache because it depends
# on the user-selected inputs.
@st.cache_data
def get_available_scans_for_day(campaign, scan_type, selected_date, start_time, end_time):
    """
    Finds all available radar scan start times for a given day within a specified
    UTC time window.
    """
    if not selected_date:
        return []
        
    st.info(f"Searching for {scan_type} scans on {selected_date.strftime('%Y-%m-%d')}...")
    
    locator = initialise_locator()
    
    try:
        all_day_files = locator.search(
            "ARGUS", scan_type, campaign=campaign,
            year=selected_date.year, mon=selected_date.month, day=selected_date.day, hour="**",
        )
    except Exception as e:
        st.warning(f"Could not search for files: {e}")
        return []

    if not all_day_files:
        st.warning(f"No radar files found for {selected_date.strftime('%Y-%m-%d')}.")
        return []

    # --- Time Filtering (UTC) ---
    utc_start = dt.datetime.combine(selected_date, start_time)
    utc_end = dt.datetime.combine(selected_date, end_time)

    scan_times = []
    for f_path in sorted(all_day_files):
        filename_base = os.path.basename(f_path)
        match = re.search(r'(\d{4}\d{2}\d{2}-\d{2}\d{2}\d{2})', filename_base)
        if match:
            try:
                # The file time is already in UTC, and it's "naive"
                scan_dt_utc = dt.datetime.strptime(match.group(1), "%Y%m%d-%H%M%S")
                
                # Check if the scan is within our UTC time bounds
                if utc_start <= scan_dt_utc <= utc_end:
                    scan_times.append(scan_dt_utc)
            except ValueError:
                continue

    st.success(f"Found {len(scan_times)} scans between {start_time.strftime('%H:%M')} and {end_time.strftime('%H:%M')} UTC.")
    return scan_times


def save_case(save_path, selected_dt, camera_name):
    """Appends the details of a 'good' case to a CSV file."""
    case_data = {
        "timestamp": [selected_dt.isoformat()],
        "camera_config": [camera_name],
        "saved_at": [dt.datetime.now().isoformat()]
    }
    df_new = pd.DataFrame(case_data)

    if Path(save_path).exists():
        df_existing = pd.read_csv(save_path)
        df_combined = pd.concat([df_existing, df_new], ignore_index=True)
    else:
        df_combined = df_new

    df_combined.to_csv(save_path, index=False)
    st.success(f"Case saved to {save_path}")


# --- Streamlit User Interface ---

# --- Default/Initial View Configuration ---
DEFAULT_YEAR = 2025
DEFAULT_MONTH = 5
DEFAULT_DAY = 31
DEFAULT_SCAN_TIME_STR = "14:53:37" # H:M:S format

st.set_page_config(layout="wide")
st.title("🛰️ Argus Radar & Camera Interface Explorer")

# --- Sidebar for Controls ---
with st.sidebar:
    st.header("Controls")

    # 1. Select Campaign
    campaign = st.selectbox("Select Campaign", ["COBALT"])

    # 2. New Date Selection Method
    st.subheader("Date Selection")
    
    # --- One-time initialization from URL params or defaults ---
    if 'init_done' not in st.session_state:
        st.session_state.init_done = True
        params = st.query_params
        try:
            year = int(params.get("year", DEFAULT_YEAR))
            month = int(params.get("month", DEFAULT_MONTH))
            day = int(params.get("day", DEFAULT_DAY))
            scan_str = params.get("scan", DEFAULT_SCAN_TIME_STR)
            
            # Store the desired initial state. We will validate it against available data later.
            st.session_state.selected_year = year
            st.session_state.selected_month = month
            st.session_state.selected_date = dt.date(year, month, day)
            st.session_state.selected_scan_dt = dt.datetime.combine(
                dt.date(year, month, day), 
                dt.time.fromisoformat(scan_str)
            )
        except (ValueError, TypeError):
            # If params are bad, clear any potentially bad state
            for key in ['selected_year', 'selected_month', 'selected_date', 'selected_scan_dt']:
                if key in st.session_state:
                    del st.session_state[key]
    
    with st.spinner("Finding available years and months..."):
        active_year_months = find_active_year_months(campaign)

    if not active_year_months:
        st.error("No data found for this campaign in the last 10 years.")
        st.stop()
    
    active_years = sorted(active_year_months.keys(), reverse=True)

    # --- Validate and Set Year, Month, Date, and Scan ---
    # Validate Year or set to default
    if 'selected_year' not in st.session_state or st.session_state.selected_year not in active_years:
        st.session_state.selected_year = active_years[0]
    
    # Validate Month or set to default
    active_months_for_year = active_year_months[st.session_state.selected_year]
    if 'selected_month' not in st.session_state or st.session_state.selected_month not in active_months_for_year:
        st.session_state.selected_month = active_months_for_year[0]

    col_date1, col_date2 = st.columns(2)
    with col_date1:
        st.selectbox("Year", active_years, key='selected_year')
    with col_date2:
        st.selectbox("Month", active_months_for_year, key='selected_month')

    with st.spinner("Finding dates with data..."):
        available_dates = find_available_dates_in_month(campaign, st.session_state.selected_year, st.session_state.selected_month)

    if not available_dates:
        st.warning("No data found for this month. Please select a different month/year.")
        st.stop()
    
    # Validate Date or set to default
    if 'selected_date' not in st.session_state or st.session_state.selected_date not in available_dates:
        st.session_state.selected_date = available_dates[0]

    st.selectbox(
        "Select an available date",
        options=available_dates,
        key='selected_date',
        format_func=lambda d: d.strftime("%Y-%m-%d") # Simplified date format
    )

    # 3. Time bounds
    st.subheader("Time Window (UTC)")
    col_time1, col_time2 = st.columns(2)
    with col_time1:
        start_time = st.time_input("Start time", value=time(6, 0))
    with col_time2:
        end_time = st.time_input("End time", value=time(18, 0))


    # 4. Find and Select Scan Time
    st.subheader("Scan Selection")
    scan_times = get_available_scans_for_day(campaign, "rhi", st.session_state.selected_date, start_time, end_time)

    if not scan_times:
        st.warning("No scans available for this date and time window. Please adjust the filters.")
        st.stop()

    # Validate Scan Time or set to default
    if 'selected_scan_dt' not in st.session_state or st.session_state.selected_scan_dt not in scan_times:
        st.session_state.selected_scan_dt = scan_times[0]

    # Buttons to navigate scans
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Next Scan", use_container_width=True):
            current_index = scan_times.index(st.session_state.selected_scan_dt)
            next_index = (current_index + 1) % len(scan_times)
            st.session_state.selected_scan_dt = scan_times[next_index]
            st.rerun()

    with col2:
        if st.button("Random Scan", use_container_width=True):
            random_year = random.choice(list(active_year_months.keys()))
            random_month = random.choice(active_year_months[random_year])
            random_available_dates = find_available_dates_in_month(campaign, random_year, random_month)
            if random_available_dates:
                random_date = random.choice(random_available_dates)
                random_scan_times = get_available_scans_for_day(campaign, "rhi", random_date, start_time, end_time)
                if random_scan_times:
                    st.session_state.selected_year = random_year
                    st.session_state.selected_month = random_month
                    st.session_state.selected_date = random_date
                    st.session_state.selected_scan_dt = random.choice(random_scan_times)
                    st.rerun()

    st.selectbox(
        "Select Scan Time (UTC)",
        options=scan_times,
        key='selected_scan_dt',
        format_func=lambda d: d.strftime("%H:%M:%S")
    )

    # 5. Select Camera
    st.subheader("Camera")
    camera_configs = get_camera_configs()
    if 'camera_choice_name' not in st.session_state:
        st.session_state.camera_choice_name = list(camera_configs.keys())[0]

    st.selectbox(
        "Select Camera Configuration",
        options=camera_configs.keys(),
        key='camera_choice_name'
    )

    # 6. Save Case
    st.divider()
    st.header("Save Case")
    save_filename = st.text_input("Save file name", "good_cases.csv")
    if st.button("Save Current Case"):
        save_case(save_filename, st.session_state.selected_scan_dt, st.session_state.camera_choice_name)


# --- Main Display Area ---
if st.session_state.get("selected_scan_dt"):
    selected_scan_dt_from_state = st.session_state.selected_scan_dt
    camera_configs = get_camera_configs()

    with st.spinner(f"Generating plot for {selected_scan_dt_from_state.strftime('%H:%M:%S')} UTC..."):
        try:
            radar = get_radar(campaign)
            
            # Create the interface object
            cri = RadarInterface(radar, camera_configs[st.session_state.camera_choice_name])
            
            # Use the show method from radar_interface.py
            ax = cri.show(selected_scan_dt_from_state, var='DBZ')
            
            fig = ax[-1].get_figure()
            # Display in Streamlit
            st.pyplot(fig)

        except ValueError as e:
            st.error(f"Error during plot generation: {e}")
            st.info("This often happens if the requested datetime does not exactly match a scan start time.")
        except FileNotFoundError as e:
            st.error(f"Data file not found: {e}")
            st.info("Please ensure the data paths in your configuration are correct and the files exist.")
        except Exception as e:
            st.error(f"An unexpected error occurred: {e}")

else:
    st.info("Please select a scan time from the sidebar.")

# Display saved cases if the file exists
if Path(save_filename).exists():
    with st.expander("View Saved Cases"):
        df_saved = pd.read_csv(save_filename)
        st.dataframe(df_saved)
