import streamlit as st
import datetime as dt
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
from pathlib import Path

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
        "Single Camera (3-7)": UndistortedCamera.from_config("COBALT", "3-7"),
        # Add other CameraArray or single Camera objects here
    }

# This function is not cached with the main resource cache because it depends
# on the date selected by the user.
@st.cache_data
def get_available_scans_for_day(campaign, scan_type, selected_date):
    """
    Finds all available radar scan start times for a given day.
    This function is based on the logic in your locator.py file.
    """
    st.info(f"Searching for {scan_type} scans on {selected_date.strftime('%Y-%m-%d')}...")
    
    # This locator is used just for searching, so it's lightweight.
    locator = initialise_locator()
    
    try:
        # Use the locator to find all files for the given day
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
        st.warning(f"Could not search for files: {e}")
        return []

    if not all_day_files:
        st.warning(f"No radar files found for {selected_date.strftime('%Y-%m-%d')}.")
        return []

    # Extract datetimes from filenames
    scan_times = []
    for f_path in sorted(all_day_files):
        filename_base = os.path.basename(f_path)
        # Regex from your locator.py to find the timestamp
        match = re.search(r'(\d{4}\d{2}\d{2}-\d{2}\d{2}\d{2})', filename_base)
        if match:
            try:
                file_start_dt = dt.datetime.strptime(match.group(1), "%Y%m%d-%H%M%S")
                scan_times.append(file_start_dt)
            except ValueError:
                continue # Skip if filename format is unexpected

    st.success(f"Found {len(scan_times)} scans.")
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

st.set_page_config(layout="wide")
st.title("🛰️ Argus Radar & Camera Interface Explorer")

# --- Sidebar for Controls ---
with st.sidebar:
    st.header("Controls")

    # 1. Select Campaign
    # Hardcoded for now, but could be dynamic
    campaign = st.selectbox("Select Campaign", ["COBALT"])

    # 2. Select Date
    selected_date = st.date_input("Select Date", dt.date(2025, 5, 31))

    # 3. Find and Select Scan Time
    scan_times = get_available_scans_for_day(campaign, "rhi", selected_date)

    if not scan_times:
        st.warning("No scans available for this date. Please select another day.")
        st.stop()

    # We display the time but the value is the full datetime object
    selected_scan_dt = st.selectbox(
        "Select Scan Time",
        options=scan_times,
        format_func=lambda d: d.strftime("%H:%M:%S")
    )

    # 4. Select Camera
    camera_configs = get_camera_configs()
    camera_choice_name = st.selectbox("Select Camera Configuration", camera_configs.keys())
    selected_camera = camera_configs[camera_choice_name]

    # 5. Generate Button
    generate_plot = st.button("Generate Visualization", type="primary")
    
    # 6. Save Case
    st.divider()
    st.header("Save Case")
    save_filename = st.text_input("Save file name", "good_cases.csv")
    if st.button("Save Current Case"):
        if selected_scan_dt:
            save_case(save_filename, selected_scan_dt, camera_choice_name)


# --- Main Display Area ---
if generate_plot and selected_scan_dt:
    with st.spinner("Generating plot... This may take a moment."):
        try:
            # Initialize the core arguslib objects
            radar = get_radar(campaign)
            
            # Create the interface object
            cri = RadarInterface(radar, selected_camera)
            
            # Use the show method from radar_interface.py
            ax = cri.show(selected_scan_dt, var='DBZ')
            
            fig = ax[-1].get_figure()
            # Display in Streamlit
            st.pyplot(fig)

        except ValueError as e:
            st.error(f"Error during plot generation: {e}")
            st.info("This often happens if the requested datetime does not exactly match a scan start time. The dropdown list should prevent this, but data inconsistencies can occur.")
        except FileNotFoundError as e:
            st.error(f"Data file not found: {e}")
            st.info("Please ensure the data paths in your configuration are correct and the files exist.")
        except Exception as e:
            st.error(f"An unexpected error occurred: {e}")

elif not selected_scan_dt:
    st.info("Please select a scan time from the sidebar and click 'Generate Visualization'.")

# Display saved cases if the file exists
if Path(save_filename).exists():
    with st.expander("View Saved Cases"):
        df_saved = pd.read_csv(save_filename)
        st.dataframe(df_saved)
