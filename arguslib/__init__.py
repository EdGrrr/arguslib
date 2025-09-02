"""
Arguslib
========

A library for bringing together and visualising different atmospheric observations,
including ground-based cameras, radar, and aircraft data.

This top-level __init__ file exposes the primary user-facing classes
for easier access.
"""

# --- Core Instrument and Geolocation Classes ---
from .instruments.instruments import Position

# --- Camera Classes ---
from .camera.camera import Camera
from .camera.undistorted_camera import UndistortedCamera
from .camera.direct_camera import DirectCamera, DirectUndistortedCamera
from .camera.camera_array import CameraArray

# --- Radar Classes ---
from .radar.radar import Radar

# --- Aircraft Classes ---
from .aircraft.fleet import Fleet

# --- High-Level Interface Classes ---
from .aircraft.aircraft_interface import (
    AircraftInterface,
    AutomaticADSBAircraftInterface,
)
from .radar.radar_interface import RadarInterface
from .camera.video_interface import VideoInterface
