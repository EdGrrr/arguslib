"""
Arguslib
========

A library for bringing together and visualising different atmospheric observations,
including ground-based cameras, radar, and aircraft data.

This top-level __init__ file exposes the primary user-facing classes
for easier access.
"""

# --- Core Instrument and Geolocation Classes ---
from .instruments import Position

# --- Camera Classes ---
from .camera import (
    Camera,
    UndistortedCamera,
    DirectCamera,
    DirectUndistortedCamera,
    CameraArray,
    VideoInterface,
)

# --- Radar Classes ---
from .radar import Radar, RadarInterface

# --- Aircraft Classes ---
from .aircraft import (
    AircraftInterface,
    AutomaticADSBAircraftInterface,
)
from .aircraft import Fleet
from .mapping import MapInstrument
