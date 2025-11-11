"""
Handles locating and loading camera data loader plugins.

This is used only to extend arguslib with custom implementations of
CameraData, e.g. for different systems which save data in formats other than
the arguslib standard of timestamped MP4 files.

File composed with input from Google Gemini.
"""

import importlib.metadata
import sys
from .locator import CameraData as DefaultCameraData

# A cache for the registry, so we don't search every time
_LOADER_REGISTRY = None


def _build_loader_registry():
    """
    Finds all installed data loader plugins and builds a registry.
    Maps the registered name (e.g., 'sirta') to the loader class.
    """
    global _LOADER_REGISTRY
    _LOADER_REGISTRY = {}

    # Use the modern (3.10+) .select() method if available
    if sys.version_info >= (3, 10):
        try:
            eps = importlib.metadata.entry_points(group="arguslib.data_loaders")
        except Exception:
            eps = []  # No plugins found
    # Use the older (3.8, 3.9) .get() method as a fallback
    else:
        try:
            eps = importlib.metadata.entry_points().get("arguslib.data_loaders", [])
        except Exception:
            eps = []

    for ep in eps:
        try:
            # We store the function to load, not the loaded class
            # This avoids importing plugins until they're needed
            _LOADER_REGISTRY[ep.name.lower()] = ep.load
        except Exception as e:
            # Handle a broken plugin installation
            print(f"Warning: Failed to load arguslib plugin '{ep.name}': {e}")


def get_data_loader_class(name: str):
    """
    Gets a data loader class by its registered name (e.g., "sirta").

    If the name isn't found in the plugin registry, it returns the
    default CameraData loader.
    """
    if _LOADER_REGISTRY is None:
        _build_loader_registry()

    # Find the loader function by its registered name (case-insensitive)
    loader_func = _LOADER_REGISTRY.get(name.lower())

    if loader_func:
        try:
            # Call the .load() function now, which does the import
            return loader_func()
        except ImportError as e:
            # Plugin is registered but not properly installed
            print(f"Warning: Could not import plugin for '{name}': {e}")
            # Fall through to default

    # Fallback to the default
    return DefaultCameraData
