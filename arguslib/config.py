"""
Handles loading and merging of configuration files from standard locations.
"""

from pathlib import Path
import yaml
from typing import Any, Dict

# Define the standard search paths in one place.
# The order is from least specific (system-wide) to most specific (user).
CONFIG_SEARCH_PATHS = [
    Path("/etc/arguslib/"),
    Path.home() / ".arguslib",
    Path.home() / ".config" / "arguslib",
]


def load_config(filename: str) -> Dict[str, Any]:
    """
    Loads and merges YAML configuration from standard locations.

    It searches for `filename` in the defined `CONFIG_SEARCH_PATHS`.
    Configurations are merged, with values from files found later in the
    search path (i.e., user-specific) overriding earlier ones (system-wide).

    Args:
        filename: The name of the YAML file (e.g., 'cameras.yml').

    Returns:
        A dictionary containing the merged configuration.

    Raises:
        FileNotFoundError: If no configuration file is found in any of the
                           search paths.
    """
    configs = []
    for config_dir in CONFIG_SEARCH_PATHS:
        config_file = config_dir / filename
        if config_file.exists():
            with open(config_file, "r") as f:
                loaded_yaml = yaml.safe_load(f)
                if loaded_yaml:  # Ensure file is not empty
                    configs.append(loaded_yaml)

    if not configs:
        raise FileNotFoundError(
            f"No configuration file named '{filename}' found in search paths."
        )

    # Merge configs. The last one found (most specific) wins.
    merged_config = {}
    for config in configs:
        merged_config.update(config)

    return merged_config


def load_path_from_config(filename: str) -> Path:
    """Loads a single path from a text file in standard config locations."""
    # Search in reverse order to find user-specific file first
    for config_dir in reversed(CONFIG_SEARCH_PATHS):
        config_file = config_dir / filename
        if config_file.exists():
            return Path(config_file.read_text().strip())

    raise FileNotFoundError(f"Path configuration file '{filename}' not found.")
