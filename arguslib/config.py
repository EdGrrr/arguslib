"""
Handles loading and merging of configuration files from standard locations.
"""

from pathlib import Path
from typing import Any, Dict, Iterable, BinaryIO, Union, Tuple, List
import yaml
import importlib.resources as ilres
import importlib.metadata as ilmeta

# Define the standard search paths in one place.
# The order is from least specific (system-wide) to most specific (user).
CONFIG_SEARCH_PATHS = [
    Path("/etc/arguslib/"),
    Path.home() / ".arguslib",
    Path.home() / ".config" / "arguslib",
]


def _iter_packaged_default_streams(filename: str) -> Iterable[BinaryIO]:
    """
    Yields open binary streams for default config files packaged with arguslib
    and any installed plugins advertising a defaults directory.
    """
    # 0) Core package defaults
    try:
        base = ilres.files("arguslib") / "defaults" / filename
        if hasattr(base, "is_file") and base.is_file():
            yield base.open("rb")
    except Exception:
        pass

    # 1) Plugin-provided defaults via entry points
    try:
        eps_all = ilmeta.entry_points()
        eps = (
            eps_all.select(group="arguslib.config_sets")
            if hasattr(eps_all, "select")
            else eps_all.get("arguslib.config_sets", [])
        )
    except Exception:
        eps = []

    for ep in eps:
        try:
            obj = ep.load()
            if callable(obj):
                obj = obj()

            if isinstance(obj, str):
                pkg, _, subdir = obj.partition(":")
                subdir = subdir or "defaults"
                p = ilres.files(pkg) / subdir / filename
            else:
                p = obj / filename  # Traversable

            if hasattr(p, "is_file") and p.is_file():
                yield p.open("rb")
        except Exception:
            continue


def _deep_copy_dict(d: dict) -> dict:
    """Lightweight deep copy for nested dicts."""
    out = {}
    for k, v in d.items():
        out[k] = _deep_copy_dict(v) if isinstance(v, dict) else v
    return out


def _deep_update(dst: dict, src: dict) -> dict:
    """
    Recursively merge src into dst.

    Special 'defaults' key logic:
    If a dict contains a 'defaults' key, those values are applied to all
    sibling dictionary entries before merging the specific values.
    """
    # 1. Apply 'defaults' from src to src's children (intra-file defaults)
    src_defaults = src.get("defaults")
    if isinstance(src_defaults, dict):
        for k, v in src.items():
            if k == "defaults":
                continue
            if isinstance(v, dict):
                # Create a new dict that is defaults + specific overrides
                src[k] = _deep_copy_dict(src_defaults) | v

    # 2. Apply 'defaults' from dst to src's children (inter-file inheritance)
    # This ensures that defaults defined in base configs apply to items in user configs
    dst_defaults = dst.get("defaults")
    if isinstance(dst_defaults, dict):
        for k, v in src.items():
            if k == "defaults":
                continue
            if isinstance(v, dict):
                # Apply inherited defaults, but let src's values (including its own defaults) win
                src[k] = _deep_copy_dict(dst_defaults) | src[k]

    # 3. Merge src into dst
    for k, v in src.items():
        if k == "defaults":
            # Persist defaults into dst so they are available for subsequent config files
            if k not in dst:
                dst[k] = {}
            if isinstance(v, dict):
                _deep_update(dst[k], v)
            continue

        if isinstance(v, dict):
            # Force recursion to ensure defaults expansion happens within v
            if k not in dst or not isinstance(dst[k], dict):
                dst[k] = {}
            _deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def load_config(filename: str) -> Dict[str, Any]:
    """
    Loads configuration from all sources (defaults, plugins, user) and performs
    a deep merge. Later sources override earlier ones.
    """
    merged: Dict[str, Any] = {}
    found_any = False

    # 1. Load from Packaged Defaults & Plugins
    for stream in _iter_packaged_default_streams(filename):
        try:
            doc = yaml.safe_load(stream)
            if isinstance(doc, dict) and doc:
                found_any = True
                _deep_update(merged, doc)
        finally:
            try:
                stream.close()
            except Exception:
                pass

    # 2. Load from System/User Paths
    for config_dir in CONFIG_SEARCH_PATHS:
        config_file = config_dir / filename
        if config_file.exists():
            with open(config_file, "r") as f:
                doc = yaml.safe_load(f)
                if isinstance(doc, dict) and doc:
                    found_any = True
                    _deep_update(merged, doc)

    if not found_any:
        raise FileNotFoundError(f"No configuration file named '{filename}' found.")

    return merged


def list_cameras() -> List[Tuple[str, str]]:
    """
    Returns a list of (campaign, camera_id) pairs available in cameras.yml.
    """
    data = load_config("cameras.yml")
    out: List[Tuple[str, str]] = []
    for campaign, cams in data.items():
        if isinstance(cams, dict):
            for cam_id in cams.keys():
                if cam_id != "defaults":
                    out.append((campaign, cam_id))
    return sorted(out)


def resolve_config_resource(spec: Union[str, Path]) -> object:
    """
    Resolve a config spec (str/Path) to a Path or Traversable resource.
    """
    p = Path(str(spec)).expanduser()

    # Absolute path or existing relative path
    if p.is_absolute() and p.exists():
        return p
    if p.exists():
        return p

    # User config dirs
    for cfg_dir in CONFIG_SEARCH_PATHS:
        candidate = cfg_dir / p.name
        if candidate.exists():
            return candidate

    # Packaged defaults
    # We reuse the stream iterator logic but just return the resource
    # (Re-implementing simplified resource finder to avoid stream opening)
    try:
        base = ilres.files("arguslib") / "defaults" / p.name
        if hasattr(base, "is_file") and base.is_file():
            return base
    except Exception:
        pass

    # Fallback
    return p
