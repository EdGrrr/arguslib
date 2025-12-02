"""
Handles loading and merging of configuration files from standard locations.
"""

from pathlib import Path
from typing import Any, Dict
from typing import Iterable, BinaryIO
from typing import Union
from typing import Iterable, Tuple, List
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

    Entry point group for plugins:
      - group: 'arguslib.config_sets'
      - each entry returns either:
          * a Traversable directory (importlib.resources.files(...)/"defaults"), or
          * a string 'package:subdir' (subdir defaults to 'defaults' if omitted)
    """
    # 0) Core package defaults: arguslib/defaults/<filename>
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
            # NEW: call provider if it is a function
            if callable(obj):
                obj = obj()
        except Exception:
            continue

        try:
            if isinstance(obj, str):
                pkg, _, subdir = obj.partition(":")
                subdir = subdir or "defaults"
                p = ilres.files(pkg) / subdir / filename
                if hasattr(p, "is_file") and p.is_file():
                    yield p.open("rb")
            else:
                p = obj / filename  # Traversable
                if hasattr(p, "is_file") and p.is_file():
                    yield p.open("rb")
        except Exception:
            continue


def _iter_packaged_default_resources(filename: str) -> Iterable[object]:
    """
    Yields Traversable resources for packaged defaults: core and plugin-provided.
    """
    # Core defaults
    try:
        base = ilres.files("arguslib") / "defaults" / filename
        if hasattr(base, "is_file") and base.is_file():
            yield base
    except Exception:
        pass

    # Plugin defaults via entry points
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
            # NEW: call provider if it is a function
            if callable(obj):
                obj = obj()
        except Exception:
            continue

        try:
            if isinstance(obj, str):
                pkg, _, subdir = obj.partition(":")
                subdir = subdir or "defaults"
                p = ilres.files(pkg) / subdir / filename
                if hasattr(p, "is_file") and p.is_file():
                    yield p
            else:
                p = obj / filename  # Traversable
                if hasattr(p, "is_file") and p.is_file():
                    yield p
        except Exception:
            continue


def load_config(filename: str) -> Dict[str, Any]:
    """
    Loads and merges YAML configuration (shallow merge; later wins).
    """
    configs: list[dict] = []
    # 0) Packaged defaults
    for stream in _iter_packaged_default_streams(filename):
        try:
            doc = yaml.safe_load(stream)
            if isinstance(doc, dict) and doc:
                configs.append(doc)
        finally:
            try:
                stream.close()
            except Exception:
                pass
    # 1) System/user overrides
    for config_dir in CONFIG_SEARCH_PATHS:
        config_file = config_dir / filename
        if config_file.exists():
            with open(config_file, "r") as f:
                doc = yaml.safe_load(f)
                if isinstance(doc, dict) and doc:
                    configs.append(doc)
    if not configs:
        raise FileNotFoundError(f"No configuration file named '{filename}' found.")
    merged: Dict[str, Any] = {}
    for c in configs:
        merged.update(c)
    return merged


def _iter_config_dicts(filename: str) -> Iterable[dict]:
    """
    Yields config dicts in precedence order:
      1) core packaged defaults, then plugin defaults
      2) system/user overrides (in CONFIG_SEARCH_PATHS order)
    """
    # 0) Packaged defaults
    for stream in _iter_packaged_default_streams(filename):
        try:
            doc = yaml.safe_load(stream)
            if isinstance(doc, dict) and doc:
                yield doc
        finally:
            try:
                stream.close()
            except Exception:
                pass
    # 1) System/user overrides
    for config_dir in CONFIG_SEARCH_PATHS:
        config_file = config_dir / filename
        if config_file.exists():
            with open(config_file, "r") as f:
                doc = yaml.safe_load(f)
                if isinstance(doc, dict) and doc:
                    yield doc


def _deep_update(dst: dict, src: dict) -> dict:
    """
    Recursively merge src into dst (maps are merged, scalars/lists replaced).
    Returns dst.
    """
    for k, v in src.items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def load_config_deep(filename: str) -> Dict[str, Any]:
    """
    Deep-merge all sources so nested maps are unioned (later wins per key).
    Useful for discovery/listing without losing default entries.
    """
    merged: Dict[str, Any] = {}
    any_found = False
    for doc in _iter_config_dicts(filename):
        any_found = True
        _deep_update(merged, doc)
    if not any_found:
        raise FileNotFoundError(f"No configuration file named '{filename}' found.")
    return merged


def list_cameras(deep_merge: bool = True) -> List[Tuple[str, str]]:
    """
    Returns a list of (campaign, camera_id) pairs available in cameras.yml.
    If deep_merge=True, unions defaults, plugins, system, and user configs.
    """
    data = load_config_deep("cameras.yml") if deep_merge else load_config("cameras.yml")
    out: List[Tuple[str, str]] = []
    for campaign, cams in (data or {}).items():
        if isinstance(cams, dict):
            out.extend((campaign, cam_id) for cam_id in cams.keys())
    return sorted(out)


def list_instrument_configs(filename: str, deep_merge: bool = True) -> Dict[str, Any]:
    """
    Generic accessor to fetch an instrument config mapping (deep-merged by default).
    Example: list_instrument_configs("cameras.yml")
    """
    return load_config_deep(filename) if deep_merge else load_config(filename)


def load_path_from_config(filename: str) -> Path:
    """Loads a single path from a text file in standard config locations."""
    # User-first search
    for config_dir in reversed(CONFIG_SEARCH_PATHS):
        config_file = config_dir / filename
        if config_file.exists():
            return Path(config_file.read_text().strip())
    # Packaged defaults
    for stream in _iter_packaged_default_streams(filename):
        try:
            content = stream.read().decode("utf-8")
            if content:
                return Path(content.strip())
        finally:
            try:
                stream.close()
            except Exception:
                pass
    raise FileNotFoundError(f"Path configuration file '{filename}' not found.")


def resolve_config_resource(spec: Union[str, Path]) -> object:
    """
    Resolve a config spec (str/Path) to either:
      - a filesystem Path, or
      - a Traversable resource (importlib.resources)
    Search order: user/system config dirs, then packaged defaults (core, plugins).
    """
    p = Path(str(spec)).expanduser()
    if p.is_absolute() and p.exists():
        return p
    if p.exists():
        return p
    for cfg_dir in CONFIG_SEARCH_PATHS:
        candidate = cfg_dir / p.name
        if candidate.exists():
            return candidate
    for res in _iter_packaged_default_resources(p.name):
        return res
    return p


def debug_config_resolution(filename: str) -> Dict[str, Any]:
    """
    Returns details about where '<filename>' would be searched/loaded from:
    - CONFIG_SEARCH_PATHS on disk
    - core packaged defaults
    - plugin-provided defaults via entry points
    """
    info: Dict[str, Any] = {
        "filename": filename,
        "CONFIG_SEARCH_PATHS": [str(p) for p in CONFIG_SEARCH_PATHS],
        "core_packaged_default": None,
        "plugin_providers": [],
        "filesystem_hits": [],
    }
    # Core packaged default
    try:
        base = ilres.files("arguslib") / "defaults" / filename
        if hasattr(base, "is_file") and base.is_file():
            info["core_packaged_default"] = str(base)
    except Exception as e:
        info["core_packaged_default"] = f"error: {e!r}"

    # Plugin defaults via entry points
    try:
        eps_all = ilmeta.entry_points()
        eps = (
            eps_all.select(group="arguslib.config_sets")
            if hasattr(eps_all, "select")
            else eps_all.get("arguslib.config_sets", [])
        )
    except Exception as e:
        eps = []
        info["plugin_ep_error"] = repr(e)

    for ep in eps:
        ep_info = {
            "name": getattr(ep, "name", None),
            "value": getattr(ep, "value", None),
            "module": getattr(ep, "module", None),
            "attr": getattr(ep, "attr", None),
        }
        try:
            obj = ep.load()
            if callable(obj):
                obj = obj()
            ep_info["loaded"] = True
            if isinstance(obj, str):
                pkg, _, subdir = obj.partition(":")
                subdir = subdir or "defaults"
                d = ilres.files(pkg) / subdir
            else:
                d = obj  # Traversable
            ep_info["defaults_dir"] = str(d)
            f = d / filename
            is_file = hasattr(f, "is_file") and f.is_file()
            ep_info["contains_file"] = is_file
            if is_file:
                ep_info["resource"] = str(f)
        except Exception as e:
            ep_info["loaded"] = False
            ep_info["error"] = repr(e)
        info["plugin_providers"].append(ep_info)

    # Filesystem hits
    for config_dir in CONFIG_SEARCH_PATHS:
        cf = config_dir / filename
        if cf.exists():
            info["filesystem_hits"].append(str(cf))

    return info


def print_config_debug(filename: str) -> None:
    """Pretty-print debug_config_resolution(filename)."""
    import pprint

    pprint.pprint(debug_config_resolution(filename))
