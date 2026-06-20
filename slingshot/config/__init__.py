"""Compatibility wrapper that preserves v3 models and fixes safe serialization."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import yaml


_SOURCE = Path(__file__).resolve().parent.parent / "config.py"
_MODULE_NAME = "slingshot._legacy_config_impl"
_SPEC = importlib.util.spec_from_file_location(_MODULE_NAME, _SOURCE)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError(f"Unable to load legacy config module from {_SOURCE}")
_LEGACY = importlib.util.module_from_spec(_SPEC)
sys.modules[_MODULE_NAME] = _LEGACY
_SPEC.loader.exec_module(_LEGACY)

for _name in dir(_LEGACY):
    if not _name.startswith("_"):
        globals()[_name] = getattr(_LEGACY, _name)


def save_config(config, output_path: str, format: str = "yaml") -> None:
    """Serialize legacy configs using JSON-compatible values and safe YAML."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(config, "model_dump"):
        data = config.model_dump(mode="json")
    else:
        data = json.loads(config.json())
    with path.open("w", encoding="utf-8") as stream:
        if format.lower() in {"yaml", "yml"}:
            yaml.safe_dump(data, stream, default_flow_style=False, sort_keys=False)
        elif format.lower() == "json":
            json.dump(data, stream, indent=2)
        else:
            raise ValueError(f"Unsupported format: {format}")


__all__ = [
    name for name in globals() if not name.startswith("_") and name != "annotations"
]
