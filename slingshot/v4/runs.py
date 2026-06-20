"""Version-aware run discovery that excludes legacy science from v4 aggregates."""

from __future__ import annotations

import json
from pathlib import Path


def classify_run(path: str | Path) -> dict:
    run_path = Path(path)
    manifest_path = run_path / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        schema_version = int(manifest.get("schema_version", 0))
        legacy = bool(manifest.get("legacy_science_model", schema_version < 4))
        return {
            "path": str(run_path),
            "schema_version": schema_version,
            "legacy_science_model": legacy,
            "manifest": manifest,
        }
    legacy_markers = any(
        (run_path / name).exists()
        for name in ("results.pkl", "summary.csv", "REPORT.md")
    )
    return {
        "path": str(run_path),
        "schema_version": 3 if legacy_markers else 0,
        "legacy_science_model": legacy_markers,
        "manifest": None,
    }


def discover_v4_runs(root: str | Path = "results") -> list[Path]:
    root_path = Path(root)
    if not root_path.exists():
        return []
    runs = []
    for path in sorted(item for item in root_path.iterdir() if item.is_dir()):
        classification = classify_run(path)
        if (
            classification["schema_version"] == 4
            and not classification["legacy_science_model"]
        ):
            runs.append(path)
    return runs
