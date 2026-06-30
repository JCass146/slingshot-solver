"""Exploratory top-candidate artifacts for v4 campaign runs.

Candidate diagnostics are deliberately separate from the primary planar-width
estimand.  They rank finite-sample examples using current v4 energy metrics
and write a reproducible table that plots and reports can consume.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable


CANDIDATE_FIELDNAMES = [
    "rank",
    "seed",
    "sample_index",
    "v_inf_kms",
    "energy_gain_dimensionless",
    "delta_specific_energy_com",
    "delta_v_inf",
    "turning_quadratic",
    "deflection_deg",
    "periapsis_planet_km",
    "periapsis_star_km",
    "impact_parameter_km",
    "impact_parameter_au",
    "incoming_direction_rad",
    "binary_mean_anomaly_rad",
    "work_star",
    "work_planet",
    "work_sum",
    "work_energy_closure_relative",
    "solver_nfev",
    "integration_time_sec",
    "outcome",
]

_REQUIRED_FINITE = (
    "energy_gain_dimensionless",
    "delta_specific_energy_com",
)


def _as_float(row: dict, key: str, default: float = math.nan) -> float:
    value = row.get(key, default)
    if value in ("", None, "nan", "NaN"):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _as_bool(value, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    return default


def _identity_value(row: dict, key: str) -> tuple[int, float | str]:
    value = row.get(key, "")
    try:
        return (0, float(value))
    except (TypeError, ValueError):
        return (1, str(value))


def is_candidate_eligible(row: dict) -> bool:
    """Return True when a sample row is valid for exploratory ranking."""
    if str(row.get("outcome", "")).strip() != "escaped":
        return False
    if not _as_bool(row.get("solver_success", True), default=True):
        return False
    return all(math.isfinite(_as_float(row, key)) for key in _REQUIRED_FINITE)


def candidate_sort_key(row: dict) -> tuple:
    """Energy-first deterministic sort key for top-candidate ranking."""
    return (
        -_as_float(row, "energy_gain_dimensionless"),
        -_as_float(row, "delta_specific_energy_com"),
        _as_float(row, "work_energy_closure_relative", default=math.inf),
        _identity_value(row, "v_inf_kms"),
        _identity_value(row, "seed"),
        _identity_value(row, "sample_index"),
    )


def select_top_candidates(rows: Iterable[dict], top_n: int = 30) -> list[dict]:
    """Select ranked exploratory candidates from sample records."""
    eligible = [row for row in rows if is_candidate_eligible(row)]
    eligible.sort(key=candidate_sort_key)
    selected = []
    for rank, row in enumerate(eligible[: int(top_n)], start=1):
        out = {"rank": rank}
        for key in CANDIDATE_FIELDNAMES:
            if key == "rank":
                continue
            out[key] = row.get(key, "")
        selected.append(out)
    return selected


def _refine_candidate_for_run(run_path: Path, row: dict) -> dict:
    """Re-integrate one selected candidate to refresh diagnostics from current code."""
    import numpy as np

    from ..constants import G_KM
    from .config import load_config
    from .dynamics import init_binary_barycentric, integrate_encounter
    from .metrics import analyze_integration
    from .sampling import state_at_inbound_boundary
    from .validation import physical_values

    config = load_config(run_path / "config.yaml")
    values = physical_values(config)
    total_mu = G_KM * (values["star_mass_kg"] + values["planet_mass_kg"])
    position, velocity = state_at_inbound_boundary(
        _as_float(row, "v_inf_kms"),
        _as_float(row, "impact_parameter_km"),
        _as_float(row, "incoming_direction_rad"),
        values["boundary_radius_km"],
        total_mu,
    )
    binary_state = init_binary_barycentric(
        values["semi_major_axis_km"],
        config.orbit.eccentricity,
        _as_float(row, "binary_mean_anomaly_rad"),
        config.orbit.argument_periapsis_rad,
        values["star_mass_kg"],
        values["planet_mass_kg"],
        config.orbit.prograde,
        config.system.bulk_velocity_x_kms,
        config.system.bulk_velocity_y_kms,
    )
    bulk_velocity = np.array(
        [config.system.bulk_velocity_x_kms, config.system.bulk_velocity_y_kms]
    )
    initial_state = np.concatenate(
        [binary_state, position, velocity + bulk_velocity, np.zeros(2)]
    )
    integration = integrate_encounter(
        initial_state,
        values["star_mass_kg"],
        values["planet_mass_kg"],
        values["star_radius_km"],
        values["planet_radius_km"],
        values["boundary_radius_km"],
        config.numerical.max_time_sec,
        method=config.numerical.method,
        rtol=config.numerical.rtol,
        atol=config.numerical.atol,
        softening_km=config.numerical.softening_km,
        max_step_sec=config.numerical.max_step_sec,
    )
    metrics = analyze_integration(
        integration,
        initial_state,
        values["star_mass_kg"],
        values["planet_mass_kg"],
        values["semi_major_axis_km"],
        config.numerical.softening_km,
    )
    refined = dict(row)
    for key, value in metrics.items():
        if key in CANDIDATE_FIELDNAMES:
            refined[key] = value
    return refined


def _refine_candidates_for_run(run_path: Path, candidates: list[dict]) -> list[dict]:
    if not candidates or not (run_path / "config.yaml").exists():
        return candidates
    refined = []
    for row in candidates:
        try:
            refined.append(_refine_candidate_for_run(run_path, row))
        except Exception:
            refined.append(row)
    return refined
def read_csv_rows(path: str | Path) -> list[dict]:
    """Read CSV rows, returning an empty list for missing or empty files."""
    csv_path = Path(path)
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return []
    with csv_path.open("r", newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_top_candidates_csv(path: str | Path, rows: Iterable[dict], top_n: int = 30) -> list[dict]:
    """Write a deterministic top-candidates table and return written rows."""
    candidates = select_top_candidates(rows, top_n=top_n)
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CANDIDATE_FIELDNAMES)
        writer.writeheader()
        writer.writerows(candidates)
    return candidates


def write_top_candidates_for_run(run_dir: str | Path, top_n: int = 30) -> list[dict]:
    """Regenerate top_candidates.csv from a completed run's samples.csv."""
    run_path = Path(run_dir)
    samples = read_csv_rows(run_path / "samples.csv")
    candidates = select_top_candidates(samples, top_n=top_n)
    candidates = _refine_candidates_for_run(run_path, candidates)
    with (run_path / "top_candidates.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CANDIDATE_FIELDNAMES)
        writer.writeheader()
        writer.writerows(candidates)
    return candidates


def load_top_candidates(run_dir: str | Path) -> list[dict]:
    """Load top candidates for plotting/reporting."""
    return read_csv_rows(Path(run_dir) / "top_candidates.csv")

