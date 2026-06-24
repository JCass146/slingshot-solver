"""COM-frame scientific metrics and moving-potential work accounting."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np

from ..constants import G_KM
from .dynamics import EncounterIntegration, center_of_mass_state


LEGACY_METRIC_ALIASES = {
    "energy_half_dv_vec_sq": "turning_quadratic",
    "half_dv_vec_sq": "turning_quadratic",
    "bary_delta_v": "delta_speed_com",
    "delta_v": "delta_speed_com",
}


def wrap_angle(angle_rad: float) -> float:
    return float((angle_rad + np.pi) % (2.0 * np.pi) - np.pi)


def specific_energy_com(
    state: np.ndarray,
    star_mass_kg: float,
    planet_mass_kg: float,
    softening_km: float = 0.0,
) -> float:
    """Instantaneous test-particle specific energy in the binary COM frame."""
    com_position, com_velocity = center_of_mass_state(
        state, star_mass_kg, planet_mass_kg
    )
    test_position = state[8:10] - com_position
    test_velocity = state[10:12] - com_velocity
    star_position = state[0:2] - com_position
    planet_position = state[4:6] - com_position
    star_distance = np.sqrt(
        float(np.dot(test_position - star_position, test_position - star_position))
        + softening_km**2
    )
    planet_distance = np.sqrt(
        float(
            np.dot(
                test_position - planet_position,
                test_position - planet_position,
            )
        )
        + softening_km**2
    )
    return float(
        0.5 * np.dot(test_velocity, test_velocity)
        - G_KM * star_mass_kg / star_distance
        - G_KM * planet_mass_kg / planet_distance
    )


def analyze_integration(
    integration: EncounterIntegration,
    initial_state: np.ndarray,
    star_mass_kg: float,
    planet_mass_kg: float,
    semi_major_axis_km: float,
    softening_km: float = 0.0,
) -> dict[str, Any]:
    """Compute boost-invariant gain, turning, work, and outcome metrics."""
    solution = integration.solution
    final_state = solution.y[:, -1]
    if initial_state.shape == (12,):
        initial_state = np.concatenate([initial_state, np.zeros(2)])

    _, initial_com_velocity = center_of_mass_state(
        initial_state, star_mass_kg, planet_mass_kg
    )
    _, final_com_velocity = center_of_mass_state(
        final_state, star_mass_kg, planet_mass_kg
    )
    initial_velocity = initial_state[10:12] - initial_com_velocity
    final_velocity = final_state[10:12] - final_com_velocity
    initial_speed = float(np.linalg.norm(initial_velocity))
    final_speed = float(np.linalg.norm(final_velocity))

    energy_initial = specific_energy_com(
        initial_state, star_mass_kg, planet_mass_kg, softening_km
    )
    energy_final = specific_energy_com(
        final_state, star_mass_kg, planet_mass_kg, softening_km
    )
    delta_energy = energy_final - energy_initial
    v_inf_initial = np.sqrt(2.0 * energy_initial) if energy_initial > 0.0 else np.nan
    v_inf_final = np.sqrt(2.0 * energy_final) if energy_final > 0.0 else np.nan
    delta_v_inf = (
        float(v_inf_final - v_inf_initial)
        if np.isfinite(v_inf_initial) and np.isfinite(v_inf_final)
        else np.nan
    )
    velocity_change = final_velocity - initial_velocity
    turning_magnitude = float(np.linalg.norm(velocity_change))
    turning_quadratic = 0.5 * turning_magnitude**2
    deflection = wrap_angle(
        np.arctan2(final_velocity[1], final_velocity[0])
        - np.arctan2(initial_velocity[1], initial_velocity[0])
    )

    work_star = float(final_state[12] - initial_state[12])
    work_planet = float(final_state[13] - initial_state[13])
    work_sum = work_star + work_planet
    closure_error = delta_energy - work_sum
    closure_scale = max(abs(delta_energy), abs(work_sum), 1.0)
    work_norm = abs(work_star) + abs(work_planet)
    planet_work_fraction_abs = abs(work_planet) / work_norm if work_norm else 0.0
    star_work_fraction_abs = abs(work_star) / work_norm if work_norm else 0.0

    total_mu = G_KM * (star_mass_kg + planet_mass_kg)
    circular_speed_squared = total_mu / semi_major_axis_km
    energy_gain_dimensionless = delta_energy / circular_speed_squared
    return {
        "outcome": integration.outcome,
        "escaped": integration.outcome == "escaped",
        "collision": integration.outcome in {"star_collision", "planet_collision"},
        "initial_specific_energy_com": energy_initial,
        "final_specific_energy_com": energy_final,
        "delta_specific_energy_com": delta_energy,
        "energy_gain_dimensionless": float(energy_gain_dimensionless),
        "initial_v_inf_kms": float(v_inf_initial) if np.isfinite(v_inf_initial) else np.nan,
        "final_v_inf_kms": float(v_inf_final) if np.isfinite(v_inf_final) else np.nan,
        "delta_v_inf": delta_v_inf,
        "initial_speed_com": initial_speed,
        "final_speed_com": final_speed,
        "delta_speed_com": final_speed - initial_speed,
        "delta_kinetic_energy_com": 0.5 * (final_speed**2 - initial_speed**2),
        "turning_magnitude": turning_magnitude,
        "turning_quadratic": turning_quadratic,
        "deflection_rad": deflection,
        "deflection_deg": float(np.degrees(deflection)),
        "periapsis_planet_km": integration.periapsis_planet_km,
        "periapsis_star_km": integration.periapsis_star_km,
        "work_star": work_star,
        "work_planet": work_planet,
        "work_sum": work_sum,
        "work_energy_closure_error": closure_error,
        "work_energy_closure_relative": abs(closure_error) / closure_scale,
        "planet_work_fraction_abs": planet_work_fraction_abs,
        "star_work_fraction_abs": star_work_fraction_abs,
        "solver_success": bool(solution.success),
        "solver_status": int(solution.status),
        "solver_nfev": int(solution.nfev),
        "solver_message": str(solution.message),
        "integration_time_sec": float(solution.t[-1]),
    }


def resolve_metric(metrics: dict[str, Any], name: str) -> Any:
    """Resolve current metrics while warning on readable legacy aliases."""
    resolved = LEGACY_METRIC_ALIASES.get(name, name)
    if resolved != name:
        warnings.warn(
            f"Metric {name!r} is legacy; use {resolved!r}.",
            DeprecationWarning,
            stacklevel=2,
        )
    return metrics[resolved]
