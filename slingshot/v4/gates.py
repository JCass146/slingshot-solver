"""Publication validation gates using a verified non-colliding reference ray."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from ..constants import AU_KM, G_KM
from .config import V4Config
from .dynamics import init_binary_barycentric, integrate_encounter, restricted_ode
from .metrics import analyze_integration
from .sampling import state_at_inbound_boundary
from .validation import (
    compare_rebound_final_state,
    jacobi_constant,
    physical_values,
    rebound_available,
    run_quick_validation,
)


def _reference(
    config: V4Config,
    method: str,
    bulk_velocity: tuple[float, float] = (0.0, 0.0),
):
    values = physical_values(config)
    boundary = max(
        2.5 * values["semi_major_axis_km"],
        min(values["boundary_radius_km"], AU_KM),
    )
    total_mu = G_KM * (values["star_mass_kg"] + values["planet_mass_kg"])
    position, velocity = state_at_inbound_boundary(
        max(60.0, config.asymptotic_sampling.v_inf_kms[0]),
        min(0.15 * AU_KM, 0.2 * values["b_max_km"]),
        0.30,
        boundary,
        total_mu,
    )
    binary = init_binary_barycentric(
        values["semi_major_axis_km"],
        config.orbit.eccentricity,
        1.2,
        config.orbit.argument_periapsis_rad,
        values["star_mass_kg"],
        values["planet_mass_kg"],
        config.orbit.prograde,
        bulk_velocity[0],
        bulk_velocity[1],
    )
    initial_state = np.concatenate(
        [binary, position, velocity + np.asarray(bulk_velocity), np.zeros(2)]
    )
    integration = integrate_encounter(
        initial_state,
        values["star_mass_kg"],
        values["planet_mass_kg"],
        values["star_radius_km"],
        values["planet_radius_km"],
        boundary,
        min(config.numerical.max_time_sec, 3.0e7),
        method=method,
        rtol=min(config.numerical.rtol, 1e-10),
        atol=min(config.numerical.atol, 1e-10),
        softening_km=0.0,
    )
    metrics = analyze_integration(
        integration,
        initial_state,
        values["star_mass_kg"],
        values["planet_mass_kg"],
        values["semi_major_axis_km"],
        0.0,
    )
    return initial_state, integration, metrics, values


def validate_work_energy(config: V4Config) -> dict[str, Any]:
    _, integration, metrics, _ = _reference(config, "DOP853")
    tolerance = config.validation.work_energy_relative_tolerance
    return {
        "name": "work_energy_closure",
        "passed": integration.outcome == "escaped"
        and metrics["work_energy_closure_relative"] <= tolerance,
        "outcome": integration.outcome,
        "relative_error": metrics["work_energy_closure_relative"],
        "tolerance": tolerance,
    }


def validate_galilean_invariance(config: V4Config) -> dict[str, Any]:
    _, base_integration, base, _ = _reference(config, "DOP853")
    _, boosted_integration, boosted, _ = _reference(
        config, "DOP853", (23.0, -17.0)
    )
    keys = (
        "delta_specific_energy_com",
        "delta_v_inf",
        "delta_speed_com",
        "turning_quadratic",
        "deflection_rad",
    )
    comparisons = {}
    for key in keys:
        first = float(base[key])
        second = float(boosted[key])
        comparisons[key] = abs(first - second) / max(
            abs(first), abs(second), 1.0
        )
    max_relative = max(comparisons.values())
    tolerance = 1e-6
    return {
        "name": "galilean_invariance",
        "passed": base_integration.outcome == boosted_integration.outcome == "escaped"
        and max_relative <= tolerance,
        "max_relative_difference": max_relative,
        "tolerance": tolerance,
        "comparisons": comparisons,
    }


def validate_integrator_agreement(config: V4Config) -> dict[str, Any]:
    _, dop_integration, dop, _ = _reference(config, "DOP853")
    _, radau_integration, radau, _ = _reference(config, "Radau")
    comparisons = {}
    for key in (
        "delta_specific_energy_com",
        "periapsis_planet_km",
        "deflection_rad",
    ):
        first = float(dop[key])
        second = float(radau[key])
        comparisons[key] = abs(first - second) / max(abs(first), abs(second), 1.0)
    max_relative = max(comparisons.values())
    tolerance = config.validation.integrator_agreement_relative_tolerance
    return {
        "name": "dop853_radau_agreement",
        "passed": dop_integration.outcome == radau_integration.outcome == "escaped"
        and max_relative <= tolerance,
        "max_relative_difference": max_relative,
        "tolerance": tolerance,
        "comparisons": comparisons,
    }


def validate_circular_jacobi(config: V4Config) -> dict[str, Any]:
    values = physical_values(config)
    binary = init_binary_barycentric(
        values["semi_major_axis_km"],
        0.0,
        0.4,
        0.0,
        values["star_mass_kg"],
        values["planet_mass_kg"],
    )
    position = np.array(
        [
            1.7 * values["semi_major_axis_km"],
            0.2 * values["semi_major_axis_km"],
        ]
    )
    circular_speed = np.sqrt(
        G_KM
        * (values["star_mass_kg"] + values["planet_mass_kg"])
        / np.linalg.norm(position)
    )
    initial_state = np.concatenate(
        [binary, position, np.array([-2.0, 0.7 * circular_speed]), np.zeros(2)]
    )
    solution = solve_ivp(
        restricted_ode,
        (0.0, 2.0e6),
        initial_state,
        args=(values["star_mass_kg"], values["planet_mass_kg"], 0.0),
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
    )
    constants = np.array(
        [
            jacobi_constant(
                solution.y[:, index],
                solution.t[index],
                values["star_mass_kg"],
                values["planet_mass_kg"],
                values["semi_major_axis_km"],
            )
            for index in range(solution.t.size)
        ]
    )
    relative_drift = float(
        (np.max(constants) - np.min(constants))
        / max(abs(float(constants[0])), 1.0)
    )
    tolerance = config.validation.jacobi_relative_tolerance
    return {
        "name": "circular_jacobi_conservation",
        "passed": solution.success and relative_drift <= tolerance,
        "relative_drift": relative_drift,
        "tolerance": tolerance,
    }


def validate_rebound(config: V4Config) -> dict[str, Any]:
    if not rebound_available():
        return {
            "name": "rebound_ias15_agreement",
            "passed": False,
            "required": False,
            "skipped": True,
            "reason": "REBOUND not installed",
        }
    initial_state, integration, _, values = _reference(config, "DOP853")
    rebound_state = compare_rebound_final_state(
        initial_state,
        values["star_mass_kg"],
        values["planet_mass_kg"],
        float(integration.solution.t[-1]),
    )
    scipy_state = integration.solution.y[:12, -1]
    relative = float(
        np.linalg.norm(rebound_state - scipy_state)
        / max(np.linalg.norm(scipy_state), 1.0)
    )
    tolerance = config.validation.integrator_agreement_relative_tolerance
    return {
        "name": "rebound_ias15_agreement",
        "passed": relative <= tolerance,
        "required": False,
        "skipped": False,
        "relative_state_difference": relative,
        "tolerance": tolerance,
    }


def run_publication_validation(config: V4Config) -> dict[str, Any]:
    quick = run_quick_validation(config)
    gates = list(quick["gates"])
    gates.extend(
        [
            validate_work_energy(config),
            validate_galilean_invariance(config),
            validate_integrator_agreement(config),
            validate_circular_jacobi(config),
            validate_rebound(config),
        ]
    )
    required = [
        gate
        for gate in gates
        if gate.get("required", True) and not gate.get("skipped", False)
    ]
    return {"passed": all(gate["passed"] for gate in required), "gates": gates}
