"""Numerical and physical validation gates for the v4 research core."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from ..constants import AU_KM, G_KM, M_JUP, M_SUN, R_JUP, R_SUN
from .config import V4Config
from .dynamics import (
    binary_elements_from_state,
    init_binary_barycentric,
    integrate_encounter,
)
from .metrics import analyze_integration
from .sampling import state_at_inbound_boundary


def physical_values(config: V4Config) -> dict[str, float]:
    return {
        "star_mass_kg": config.system.star_mass_msun * M_SUN,
        "planet_mass_kg": config.system.planet_mass_mjup * M_JUP,
        "star_radius_km": config.system.star_radius_rsun * R_SUN,
        "planet_radius_km": config.system.planet_radius_rjup * R_JUP,
        "semi_major_axis_km": config.orbit.semi_major_axis_au * AU_KM,
        "boundary_radius_km": config.asymptotic_sampling.boundary_radius_au
        * AU_KM,
        "b_max_km": config.asymptotic_sampling.b_max_au * AU_KM,
    }


def validate_binary_elements(config: V4Config) -> dict[str, Any]:
    values = physical_values(config)
    state = init_binary_barycentric(
        semi_major_axis_km=values["semi_major_axis_km"],
        eccentricity=config.orbit.eccentricity,
        mean_anomaly_rad=config.orbit.mean_anomaly_rad,
        argument_periapsis_rad=config.orbit.argument_periapsis_rad,
        star_mass_kg=values["star_mass_kg"],
        planet_mass_kg=values["planet_mass_kg"],
        prograde=config.orbit.prograde,
        bulk_velocity_x_kms=config.system.bulk_velocity_x_kms,
        bulk_velocity_y_kms=config.system.bulk_velocity_y_kms,
    )
    recovered = binary_elements_from_state(
        state, values["star_mass_kg"], values["planet_mass_kg"]
    )
    semi_major_error = abs(
        recovered["semi_major_axis_km"] - values["semi_major_axis_km"]
    ) / values["semi_major_axis_km"]
    eccentricity_scale = max(config.orbit.eccentricity, 1.0)
    eccentricity_error = abs(
        recovered["eccentricity"] - config.orbit.eccentricity
    ) / eccentricity_scale
    tolerance = config.validation.binary_elements_relative_tolerance
    return {
        "name": "binary_elements",
        "passed": semi_major_error <= tolerance and eccentricity_error <= tolerance,
        "semi_major_axis_relative_error": semi_major_error,
        "eccentricity_relative_error": eccentricity_error,
        "tolerance": tolerance,
    }


def numerical_two_body_deflection(
    v_inf_kms: float,
    impact_parameter_km: float,
    boundary_radius_km: float,
    mu_km3_s2: float,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> dict[str, float]:
    """Compare a numerical central-force encounter with analytic scattering."""
    initial_position, initial_velocity = state_at_inbound_boundary(
        v_inf_kms,
        impact_parameter_km,
        0.0,
        boundary_radius_km,
        mu_km3_s2,
    )
    initial_state = np.concatenate([initial_position, initial_velocity])

    def ode(time, state):
        del time
        radius = np.linalg.norm(state[:2])
        acceleration = -mu_km3_s2 * state[:2] / radius**3
        return np.array([state[2], state[3], acceleration[0], acceleration[1]])

    def outbound(time, state):
        if time <= 0.0:
            return -1.0
        return np.linalg.norm(state[:2]) - boundary_radius_km

    outbound.terminal = True
    outbound.direction = 1.0
    travel_time = 4.0 * boundary_radius_km / v_inf_kms
    solution = solve_ivp(
        ode,
        (0.0, travel_time),
        initial_state,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=outbound,
    )
    if not len(solution.t_events[0]):
        raise RuntimeError("Two-body validation trajectory did not exit the boundary")

    final_position = solution.y[:2, -1]
    final_velocity = solution.y[2:4, -1]
    initial_energy = (
        0.5 * np.dot(initial_velocity, initial_velocity)
        - mu_km3_s2 / np.linalg.norm(initial_position)
    )
    final_energy = (
        0.5 * np.dot(final_velocity, final_velocity)
        - mu_km3_s2 / np.linalg.norm(final_position)
    )
    initial_angular_momentum = (
        initial_position[0] * initial_velocity[1]
        - initial_position[1] * initial_velocity[0]
    )
    final_angular_momentum = (
        final_position[0] * final_velocity[1]
        - final_position[1] * final_velocity[0]
    )
    eccentricity = np.sqrt(
        1.0
        + impact_parameter_km**2 * v_inf_kms**4 / mu_km3_s2**2
    )
    analytic_deflection = 2.0 * np.arcsin(1.0 / eccentricity)

    # Compute asymptotic deflection from the outgoing orbit elements.
    # For a pure-Keplerian hyperbola the eccentricity magnitude is constant;
    # we recover it from the final state and use 2*arcsin(1/e_num) so that
    # the comparison is asymptote-to-asymptote rather than
    # finite-boundary-velocity-to-finite-boundary-velocity.
    r_final = np.linalg.norm(final_position)
    v_sq_final = float(np.dot(final_velocity, final_velocity))
    rdotv_final = float(np.dot(final_position, final_velocity))
    e_vec_final = (
        (1.0 / mu_km3_s2)
        * (
            (v_sq_final - mu_km3_s2 / r_final) * final_position
            - rdotv_final * final_velocity
        )
    )
    e_numerical = float(np.linalg.norm(e_vec_final))
    # Guard against tiny rounding that makes arcsin argument exceed 1
    numerical_asymptotic_deflection = 2.0 * np.arcsin(
        np.clip(1.0 / max(e_numerical, 1.0 + 1e-15), 0.0, 1.0)
    )

    deflection_relative_error = float(
        abs(numerical_asymptotic_deflection - analytic_deflection)
        / max(abs(analytic_deflection), 1e-15)
    )
    return {
        "analytic_deflection_rad": float(analytic_deflection),
        "numerical_asymptotic_deflection_rad": float(numerical_asymptotic_deflection),
        "deflection_relative_error": deflection_relative_error,
        "energy_relative_error": float(
            abs(final_energy - initial_energy) / max(abs(initial_energy), 1.0)
        ),
        "angular_momentum_relative_error": float(
            abs(final_angular_momentum - initial_angular_momentum)
            / max(abs(initial_angular_momentum), 1.0)
        ),
    }


def jacobi_constant(
    state: np.ndarray,
    time_sec: float,
    star_mass_kg: float,
    planet_mass_kg: float,
    semi_major_axis_km: float,
) -> float:
    """Dimensional Jacobi constant for a circular restricted binary."""
    del time_sec
    com_position = (
        star_mass_kg * state[0:2] + planet_mass_kg * state[4:6]
    ) / (star_mass_kg + planet_mass_kg)
    com_velocity = (
        star_mass_kg * state[2:4] + planet_mass_kg * state[6:8]
    ) / (star_mass_kg + planet_mass_kg)
    position = state[8:10] - com_position
    velocity = state[10:12] - com_velocity
    mean_motion = np.sqrt(
        G_KM * (star_mass_kg + planet_mass_kg) / semi_major_axis_km**3
    )
    rotating_velocity = velocity - mean_motion * np.array(
        [-position[1], position[0]]
    )
    star_distance = np.linalg.norm(state[8:10] - state[0:2])
    planet_distance = np.linalg.norm(state[8:10] - state[4:6])
    effective_potential = (
        G_KM * star_mass_kg / star_distance
        + G_KM * planet_mass_kg / planet_distance
        + 0.5 * mean_motion**2 * np.dot(position, position)
    )
    return float(2.0 * effective_potential - np.dot(rotating_velocity, rotating_velocity))


def rebound_available() -> bool:
    try:
        import rebound  # noqa: F401
    except ImportError:
        return False
    return True


def compare_rebound_final_state(
    initial_state: np.ndarray,
    star_mass_kg: float,
    planet_mass_kg: float,
    final_time_sec: float,
) -> np.ndarray:
    """Integrate the same three-body state with optional REBOUND IAS15."""
    try:
        import rebound
    except ImportError as error:
        raise RuntimeError("REBOUND is not installed") from error
    simulation = rebound.Simulation()
    simulation.G = G_KM
    simulation.integrator = "ias15"
    simulation.add(
        m=star_mass_kg,
        x=initial_state[0],
        y=initial_state[1],
        vx=initial_state[2],
        vy=initial_state[3],
    )
    simulation.add(
        m=planet_mass_kg,
        x=initial_state[4],
        y=initial_state[5],
        vx=initial_state[6],
        vy=initial_state[7],
    )
    simulation.add(
        m=0.0,
        x=initial_state[8],
        y=initial_state[9],
        vx=initial_state[10],
        vy=initial_state[11],
    )
    simulation.integrate(final_time_sec, exact_finish_time=1)
    particles = simulation.particles
    return np.array(
        [
            particles[0].x,
            particles[0].y,
            particles[0].vx,
            particles[0].vy,
            particles[1].x,
            particles[1].y,
            particles[1].vx,
            particles[1].vy,
            particles[2].x,
            particles[2].y,
            particles[2].vx,
            particles[2].vy,
        ]
    )


def run_quick_validation(config: V4Config) -> dict[str, Any]:
    """Run deterministic gates that do not require a full campaign."""
    values = physical_values(config)
    binary_gate = validate_binary_elements(config)
    deflection = numerical_two_body_deflection(
        v_inf_kms=config.asymptotic_sampling.v_inf_kms[0],
        impact_parameter_km=min(
            0.25 * values["b_max_km"], 0.25 * values["boundary_radius_km"]
        ),
        boundary_radius_km=values["boundary_radius_km"],
        mu_km3_s2=G_KM * (values["star_mass_kg"] + values["planet_mass_kg"]),
    )
    tol = config.validation.two_body_deflection_relative_tolerance
    deflection_passed = (
        deflection["energy_relative_error"] <= tol
        and deflection["angular_momentum_relative_error"] <= tol
        and deflection["deflection_relative_error"] <= tol
    )
    gates = [
        binary_gate,
        {
            "name": "two_body_invariants",
            "passed": deflection_passed,
            "tolerance": tol,
            **deflection,
        },
        {
            "name": "newtonian_default",
            "passed": config.numerical.softening_km == 0.0,
            "softening_km": config.numerical.softening_km,
        },
        {
            "name": "rebound_available",
            "passed": rebound_available(),
            "required": False,
        },
    ]
    return {
        "passed": all(gate["passed"] for gate in gates if gate["name"] != "rebound_available"),
        "gates": gates,
    }
