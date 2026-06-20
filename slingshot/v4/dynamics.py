"""Keplerian binary initialization and event-driven planar integration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp

from ..constants import G_KM


def solve_kepler(mean_anomaly: float, eccentricity: float, tolerance: float = 1e-14) -> float:
    """Solve M = E - e sin(E) for an elliptic orbit."""
    mean_anomaly = float(np.mod(mean_anomaly, 2.0 * np.pi))
    if eccentricity == 0.0:
        return mean_anomaly
    eccentric_anomaly = mean_anomaly if eccentricity < 0.8 else np.pi
    for _ in range(64):
        residual = eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly) - mean_anomaly
        derivative = 1.0 - eccentricity * np.cos(eccentric_anomaly)
        correction = residual / derivative
        eccentric_anomaly -= correction
        if abs(correction) <= tolerance:
            return float(eccentric_anomaly)
    raise RuntimeError("Kepler equation did not converge")


def init_binary_barycentric(
    semi_major_axis_km: float,
    eccentricity: float,
    mean_anomaly_rad: float,
    argument_periapsis_rad: float,
    star_mass_kg: float,
    planet_mass_kg: float,
    prograde: bool = True,
    bulk_velocity_x_kms: float = 0.0,
    bulk_velocity_y_kms: float = 0.0,
) -> np.ndarray:
    """Return star/planet barycentric state for an elliptic relative orbit."""
    total_mass = star_mass_kg + planet_mass_kg
    mu = G_KM * total_mass
    eccentric_anomaly = solve_kepler(mean_anomaly_rad, eccentricity)
    root = np.sqrt(1.0 - eccentricity * eccentricity)
    denominator = 1.0 - eccentricity * np.cos(eccentric_anomaly)
    mean_motion = np.sqrt(mu / semi_major_axis_km**3)

    x_relative = semi_major_axis_km * (np.cos(eccentric_anomaly) - eccentricity)
    y_relative = semi_major_axis_km * root * np.sin(eccentric_anomaly)
    vx_relative = (
        -semi_major_axis_km * mean_motion * np.sin(eccentric_anomaly) / denominator
    )
    vy_relative = (
        semi_major_axis_km
        * mean_motion
        * root
        * np.cos(eccentric_anomaly)
        / denominator
    )
    if not prograde:
        y_relative = -y_relative
        vy_relative = -vy_relative

    cosine = np.cos(argument_periapsis_rad)
    sine = np.sin(argument_periapsis_rad)
    relative_position = np.array(
        [
            cosine * x_relative - sine * y_relative,
            sine * x_relative + cosine * y_relative,
        ]
    )
    relative_velocity = np.array(
        [
            cosine * vx_relative - sine * vy_relative,
            sine * vx_relative + cosine * vy_relative,
        ]
    )

    star_fraction = planet_mass_kg / total_mass
    planet_fraction = star_mass_kg / total_mass
    bulk_velocity = np.array([bulk_velocity_x_kms, bulk_velocity_y_kms])
    star_position = -star_fraction * relative_position
    planet_position = planet_fraction * relative_position
    star_velocity = -star_fraction * relative_velocity + bulk_velocity
    planet_velocity = planet_fraction * relative_velocity + bulk_velocity
    return np.array(
        [
            star_position[0],
            star_position[1],
            star_velocity[0],
            star_velocity[1],
            planet_position[0],
            planet_position[1],
            planet_velocity[0],
            planet_velocity[1],
        ],
        dtype=float,
    )


def binary_elements_from_state(
    state: np.ndarray, star_mass_kg: float, planet_mass_kg: float
) -> dict:
    """Recover relative semi-major axis and eccentricity from a binary state."""
    relative_position = state[4:6] - state[0:2]
    relative_velocity = state[6:8] - state[2:4]
    radius = np.linalg.norm(relative_position)
    speed_squared = float(np.dot(relative_velocity, relative_velocity))
    mu = G_KM * (star_mass_kg + planet_mass_kg)
    specific_energy = 0.5 * speed_squared - mu / radius
    semi_major_axis = -mu / (2.0 * specific_energy)
    angular_momentum = (
        relative_position[0] * relative_velocity[1]
        - relative_position[1] * relative_velocity[0]
    )
    eccentricity = np.sqrt(
        max(0.0, 1.0 + 2.0 * specific_energy * angular_momentum**2 / mu**2)
    )
    return {
        "semi_major_axis_km": float(semi_major_axis),
        "eccentricity": float(eccentricity),
        "specific_energy": float(specific_energy),
        "angular_momentum": float(angular_momentum),
    }


def center_of_mass_state(
    state: np.ndarray, star_mass_kg: float, planet_mass_kg: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return binary center-of-mass position and velocity."""
    total_mass = star_mass_kg + planet_mass_kg
    position = (
        star_mass_kg * state[0:2] + planet_mass_kg * state[4:6]
    ) / total_mass
    velocity = (
        star_mass_kg * state[2:4] + planet_mass_kg * state[6:8]
    ) / total_mass
    return position, velocity


def restricted_ode(
    time: float,
    state: np.ndarray,
    star_mass_kg: float,
    planet_mass_kg: float,
    softening_squared_km2: float,
) -> np.ndarray:
    """Restricted binary plus massless test-particle ODE with work integrals."""
    del time
    xs, ys, vxs, vys, xp, yp, vxp, vyp, xt, yt, vxt, vyt = state[:12]
    dx_binary = xp - xs
    dy_binary = yp - ys
    binary_r2 = dx_binary**2 + dy_binary**2
    binary_r3 = binary_r2 * np.sqrt(binary_r2)

    q_star = np.array([xt - xs, yt - ys])
    q_planet = np.array([xt - xp, yt - yp])
    star_r3 = (float(np.dot(q_star, q_star)) + softening_squared_km2) ** 1.5
    planet_r3 = (float(np.dot(q_planet, q_planet)) + softening_squared_km2) ** 1.5

    star_ax = G_KM * planet_mass_kg * dx_binary / binary_r3
    star_ay = G_KM * planet_mass_kg * dy_binary / binary_r3
    planet_ax = -G_KM * star_mass_kg * dx_binary / binary_r3
    planet_ay = -G_KM * star_mass_kg * dy_binary / binary_r3
    test_acceleration = (
        -G_KM * star_mass_kg * q_star / star_r3
        - G_KM * planet_mass_kg * q_planet / planet_r3
    )

    com_velocity = (
        star_mass_kg * np.array([vxs, vys])
        + planet_mass_kg * np.array([vxp, vyp])
    ) / (star_mass_kg + planet_mass_kg)
    star_velocity_com = np.array([vxs, vys]) - com_velocity
    planet_velocity_com = np.array([vxp, vyp]) - com_velocity
    work_star_rate = -G_KM * star_mass_kg * float(
        np.dot(q_star, star_velocity_com)
    ) / star_r3
    work_planet_rate = -G_KM * planet_mass_kg * float(
        np.dot(q_planet, planet_velocity_com)
    ) / planet_r3

    derivative = np.array(
        [
            vxs,
            vys,
            star_ax,
            star_ay,
            vxp,
            vyp,
            planet_ax,
            planet_ay,
            vxt,
            vyt,
            test_acceleration[0],
            test_acceleration[1],
            work_star_rate,
            work_planet_rate,
        ]
    )
    return derivative


@dataclass
class EncounterIntegration:
    solution: object
    outcome: str
    periapsis_planet_km: float
    periapsis_star_km: float
    periapsis_time_sec: Optional[float]
    star_periapsis_time_sec: Optional[float]


def integrate_encounter(
    initial_state: np.ndarray,
    star_mass_kg: float,
    planet_mass_kg: float,
    star_radius_km: float,
    planet_radius_km: float,
    boundary_radius_km: float,
    max_time_sec: float,
    method: str = "DOP853",
    rtol: float = 1e-10,
    atol: float = 1e-10,
    softening_km: float = 0.0,
    max_step_sec: Optional[float] = None,
) -> EncounterIntegration:
    """Integrate until outbound boundary crossing or a physical collision."""

    if initial_state.shape == (12,):
        initial_state = np.concatenate([initial_state, np.zeros(2)])
    if initial_state.shape != (14,):
        raise ValueError("initial_state must contain 12 phase variables and 2 work variables")

    def outbound_event(time, state, *args):
        del args
        if time <= 0.0:
            return -1.0
        com_position, _ = center_of_mass_state(
            state, star_mass_kg, planet_mass_kg
        )
        return np.linalg.norm(state[8:10] - com_position) - boundary_radius_km

    outbound_event.terminal = True
    outbound_event.direction = 1.0

    def star_collision_event(time, state, *args):
        del time, args
        return np.linalg.norm(state[8:10] - state[0:2]) - star_radius_km

    star_collision_event.terminal = True
    star_collision_event.direction = -1.0

    def planet_collision_event(time, state, *args):
        del time, args
        return np.linalg.norm(state[8:10] - state[4:6]) - planet_radius_km

    planet_collision_event.terminal = True
    planet_collision_event.direction = -1.0

    def planet_periapsis_event(time, state, *args):
        del time, args
        relative_position = state[8:10] - state[4:6]
        relative_velocity = state[10:12] - state[6:8]
        return float(np.dot(relative_position, relative_velocity))

    planet_periapsis_event.terminal = False
    planet_periapsis_event.direction = 1.0

    def star_periapsis_event(time, state, *args):
        del time, args
        relative_position = state[8:10] - state[0:2]
        relative_velocity = state[10:12] - state[2:4]
        return float(np.dot(relative_position, relative_velocity))

    star_periapsis_event.terminal = False
    star_periapsis_event.direction = 1.0

    solve_kwargs = {}
    if max_step_sec is not None:
        solve_kwargs["max_step"] = max_step_sec
    solution = solve_ivp(
        restricted_ode,
        (0.0, max_time_sec),
        initial_state,
        args=(star_mass_kg, planet_mass_kg, softening_km**2),
        method=method,
        rtol=rtol,
        atol=atol,
        dense_output=True,
        events=[
            outbound_event,
            star_collision_event,
            planet_collision_event,
            planet_periapsis_event,
            star_periapsis_event,
        ],
        **solve_kwargs,
    )
    if not solution.success:
        outcome = "integration_failed"
    elif len(solution.t_events[1]):
        outcome = "star_collision"
    elif len(solution.t_events[2]):
        outcome = "planet_collision"
    elif len(solution.t_events[0]):
        outcome = "escaped"
    else:
        outcome = "time_limit"

    star_distance = np.linalg.norm(solution.y[8:10] - solution.y[0:2], axis=0)
    planet_distance = np.linalg.norm(solution.y[8:10] - solution.y[4:6], axis=0)
    planet_periapsis_time = (
        float(solution.t_events[3][0]) if len(solution.t_events[3]) else None
    )
    star_periapsis_time = (
        float(solution.t_events[4][0]) if len(solution.t_events[4]) else None
    )
    planet_periapsis = float(np.min(planet_distance))
    star_periapsis = float(np.min(star_distance))
    if planet_periapsis_time is not None and solution.sol is not None:
        event_state = solution.sol(planet_periapsis_time)
        planet_periapsis = float(
            np.linalg.norm(event_state[8:10] - event_state[4:6])
        )
    if star_periapsis_time is not None and solution.sol is not None:
        event_state = solution.sol(star_periapsis_time)
        star_periapsis = float(
            np.linalg.norm(event_state[8:10] - event_state[0:2])
        )
    return EncounterIntegration(
        solution=solution,
        outcome=outcome,
        periapsis_planet_km=planet_periapsis,
        periapsis_star_km=star_periapsis,
        periapsis_time_sec=planet_periapsis_time,
        star_periapsis_time_sec=star_periapsis_time,
    )
