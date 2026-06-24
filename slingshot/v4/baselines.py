"""Matched moving-body two-body controls for research comparisons."""

from __future__ import annotations

from typing import Literal

import numpy as np

from ..constants import G_KM
from ..core.twobody_scatter import gravity_assist_no_burn
from .dynamics import center_of_mass_state


def matched_moving_body_control(
    initial_state: np.ndarray,
    star_mass_kg: float,
    planet_mass_kg: float,
    body: Literal["star", "planet"],
    body_radius_km: float,
) -> dict:
    """Evaluate an isolated moving-body flyby from the same finite entry state.

    Both energy endpoints are asymptotic. The finite inbound state is used to
    infer the relative hyperbola, then its incoming and outgoing v-infinity
    vectors are transformed into the binary COM frame.
    """
    com_position, com_velocity = center_of_mass_state(
        initial_state, star_mass_kg, planet_mass_kg
    )
    if body == "star":
        body_position = initial_state[0:2] - com_position
        body_velocity = initial_state[2:4] - com_velocity
        body_mass = star_mass_kg
    elif body == "planet":
        body_position = initial_state[4:6] - com_position
        body_velocity = initial_state[6:8] - com_velocity
        body_mass = planet_mass_kg
    else:
        raise ValueError("body must be 'star' or 'planet'")

    test_position = initial_state[8:10] - com_position
    finite_test_velocity = initial_state[10:12] - com_velocity
    relative_position = test_position - body_position
    result = gravity_assist_no_burn(
        float(relative_position[0]),
        float(relative_position[1]),
        float(finite_test_velocity[0]),
        float(finite_test_velocity[1]),
        (float(body_velocity[0]), float(body_velocity[1])),
        G_KM * body_mass,
    )
    collision = result.rp <= body_radius_km
    finite_relative_velocity = finite_test_velocity - body_velocity
    finite_relative_speed = float(np.linalg.norm(finite_relative_velocity))
    incoming_direction = finite_relative_velocity / finite_relative_speed
    initial_velocity = body_velocity + result.vinf * incoming_direction
    final_velocity = np.array([result.umF, result.vmF])
    velocity_change = final_velocity - initial_velocity
    initial_speed = float(np.linalg.norm(initial_velocity))
    final_speed = float(np.linalg.norm(final_velocity))
    deflection = float(
        (
            np.arctan2(final_velocity[1], final_velocity[0])
            - np.arctan2(initial_velocity[1], initial_velocity[0])
            + np.pi
        )
        % (2.0 * np.pi)
        - np.pi
    )
    return {
        "body": body,
        "valid": bool(result.e > 1.0 and result.epsilon > 0.0 and not collision),
        "collision": bool(collision),
        "relative_v_inf_kms": float(result.vinf),
        "relative_periapsis_km": float(result.rp),
        "delta_specific_energy_com": float(
            0.5 * (final_speed**2 - initial_speed**2)
        ),
        "delta_speed_com": final_speed - initial_speed,
        "turning_magnitude": float(np.linalg.norm(velocity_change)),
        "turning_quadratic": float(0.5 * np.dot(velocity_change, velocity_change)),
        "deflection_rad": deflection,
        "deflection_deg": float(np.degrees(deflection)),
    }
