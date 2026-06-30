"""Orbital-element and event-integration tests for the research core."""

from types import SimpleNamespace

import numpy as np
import pytest


def test_kepler_equation_residual():
    from slingshot.v4.dynamics import solve_kepler

    mean_anomaly = 2.1
    eccentricity = 0.62
    eccentric_anomaly = solve_kepler(mean_anomaly, eccentricity)
    residual = eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly) - mean_anomaly
    assert residual == pytest.approx(0.0, abs=1e-13)


@pytest.mark.parametrize("mean_anomaly", [0.0, 0.7, 2.4, 5.8])
def test_keplerian_initializer_recovers_elements(mean_anomaly):
    from slingshot.constants import AU_KM, M_JUP, M_SUN
    from slingshot.v4.dynamics import (
        binary_elements_from_state,
        init_binary_barycentric,
    )

    semi_major_axis = 0.303 * AU_KM
    eccentricity = 0.478
    state = init_binary_barycentric(
        semi_major_axis,
        eccentricity,
        mean_anomaly,
        0.31,
        1.35 * M_SUN,
        5.84 * M_JUP,
        bulk_velocity_x_kms=17.0,
        bulk_velocity_y_kms=-4.0,
    )
    recovered = binary_elements_from_state(
        state, 1.35 * M_SUN, 5.84 * M_JUP
    )
    assert recovered["semi_major_axis_km"] == pytest.approx(
        semi_major_axis, rel=1e-12
    )
    assert recovered["eccentricity"] == pytest.approx(eccentricity, abs=1e-12)


def test_barycenter_and_bulk_velocity():
    from slingshot.constants import AU_KM, M_JUP, M_SUN
    from slingshot.v4.dynamics import center_of_mass_state, init_binary_barycentric

    star_mass = 1.2 * M_SUN
    planet_mass = 4.0 * M_JUP
    state = init_binary_barycentric(
        0.2 * AU_KM,
        0.4,
        1.0,
        0.0,
        star_mass,
        planet_mass,
        bulk_velocity_x_kms=12.0,
        bulk_velocity_y_kms=-3.0,
    )
    position, velocity = center_of_mass_state(state, star_mass, planet_mass)
    assert position == pytest.approx([0.0, 0.0], abs=1e-6)
    assert velocity == pytest.approx([12.0, -3.0], abs=1e-12)

def test_minimum_separation_uses_closest_event_not_first():
    from slingshot.v4.dynamics import _minimum_separation_km

    y = np.zeros((14, 3), dtype=float)
    y[8, :] = [50.0, 40.0, 30.0]

    def dense_state(time):
        state = np.zeros(14, dtype=float)
        state[8] = 20.0 if time == 1.0 else 3.0
        return state

    solution = SimpleNamespace(
        y=y,
        t=np.array([0.0, 1.0, 2.0]),
        sol=dense_state,
    )
    distance, time = _minimum_separation_km(
        solution,
        slice(8, 10),
        slice(4, 6),
        np.array([1.0, 2.0]),
    )
    assert distance == pytest.approx(3.0)
    assert time == pytest.approx(2.0)