"""Asymptotic-state reconstruction tests."""

import numpy as np
import pytest


@pytest.mark.parametrize("impact_parameter", [-2e7, 0.0, 2e7])
def test_boundary_state_preserves_asymptotic_invariants(impact_parameter):
    from slingshot.v4.sampling import state_at_inbound_boundary

    v_inf = 35.0
    mu = 1.7e11
    boundary = 1.2e9
    position, velocity = state_at_inbound_boundary(
        v_inf,
        impact_parameter,
        0.73,
        boundary,
        mu,
    )
    energy = 0.5 * np.dot(velocity, velocity) - mu / np.linalg.norm(position)
    angular_momentum = position[0] * velocity[1] - position[1] * velocity[0]
    assert np.linalg.norm(position) == pytest.approx(boundary, rel=1e-12)
    assert energy == pytest.approx(0.5 * v_inf**2, rel=1e-12)
    assert angular_momentum == pytest.approx(
        impact_parameter * v_inf, abs=1e-4
    )
    assert np.dot(position, velocity) < 0.0


def test_draw_samples_are_reproducible():
    from slingshot.v4.sampling import draw_samples

    kwargs = dict(
        count=4,
        v_inf_kms=40.0,
        b_max_km=1e8,
        boundary_radius_km=8e8,
        total_mu_km3_s2=1.6e11,
    )
    first = draw_samples(np.random.default_rng(42), **kwargs)
    second = draw_samples(np.random.default_rng(42), **kwargs)
    assert [sample.impact_parameter_km for sample in first] == pytest.approx(
        [sample.impact_parameter_km for sample in second]
    )
    assert [sample.binary_mean_anomaly_rad for sample in first] == pytest.approx(
        [sample.binary_mean_anomaly_rad for sample in second]
    )
