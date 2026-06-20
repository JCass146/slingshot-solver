"""Matched two-body control tests."""

import numpy as np
import pytest


def test_stationary_body_conserves_asymptotic_speed():
    from slingshot.constants import M_JUP, M_SUN, R_JUP
    from slingshot.v4.baselines import matched_moving_body_control

    initial = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            1e9,
            0.0,
            0.0,
            0.0,
            -1e10,
            1e7,
            50.0,
            0.0,
        ]
    )
    result = matched_moving_body_control(
        initial,
        M_SUN,
        M_JUP,
        "planet",
        R_JUP,
    )
    assert result["valid"]
    assert result["delta_specific_energy_com"] == pytest.approx(0.0, abs=1e-9)
    assert result["turning_quadratic"] > 0.0


def test_unknown_body_is_rejected():
    from slingshot.constants import M_JUP, M_SUN, R_JUP
    from slingshot.v4.baselines import matched_moving_body_control

    with pytest.raises(ValueError):
        matched_moving_body_control(
            np.zeros(12), M_SUN, M_JUP, "moon", R_JUP
        )
