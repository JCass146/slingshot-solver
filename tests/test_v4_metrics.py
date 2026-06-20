"""Scientific metric semantics and boost-invariance tests."""

from types import SimpleNamespace

import numpy as np
import pytest


def _integration(initial_state, final_state):
    solution = SimpleNamespace(
        y=np.column_stack([initial_state, final_state]),
        success=True,
        status=1,
        nfev=20,
        message="ok",
        t=np.array([0.0, 1.0]),
    )
    return SimpleNamespace(
        solution=solution,
        outcome="escaped",
        periapsis_planet_km=1e6,
        periapsis_star_km=2e7,
    )


def test_turning_quadratic_is_not_energy_gain():
    from slingshot.constants import AU_KM, M_JUP, M_SUN
    from slingshot.v4.metrics import analyze_integration

    initial = np.array(
        [0, 0, 0, 0, 1e7, 0, 0, 1, -2e8, 3e7, 10, 0, 0, 0],
        dtype=float,
    )
    final = initial.copy()
    final[10:12] = [-10.0, 0.0]
    metrics = analyze_integration(
        _integration(initial, final),
        initial,
        M_SUN,
        M_JUP,
        AU_KM,
    )
    assert metrics["delta_kinetic_energy_com"] == pytest.approx(0.0, abs=1e-6)
    assert metrics["turning_quadratic"] > 100.0


def test_com_metrics_are_galilean_invariant():
    from slingshot.constants import AU_KM, M_JUP, M_SUN
    from slingshot.v4.metrics import analyze_integration

    initial = np.array(
        [-1e5, 0, 0, -1, 1e7, 0, 0, 100, -2e8, 4e7, 50, -5, 0, 0],
        dtype=float,
    )
    final = initial.copy()
    final[8:10] = [2e8, -3e7]
    final[10:12] = [65, 12]
    base = analyze_integration(
        _integration(initial, final), initial, M_SUN, M_JUP, AU_KM
    )
    boosted_initial = initial.copy()
    boosted_final = final.copy()
    for indices in ((2, 3), (6, 7), (10, 11)):
        boosted_initial[list(indices)] += [23.0, -17.0]
        boosted_final[list(indices)] += [23.0, -17.0]
    boosted = analyze_integration(
        _integration(boosted_initial, boosted_final),
        boosted_initial,
        M_SUN,
        M_JUP,
        AU_KM,
    )
    for key in (
        "delta_specific_energy_com",
        "delta_speed_com",
        "delta_v_inf",
        "turning_quadratic",
        "deflection_rad",
    ):
        assert boosted[key] == pytest.approx(base[key], rel=1e-12, abs=1e-12)


def test_legacy_metric_alias_warns():
    from slingshot.v4.metrics import resolve_metric

    with pytest.warns(DeprecationWarning):
        value = resolve_metric({"turning_quadratic": 12.0}, "half_dv_vec_sq")
    assert value == 12.0
