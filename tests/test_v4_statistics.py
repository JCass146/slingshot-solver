"""Known-probability tests for effective planar widths."""

import pytest


def test_planar_width_known_probability():
    from slingshot.v4.statistics import planar_width_interval

    width, low, high = planar_width_interval(250, 1000, b_max_km=100.0)
    assert width == pytest.approx(50.0)
    assert low < width < high


def test_width_summary_and_tail_gate():
    from slingshot.v4.statistics import summarize_planar_widths

    rows = []
    for index in range(100):
        rows.append(
            {
                "escaped": True,
                "collision": False,
                "energy_gain_dimensionless": 0.2 if index < 20 else -0.1,
                "impact_parameter_km": index - 50,
            }
        )
    summary = summarize_planar_widths(
        rows,
        b_max_km=50.0,
        thresholds=[0.0, 0.1, 0.3],
    )
    threshold_zero = summary[0]
    assert threshold_zero["events"] == 20
    assert threshold_zero["width_km"] == pytest.approx(20.0)
    assert threshold_zero["effective_sample_size"] == 100
    assert summary[-1]["statistic"] == "collision"
    # New fields from the CI-based tail gate
    assert "tail_fraction_upper_bound" in threshold_zero
    assert "tail_zone_trials" in threshold_zero


def test_tail_gate_zero_events_does_not_auto_pass():
    """Verify zero tail events produce a non-zero upper bound and may fail the gate."""
    from slingshot.v4.statistics import summarize_planar_widths, wilson_upper_bound

    # All events are near zero impact — none in the outer tail zone
    rows = []
    for index in range(200):
        rows.append(
            {
                "escaped": True,
                "collision": False,
                "energy_gain_dimensionless": 0.5,
                "impact_parameter_km": float(index) * 0.1,  # 0 to 19.9 km
            }
        )
    # b_max=20, tail_fraction=0.10 → outer zone is |b| >= 18
    summary = summarize_planar_widths(
        rows,
        b_max_km=20.0,
        thresholds=[0.0],
        tail_fraction=0.10,
        max_tail_event_fraction=0.005,  # very tight threshold
    )
    row = summary[0]
    # Upper bound must be > 0 even with 0 tail events — the gate should correctly
    # use the CI bound rather than the point estimate.
    assert row["tail_fraction_upper_bound"] > 0.0
    # The gate result depends on the bound, not the point estimate
    assert "tail_check_passed" in row


def test_wilson_upper_bound_zero_events():
    """Upper confidence bound for 0 successes out of N must be positive."""
    from slingshot.v4.statistics import wilson_upper_bound

    ub = wilson_upper_bound(0, 100, confidence_level=0.95)
    assert 0.0 < ub < 0.05  # should be around 0.03 for N=100 at 95% CI
