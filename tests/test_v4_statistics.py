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
