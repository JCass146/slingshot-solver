"""Effective planar-width estimators and confidence intervals."""

from __future__ import annotations

from statistics import NormalDist
from typing import Iterable

import numpy as np


def wilson_interval(
    successes: int, trials: int, confidence_level: float = 0.95
) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if trials <= 0:
        return np.nan, np.nan
    z_score = NormalDist().inv_cdf(0.5 + confidence_level / 2.0)
    proportion = successes / trials
    denominator = 1.0 + z_score**2 / trials
    center = (proportion + z_score**2 / (2.0 * trials)) / denominator
    half_width = (
        z_score
        * np.sqrt(
            proportion * (1.0 - proportion) / trials
            + z_score**2 / (4.0 * trials**2)
        )
        / denominator
    )
    return max(0.0, center - half_width), min(1.0, center + half_width)


def planar_width_interval(
    successes: int,
    trials: int,
    b_max_km: float,
    confidence_level: float = 0.95,
) -> tuple[float, float, float]:
    """Estimate W = 2 b_max p and transform its Wilson interval."""
    probability = successes / trials if trials else np.nan
    probability_low, probability_high = wilson_interval(
        successes, trials, confidence_level
    )
    scale = 2.0 * b_max_km
    return (
        scale * probability,
        scale * probability_low,
        scale * probability_high,
    )


def summarize_planar_widths(
    records: Iterable[dict],
    b_max_km: float,
    thresholds: list[float],
    confidence_level: float = 0.95,
    tail_fraction: float = 0.10,
    max_tail_event_fraction: float = 0.01,
) -> list[dict]:
    """Summarize threshold widths, collision width, tails, and gain quantiles."""
    rows = list(records)
    trials = len(rows)
    escaped_gains = np.array(
        [
            row["energy_gain_dimensionless"]
            for row in rows
            if row.get("escaped") and np.isfinite(row.get("energy_gain_dimensionless", np.nan))
        ],
        dtype=float,
    )
    output = []
    tail_limit = (1.0 - tail_fraction) * b_max_km
    for threshold in thresholds:
        event_mask = [
            bool(row.get("escaped"))
            and float(row.get("energy_gain_dimensionless", -np.inf)) > threshold
            for row in rows
        ]
        event_count = int(sum(event_mask))
        width, width_low, width_high = planar_width_interval(
            event_count, trials, b_max_km, confidence_level
        )
        tail_events = sum(
            event
            and abs(float(row["impact_parameter_km"])) >= tail_limit
            for row, event in zip(rows, event_mask)
        )
        tail_event_fraction = tail_events / event_count if event_count else 0.0
        output.append(
            {
                "statistic": "energy_threshold",
                "threshold": float(threshold),
                "trials": trials,
                "events": event_count,
                "effective_sample_size": float(trials),
                "width_km": width,
                "width_low_km": width_low,
                "width_high_km": width_high,
                "tail_events": int(tail_events),
                "tail_event_fraction": float(tail_event_fraction),
                "tail_check_passed": tail_event_fraction <= max_tail_event_fraction,
                "median_gain": (
                    float(np.median(escaped_gains)) if escaped_gains.size else np.nan
                ),
                "gain_q90": (
                    float(np.quantile(escaped_gains, 0.90))
                    if escaped_gains.size
                    else np.nan
                ),
                "gain_q95": (
                    float(np.quantile(escaped_gains, 0.95))
                    if escaped_gains.size
                    else np.nan
                ),
                "gain_q99": (
                    float(np.quantile(escaped_gains, 0.99))
                    if escaped_gains.size
                    else np.nan
                ),
            }
        )

    collision_count = sum(bool(row.get("collision")) for row in rows)
    width, width_low, width_high = planar_width_interval(
        collision_count, trials, b_max_km, confidence_level
    )
    output.append(
        {
            "statistic": "collision",
            "threshold": np.nan,
            "trials": trials,
            "events": int(collision_count),
            "effective_sample_size": float(trials),
            "width_km": width,
            "width_low_km": width_low,
            "width_high_km": width_high,
            "tail_events": 0,
            "tail_event_fraction": 0.0,
            "tail_check_passed": True,
            "median_gain": np.nan,
            "gain_q90": np.nan,
            "gain_q95": np.nan,
            "gain_q99": np.nan,
        }
    )
    return output
