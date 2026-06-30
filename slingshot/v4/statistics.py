"""Effective planar-width estimators and confidence intervals."""

from __future__ import annotations

from statistics import NormalDist
from typing import Iterable

import numpy as np


ABILITY_CLAIM_THRESHOLD_MIN = 0.01


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


def wilson_upper_bound(
    successes: int, trials: int, confidence_level: float = 0.95
) -> float:
    """One-sided Wilson upper confidence bound for a binomial proportion.

    Used for the tail-support gate when event counts are small or zero: a zero
    observation does NOT mean the tail contribution is negligible — we report
    the upper bound so the gate can require that even the upper bound is below
    the declared threshold.
    """
    if trials <= 0:
        return np.nan
    z_score = NormalDist().inv_cdf(confidence_level)
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
    return float(min(1.0, center + half_width))


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
        # Outer-strip trials: samples whose |b| falls in the outer tail zone
        tail_zone_trials = sum(
            1 for row in rows
            if abs(float(row["impact_parameter_km"])) >= tail_limit
        )
        # Use one-sided upper Wilson bound so that zero observed tail events
        # does NOT automatically pass — require the upper bound to be below
        # the declared threshold.
        tail_fraction_upper_bound = wilson_upper_bound(
            int(tail_events),
            max(tail_zone_trials, 1),
            confidence_level,
        )
        tail_event_fraction = tail_events / event_count if event_count else 0.0
        # Gate passes only when the upper CI bound is below the threshold AND
        # there are enough outer-strip trials to be informative.
        tail_check_passed = (
            tail_zone_trials > 0
            and tail_fraction_upper_bound <= max_tail_event_fraction
        )
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
                "tail_zone_trials": int(tail_zone_trials),
                "tail_event_fraction": float(tail_event_fraction),
                "tail_fraction_upper_bound": float(tail_fraction_upper_bound),
                "tail_check_passed": tail_check_passed,
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


def _tail_gate_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def summarize_tail_gate_status(
    summary_rows: Iterable[dict],
    claim_threshold_min: float = ABILITY_CLAIM_THRESHOLD_MIN,
) -> dict:
    """Summarize formal ability-tail gates separately from q=0 diagnostics.

    The q=0 row is useful for diagnosing weak distant perturbations, but the
    defensible v4 ability claim starts at the first declared positive threshold
    at or above claim_threshold_min.
    """
    combined = [
        row for row in summary_rows
        if row.get("scope") == "combined"
        and row.get("statistic") == "energy_threshold"
    ]
    q0_rows = []
    claim_rows = []
    for row in combined:
        try:
            threshold = float(row.get("threshold", "nan"))
        except (TypeError, ValueError):
            continue
        if abs(threshold) <= 1e-12:
            q0_rows.append(row)
        if threshold >= claim_threshold_min - 1e-12:
            claim_rows.append(row)

    failed_claim_rows = [
        row for row in claim_rows if not _tail_gate_bool(row.get("tail_check_passed"))
    ]
    failed_q0_rows = [
        row for row in q0_rows if not _tail_gate_bool(row.get("tail_check_passed"))
    ]

    def _failures(rows: list[dict]) -> list[dict]:
        failures = []
        for row in rows:
            failures.append(
                {
                    "v_inf_kms": float(row.get("v_inf_kms", np.nan)),
                    "threshold": float(row.get("threshold", np.nan)),
                    "tail_fraction_upper_bound": float(
                        row.get("tail_fraction_upper_bound", np.nan)
                    ),
                }
            )
        return failures

    return {
        "tail_checks_passed": bool(claim_rows) and not failed_claim_rows,
        "ability_tail_checks_passed": bool(claim_rows) and not failed_claim_rows,
        "ability_tail_threshold_min": float(claim_threshold_min),
        "ability_tail_thresholds": sorted(
            {float(row.get("threshold", np.nan)) for row in claim_rows}
        ),
        "ability_tail_failures": _failures(failed_claim_rows),
        "q0_tail_diagnostic_passed": bool(q0_rows) and not failed_q0_rows,
        "q0_tail_diagnostic_failures": _failures(failed_q0_rows),
    }