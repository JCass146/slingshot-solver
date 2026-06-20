"""Reproducible planar-width campaign runner and artifact writer."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from ..constants import AU_KM, G_KM
from .config import V4Config, load_config, save_config
from .dynamics import init_binary_barycentric, integrate_encounter
from .metrics import analyze_integration
from .report import generate_report
from .sampling import draw_samples
from .statistics import summarize_planar_widths
from .validation import physical_values, run_quick_validation


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _config_hash(config: V4Config) -> str:
    payload = json.dumps(
        config.model_dump(mode="json"), sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _evaluate_sample(
    config: V4Config,
    sample,
    seed: int,
    sample_index: int,
    values: dict[str, float],
) -> dict:
    binary_state = init_binary_barycentric(
        semi_major_axis_km=values["semi_major_axis_km"],
        eccentricity=config.orbit.eccentricity,
        mean_anomaly_rad=sample.binary_mean_anomaly_rad,
        argument_periapsis_rad=config.orbit.argument_periapsis_rad,
        star_mass_kg=values["star_mass_kg"],
        planet_mass_kg=values["planet_mass_kg"],
        prograde=config.orbit.prograde,
        bulk_velocity_x_kms=config.system.bulk_velocity_x_kms,
        bulk_velocity_y_kms=config.system.bulk_velocity_y_kms,
    )
    bulk_velocity = np.array(
        [
            config.system.bulk_velocity_x_kms,
            config.system.bulk_velocity_y_kms,
        ]
    )
    initial_state = np.concatenate(
        [
            binary_state,
            sample.position_km,
            sample.velocity_kms + bulk_velocity,
            np.zeros(2),
        ]
    )
    integration = integrate_encounter(
        initial_state=initial_state,
        star_mass_kg=values["star_mass_kg"],
        planet_mass_kg=values["planet_mass_kg"],
        star_radius_km=values["star_radius_km"],
        planet_radius_km=values["planet_radius_km"],
        boundary_radius_km=values["boundary_radius_km"],
        max_time_sec=config.numerical.max_time_sec,
        method=config.numerical.method,
        rtol=config.numerical.rtol,
        atol=config.numerical.atol,
        softening_km=config.numerical.softening_km,
        max_step_sec=config.numerical.max_step_sec,
    )
    metrics = analyze_integration(
        integration=integration,
        initial_state=initial_state,
        star_mass_kg=values["star_mass_kg"],
        planet_mass_kg=values["planet_mass_kg"],
        semi_major_axis_km=values["semi_major_axis_km"],
        softening_km=config.numerical.softening_km,
    )
    return {
        "sample_index": sample_index,
        "seed": seed,
        "v_inf_kms": sample.v_inf_kms,
        "impact_parameter_km": sample.impact_parameter_km,
        "impact_parameter_au": sample.impact_parameter_km / AU_KM,
        "incoming_direction_rad": sample.incoming_direction_rad,
        "binary_mean_anomaly_rad": sample.binary_mean_anomaly_rad,
        "proposal_specific_energy_km2_s2": sample.specific_energy_km2_s2,
        "proposal_angular_momentum_km2_s": sample.angular_momentum_km2_s,
        **metrics,
    }


def run_campaign(
    config_path: str | Path,
    output_dir: Optional[str | Path] = None,
    samples_per_bin: Optional[int] = None,
    seeds: Optional[Sequence[int]] = None,
    verbose: bool = True,
) -> dict:
    """Run all configured v-infinity bins and write compact v4 artifacts."""
    config = load_config(config_path)
    if samples_per_bin is not None:
        config.asymptotic_sampling.samples_per_bin = int(samples_per_bin)
    if seeds is not None:
        config.asymptotic_sampling.seeds = [int(seed) for seed in seeds]

    started = datetime.now(timezone.utc)
    if output_dir is None:
        stamp = started.strftime("%Y%m%d_%H%M%S")
        output_path = Path(
            f"results/v4_{config.metadata.case_name}_{stamp}"
        )
    else:
        output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    save_config(config, output_path / "config.yaml")

    validation = run_quick_validation(config)
    values = physical_values(config)
    total_mu = G_KM * (values["star_mass_kg"] + values["planet_mass_kg"])
    sample_records: list[dict] = []
    summary_records: list[dict] = []

    for v_inf in config.asymptotic_sampling.v_inf_kms:
        combined_records = []
        for seed in config.asymptotic_sampling.seeds:
            if verbose:
                print(
                    f"v∞={v_inf:g} km/s seed={seed} "
                    f"N={config.asymptotic_sampling.samples_per_bin}"
                )
            rng = np.random.default_rng(seed)
            proposals = draw_samples(
                rng=rng,
                count=config.asymptotic_sampling.samples_per_bin,
                v_inf_kms=v_inf,
                b_max_km=values["b_max_km"],
                boundary_radius_km=values["boundary_radius_km"],
                total_mu_km3_s2=total_mu,
                sample_incoming_direction=config.asymptotic_sampling.sample_incoming_direction,
                sample_binary_mean_anomaly=config.asymptotic_sampling.sample_binary_mean_anomaly,
                fixed_mean_anomaly_rad=config.orbit.mean_anomaly_rad,
            )
            seed_records = []
            for sample_index, sample in enumerate(proposals):
                record = _evaluate_sample(
                    config, sample, seed, sample_index, values
                )
                seed_records.append(record)
                sample_records.append(record)
                if verbose and (sample_index + 1) % max(1, len(proposals) // 10) == 0:
                    print(f"  completed {sample_index + 1}/{len(proposals)}")
            combined_records.extend(seed_records)
            seed_summary = summarize_planar_widths(
                seed_records,
                b_max_km=values["b_max_km"],
                thresholds=config.planar_width.dimensionless_energy_thresholds,
                confidence_level=config.planar_width.confidence_level,
                tail_fraction=config.planar_width.tail_fraction,
                max_tail_event_fraction=config.planar_width.max_tail_event_fraction,
            )
            for row in seed_summary:
                row.update(
                    {
                        "scope": "seed",
                        "seed": seed,
                        "v_inf_kms": v_inf,
                        "width_au": row["width_km"] / AU_KM,
                        "width_low_au": row["width_low_km"] / AU_KM,
                        "width_high_au": row["width_high_km"] / AU_KM,
                    }
                )
                summary_records.append(row)

        combined_summary = summarize_planar_widths(
            combined_records,
            b_max_km=values["b_max_km"],
            thresholds=config.planar_width.dimensionless_energy_thresholds,
            confidence_level=config.planar_width.confidence_level,
            tail_fraction=config.planar_width.tail_fraction,
            max_tail_event_fraction=config.planar_width.max_tail_event_fraction,
        )
        for row in combined_summary:
            row.update(
                {
                    "scope": "combined",
                    "seed": "",
                    "v_inf_kms": v_inf,
                    "width_au": row["width_km"] / AU_KM,
                    "width_low_au": row["width_low_km"] / AU_KM,
                    "width_high_au": row["width_high_km"] / AU_KM,
                }
            )
            summary_records.append(row)

    _write_csv(output_path / "samples.csv", sample_records)
    _write_csv(output_path / "width_summary.csv", summary_records)
    completed = datetime.now(timezone.utc)
    closure_values = [
        float(row["work_energy_closure_relative"])
        for row in sample_records
        if np.isfinite(row.get("work_energy_closure_relative", np.nan))
    ]
    tail_passed = all(
        bool(row["tail_check_passed"])
        for row in summary_records
        if row["scope"] == "combined" and row["statistic"] == "energy_threshold"
    )
    campaign_validation = {
        "quick": validation,
        "work_energy_max_relative": max(closure_values) if closure_values else np.nan,
        "work_energy_passed": (
            bool(closure_values)
            and max(closure_values)
            <= config.validation.work_energy_relative_tolerance
        ),
        "tail_checks_passed": tail_passed,
    }
    campaign_validation["passed"] = (
        validation["passed"]
        and campaign_validation["work_energy_passed"]
        and campaign_validation["tail_checks_passed"]
    )
    manifest = {
        "schema_version": 4,
        "science_model": "defensible_planar_width",
        "legacy_science_model": False,
        "package_version": "4.0.0",
        "git_commit": _git_commit(),
        "config_sha256": _config_hash(config),
        "started_utc": started.isoformat(),
        "completed_utc": completed.isoformat(),
        "duration_sec": (completed - started).total_seconds(),
        "sample_count": len(sample_records),
        "samples_per_bin": config.asymptotic_sampling.samples_per_bin,
        "seeds": list(config.asymptotic_sampling.seeds),
        "v_inf_kms": list(config.asymptotic_sampling.v_inf_kms),
        "integrator": config.numerical.method,
        "softening_km": config.numerical.softening_km,
        "observational_source": config.metadata.parameter_source,
        "observational_metadata": config.metadata.model_dump(mode="json"),
        "citation": config.metadata.citation,
        "validation_status": (
            "passed" if campaign_validation["passed"] else "failed"
        ),
        "validation": campaign_validation,
        "artifacts": [
            "config.yaml",
            "samples.csv",
            "width_summary.csv",
            "manifest.json",
            "REPORT.md",
        ],
    }
    with (output_path / "manifest.json").open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, allow_nan=False)
    report = generate_report(output_path, config, summary_records, manifest)
    return {
        "config": config,
        "output_dir": output_path,
        "samples": sample_records,
        "summary": summary_records,
        "manifest": manifest,
        "report": report,
    }
