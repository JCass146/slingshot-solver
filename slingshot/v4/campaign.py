"""Reproducible planar-width campaign runner and artifact writer."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version as pkg_version
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from ..constants import AU_KM, G_KM
from .candidates import write_top_candidates_for_run
from .config import V4Config, load_config, save_config
from .dynamics import init_binary_barycentric, integrate_encounter
from .metrics import analyze_integration
from .report import generate_report
from .sampling import draw_samples
from .statistics import summarize_planar_widths, summarize_tail_gate_status
from .validation import physical_values, run_quick_validation


def _package_version() -> str:
    try:
        return pkg_version("slingshot-solver")
    except (PackageNotFoundError, Exception):
        # Fall back to reading pyproject.toml when the package is not installed
        try:
            import re
            _root = Path(__file__).parent.parent.parent
            text = (_root / "pyproject.toml").read_text(encoding="utf-8")
            match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
            return match.group(1) if match else "unknown"
        except Exception:
            return "unknown"


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


def _append_seed_variance_rows(
    summary_records: list[dict],
    v_inf: float,
    seeds: list[int],
    thresholds: list[float],
) -> None:
    """Compute between-seed variance and heterogeneity for each threshold.

    For each threshold, collect the per-seed width estimates already written
    to summary_records, compute the between-seed standard deviation and a
    simple heterogeneity ratio, and append a 'seed_variance' scope row.
    """
    for threshold in thresholds:
        seed_widths = []
        for row in summary_records:
            if (
                row.get("scope") == "seed"
                and row.get("v_inf_kms") == v_inf
                and row.get("statistic") == "energy_threshold"
                and float(row.get("threshold", -1.0)) == float(threshold)
            ):
                w = row.get("width_km")
                if w is not None and np.isfinite(float(w)):
                    seed_widths.append(float(w))
        if len(seed_widths) < 2:
            continue
        arr = np.array(seed_widths)
        mean_w = float(np.mean(arr))
        std_w = float(np.std(arr, ddof=1))
        # Heterogeneity ratio: between-seed std / pooled-width-estimate.
        # A ratio near 0 indicates low between-seed variability.
        combined_rows = [
            row for row in summary_records
            if (
                row.get("scope") == "combined"
                and row.get("v_inf_kms") == v_inf
                and row.get("statistic") == "energy_threshold"
                and float(row.get("threshold", -1.0)) == float(threshold)
            )
        ]
        pooled_width = float(combined_rows[0]["width_km"]) if combined_rows else mean_w
        heterogeneity = std_w / max(abs(pooled_width), 1.0)
        summary_records.append(
            {
                "scope": "seed_variance",
                "seed": "",
                "v_inf_kms": v_inf,
                "statistic": "energy_threshold",
                "threshold": float(threshold),
                "n_seeds": len(seed_widths),
                "seed_mean_width_km": mean_w,
                "seed_std_width_km": std_w,
                "seed_heterogeneity": heterogeneity,
                "width_km": mean_w,
                "width_low_km": mean_w - std_w,
                "width_high_km": mean_w + std_w,
                "width_au": mean_w / AU_KM,
                "width_low_au": (mean_w - std_w) / AU_KM,
                "width_high_au": (mean_w + std_w) / AU_KM,
            }
        )


def run_campaign(
    config_path: str | Path,
    output_dir: Optional[str | Path] = None,
    samples_per_bin: Optional[int] = None,
    seeds: Optional[Sequence[int]] = None,
    verbose: bool = True,
) -> dict:
    """Run all configured v-infinity bins and write compact artifacts."""
    config = load_config(config_path)
    if samples_per_bin is not None:
        config.asymptotic_sampling.samples_per_bin = int(samples_per_bin)
    if seeds is not None:
        config.asymptotic_sampling.seeds = [int(seed) for seed in seeds]

    started = datetime.now(timezone.utc)
    if output_dir is None:
        stamp = started.strftime("%Y%m%d_%H%M%S")
        output_path = Path(
            f"results/{config.metadata.case_name}_{stamp}"
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
            # Derive an independent child RNG for this (v_inf, seed) pair so
            # that different speed bins do not share common random numbers
            # unless a paired CRN design is explicitly declared (P1.3).
            # Mix the v_inf value into the stream via a secondary seed
            v_inf_bits = abs(hash(float(v_inf))) & 0xFFFFFFFFFFFFFFFF
            child_seq = np.random.SeedSequence([seed, int(v_inf_bits)])
            rng = np.random.default_rng(child_seq)
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

        # Seed-level variance and heterogeneity diagnostics (P1.2)
        _append_seed_variance_rows(
            summary_records,
            v_inf,
            config.asymptotic_sampling.seeds,
            config.planar_width.dimensionless_energy_thresholds,
        )

    _write_csv(output_path / "samples.csv", sample_records)
    _write_csv(output_path / "width_summary.csv", summary_records)
    candidate_records = []
    if config.candidate_diagnostics.enabled:
        candidate_records = write_top_candidates_for_run(
            output_path, top_n=config.candidate_diagnostics.top_n
        )
    completed = datetime.now(timezone.utc)
    closure_values = [
        float(row["work_energy_closure_relative"])
        for row in sample_records
        if np.isfinite(row.get("work_energy_closure_relative", np.nan))
    ]
    tail_status = summarize_tail_gate_status(summary_records)
    # Campaign-level failure gates (P0.2)
    total_samples = len(sample_records)
    time_limit_count = sum(
        1 for row in sample_records if row.get("outcome") == "time_limit"
    )
    numerical_failure_count = sum(
        1 for row in sample_records if not row.get("solver_success", True)
        and row.get("outcome") not in {"escaped", "star_collision", "planet_collision", "time_limit"}
    )
    time_limit_fraction = time_limit_count / total_samples if total_samples else 0.0
    numerical_failure_fraction = numerical_failure_count / total_samples if total_samples else 0.0
    time_limit_passed = time_limit_fraction <= config.validation.max_time_limit_fraction
    numerical_failure_passed = (
        numerical_failure_fraction <= config.validation.max_numerical_failure_fraction
    )
    campaign_validation = {
        "quick": validation,
        "work_energy_max_relative": max(closure_values) if closure_values else np.nan,
        "work_energy_passed": (
            bool(closure_values)
            and max(closure_values)
            <= config.validation.work_energy_relative_tolerance
        ),
        **tail_status,
        "time_limit_count": time_limit_count,
        "time_limit_fraction": time_limit_fraction,
        "time_limit_passed": time_limit_passed,
        "numerical_failure_count": numerical_failure_count,
        "numerical_failure_fraction": numerical_failure_fraction,
        "numerical_failure_passed": numerical_failure_passed,
    }
    campaign_validation["passed"] = (
        validation["passed"]
        and campaign_validation["work_energy_passed"]
        and campaign_validation["tail_checks_passed"]
        and campaign_validation["time_limit_passed"]
        and campaign_validation["numerical_failure_passed"]
    )
    manifest = {
        "schema_version": 4,
        "science_model": "defensible_planar_width",
        "legacy_science_model": False,
        "package_version": _package_version(),
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
        "candidate_diagnostics": config.candidate_diagnostics.model_dump(mode="json"),
        "candidate_count": len(candidate_records),
        "best_observed_gain": (
            float(candidate_records[0]["energy_gain_dimensionless"])
            if candidate_records else None
        ),
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
            "top_candidates.csv",
            "manifest.json",
            "REPORT.md",
        ],
    }
    with (output_path / "manifest.json").open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, allow_nan=False)

    # Generate diagnostic figures before writing the report so the report can
    # include links to any figures that were successfully created.
    generated_plots = []
    try:
        from .plotting import generate_all_plots
        if verbose:
            print("Generating diagnostic figures...")
        generated_plots = generate_all_plots(output_path, verbose=verbose)
    except Exception as _plot_exc:
        if verbose:
            print(f"  Plotting skipped: {_plot_exc}")

    for artifact in [Path(path).name for path in generated_plots]:
        if artifact not in manifest["artifacts"]:
            manifest["artifacts"].append(artifact)
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
