"""Concise scientific report for planar-width campaigns."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from ..constants import AU_KM
from .config import V4Config


def generate_report(
    output_dir: str | Path,
    config: V4Config,
    summary_rows: Iterable[dict],
    manifest: dict,
) -> str:
    output_path = Path(output_dir)
    rows = list(summary_rows)
    combined = [row for row in rows if row.get("scope") == "combined"]
    lines = [
        f"# Slingshot Solver v4 Planar Research Report — {config.system.name}",
        "",
        "## Interpretation Boundary",
        "",
        "This run estimates **effective planar encounter widths** under a declared "
        "uniform signed-impact proposal. It does not estimate a three-dimensional "
        "area cross-section, an astrophysical event probability, or an occurrence rate.",
        "",
        "## Observational Model",
        "",
        f"- Case: `{config.metadata.case_name}`",
        f"- Source: {config.metadata.parameter_source}",
        f"- Citation: {config.metadata.citation}",
        f"- Binary orbit: a={config.orbit.semi_major_axis_au:.6g} AU, "
        f"e={config.orbit.eccentricity:.6g}",
        f"- Proposal boundary: {config.asymptotic_sampling.boundary_radius_au:.6g} AU",
        f"- Signed impact interval: ±{config.asymptotic_sampling.b_max_au:.6g} AU",
        f"- Samples per speed/seed: {config.asymptotic_sampling.samples_per_bin}",
        f"- Seeds: {config.asymptotic_sampling.seeds}",
        "",
        "## Primary Estimand",
        "",
        "The primary gain is the binary-COM specific-energy change "
        "`delta_specific_energy_com`. Width thresholds use "
        "`delta_specific_energy_com / v_c²`, where `v_c² = G(M★+Mp)/a`.",
        "",
        "| v∞ (km/s) | Threshold Δε/vc² | Events / N | Width (AU) | 95% interval (AU) | Tail pass |",
        "|---:|---:|---:|---:|---:|:---:|",
    ]
    for row in combined:
        if row["statistic"] != "energy_threshold":
            continue
        lines.append(
            f"| {row['v_inf_kms']:.3g} | {row['threshold']:.3g} | "
            f"{row['events']} / {row['trials']} | "
            f"{row['width_km'] / AU_KM:.6g} | "
            f"[{row['width_low_km'] / AU_KM:.6g}, "
            f"{row['width_high_km'] / AU_KM:.6g}] | "
            f"{'yes' if row['tail_check_passed'] else 'no'} |"
        )
    lines.extend(
        [
            "",
            "## Metric Semantics",
            "",
            "- `turning_quadratic = ½|vout−vin|²` is retained only as a turning diagnostic.",
            "- It is not kinetic-energy gain and is not used for width thresholds.",
            "- Stellar and planetary work are integrated separately from the moving potentials.",
            "- Work-energy closure is recorded for every trajectory.",
            "",
            "## Validation",
            "",
            f"- Quick validation passed: {manifest['validation']['passed']}",
            f"- Newtonian scientific run: {config.numerical.softening_km == 0.0}",
            f"- Git commit: `{manifest.get('git_commit', 'unknown')}`",
            "",
            "## Limitations",
            "",
            "- The model is planar; inclination effects require the later 3D milestone.",
            "- Quinn and Ortiz parameter sets are treated as separate observational models.",
            "- Reported extrema are diagnostics only, not converged physical limits.",
            "- Physical rates require a 3D cross-section and an external ISO velocity distribution.",
            "",
        ]
    )
    report = "\n".join(lines)
    (output_path / "REPORT.md").write_text(report, encoding="utf-8")
    return report
