"""Comprehensive scientific report generator for planar-width campaigns."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

from ..constants import AU_KM
from .config import V4Config

# Figures produced by plotting.generate_all_plots, in display order
_FIGURE_CATALOGUE = [
    ("width_vs_vinf.png",
     "**Figure 1 — Planar-width estimates vs v∞.**  "
     "Combined Wilson 95% CI bands for every energy threshold, with individual "
     "per-seed points overlaid at Δε/vc²>0. The gray dashed line marks 2b_max "
     "(the sampling ceiling). All widths are lower bounds until the tail gate passes."),
    ("outcome_fractions.png",
     "**Figure 2 — Outcome fractions vs v∞.**  "
     "Stacked bar chart showing the fraction of all samples that escape, collide "
     "with the star, collide with the planet, hit the time limit, or fail numerically. "
     "Star-collision fraction decreases strongly with speed; escape fraction increases."),
    ("collision_vs_escape.png",
     "**Figure 3 — Escape vs collision widths vs v∞.**  "
     "Effective planar widths for escape (Δε/vc²>0) and any collision, with 95% "
     "Wilson CI bands. At low v∞ the collision cross-section exceeds the escape "
     "width; at high v∞ the ordering reverses."),
    ("tail_support.png",
     "**Figure 4 — Tail-support diagnostic (middle speed bin).**  "
     "Escape event probability binned by |b|/b_max with Wilson CI. If the event "
     "rate in the outer 10% strip does not fall toward zero, b_max is too small "
     "and the width estimates are underestimates."),
    ("seed_stability.png",
     "**Figure 5 — Seed stability (Δε/vc²>0).**  "
     "Per-seed planar-width curves (dashed) vs the pooled Wilson CI (solid black). "
     "Low between-seed scatter indicates stable estimation."),
    ("sampling_distributions.png",
     "**Figure 6 — Proposal and acceptance distributions.**  "
     "Marginal histograms of signed impact parameter, incoming direction, binary "
     "mean anomaly, and speed bin. Departure from uniform acceptance reveals "
     "parameter regions that drive events."),
    ("gain_ecdf.png",
     "**Figure 7 — COM energy-gain ECDF by v∞ (escaped trajectories only).**  "
     "Empirical CDFs of Δε/vc² conditioned on escape. Higher speed bins shift "
     "the distribution rightward."),
    ("deflection_distribution.png",
     "**Figure 8 — Deflection-angle distributions by v∞ (escaped).**  "
     "COM-frame deflection angles for escaped trajectories."),
    ("velocity_phase_space.png",
     "**Figure 9 — COM-frame velocity phase space.**  "
     "Initial vs final COM speed coloured by energy gain, and Δ|v|_COM distributions. "
     "The radial-normal decomposition is omitted — the legacy radial basis was "
     "computed relative to the star-planet axis rather than the particle-planet vector."),
    ("phase_map.png",
     "**Figure 10 — Phase map: conditional mean energy gain (middle speed bin).**  "
     "Cell-averaged Δε/vc² over (impact parameter, binary mean anomaly) for escaped "
     "trajectories, with a support-count panel. Uses conditional means, not per-cell "
     "maxima over a hidden velocity dimension."),
    ("periapsis_distributions.png",
     "**Figure 11 — Periapsis distributions.**  "
     "Minimum distance to planet and star in units of the body radius. "
     "Stellar-surface-crossing trajectories are tracked as collisions."),
    ("work_energy_diagnostics.png",
     "**Figure 12 — Work-energy diagnostics.**  "
     "Closure residual distribution, signed work fractions, and planet work vs "
     "energy gain for escaped trajectories."),
    ("parameter_correlations.png",
     "**Figure 13 — Parameter correlations (escaped trajectories).**  "
     "Energy gain vs planet periapsis, energy gain vs deflection, and deflection "
     "vs star proximity, all coloured by v∞."),
    ("candidate_ranking.png",
     "**Figure 14 — Top-30 candidate ranking.**  "
     "Top 30 escaped trajectories by Δε/vc²: panels show gain, deflection, and "
     "planet periapsis. These are illustrations — not converged optima."),
    ("pareto_front.png",
     "**Figure 15 — Pareto fronts (escaped trajectories).**  "
     "Left: (maximise gain, minimise periapsis). Right: (maximise gain, maximise "
     "|deflection|). Both panels use the current scientific metrics."),
    ("trajectory_tracks.png",
     "**Figure 16 — Top-10 trajectory tracks in the planet frame.**  "
     "Re-integrated from asymptotic parameters in samples.csv; coloured by Δε/vc². "
     "These are examples, not a probability density."),
]


def _fmt(value, fmt: str = ".4g") -> str:
    try:
        v = float(value)
        if math.isnan(v):
            return "—"
        return format(v, fmt)
    except (TypeError, ValueError):
        return str(value) if value not in (None, "", "nan") else "—"


def _pct(value) -> str:
    try:
        return f"{float(value)*100:.2f}%"
    except (TypeError, ValueError):
        return "—"


def generate_report(
    output_dir: str | Path,
    config: V4Config,
    summary_rows: Iterable[dict],
    manifest: dict,
) -> str:
    output_path = Path(output_dir)
    rows = list(summary_rows)
    combined = [r for r in rows if r.get("scope") == "combined"]
    svar_rows = [r for r in rows if r.get("scope") == "seed_variance"]
    val = manifest.get("validation", {})
    quick_gates = val.get("quick", {}).get("gates", [])

    lines = [
        f"# Slingshot Solver — Planar Research Report",
        f"## {config.system.name}",
        "",
        "| | |",
        "|---|---|",
        f"| **Started** | {manifest.get('started_utc', '—')} |",
        f"| **Duration** | {manifest.get('duration_sec', 0)/3600:.2f} h |",
        f"| **Validation** | `{manifest.get('validation_status', '—').upper()}` |",
        "",
        "---",
        "",
        "## Interpretation Boundary",
        "",
        "This report presents **effective planar encounter widths** estimated under "
        "a declared uniform signed-impact-parameter proposal over a fixed speed grid. "
        "The reported widths are **not** three-dimensional area cross-sections, "
        "astrophysical event probabilities, or occurrence rates. "
        "Sample extrema (e.g. highest energy-gain trajectories) are exploratory "
        "illustrations only — they are not converged physical limits.",
        "",
        "---",
        "",
        "## Observational Model",
        "",
        "| Parameter | Value |",
        "|---|---|",
        f"| Case | `{config.metadata.case_name}` |",
        f"| Source | {config.metadata.parameter_source} |",
        f"| Citation | {config.metadata.citation} |",
        f"| Star mass | {config.system.star_mass_msun:.4g} M☉ |",
        f"| Star radius | {config.system.star_radius_rsun:.4g} R☉ |",
        f"| Planet mass | {config.system.planet_mass_mjup:.4g} M_Jup |",
        f"| Planet radius | {config.system.planet_radius_rjup:.4g} R_Jup |",
        f"| Orbit model | {config.orbit.model} |",
        f"| Semi-major axis | {config.orbit.semi_major_axis_au:.6g} AU |",
        f"| Eccentricity | {config.orbit.eccentricity:.6g} |",
        f"| Bulk velocity | ({config.system.bulk_velocity_x_kms:.3g},"
        f" {config.system.bulk_velocity_y_kms:.3g}) km/s |",
    ]

    uncert = config.metadata.uncertainties
    if uncert:
        lines += [
            "",
            "### Observational Uncertainties",
            "",
            "| Parameter | -sigma | +sigma |",
            "|---|---|---|",
        ]
        for param, bounds in uncert.items():
            lines.append(
                f"| {param} | {_fmt(bounds.get('minus', 'nan'))} "
                f"| {_fmt(bounds.get('plus', 'nan'))} |"
            )

    lines += [
        "",
        "---",
        "",
        "## Sampling and Integration",
        "",
        "| Parameter | Value |",
        "|---|---|",
        f"| Speed bins v_inf (km/s) | {config.asymptotic_sampling.v_inf_kms} |",
        f"| b_max | {config.asymptotic_sampling.b_max_au:.4g} AU |",
        f"| Boundary radius | {config.asymptotic_sampling.boundary_radius_au:.4g} AU |",
        f"| Samples per bin | {config.asymptotic_sampling.samples_per_bin} |",
        f"| Seeds | {config.asymptotic_sampling.seeds} |",
        f"| Total samples | {manifest.get('sample_count', '—')} |",
        f"| Integrator | {config.numerical.method} |",
        f"| rtol / atol | {config.numerical.rtol:.2e} / {config.numerical.atol:.2e} |",
        f"| Softening | {config.numerical.softening_km:.3g} km |",
        f"| Max integration time | {config.numerical.max_time_sec:.2e} s |",
        "",
        "---",
        "",
        "## Primary Estimand",
        "",
        "W(delta_epsilon/vc^2 > q | v_inf) = 2 b_max * N_event / N",
        "",
        "- v_inf: incoming asymptotic speed",
        "- b_max: half-width of the uniform signed-impact proposal",
        f"- delta_epsilon = epsilon_out - epsilon_in in the binary COM frame",
        f"- vc^2 = G(M_star + M_planet)/a",
        f"- Wilson confidence intervals at {config.planar_width.confidence_level*100:.0f}%",
        "",
        "---",
        "",
        "## Planar-Width Results",
        "",
        "### All thresholds, all speed bins (combined)",
        "",
        "| v_inf (km/s) | Threshold | Events | N | Width (AU) | CI low (AU) | CI high (AU) | Tail UB | Tail gate |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|:---:|",
    ]

    for r in sorted(
        [r for r in combined if r["statistic"] == "energy_threshold"],
        key=lambda r: (float(r["v_inf_kms"]), float(r["threshold"])),
    ):
        tail_ok = str(r.get("tail_check_passed", "")) == "True"
        lines.append(
            f"| {float(r['v_inf_kms']):.0f} | {float(r['threshold']):.3g} | "
            f"{r['events']} | {r['trials']} | "
            f"{float(r['width_km'])/AU_KM:.5f} | "
            f"{float(r['width_low_km'])/AU_KM:.5f} | "
            f"{float(r['width_high_km'])/AU_KM:.5f} | "
            f"{_fmt(r.get('tail_fraction_upper_bound', 'nan'), '.4f')} | "
            f"{'yes' if tail_ok else 'no'} |"
        )

    col_rows = sorted(
        [r for r in combined if r["statistic"] == "collision"],
        key=lambda r: float(r["v_inf_kms"]),
    )
    if col_rows:
        lines += [
            "",
            "### Collision widths (combined)",
            "",
            "| v_inf (km/s) | Events | N | Width (AU) | CI low (AU) | CI high (AU) |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
        for r in col_rows:
            lines.append(
                f"| {float(r['v_inf_kms']):.0f} | {r['events']} | {r['trials']} | "
                f"{float(r['width_km'])/AU_KM:.5f} | "
                f"{float(r['width_low_km'])/AU_KM:.5f} | "
                f"{float(r['width_high_km'])/AU_KM:.5f} |"
            )

    lines += [
        "",
        "### Gain quantiles at threshold=0 (escaped, combined)",
        "",
        "| v_inf (km/s) | Median gain | Q90 | Q95 | Q99 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for r in sorted(
        [r for r in combined
         if r["statistic"] == "energy_threshold" and float(r["threshold"]) == 0.0],
        key=lambda r: float(r["v_inf_kms"]),
    ):
        lines.append(
            f"| {float(r['v_inf_kms']):.0f} | "
            f"{_fmt(r.get('median_gain', 'nan'))} | "
            f"{_fmt(r.get('gain_q90', 'nan'))} | "
            f"{_fmt(r.get('gain_q95', 'nan'))} | "
            f"{_fmt(r.get('gain_q99', 'nan'))} |"
        )

    if svar_rows:
        lines += [
            "",
            "---",
            "",
            "## Seed-Level Variance",
            "",
            "Heterogeneity = sigma_seed / W_pooled. "
            "Values below 0.05 indicate stable estimation.",
            "",
            "| v_inf (km/s) | Threshold | n seeds | Mean W (AU) | Std W (AU) | Heterogeneity |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
        for r in sorted(
            svar_rows,
            key=lambda r: (float(r["v_inf_kms"]), float(r.get("threshold", 0))),
        ):
            lines.append(
                f"| {float(r['v_inf_kms']):.0f} | "
                f"{float(r.get('threshold', 0)):.3g} | "
                f"{r.get('n_seeds', '—')} | "
                f"{float(r['seed_mean_width_km'])/AU_KM:.5f} | "
                f"{float(r['seed_std_width_km'])/AU_KM:.5f} | "
                f"{_fmt(r['seed_heterogeneity'], '.4f')} |"
            )

    lines += [
        "",
        "---",
        "",
        "## Validation",
        "",
        "### Campaign gates",
        "",
        "| Gate | Status | Value | Threshold |",
        "|---|:---:|---:|---:|",
        f"| Work-energy closure | {'PASS' if val.get('work_energy_passed') else 'FAIL'} | "
        f"{_fmt(val.get('work_energy_max_relative', 'nan'), '.3e')} | "
        f"{config.validation.work_energy_relative_tolerance:.2e} |",
        f"| Tail CI checks | {'PASS' if val.get('tail_checks_passed') else 'FAIL'} | — | "
        f"UB <= {config.planar_width.max_tail_event_fraction:.3g} |",
        f"| Time-limit fraction | {'PASS' if val.get('time_limit_passed', True) else 'FAIL'} | "
        f"{_pct(val.get('time_limit_fraction', 0))} "
        f"({val.get('time_limit_count', 0)} samples) | "
        f"{config.validation.max_time_limit_fraction:.0%} |",
        f"| Numerical failure | {'PASS' if val.get('numerical_failure_passed', True) else 'FAIL'} | "
        f"{_pct(val.get('numerical_failure_fraction', 0))} | "
        f"{config.validation.max_numerical_failure_fraction:.0%} |",
        "",
        "### Quick gates",
        "",
        "| Gate | Status | Notes |",
        "|---|:---:|---|",
    ]
    for gate in quick_gates:
        status = "PASS" if gate["passed"] else "FAIL"
        req = "" if gate.get("required", True) else " (diagnostic)"
        tol = f"tol={gate.get('tolerance', '')}" if "tolerance" in gate else ""
        lines.append(f"| {gate['name']}{req} | {status} | {tol} |")

    lines += [
        "",
        "---",
        "",
        "## Metric Semantics",
        "",
        "| Metric | Definition | Notes |",
        "|---|---|---|",
        "| `delta_specific_energy_com` | epsilon_out - epsilon_in in binary COM frame | Primary gain metric |",
        "| `energy_gain_dimensionless` | delta_epsilon / vc^2 | Threshold variable for width |",
        "| `delta_v_inf` | v_inf_out - v_inf_in | Defined only when both endpoints are unbound |",
        "| `turning_quadratic` | 0.5 * norm(v_out - v_in)^2 | Turning diagnostic only — not energy gain |",
        "| `work_star` | integral F_star dot v dt | Work by moving stellar potential |",
        "| `work_planet` | integral F_planet dot v dt | Work by moving planetary potential |",
        "| `work_energy_closure_relative` | abs(delta_epsilon - W_total) / scale | Integration quality |",
        "| `periapsis_planet_km` | Min test-particle–planet distance | Physical encounter depth |",
        "| `deflection_rad` | COM-frame velocity deflection angle | Boost-invariant |",
        "",
        "> `turning_quadratic` can be large for a constant-speed reversal with zero energy gain.",
        "",
        "---",
        "",
        "## Diagnostic Figures",
        "",
        "All figures read from `samples.csv` and `width_summary.csv`. "
        "Regenerate at any time with: `slingshot plot <run_dir>`",
        "",
    ]

    included = 0
    for filename, caption in _FIGURE_CATALOGUE:
        if (output_path / filename).exists():
            section = filename.replace(".png", "").replace("_", " ").title()
            lines += [
                f"### {section}",
                "",
                f"![{filename}]({filename})",
                "",
                caption,
                "",
            ]
            included += 1

    if included == 0:
        lines += [
            "*Figures not yet generated. Run:*",
            "```",
            "slingshot plot <run_dir>",
            "```",
            "",
        ]

    lines += [
        "---",
        "",
        "## Limitations",
        "",
        "1. **Planar model.** All trajectories are confined to the orbital plane. "
        "3D cross-sections require the future 3D milestone.",
        "",
        "2. **Proposal dependence.** W is defined under a specific uniform signed-impact "
        "proposal. Astrophysical rates require an external velocity distribution and a 3D cross-section.",
        "",
        "3. **Two discrete models.** Quinn et al. and Ortiz et al. are treated as separate "
        "observational models and should not be averaged without a hierarchical justification.",
        "",
        "4. **Sample extrema.** The highest-gain trajectories in any finite sample are "
        "dominated by extreme-value instability and are not converged optima.",
        "",
        "5. **Tail gate.** If tail_check_passed is False, the width estimate is a lower bound "
        "and b_max should be increased.",
        "",
        "6. **Time-limited trajectories.** Trajectories that hit max_time_sec are counted "
        "as non-events, reducing width estimates.",
        "",
        "---",
        "",
        "## Run Reference",
        "",
        "| Field | Value |",
        "|---|---|",
        f"| Package version | {manifest.get('package_version', '—')} |",
        f"| Schema version | {manifest.get('schema_version', '—')} |",
        f"| Science model | {manifest.get('science_model', '—')} |",
        f"| Git commit | `{manifest.get('git_commit', '—')[:12]}` |",
        "",
    ]

    report = "\n".join(lines)
    (output_path / "REPORT.md").write_text(report, encoding="utf-8")
    return report
