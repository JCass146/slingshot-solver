# Hybrid v4 Methodology

## Purpose

This document defines the hybrid v4 methodology for Slingshot Solver. The goal
is to keep effective planar widths as the primary defensible scientific
claim while restoring the candidate-level intuition that made the exploratory
workflow useful.

The hybrid method separates three kinds of statements:

- population-level width estimates that can support quantitative claims;
- energy-gain ability summaries that describe the high-gain tail;
- finite-sample candidate examples that illustrate mechanisms but do not
  establish physical limits.

This document records the hybrid methodology implemented for v4 calibration
runs. It extends, but does not replace, the v4 planar-width estimand.

## Claim Hierarchy

### Primary Claim: Effective Planar Width

The primary claim remains:

```text
W(Delta epsilon / vc^2 > q | v_inf) = 2 * b_max * N_event / N
```

where:

- `v_inf` is the incoming asymptotic speed bin;
- `b_max` is the declared signed-impact-parameter sampling bound;
- `N` is the number of samples in the speed bin;
- `N_event` is the number of samples that escape and exceed threshold `q`;
- `Delta epsilon` is the COM-frame specific-energy change;
- `vc^2 = G(M_star + M_planet) / a` is the normalization scale.

This is a one-dimensional planar width under a declared proposal. It is not a
three-dimensional area cross-section, astrophysical event probability, or
occurrence rate.

The width estimate must be reported with Wilson confidence intervals and the
outer-tail gate. Production claim validation applies to declared ability
thresholds `q >= 0.01`; the `q = 0` row is retained as a diagnostic sign-width
for weak distant perturbations. If a claim-threshold tail gate fails, the width
is treated as a lower bound and `b_max` should be increased.

### Secondary Claim: Energy-Gain Ability

The energy-gain ability of an encounter should be quantified by the high-gain
distribution, not by a single best sample.

For each `v_inf`, report:

- widths above declared gain thresholds;
- median gain among escaped samples;
- high-gain quantiles such as Q90, Q95, and Q99;
- seed-level stability and heterogeneity;
- the strongest observed sample as an exploratory companion statistic.

This supports language such as:

```text
For v_inf = X, the effective planar width for gain above q is W, and the
escaped-sample Q99 gain is Y.
```

### Exploratory Claim: Top Candidates

Candidate diagnostics answer a different question:

```text
Which sampled trajectories best illustrate the high-gain tail?
```

They may identify `best_observed_gain`, ranked candidates, Pareto-front
members, and representative trajectories. These are finite-sample examples,
not converged optima or physical ceilings.

Allowed language:

```text
The strongest observed candidate achieved X.
```

Disallowed language without a separate convergence or optimization study:

```text
The maximum physical gain is X.
```

## Planar Width Estimand

The planar-width estimator is retained unchanged:

```text
width_km = 2 * b_max_km * event_count / trial_count
```

For each speed bin and threshold, the workflow should report:

- `events`;
- `trials`;
- `width_km` and `width_au`;
- Wilson low/high confidence bounds;
- `tail_fraction_upper_bound`;
- `tail_check_passed`;
- median and high-quantile escaped gains.

The event condition for energy-gain widths is:

```text
escaped == true
and energy_gain_dimensionless > threshold
```

Numerical failures and time-limited trajectories must be reported separately.
If their fractions exceed configured tolerances, the campaign is not a clean
publication-quality run.

## Energy-Gain Ability

The preferred ability summary is threshold-and-quantile based:

```text
ability(v_inf) =
  {
    W(gain > q_0),
    W(gain > q_1),
    ...,
    Q90(gain | escaped),
    Q95(gain | escaped),
    Q99(gain | escaped),
    best_observed_gain
  }
```

where:

```text
gain = energy_gain_dimensionless
     = delta_specific_energy_com / vc^2
```

This keeps the statistically meaningful population estimate central while
still letting reports discuss the high-gain edge of the sampled distribution.

The sample maximum should always be labeled as:

- `best_observed_gain`;
- `sample_max_gain`;
- or `strongest_observed_candidate`.

It should not be called:

- `maximum_gain`;
- `physical_limit`;
- `energy_ceiling`;
- or `converged optimum`.

## Valid v4 Metrics

The hybrid methodology uses only current v4 scientific metrics for claims
about energy gain:

| Metric | Role |
|---|---|
| `delta_specific_energy_com` | Primary COM-frame specific-energy change |
| `energy_gain_dimensionless` | Primary gain variable for widths, quantiles, and ranking |
| `delta_v_inf` | Asymptotic speed change when both endpoints are unbound |
| `deflection_rad` / `deflection_deg` | Boost-invariant turning angle diagnostic |
| `periapsis_planet_km` | Encounter-depth and risk diagnostic |
| `periapsis_star_km` | Physical-validity and collision-risk diagnostic |
| `work_star`, `work_planet`, `work_sum` | Work-energy mechanism and closure diagnostics |
| `work_energy_closure_relative` | Per-trajectory integration quality |
| `turning_quadratic` | Turning diagnostic only |

The following legacy quantities must not be used as scientific energy-gain
estimators:

- `energy_from_planet_orbit`;
- `bary_delta_v_pct`;
- scalar speed-change maxima;
- `0.5 * |Delta V|^2` as an energy-gain proxy.

`turning_quadratic` may be plotted beside true gain because it helps explain
geometry, but captions must state that it is not energy gain.

## Top-Candidate Diagnostics

The implementation generates:

```text
top_candidates.csv
```

from `samples.csv`.

Eligible candidates:

- `escaped == true`;
- `solver_success == true`;
- finite `energy_gain_dimensionless`;
- finite `delta_specific_energy_com`;
- finite `work_energy_closure_relative`.

Default ranking:

1. descending `energy_gain_dimensionless`;
2. descending `delta_specific_energy_com`;
3. ascending `work_energy_closure_relative`;
4. stable sample identity, such as `v_inf_kms`, `seed`, and `sample_index`.

Recommended columns:

- `rank`;
- `seed`;
- `sample_index`;
- `v_inf_kms`;
- `energy_gain_dimensionless`;
- `delta_specific_energy_com`;
- `delta_v_inf`;
- `turning_quadratic`;
- `deflection_deg`;
- `periapsis_planet_km`;
- `periapsis_star_km`;
- `impact_parameter_km`;
- `incoming_direction_rad`;
- `binary_mean_anomaly_rad`;
- `work_star`;
- `work_planet`;
- `work_sum`;
- `work_energy_closure_relative`;
- `solver_nfev`;
- `integration_time_sec`.

The candidate table is exploratory. It must not alter `width_summary.csv` or
replace population-level width estimates.

## Best-Candidate Example Plots

The implementation writes:

```text
best_candidate.png
```

This figure should plot the rank-1 candidate trajectory reconstructed from
stored asymptotic proposal parameters. The plot is an example visualization,
not an estimator.

Recommended figure content:

- trajectory in the binary barycentric frame;
- planet and star markers or relative locations where useful;
- inbound and outbound directions;
- title or caption with `v_inf`, gain, deflection, and periapsis;
- visible warning in the report caption that this is a finite-sample example.

Re-integration policy:

- Re-integrate selected candidates when regenerating `top_candidates.csv` so
  trajectory diagnostics use the current metric code.
- Store refreshed selected-candidate diagnostics in `top_candidates.csv`; do not
  overwrite population-level `samples.csv` unless a full repair/replay pass is run.
- If a re-integrated plot fails, keep the statistical run valid and report the
  plotting failure separately.

The existing `trajectory_tracks.png` should remain a top-N illustration in the
binary barycentric frame, colored by `energy_gain_dimensionless`.

## Pareto and Tradeoff Views

Pareto plots are useful for showing tradeoffs, but they are not primary
population statistics.

Recommended v4 Pareto objectives:

- maximize `energy_gain_dimensionless`;
- maximize `abs(deflection_deg)`;
- minimize `periapsis_planet_km` only when explicitly framed as encounter
  depth or risk.

A safer default two-panel view is:

- gain versus planet periapsis;
- gain versus absolute deflection.

Captions should explain that Pareto-front members are selected examples from
the sampled population, not a proof of a global frontier.

## Report Wording Rules

Use these phrases:

- "effective planar width";
- "width for gain above threshold";
- "escaped-sample high-gain quantile";
- "strongest observed candidate";
- "best observed sample";
- "finite-sample candidate example";
- "illustrative trajectory";
- "turning diagnostic".

Avoid these phrases unless backed by a dedicated optimization or convergence
study:

- "maximum physical gain";
- "true maximum";
- "energy ceiling";
- "optimal trajectory";
- "global Pareto frontier";
- "event probability" without the proposal qualification;
- "3D cross-section";
- "astrophysical rate".

Recommended wording:

```text
The strongest observed candidate in this run achieved gain X. This is an
exploratory finite-sample extremum. The defensible population result is the
width above declared thresholds and the high-gain quantile table.
```

## Required Artifacts

The hybrid methodology should preserve all current v4 artifacts:

- `config.yaml`;
- `samples.csv`;
- `width_summary.csv`;
- `manifest.json`;
- `REPORT.md`;
- `width_vs_vinf.png`;
- `gain_ecdf.png`;
- `tail_support.png`;
- `seed_stability.png`;
- `work_energy_diagnostics.png`;
- `candidate_ranking.png`;
- standalone `candidate_ranking_*.png` panels;
- `pareto_front.png`;
- `trajectory_tracks.png`.

Candidate artifacts:

- `top_candidates.csv`;
- `best_candidate.png`;
- optional `candidate_tradeoffs.png` if the existing ranking plot becomes too
  crowded.

The manifest should list generated candidate artifacts and the report should
embed them when present.

## Validation and Acceptance Criteria

A run may make the primary planar-width claim only when:

- quick validation passes;
- work-energy closure passes the configured tolerance;
- tail gates pass for the reported thresholds or failed tails are clearly
  labeled as lower bounds;
- time-limit and numerical-failure fractions pass configured campaign gates;
- the frozen config, seeds, package version, and Git commit are recorded.

Candidate diagnostics may be shown when:

- candidates are drawn from valid escaped samples;
- per-candidate closure errors are shown or summarized;
- collision and periapsis information are included;
- captions state that examples are illustrative.

Candidate diagnostics should be suppressed or marked unavailable when:

- no escaped finite-gain samples exist;
- candidate rows have non-finite required metrics;
- plotting re-integration fails for all selected examples.

## Non-Claims and Failure Modes

The hybrid methodology still does not claim:

- a 3D isotropic cross-section;
- an astrophysical event rate;
- a physical maximum energy gain;
- a globally optimized trajectory;
- a mechanism attribution from legacy heuristic energy metrics.

Known failure modes:

- Tail-gate failure at `q >= 0.01` means `b_max` is too small for the claimed threshold.
- Tail-gate failure at `q = 0` is reported as diagnostic, not as production validation failure.
- Sparse event counts make high thresholds noisy even when point estimates are
  nonzero.
- Sample maxima are unstable under larger sample sizes and new seeds.
- Candidate plots can be visually compelling while remaining statistically
  non-representative.
- Turning can be large even when true energy gain is small.

The hybrid methodology is therefore intentionally two-voiced: widths and
quantiles carry the statistical claim; candidates and trajectories restore
physical intuition.
