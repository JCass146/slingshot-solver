# Hybrid v4 Parameter Exploration Plan

**Status:** Revised for the hybrid v4 methodology.

**Research question:** Under a declared planar asymptotic proposal, when do
binary encounters provide measurable energy-gain ability beyond appropriate
single-body or single-scatterer baselines?

**Core answer format:** The defensible result is an effective planar width,

```text
W(Delta epsilon / v_c^2 > q | v_inf),
```

reported with Wilson intervals, seed diagnostics, and tail-gate status. Top
candidates, best-observed gains, Pareto fronts, and barycentric trajectory
plots are retained as exploratory calibration tools, not as maximum-gain claims.

This document is the canonical root-level plan for parameter exploration. New
campaign reports may live under `results/` or analysis-specific directories,
but new root Markdown files should be avoided unless they become canonical
project documents.

---

## 0. Methodology Guardrails

Before changing parameters, keep these rules fixed across all comparisons:

- Primary claim thresholds start at `q >= 0.01`.
- The `q = 0` row is diagnostic only. It is useful for weak sign-width support,
  but it is not a production ability claim and should not decide grid winners.
- The scientific gain variable is `energy_gain_dimensionless =
  delta_specific_energy_com / v_c^2`.
- `turning_quadratic` is a turning diagnostic only.
- Legacy scalar speed changes, `bary_delta_v_pct`, and
  `energy_from_planet_orbit` must not be used as gain estimators.
- Candidate maxima are finite-sample observations. Use
  `best_observed_gain`, not `maximum physical gain`.
- Candidate trajectories should be plotted in the barycentric frame and labeled
  illustrative.
- Reports may coalesce multi-panel candidate views, but standalone plots should
  also be emitted for documentation and inspection.

The main scouting metrics are:

| Metric | Role |
|---|---|
| `W(gain > 0.01 | v_inf)` | Primary low-tail ability width |
| `W(gain > 0.03 | v_inf)` | Primary stronger-tail ability width |
| `W(gain > 0.10 | v_inf)` | Secondary rare/high-gain width |
| `Q95`, `Q99` among escaped samples | High-gain ability summaries |
| seed heterogeneity | Stability and calibration diagnostic |
| `best_observed_gain` | Exploratory strongest observed sample |
| star/planet periapsis cuts | Physical-risk and mission-plausibility diagnostic |

Median gain is allowed as a descriptive statistic, but it should not drive
parameter selection unless the distribution is clearly shifted away from zero.
For these slingshot distributions, the tail usually matters more than the
median.

---

## Phase 0: Baseline Design and Sanity Runs

**Objective:** Define the comparison baseline before claiming that a binary
system improves slingshot ability.

Do not assume that `planet_mass_mjup: 0` is automatically a valid v4 baseline.
A degenerate binary can break element recovery, collision logic, or the meaning
of `v_c^2`. Use explicit baseline support or a clearly marked diagnostic config.

### Baseline A: Isolated-Star Sanity Check

**Purpose:** Confirm that a stationary single-body potential does not create
spurious COM-frame energy gain beyond numerical and finite-boundary residuals.

Expected behavior:

- `delta_specific_energy_com` should center near zero.
- Work-energy closure should pass.
- Widths above `q >= 0.01` should be consistent with zero after uncertainty.
- Any nonzero result is a numerical or boundary diagnostic, not binary ability.

This is a sanity baseline, not necessarily the physical comparator for a moving
gravity-assist encounter.

### Baseline B: Single Moving Scatterer / Analytic Gravity Assist

**Purpose:** Compare the binary result with the best single-scatterer intuition:
a moving massive body can exchange energy with a probe, but lacks the full
star-companion time-dependent binary field.

Acceptable implementations:

- an analytic two-body gravity-assist envelope around a moving companion;
- a dedicated single-moving-scatterer numerical mode;
- a carefully documented approximation that uses the same COM-frame gain
  variable and the same proposal definition.

This is the better baseline for the question, "Does the binary add ability
beyond a single moving encounter?"

### Baseline C: Current Kepler-432 Binary Models

Use Quinn and Ortiz as separate observational models. Do not combine them into
a single posterior without a separate hierarchical model.

Current priority:

1. Keep the recent Ortiz full calibration as the current binary reference.
2. Run or refresh a matched Quinn campaign with the same methodology.
3. Compare each against the selected baseline(s), not against legacy v3 runs.

### Phase 0 Deliverables

Keep deliverables out of the root unless they become canonical docs:

- `configs/` baseline config(s), if the code supports them cleanly;
- run directories under `results/`;
- baseline comparison tables in the relevant run report or an analysis output
  directory;
- README/CHANGELOG updates only when methodology or user-facing workflow
  changes.

Acceptance criteria:

- Baseline validation gates are meaningful for the selected mode.
- Numerical/time-limit failure fractions pass configured gates.
- Work-energy closure passes.
- The report clearly states what the baseline does and does not test.

---

## Phase 1: Kepler-432 Comparative Results

**Objective:** Quantify Kepler-432 binary ability relative to the baseline(s)
using the hybrid v4 metrics.

Run set:

- Kepler-432 b, Ortiz et al. model;
- Kepler-432 b, Quinn et al. model;
- selected baseline A and/or B from Phase 0.

Analysis for each `v_inf` and claim threshold:

- width ratio `R = W_binary / W_baseline` for `q = 0.01` and `q = 0.03`;
- optional ratio for `q = 0.10` when event counts are sufficient;
- Wilson intervals on widths and propagated/bootstrap uncertainty on ratios;
- `Q95` and `Q99` escaped-gain comparison;
- seed heterogeneity comparison;
- top-candidate tables and barycentric best-candidate plots for intuition.

Report `q = 0` separately as a diagnostic sign-width. Do not use it for the
binary-advantage claim.

Add optional physical-risk views:

```text
periapsis_star > 3 R_star
periapsis_star > 5 R_star
periapsis_star > 10 R_star
periapsis_planet > 1 R_planet
```

These are not replacements for the primary width, but they help calibrate which
high-gain candidates are dynamically interesting versus physically hazardous.

Recommended wording:

```text
For the declared planar proposal, the Ortiz model has an effective width of W
for gain above q = 0.03 at v_inf = X. The strongest observed candidate achieved
Y, but this is an exploratory sample extremum rather than a physical limit.
```

---

## Phase 2: Scaling Hypotheses Before the Grid

**Objective:** Write down falsifiable expectations before inspecting a broad
parameter grid.

Keep the hypothesis register in this file unless it grows enough to justify a
new canonical document. Avoid creating one-off root Markdown files.

### Dimensionless Parameters

| Parameter | Symbol | Notes |
|---|---|---|
| Mass ratio | `mu = M_companion / (M_star + M_companion)` | For large `mu`, this is a binary companion, not a planet. |
| Eccentricity | `e` | Controls orbital-speed variation and phase-dependent encounter geometry. |
| Separation scale | `a` | Sets `v_c^2 = G M_tot / a`; compare using `v_inf / v_c`. |
| Speed ratio | `s = v_inf / v_c` | Preferred speed coordinate for scaling. |
| Encounter width scale | `b_max / a` | Tail gate should determine whether the sampled strip is wide enough. |

### Initial Hypotheses

**H1: Width increases with companion mass ratio until collision/risk effects
or phase-space narrowing dominate.**

Prediction: `W(gain > 0.03)` should generally increase with `mu`, but high `mu`
may change the physical interpretation from star-planet to stellar binary.
Candidate plots must label this honestly.

**H2: Eccentricity changes the high-gain tail, but circular binaries are not a
zero-gain limit.**

A circular binary is still time-dependent in the inertial barycentric frame.
Do not claim `e -> 0` implies zero gain. The safer prediction is that
eccentricity changes phase concentration, close-approach geometry, and the
relative duration of high-speed orbital phases.

**H3: Semi-major-axis trends must be interpreted through `v_inf / v_c`.**

Changing `a` changes both the physical encounter scale and the normalization
`v_c^2`. A fixed physical `v_inf` grid mixes those effects. For scaling runs,
prefer `v_inf / v_c` bins, then translate to km/s for selected physical cases.

**H4: The best raw candidates may be star-grazing.**

Prediction: imposing star-clearance cuts will reduce top observed gains and may
change the apparent best parameter region. Report both raw and clearance-filtered
candidate rankings.

---

## Phase 3A: Small Pilot Grid

**Objective:** Validate the grid machinery and identify whether the proposed
parameter ranges are informative before paying for the full 27-configuration
sweep.

Recommended pilot grid:

| Parameter | Values |
|---|---|
| `mu` | Kepler-like low value, `0.10` |
| `e` | `0.20`, `0.80` |
| `a_au` | `0.303`, `0.80` |

This is `2 x 2 x 2 = 8` configurations.

Sampling:

- speeds: preferably `v_inf / v_c = [0.25, 0.5, 1.0, 2.0]`, converted into
  config-level `v_inf_kms`;
- optional physical sanity speeds: `[10, 20, 40, 80]` km/s for Kepler-432-like
  cases;
- `samples_per_bin: 500`;
- `seeds: [42, 43]`;
- thresholds: `[0, 0.01, 0.03, 0.1]`;
- start with `b_max / a` near the successful Ortiz calibration scale and let
  the tail gate decide whether to expand it;
- keep `boundary_radius > b_max`, preferably scaled with `a`.

Trajectory count estimate:

```text
8 configs x 4 speeds x 500 samples/bin x 2 seeds = 32,000 trajectories
```

Decision metrics:

- tail-gated `W(gain > 0.01)`;
- tail-gated `W(gain > 0.03)`;
- event counts for `gain > 0.10`;
- `Q99(gain | escaped)`;
- time-limit and numerical-failure fractions;
- seed heterogeneity;
- raw and clearance-filtered top candidates.

Acceptance criteria before Phase 3B:

- generated configs validate cleanly;
- reports include the hybrid candidate plots and standalone ranking panels;
- `q >= 0.01` tail gates pass or clearly identify which configs need larger
  `b_max`;
- no parameter cell is dominated by time limits or numerical failures;
- the pilot reveals enough structure to justify the full grid.

---

## Phase 3B: Full Coarse Grid

**Objective:** Explore a coarse response surface after the pilot confirms that
the grid design is meaningful.

Recommended grid:

| Parameter | Values |
|---|---|
| `mu` | `0.01`, `0.10`, `0.30` |
| `e` | `0.20`, `0.50`, `0.80` |
| `a_au` | `0.20`, `0.50`, `1.00` |

Total configurations:

```text
27 configs x 4 speeds x 500 samples/bin x 2 seeds = 108,000 trajectories
```

Important correction: this is the total trajectory count for the full grid, not
54k per configuration.

Use the same thresholds as the pilot:

```text
[0, 0.01, 0.03, 0.1]
```

Keep `q = 0.3` and `q = 1.0` only as optional rare-event diagnostics unless the
pilot shows those thresholds have enough support.

For `mu = 0.30`, call the object a companion or binary component rather than a
planet. Use physically plausible radii and collision labels for the chosen
interpretation.

Grid output should be tabular and reproducible:

| config_id | mu | e | a_au | v_inf_kms | v_inf_over_vc | W_q0p01 | W_q0p03 | Q99 | tail_pass | time_limit_frac | seed_heterogeneity |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|

Plotting priorities:

- heatmaps of `W(gain > 0.01)` and `W(gain > 0.03)`;
- `Q99` heatmaps;
- event-count/support heatmaps for `q = 0.10`;
- seed-heterogeneity heatmaps;
- clearance-filtered top-candidate comparison plots;
- representative barycentric best-candidate examples from promising cells.

---

## Phase 4: Focused High-Precision Runs

**Objective:** Re-run only the promising region with production precision.

Selection rule:

- choose 1-3 regions with high `W(gain > 0.03)`, acceptable seed stability, and
  physically interpretable top candidates;
- include at least one control region that was predicted to be weak;
- avoid selecting solely on `best_observed_gain`.

Production settings:

- `samples_per_bin: 2000` or higher if pilot event counts require it;
- `seeds: [42, 43, 44, 45, 46]`;
- full speed grid in either physical `v_inf_kms` or dimensionless
  `v_inf / v_c`, declared before the run;
- thresholds: `[0, 0.01, 0.03, 0.1, 0.3, 1.0]`;
- full validation and report regeneration.

Deliverables:

- full run directories with `REPORT.md`;
- comparison tables for claim thresholds;
- candidate diagnostics for calibration;
- final decision on whether the parameter region is worth a 3D extension.

---

## Phase 5: 3D Extension Placeholder

The current project estimates a planar width. A 3D extension would require a
new proposal, new geometry, and a true area cross-section or solid-angle
marginalization. Do not describe planar widths as 3D cross-sections.

---

## Required Artifacts Per Run

Each production or calibration run should include:

- `config.yaml`;
- `samples.csv`;
- `width_summary.csv`;
- `top_candidates.csv`;
- `manifest.json`;
- `REPORT.md`;
- `width_vs_vinf.png`;
- `gain_ecdf.png`;
- `tail_support.png`;
- `seed_stability.png`;
- `work_energy_diagnostics.png`;
- `best_candidate.png`;
- `candidate_ranking.png`;
- standalone `candidate_ranking_*.png` files;
- `pareto_front.png`;
- `trajectory_tracks.png`.

Reports may coalesce figures for readability, but standalone plots should stay
available for documentation and simulator calibration.

---

## Failure Modes and Stopping Rules

### Sample-Maximum Inflation

Symptom: `best_observed_gain` increases substantially when sample size or seed
count increases.

Response: report the new observation as a sample maximum, not a physical limit.
Use quantiles and threshold widths for ability claims.

### Sparse High Thresholds

Symptom: `q = 0.3` or `q = 1.0` has very few events.

Response: keep the rows, report the wide intervals, and treat them as
exploratory until sample support improves.

### Tail-Gate Failure

Symptom: `q >= 0.01` tail gate fails.

Response: increase `b_max`, rerun, or label the width as a lower bound. Do not
remove the tail gate to make the result look cleaner.

### Star-Grazing Dominance

Symptom: top candidates all pass within a few stellar radii.

Response: add clearance-filtered rankings and discuss mission plausibility
separately from dynamical possibility.

### Methodology Breach

Stop and revise the claim if a report begins to:

- call a finite-sample maximum a physical ceiling;
- use `q = 0` as the main ability claim;
- compare legacy v3 gains against v4 COM energy gain;
- pool Quinn and Ortiz as one posterior without a model;
- describe planar width as an astrophysical rate or 3D cross-section.

---

## Publication-Oriented Result Language

Allowed:

```text
The strongest observed candidate achieved X in this finite sample.
```

Allowed:

```text
The effective planar width for gain above q = 0.03 is W with a pointwise
Wilson interval of [lo, hi].
```

Allowed:

```text
The largest measured width in the explored grid occurred at these parameters.
```

Not allowed without a separate optimization or convergence study:

```text
The maximum physical gain is X.
```

Not allowed for the planar method:

```text
The astrophysical event probability is W.
```

---

## Phase Checklist

Before Phase 1:

- baseline mode is defined and validated;
- baseline report explains what comparison it supports;
- Ortiz and Quinn runs use comparable settings.

Before Phase 2:

- binary/baseline comparison has claim-threshold widths;
- `q = 0` is clearly diagnostic;
- candidate examples are illustrative only.

Before Phase 3A:

- scaling hypotheses in this file are reviewed;
- config-generation script computes derived masses, speeds, and scales;
- `b_max / a` and `boundary_radius / a` policy is declared.

Before Phase 3B:

- pilot configs pass validation or reveal fixable settings;
- tail gates pass for claim thresholds or rerun settings are chosen;
- grid output schema is stable.

Before Phase 4:

- selected regions are based on widths, quantiles, and stability;
- high-gain candidates survive basic physical-risk review;
- sample-size needs are estimated from pilot event counts.

Before publication-facing claims:

- reports have current hybrid plots and candidate tables;
- seed heterogeneity and tail support are visible;
- no language claims a global optimum, physical maximum, 3D cross-section, or
  astrophysical rate.

---

## Canonical References

- [README.md](README.md) - user workflow, artifacts, validation, and docs map.
- [CHANGELOG.md](CHANGELOG.md) - implemented methodology and cleanup history.
- [HYBRID_V4_METHODOLOGY.md](HYBRID_V4_METHODOLOGY.md) - claim hierarchy and wording rules.
- [docs/slingshot_derivations.md](docs/slingshot_derivations.md) - mathematical foundations.
- [docs/archive/METHODOLOGY_AUDIT_2026-06-20.md](docs/archive/METHODOLOGY_AUDIT_2026-06-20.md) - historical audit only.