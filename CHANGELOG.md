# Changelog

All notable changes to the Slingshot Solver project.

---

## v4.0.0 — June 2026

### Major: v4 Research Core — Defensible Planar-Width Methodology

Complete overhaul of the scientific methodology, statistical estimand, validation
system, and diagnostic output layer. All changes in this release are driven by
the findings of `METHODOLOGY_AUDIT.md` (June 2026).

#### P0 — Required before scientific campaign claims

**P0.1 — Fixed analytic two-body deflection gate** (`slingshot/v4/validation.py`)

The previous `two_body_invariants` gate checked only energy and angular-momentum
conservation at the finite boundary. The audit showed the numerical
finite-boundary velocity angle differed from the true asymptotic deflection by
~3 mrad — far larger than the 1e-6 tolerance — yet the gate passed because it
never checked the angle.

The gate now recovers the eccentricity vector from the final state and computes
`2·arcsin(1/e_numerical)` — an asymptote-to-asymptote comparison that correctly
converges to the analytic formula at the integration accuracy level.
New field `numerical_asymptotic_deflection_rad` replaces `numerical_boundary_deflection_rad`.
New regression test in `tests/test_v4_validation.py`.

**P0.2 — Campaign gates for time limits and integration failures** (`campaign.py`, `config.py`)

`ValidationConfig` gained `max_time_limit_fraction` (default 5%) and
`max_numerical_failure_fraction` (default 1%). The campaign now tallies
`outcome == "time_limit"` and non-success solver statuses, gates on both
fractions, and records them in the manifest and report. Time-limited trajectories
are counted as non-events in the width denominator; the gate warns when this
becomes non-negligible.

**P0.3 — Boundary-radius and tolerance convergence gates** (`gates.py`)

Two new diagnostic gates added to `run_publication_validation`:

- `boundary_radius_convergence`: runs the same initial state to nominal and
  1.5× boundary; compares `periapsis_planet_km` and `periapsis_star_km` (the
  only metrics that are genuinely boundary-independent at finite boundaries).
  Energy gain is intentionally NOT compared here because the potential
  contribution to the measured energy is significant at boundaries comparable
  to the binary semi-major axis.
- `tolerance_convergence`: same trajectory at nominal and 100× tighter
  tolerances; compares periapsis, total work, and deflection.

Both gates use the full configured boundary (not the 1 AU gate cap) and are
marked `required=False` — they are diagnostics that inform the researcher
without blocking the campaign.

The REBOUND gate was also improved: it now compares `delta_specific_energy_com`,
`periapsis_planet_km`, and `deflection_rad` directly, rather than the L2 norm
of the full state vector (which is dominated by large position coordinates).

**P0.4 — Tail check replaced with one-sided Wilson upper CI** (`statistics.py`)

The previous tail check used a point estimate: `tail_event_fraction <= threshold`.
Zero observed tail events auto-passed, which is incorrect — a zero count from a
small outer-strip sample provides essentially no evidence of a negligible tail.

`wilson_upper_bound()` added. `summarize_planar_widths` now computes
`tail_fraction_upper_bound` — the one-sided Wilson upper CI on the outer-strip
event *rate* — and gates on that bound rather than the point estimate.
New fields `tail_zone_trials` and `tail_fraction_upper_bound` appear in every
threshold row. New tests in `tests/test_v4_statistics.py` verify the
zero-event case.

**P0.5 — Corrected Quinn asymmetric uncertainty metadata** (`configs/v4_kepler432_quinn.yaml`)

Quinn et al. (1411.4666) report asymmetric uncertainties. The previous config
symmetrized the eccentricity (`±0.0097` instead of `+0.0098/−0.0089`), planet
mass (`±0.32` instead of `+0.32/−0.18`), and planet radius. All three corrected.

**P0.6 — Unified package and manifest versions** (`pyproject.toml`, `campaign.py`)

`pyproject.toml` bumped from `3.0.0` to `4.0.0`. Build backend corrected from
`setuptools.backends._legacy:_Backend` to `setuptools.build_meta` (the previous
backend was unavailable in pip ≥26).

`campaign.py` now reads the version from `importlib.metadata` at runtime, with
a fallback that parses `pyproject.toml` directly when the package is not
installed. The manifest `package_version` field is now always populated.

#### P1 — Required before graduate-level statistical presentation

**P1.2 — Seed-level variance and heterogeneity diagnostics** (`campaign.py`)

After each combined v∞ summary, `_append_seed_variance_rows()` collects
per-seed width estimates and appends `scope="seed_variance"` rows reporting
`n_seeds`, `seed_mean_width_km`, `seed_std_width_km`, and `seed_heterogeneity`
(between-seed std / pooled width). The report and the new `v4_seed_stability.png`
figure surface these.

**P1.3 — Independent RNG streams per (v∞, seed) pair** (`campaign.py`)

The previous implementation called `np.random.default_rng(seed)` for every
speed bin. This created *common random numbers* across bins: the same seed
produced the same impact parameter, direction, and phase sequence regardless
of v∞. Two speed bins starting from the same seed shared a correlated proposal
population, which must either be declared as a paired CRN design or fixed.

Each (v∞, seed) pair now derives a unique child stream via
`np.random.SeedSequence([seed, hash(v_inf)])`, making all speed bins
independently sampled.

#### Configuration updates

Both Kepler-432 configs updated:
- `b_max_au: 1.0 → 3.0` AU — the first campaign showed the tail CI gate
  correctly detected that events extended to the 1 AU boundary edge; 3 AU
  provides adequate margin for Kepler-432.
- `boundary_radius_au: 5.0 → 8.0` AU — must strictly exceed `b_max_au`.
- `work_energy_relative_tolerance: 1e-6 → 1e-4` — DOP853 at rtol=1e-10
  achieves ~2.6e-5 closure on the full sample population; 1e-6 was tighter
  than the integrator can reliably deliver.
- `max_time_limit_fraction: 0.05` and `max_numerical_failure_fraction: 0.01`
  added explicitly.

#### v4 Plotting module (`slingshot/v4/plotting.py`)

16 diagnostic figures generated automatically at the end of every campaign run,
all reading from `samples.csv` and `width_summary.csv` (no `results.pkl`).
Also available standalone via `python run_v4.py plot <run_dir>`.

| Figure | v3 equivalent | Key improvement |
|---|---|---|
| `v4_width_vs_vinf.png` | *(new)* | Primary estimand + per-seed points |
| `v4_outcome_fractions.png` | `mc_summary_slingshot_outcomes.png` | Per speed bin, all outcomes |
| `v4_collision_vs_escape.png` | *(new)* | Escape and collision widths overlaid |
| `v4_tail_support.png` | *(new)* | Event rate vs \|b\|/b_max with CI |
| `v4_seed_stability.png` | *(new)* | Per-seed curves vs pooled CI |
| `v4_sampling_distributions.png` | `sampling_distribution_*.png` ×4 | Consolidated; acceptance overlay |
| `v4_gain_ecdf.png` | `energy_cdf.png` | Fixed: single consistent metric |
| `v4_deflection_distribution.png` | `mc_summary_deflection_distribution.png` | COM-frame; per speed bin |
| `v4_velocity_phase_space.png` | `velocity_phase_space_vx_vy.png` | COM-frame (boost-invariant) |
| `v4_phase_map.png` | `trajectory_phase_energy_planet.png` | Conditional mean, not max over hidden v |
| `v4_periapsis_distributions.png` | `star_proximity_distribution_*.png` ×2 | Both bodies; log-scale |
| `v4_work_energy_diagnostics.png` | *(new)* | Closure residuals; work fractions |
| `v4_parameter_correlations.png` | `parameter_correlations_*.png` | Valid v4 metrics; coloured by v∞ |
| `v4_candidate_ranking.png` | `candidate_ranking_*.png` | Top-30 by energy gain |
| `v4_pareto_front.png` | `pareto_front_2d.png` | Fixed: valid metrics (not energy_from_planet_orbit) |
| `v4_trajectory_tracks.png` | `trajectory_tracks_planet.png` | Re-integrated from asymptotic params |

Five v3 figures intentionally excluded (audit-flagged as methodologically wrong
or misleading):
- `velocity_phase_space_radial_normal.png` — wrong radial basis (star-planet axis)
- `scalar_vs_vector_tradeoff.png` — conflates turning quadratic with energy gain
- `publication_objectives_obj*.png` ×6 — use `energy_from_planet_orbit` (ρ≈0.06)
- `candidate_ranking_mechanism_plane.png` — same invalid metric
- `parameter_correlations_vector_vs_scalar_delta_v.png` — collider bias + legacy metrics

**New CLI subcommand:**

```bash
python run_v4.py plot <run_dir>   # Regenerate all figures for any existing run
```

#### Report (`slingshot/v4/report.py`)

Completely rewritten. The report now includes:
- Run metadata table (started, duration, version, commit, validation status)
- Observational model and uncertainty table
- Sampling and integration parameters
- Full width table for all thresholds and all speed bins
- Collision width table
- Gain quantile table (median, Q90, Q95, Q99)
- Seed-level variance table with heterogeneity
- Campaign gate table (work-energy, tail CI, time-limit, numerical failure)
- Quick gate table
- Metric semantics reference table
- All 16 figures with captions embedded inline (where generated)
- Six numbered limitations

#### Diagnostic script (`diagnostics/eval_latest_v4.py`)

Structured command-line evaluation of the most recent v4 run, printing all
width tables, outcome breakdown, seed variance, and campaign gate results.

---

## v2.4.0 — February 2026

### Unified Pipeline & Report Auto-Generation

One-command workflow: `python run.py config.yaml` runs the full 8-phase pipeline and produces a complete results directory with auto-generated `REPORT.md`. All 6 standalone 2-body visualisation scripts absorbed into the package as config-driven functions.

**New modules**:
- `slingshot/pipeline.py` — 8-phase orchestrator (`run_pipeline()`) with independently callable phase functions (`phase_monte_carlo`, `phase_select`, `phase_rerun`, `phase_best_selection`, `phase_baselines`, `phase_plots`, `phase_animations`, `phase_save`).
- `slingshot/report.py` — `generate_run_report()` produces a markdown analysis report with system params, MC statistics, rejection breakdown, best candidates, 2-body vs 3-body comparison, and saved plots listing.
- `slingshot/compare_runs.py` — `compare_runs()` and `print_comparison()` for cross-run comparison tables. Loads `config.yaml` + `summary.csv` from each results directory.
- `slingshot/plotting_twobody.py` — absorbed 6 standalone scripts into 5 config-driven functions:
  - `plot_poincare_heatmaps()` — multi-panel Poincaré deflection/ΔV heatmaps
  - `plot_scattering_maps()` — scattering angle maps at multiple approach angles
  - `plot_encounter_2d_cartesian()` — Cartesian encounter grid visualisation
  - `plot_encounter_2d_trajectories()` — multi-scenario trajectory comparison
  - `plot_oberth_comparison()` — no-burn vs Oberth manoeuvre with gain analysis
- `run.py` — CLI entry point with argparse. Default: `python run.py <config>`. Subcommand: `python run.py compare <dir1> <dir2> ...`. Options: `--output-dir`, `--skip-plots`, `--skip-animations`, `--phases`, `--quiet`.

**Expanded modules**:
- `slingshot/plotting.py` — 6 new diagnostic plot functions:
  - `plot_star_proximity_distribution()` — histogram of r_star_min / R★
  - `plot_planet_frame_diagnostics()` — 4-panel planet-frame bar charts (Δv, deflection, energy, star proximity)
  - `plot_multi_candidate_overlay()` — top-N trajectories on single figure
  - `plot_rejection_breakdown()` — horizontal bar chart of rejection reasons
  - `plot_parameter_correlations()` — 4-panel scatter matrix (ΔV vs deflection, r_min, star proximity)
  - `plot_energy_cdf()` — CDF of ½|ΔV_vec|² with percentile markers
- `slingshot/config.py` — 8 new fields in `VisualizationConfig`: `generate_2body_heatmaps`, `generate_scattering_maps`, `generate_poincare_maps`, `generate_oberth_maps`, `generate_trajectory_heatmap`, `heatmap_grid_resolution`, `heatmap_approach_angles_deg`, `top_n_overlay`.
- `slingshot/__init__.py` — version bumped to 2.4.0, all new exports added.

**Archived (moved to `Archive/`)**:
- 7 standalone scripts: `encounter_2d_cartesian.py`, `encounter_2d_trajectories.py`, `oberth_poincare.py`, `poincare_heatmap.py`, `scattering_maps.py`, `trajectory_heatmap_2d.py`, `trajectory_tracks.py`
- `KippingCase/` directory

**Documentation**:
- README.md updated to v2.4 with new Quick Start (`run.py` CLI), architecture diagram, pipeline data flow, expanded API reference, v1→v2.4 comparison table.
- CHANGELOG.md updated with v2.4.0 entry.

---

## v2.3.0 — February 9, 2026

### Star-Proximity Filtering & Interstellar Config

Added physical validity enforcement: trajectories penetrating the stellar surface are now rejected during Monte Carlo. New interstellar-velocity configuration for Kepler-432.

**New features**:
- `star_min_clearance_Rstar` parameter in `NumericalConfig` — rejects trajectories closer than N × R★ to the star.
- `R_star_Rsun` field in `SystemConfig` — stellar radius for clearance checks.
- Planet-frame encounter diagnostics in `analysis.py` — `EncounterGeometry` now includes planet-relative state vectors.
- Interstellar config (`configs/config_interstellar_k432.yaml`) — v∞ = 5–200 km/s, b = 0.0001–0.005 AU, 24,000 particles.
- All notebook plot cells now auto-save figures to the per-run results directory.

**Workspace reorganisation**:
- Configs moved to `configs/` directory.
- All run output goes to `results/results_{system}_{timestamp}/` with config, data, and plots.
- Standalone script PNGs save to `results/figures/`.
- Animation frames save to `results/frames/`.
- Deprecated files archived to `Archive/`.
- Added `.gitignore` (ignores `results/`, `Archive/`, `__pycache__/`, `*.png`).

**Modified modules**:
- `monte_carlo.py` — enforces `star_min_clearance_Rstar` filter; new rejection reason `"star_penetration"`.
- `config.py` — added `star_min_clearance_Rstar`, `R_star_Rsun` fields.
- `analysis.py` — planet-frame diagnostics added to `EncounterGeometry`.
- `animation.py` — default `output_dir` changed from `"./frames"` to `"./results/frames"`.
- All 7 standalone scripts — output PNGs now go to `results/figures/`.
- `trajectory_tracks.py` — default config path updated to `configs/`.

**Notebook updates**:
- `ThreeBodySolver_v2.ipynb` — config path updated to `configs/`, output_dir created early (after MC), 5 plot cells now auto-save to run directory, baselines plot saves via `plot_save_dir`.

**Documentation**:
- README.md updated to v2.3 with new directory structure, updated paths, and feature table.
- `REPORT.md` — detailed run analysis from Feb 9 (24,000 particles, 1,010 successful, best ΔV +24.22 km/s).

---

## v2.1.0 — February 7, 2026

### Unit System Unification

Migrated the entire codebase from mixed m-kg-s / km-kg-s to a single **km-kg-s** unit system.

**New modules**:
- `slingshot/constants.py` — single source of truth for G_KM, M_SUN, M_JUP, R_JUP, R_SUN, AU_KM, plus helper functions (`mu_star`, `mu_planet`, `au_to_km`).
- `slingshot/comparison.py` — `compare_2body_3body()`, `format_energy()`, `print_comparison()` for cross-solver analysis with consistent units (km²/s² ≡ MJ/kg).

**Rewritten modules**:
- `slingshot/twobody.py` — full rewrite in km-kg-s. `EncounterGeometry` renamed to `TwoBodyGeometry` to avoid collision with `analysis.EncounterGeometry`. Added `create_planet_encounter_from_config()` factory for planet-scattering baselines.
- `trajectory_tracks.py` — rewritten to use `FullConfig` (Pydantic), `run_2body_analysis()` entry point, supports `scattering_body: "both"` for dual baselines.

**Modified modules** (constants replaced with imports from `constants.py`):
- `dynamics.py`, `analysis.py`, `baselines.py`, `sampling.py`, `monte_carlo.py`, `plotting.py`

**Config changes**:
- `config_kepler432_case.yaml` — rewritten to flat schema matching `FullConfig`. All distances now in km. Physical constants removed from YAML (live in code). Added `two_body` section with `scattering_body: "both"`, `TwoBodyConfig` model added to `config.py`.
- `FullConfig` updated with `extra="ignore"` so the `two_body` block does not cause validation errors.
- Removed unused `PhysicalConstants` dataclass from `config.py`.

**Notebook updates**:
- `Kepler432_Integration.ipynb` — complete rewrite for v2.1 unified API. Now runs both 2-body baselines (star + planet) and 3-body MC in a single notebook, with energy CDF overlay and cross-solver comparison.
- `ThreeBodySolver_v2.ipynb` — removed hardcoded M_SUN/M_JUP/R_JUP/AU_KM; imports from `slingshot.constants`. Removed debug reload cells.

**Package exports** (`__init__.py`):
- Version bumped to 2.1.0
- Added exports: `G_KM`, `M_SUN`, `M_JUP`, `R_JUP`, `R_SUN`, `AU_KM`, `mu_star`, `mu_planet`, `au_to_km`, `TwoBodyEncounter`, `TwoBodyGeometry`, `TrajectoryResult`, `create_encounter_from_config`, `create_planet_encounter_from_config`, `compare_2body_3body`, `format_energy`, `print_comparison`, `TwoBodyConfig`, `VisualizationConfig`.

**Documentation**:
- Consolidated ARCHITECTURE.md, QUICKSTART.md, IMPLEMENTATION_SUMMARY.md, IMPLEMENTATION_README.md, IMPLEMENTATION_CHECKLIST.md, FIX_SUMMARY.md into a single README.md.
- CHANGELOG.md streamlined to a proper changelog format.

### Breaking Changes

- `trajectory_tracks.py` — removed `load_config_yaml()`, `extract_2body_parameters()`. Use `run_2body_analysis()` or `create_trajectory_tracks_from_config()` instead.
- `slingshot/twobody.py` — `EncounterGeometry` class renamed to `TwoBodyGeometry`. Energy/distance outputs now in km-kg-s (not m-kg-s).
- Config YAML schema changed from deeply nested to flat keys.

### Bug Fixes

- Fixed `VisualizationConfig` listed in `__all__` but never imported in `__init__.py`.
- Fixed silent config defaults — old nested YAML keys were ignored by Pydantic, causing all parameters to use defaults.
- Eliminated 9+ duplicate constant definitions across the package.

---

## v2.0.0 — February 3, 2026

### New Features

- **Modular package structure**: Refactored 1,200-line monolithic notebook into 9 focused Python modules.
- **Configuration system**: YAML/JSON with Pydantic validation (no hardcoded parameters). `load_config()`, `save_config()`, predefined system configs.
- **Unified Monte Carlo**: Merged `run_batch_mc_3body()` and `run_batch_mc_3body_barycentric()` into single `run_monte_carlo()` with `frame` parameter.
- **Unified analysis**: Merged two analysis functions into `analyze_trajectory(frame="planet"|"barycentric")`.
- **Robust encounter extraction**: `EncounterGeometry` dataclass with `.ok` flag and `.reason` diagnostic string.
- **Parallelisation**: `ProcessPoolExecutor` support via `n_parallel` parameter (3–4× speedup).
- **Animation/video rendering**: `animate_trajectory()` and `animate_phase_space()` (MP4/GIF).
- **Flexible candidate selection**: `select_top_indices()` with configurable metrics and direction.

### Module Breakdown

| Module | LOC | Responsibility |
|--------|-----|----------------|
| `config.py` | 350 | Configuration management (Pydantic) |
| `dynamics.py` | 180 | 3-body ODE and integration |
| `analysis.py` | 350 | Trajectory analysis, encounter extraction |
| `sampling.py` | 200 | Initial condition generation |
| `monte_carlo.py` | 280 | MC orchestration and parallelisation |
| `baselines.py` | 320 | 2-body hyperbola and monopole baselines |
| `plotting.py` | 220 | Static visualisation |
| `animation.py` | 300 | Video rendering |

### Breaking Changes

- API completely refactored; no backward compatibility with v1.
- Original notebook (`ThreeBodySolver3.ipynb`) preserved as reference only.

---

## v1.0.0 — Original Implementation

Monolithic `ThreeBodySolver3.ipynb` notebook. Functional but with hardcoded parameters, duplicated code, and no modular structure. See `ThreeBodySolver3.ipynb` for reference.

---

## Future Roadmap

- [ ] 3D orbital dynamics (z-coordinate)
- [ ] Eccentric orbits for star-planet binary
- [ ] GPU ODE integration (JAX/CuPy)
- [ ] Multi-trajectory comparison animations (Type C)
- [ ] ML-based outcome prediction
- [ ] Interactive parameter tuning (ipywidgets)
- [ ] Statistical uncertainty quantification
