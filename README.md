# Slingshot Solver

Slingshot Solver is a Python research project for studying planar gravitational
scattering in a restricted three-body star-planet system. The repository
contains two versioned workflows:

- **v4 research core** — the current scientific workflow for estimating
  effective planar encounter widths as a function of asymptotic speed.
- **v3 legacy pipeline** — the earlier exploratory Monte Carlo, plotting,
  animation, and candidate-ranking workflow.

The v4 core corrects the orbital initialization, asymptotic encounter proposal,
energy definitions, event handling, statistical estimand, and validation
strategy. Historical v3 runs are readable but not scientifically comparable
with v4 results and are excluded from v4 aggregates.

**Units:** km, kg, and s throughout. Specific energies in km²/s² (= MJ/kg).
**Package version:** 4.0.0

---

## Scientific Scope

The primary v4 result is an **effective planar encounter width**:

$$W\!\left(\Delta\epsilon / v_c^2 > q \mid v_\infty\right) = 2b_{\max}\frac{N_{\mathrm{event}}}{N}$$

where:

- $v_\infty$ is the incoming asymptotic speed;
- $b_{\max}$ defines the sampled signed-impact interval;
- $\Delta\epsilon$ is the test particle's specific-energy change in the binary
  center-of-mass (COM) frame;
- $v_c^2 = G(M_\star + M_p)/a$ sets the dimensionless energy scale.

This quantity is a one-dimensional width under a declared planar sampling
proposal. It is **not** a three-dimensional area cross-section, astrophysical
event probability, or occurrence rate.

---

## What v4 Adds

- Eccentric Keplerian binary initialization from semi-major axis,
  eccentricity, mean anomaly, argument of periapsis, masses, and direction.
- Proper inbound hyperbolic states defined by $v_\infty$, signed impact
  parameter, incoming direction, and a finite boundary radius.
- Uniform sampling in signed impact parameter, incoming direction, and binary
  mean anomaly (mean anomaly = uniform observation time).
- Center-of-mass scientific metrics invariant under Galilean boosts.
- Root-event detection for stellar collision, planetary collision, outbound
  escape, and periapsis.
- Separate work integrals for the moving stellar and planetary potentials,
  with work-energy closure checks per trajectory.
- Wilson confidence intervals, collision widths, gain quantiles, and
  one-sided confidence-bound tail gate.
- DOP853 campaign integration; Radau cross-checks; optional REBOUND IAS15.
- Versioned configurations, manifests, observational provenance, and compact
  per-sample records.
- Separate Kepler-432 models for Quinn et al. and Ortiz et al., treated as
  discrete observational-model uncertainty.
- Per-(v∞, seed) independent RNG streams for reproducibility.
- Seed-level variance and heterogeneity diagnostics.
- 16 automatically generated diagnostic figures per campaign run.
- `python run_v4.py plot <run_dir>` to regenerate figures from any archived run.

---

## Quick Start

### Install

Python 3.9 or newer is required.

```bash
git clone <repository-url>
cd slingshot-solver
python -m pip install -e .
```

Install development tools for testing:

```bash
python -m pip install -e ".[dev]"
```

Optional independent-integrator validation:

```bash
python -m pip install rebound
```

### Validate a v4 Preset

Run the fast deterministic checks (binary elements, deflection gate,
Jacobi conservation, etc.):

```bash
python run_v4.py validate configs/v4_kepler432_quinn.yaml
```

Run the full publication-gate set (includes convergence diagnostics,
Galilean invariance, integrator agreement):

```bash
python validate_v4.py configs/v4_kepler432_quinn.yaml
```

### Run a Small Pilot

```bash
python run_v4.py run configs/v4_kepler432_ortiz.yaml \
  --samples-per-bin 100 \
  --seeds 42,43 \
  --output-dir results/v4_ortiz_pilot
```

### Run the Configured Campaign

```bash
python run_v4.py run configs/v4_kepler432_ortiz.yaml
python run_v4.py run configs/v4_kepler432_quinn.yaml
```

Each preset defaults to 7 speed bins × 5 seeds × 2,000 trajectories
= 70,000 trajectories per observational model. Runtime is approximately
2–3 hours on a modern laptop.

### Regenerate Figures for an Existing Run

```bash
python run_v4.py plot results/v4_kepler432_ortiz2015_<timestamp>
```

---

## v4 Presets

| Preset | Purpose |
|---|---|
| `configs/v4_dimensionless_reference.yaml` | Dimensionless Jupiter-analog reference for method development |
| `configs/v4_kepler432_quinn.yaml` | Kepler-432 b from Quinn et al. 2015 |
| `configs/v4_kepler432_ortiz.yaml` | Kepler-432 b from Ortiz et al. 2015 |

Quinn and Ortiz results should be compared side by side. They should not be
combined into a single posterior without a hierarchical model that justifies
doing so.

---

## v4 Configuration

Every v4 configuration declares `schema_version: 4`:

```yaml
schema_version: 4

system:
  name: Kepler-432 — Ortiz et al.
  star_mass_msun: 1.35
  star_radius_rsun: 4.16
  planet_mass_mjup: 5.84
  planet_radius_rjup: 1.102
  bulk_velocity_x_kms: 10.0   # Galilean invariance test
  bulk_velocity_y_kms: 0.0

orbit:
  model: keplerian
  semi_major_axis_au: 0.303
  eccentricity: 0.478
  mean_anomaly_rad: 0.0
  argument_periapsis_rad: 0.0
  prograde: true

asymptotic_sampling:
  v_inf_kms: [10, 20, 30, 40, 60, 80, 120]
  b_max_au: 3.0               # Must be large enough for tail gate to pass
  boundary_radius_au: 8.0     # Must exceed b_max
  samples_per_bin: 2000
  seeds: [42, 43, 44, 45, 46]

planar_width:
  dimensionless_energy_thresholds: [0, 0.01, 0.03, 0.1, 0.3, 1.0]
  confidence_level: 0.95
  tail_fraction: 0.10
  max_tail_event_fraction: 0.01

numerical:
  method: DOP853
  rtol: 1.0e-10
  atol: 1.0e-10
  softening_km: 0.0           # Newtonian; non-zero only for diagnostics

validation:
  work_energy_relative_tolerance: 1.0e-4
  max_time_limit_fraction: 0.05
  max_numerical_failure_fraction: 0.01
```

**Key configuration notes:**

- `b_max_au` must be large enough that the outer-10%-strip event rate is
  near zero. The tail CI gate enforces this; increase `b_max_au` if it fails.
- `boundary_radius_au` must strictly exceed `b_max_au`.
- `work_energy_relative_tolerance: 1e-4` is appropriate for DOP853 at
  `rtol=atol=1e-10`. The tolerance `1e-6` is tighter than the integrator
  can deliver across the full sample population.
- `softening_km: 0.0` is mandatory for scientific runs.

---

## Methodology

### 1. Binary State

The star and planet are initialized in barycentric coordinates from the
configured Keplerian elements. Scientific quantities are measured in the
binary center-of-mass frame.

### 2. Asymptotic Proposal

For each fixed $v_\infty$ bin, v4 samples:

- signed impact parameter uniformly over $[-b_{\max}, b_{\max}]$;
- incoming direction uniformly over $[0, 2\pi)$;
- binary mean anomaly uniformly over $[0, 2\pi)$.

Each (v_inf, seed) pair gets an independent child RNG stream derived via
`numpy.random.SeedSequence`, so speed bins do not share common random numbers.

The asymptotic energy and angular momentum are mapped analytically onto a
finite inbound boundary using the total-mass monopole.

### 3. Integration and Outcomes

Campaigns use SciPy DOP853 (default) or Radau. Root events classify:

- outbound escape (primary event);
- stellar collision;
- planetary collision;
- close approach and periapsis (non-terminal; tracked).

### 4. Scientific Metrics

The primary gain metric is:

```
delta_specific_energy_com = epsilon_out - epsilon_in  [km²/s²]
energy_gain_dimensionless  = delta_specific_energy_com / vc²
```

`turning_quadratic = 0.5 * ||v_out - v_in||²` is retained only as a turning
diagnostic. It is **not** interpreted as energy gain.

### 5. Statistical Summary

For every (speed, threshold) pair, v4 reports:

- effective planar width W with Wilson CI;
- event count and effective sample size;
- median escaped gain and gain quantiles (Q90, Q95, Q99);
- collision width;
- one-sided Wilson upper CI on the outer-10%-strip event rate (tail gate).

Between-seed variance and heterogeneity are reported in `width_summary.csv`
under `scope=seed_variance`.

---

## Run Artifacts

Each v4 campaign writes:

```text
results/v4_<case>_<timestamp>/
├── config.yaml              # Frozen schema-v4 configuration
├── samples.csv              # Proposal vars, outcomes, metrics, work, solver diagnostics
├── width_summary.csv        # Per-seed, combined, and seed-variance width rows
├── manifest.json            # Schema, version, commit, seeds, gates, provenance
├── REPORT.md                # Comprehensive report with all tables and figures
├── v4_width_vs_vinf.png     # Primary estimand: W(v∞) with CI and seed points
├── v4_outcome_fractions.png # Stacked outcome bars per speed bin
├── v4_collision_vs_escape.png
├── v4_tail_support.png
├── v4_seed_stability.png
├── v4_sampling_distributions.png
├── v4_gain_ecdf.png
├── v4_deflection_distribution.png
├── v4_velocity_phase_space.png
├── v4_phase_map.png         # Conditional mean gain over (b, binary phase)
├── v4_periapsis_distributions.png
├── v4_work_energy_diagnostics.png
├── v4_parameter_correlations.png
├── v4_candidate_ranking.png
├── v4_pareto_front.png      # Two Pareto fronts using valid v4 metrics
└── v4_trajectory_tracks.png # Top-10 re-integrated trajectories in planet frame
```

---

## Validation System

The v4 validation system checks:

| Gate | Type | Description |
|---|---|---|
| `binary_elements` | Required | Recovery of configured semi-major axis and eccentricity |
| `two_body_invariants` | Required | Asymptotic deflection, energy, and angular-momentum conservation |
| `newtonian_default` | Required | Zero softening enforced |
| `work_energy_closure` | Required | Moving-potential work matches energy change |
| `galilean_invariance` | Required | COM metrics unchanged under bulk-velocity boost |
| `dop853_radau_agreement` | Required | DOP853 and Radau give consistent results |
| `circular_jacobi_conservation` | Required | Jacobi constant conserved in circular case |
| `boundary_radius_convergence` | Diagnostic | Periapsis stable under 1.5× boundary extension |
| `tolerance_convergence` | Diagnostic | Periapsis and work stable under 100× tighter tolerances |
| `rebound_ias15_agreement` | Diagnostic | Optional REBOUND IAS15 cross-check |

Campaign-level gates:

| Gate | Description |
|---|---|
| Work-energy closure | Max relative error across all samples ≤ configured tolerance |
| Tail CI | One-sided Wilson upper bound on outer-strip event rate ≤ `max_tail_event_fraction` |
| Time-limit fraction | Fraction of time-limited trajectories ≤ `max_time_limit_fraction` |
| Numerical failure fraction | Fraction of failed integrations ≤ `max_numerical_failure_fraction` |

Run the test suite:

```bash
pytest -q
```

---

## Project Structure

```text
slingshot-solver/
├── slingshot/
│   ├── v4/
│   │   ├── config.py       # Schema-v4 Pydantic models and serialization
│   │   ├── dynamics.py     # Keplerian binary init and event-driven integration
│   │   ├── sampling.py     # Asymptotic hyperbolic proposal
│   │   ├── metrics.py      # COM metrics and moving-potential work accounting
│   │   ├── statistics.py   # Wilson intervals, planar widths, tail CI gate
│   │   ├── validation.py   # Fast deterministic checks (binary, two-body, Jacobi)
│   │   ├── gates.py        # Publication gates (work-energy, Galilean, convergence)
│   │   ├── campaign.py     # Campaign runner, artifacts, seed variance
│   │   ├── report.py       # Comprehensive report with all tables and figure refs
│   │   ├── plotting.py     # 16 diagnostic figures from CSV artifacts
│   │   ├── runs.py         # Version-aware run discovery
│   │   └── cli.py          # run / validate / plot subcommands
│   ├── core/               # Legacy v3 dynamics and two-body tools
│   ├── analysis/           # Legacy v3 Monte Carlo and diagnostics
│   └── output/             # Legacy v3 plots, reports, and animation
├── configs/
│   ├── v4_kepler432_ortiz.yaml
│   ├── v4_kepler432_quinn.yaml
│   └── v4_dimensionless_reference.yaml
├── diagnostics/
│   └── eval_latest_v4.py   # Structured evaluation of most recent v4 run
├── tests/
├── run_v4.py               # v4 entry point (run / validate / plot)
├── validate_v4.py          # Full publication-gate runner
└── run.py                  # Legacy v3 CLI wrapper
```

---

## Legacy v3 Workflow

The v3 pipeline remains available for historical analysis:

```bash
python run.py configs/config_kepler432_case.yaml
```

Existing v3 directories are recognized by `slingshot.v4.runs.classify_run()`
and marked `legacy_science_model: true`. `discover_v4_runs()` returns only
eligible v4 runs.

**v3 results should not be presented as scientific conclusions about Kepler-432 b.
All known limitations are documented in `METHODOLOGY_AUDIT.md`.**

---

## Known Limits and Roadmap

- The current core is planar. 3D isotropic cross-sections are the next major milestone.
- Reported widths are proposal-dependent and must include their sampling definition.
- Astrophysical rates require an external velocity distribution and are deferred.
- Quinn and Ortiz represent discrete observational-model uncertainty, not a posterior.
- REBOUND IAS15 validation requires the optional REBOUND package.

---

## License

MIT

**Repository status:** v4 scientific research core (4.0.0) with retained v3 legacy pipeline.
