# Slingshot Solver v4 — Defensible Planar Research Core

The v4 core is a versioned scientific workflow beside the legacy v3 package.
It corrects the orbital model, inbound proposal, energy definitions, numerical
validation, and statistical estimand without rewriting or silently
reinterpreting historical runs.

## Scientific Claim

v4 estimates **effective planar encounter width** as a function of asymptotic
speed:

\[
W(\Delta\epsilon/v_c^2 > q \mid v_\infty)
= 2 b_{\max}\frac{N_{\rm event}}{N}.
\]

This is not a three-dimensional area cross-section, event probability, or
astrophysical rate. Those require the later 3D milestone.

## Key Corrections

- Eccentric Keplerian star-planet initial states use the configured
  semi-major axis, eccentricity, mean anomaly, and periapsis orientation.
- Incoming particles are defined by \(v_\infty\), signed impact parameter,
  incoming direction, and binary mean anomaly.
- Scientific gain is
  `delta_specific_energy_com = epsilon_out - epsilon_in`.
- `turning_quadratic = 0.5 * |v_out - v_in|^2` is explicitly a turning
  diagnostic, not an energy gain.
- Stellar and planetary moving-potential work are integrated independently
  and checked against total COM energy change.
- Periapses and collisions are located by root events.
- Newtonian gravity (`softening_km: 0`) is the scientific default.
- Wilson intervals and impact-boundary tail checks accompany every width.

## Presets

- `configs/v4_dimensionless_reference.yaml`
- `configs/v4_kepler432_quinn.yaml`
- `configs/v4_kepler432_ortiz.yaml`

Quinn and Ortiz are intentionally separate observational models. Their
results should be compared as discrete model uncertainty, not pooled as if
they were draws from one Gaussian posterior.

## Commands

Run deterministic configuration checks:

```bash
python run_v4.py validate configs/v4_kepler432_quinn.yaml
```

Run all publication gates, including DOP853/Radau agreement and optional
REBOUND IAS15:

```bash
python validate_v4.py configs/v4_kepler432_quinn.yaml
```

Run a small pilot before committing to the full 70,000-trajectory preset:

```bash
python run_v4.py run configs/v4_kepler432_quinn.yaml \
  --samples-per-bin 100 --seeds 42 \
  --output-dir results/v4_quinn_pilot
```

Run the configured campaign:

```bash
python run_v4.py run configs/v4_kepler432_quinn.yaml
```

## Artifacts

Each run writes:

- `config.yaml` — frozen schema-v4 configuration
- `samples.csv` — proposal variables, outcomes, physical metrics, and solver diagnostics
- `width_summary.csv` — per-seed and pooled effective widths
- `manifest.json` — schema, provenance, Git commit, integrator, and validation status
- `REPORT.md` — interpretation-safe scientific summary

The compact format deliberately does not retain every coarse ODE trajectory.

## Validation Gates

- Keplerian element reconstruction
- Central-force energy and angular-momentum conservation
- Circular Jacobi-constant conservation
- Galilean invariance of COM metrics
- Moving-potential work-energy closure
- DOP853 versus Radau agreement
- Optional REBOUND IAS15 agreement
- Zero-softening and proposal-boundary declarations in manifests

REBOUND is optional:

```bash
pip install rebound
```

## Legacy Runs

`slingshot.v4.runs.classify_run()` marks existing v3 result directories as
`legacy_science_model: true`. `discover_v4_runs()` excludes them from v4
aggregates. Historical results remain readable but are not mixed with the
corrected estimand.
