# Statistical Methodology and Evidence Audit

**Audit date:** June 20, 2026  
**Scope:** Existing result artifacts, legacy v3 data collection and plots, and
the implemented v4 planar-width methodology.

## Executive Verdict

The project now contains the foundation of a meaningful research method in
v4, but the data currently stored under `results/` do not support the v4
scientific claims.

- Every existing result directory is classified as schema-v3 legacy science.
- No completed v4 pilot or campaign is currently discoverable under
  `results/`.
- The newest Ortiz-labeled run uses a circular orbit at 0.0896 AU, not the
  cited eccentric Ortiz orbit at 0.303 AU and \(e=0.478\).
- Only 21 of 161 accepted trajectories in that run reached the configured
  outbound escape event. The other 140 reached the time limit.
- Twenty-five percent of accepted trajectories crossed inside the stellar
  photosphere because the stellar-clearance filter was disabled.
- The phase-gradient plots show a maximum over a hidden velocity dimension,
  not a resolved conditional mean or event probability.
- Legacy energy and mechanism plots mix turning diagnostics, scalar speed
  changes, and heuristic attribution in ways that cannot support energy-gain
  claims.

The defensible position is therefore:

> Existing v3 runs are useful for software diagnostics and hypothesis
> generation, but not for quantitative scientific inference. The v4 estimand
> and asymptotic sampling design are substantially better, but several
> statistical and validation gaps should be fixed before generating
> publication-facing campaigns.

## Available Evidence

Four legacy result directories are present:

| Run | Attempted | Accepted | Acceptance | Selected |
|---|---:|---:|---:|---:|
| Ortiz-labeled, June 20, 2026 | 3,000 | 161 | 5.37% | 17 |
| Kepler-432, June 19, 2026 | 10,000 | 505 | 5.05% | 51 |
| Kepler-432, March 8, 2026 21:16 | 1,000 | 39 | 3.90% | 5 |
| Kepler-432, March 8, 2026 19:32 | 40,000 | 1,920 | 4.80% | 192 |

These runs are not replicates of one statistical experiment. Their proposal
ranges, filters, bulk velocities, softening, selection modes, and sample sizes
differ. They should not be pooled or treated as convergence evidence.

`slingshot.v4.runs.discover_v4_runs("results")` currently returns no runs.

## Legacy v3 Data Collection

### 1. The proposal is not asymptotic sampling

The legacy sampler draws:

- finite initial speed;
- a positive quantity labeled impact parameter;
- an incoming angle;
- an azimuth;
- an initial radius set to twice the sampled impact parameter when no fixed
  radius is supplied.

The finite state is not reconstructed from asymptotic energy and angular
momentum. Consequently:

- the configured speed is not \(v_\infty\);
- the configured impact parameter is not generally the asymptotic impact
  parameter;
- the initial radius is perfectly correlated with the configured impact
  parameter;
- the angle and offset jointly alter the actual angular momentum in a way
  that is difficult to interpret as a known proposal density.

For the Ortiz-labeled run:

- 14.1% of proposals were initially bound under a total-mass monopole
  reconstruction;
- configured impact parameters of 0.5–3 AU mapped to approximate monopole
  asymptotic impact parameters ranging from 0.004 to 26.3 AU;
- accepted configurations occupied only a narrow portion of the nominal
  incoming-angle range.

This means the acceptance fraction is a property of an arbitrary finite-state
proposal, not a physical cross-section or encounter probability.

### 2. The binary model does not match Kepler-432

The legacy initializer is circular and defaults to 0.0896 AU. The pipeline does
not pass `cfg.system.a_planet_AU` into that initializer.

The newest Ortiz-labeled run also freezes:

```yaml
a_planet_AU: 0.0896
orbital_phase_rad: 0.0
```

Ortiz et al. report approximately:

- \(a=0.303\pm0.007\) AU;
- \(e=0.478\pm0.004\);
- \(M_p=5.84\pm0.05\,M_{\rm Jup}\).

Quinn et al. report \(e=0.5134^{+0.0098}_{-0.0089}\), with a distinct stellar
and planetary solution.

Therefore the stored v3 data do not represent either implemented v4
observational preset.

### 3. Binary phase is fixed

The v3 run uses `orbital_phase_rad: 0.0` for every particle. An eccentric
binary cannot be represented, and even the circular calculation samples only
one binary orientation relative to the proposal construction.

This prevents interpretation as an average over observation time.

### 4. Acceptance is dominated by post-sampling filters

The Ortiz-labeled run retained 161 of 3,000 proposals. Rejections were:

| Outcome | Count | Fraction |
|---|---:|---:|
| Flyby too distant | 2,216 | 73.9% |
| Unbound requirement failed | 423 | 14.1% |
| Flyby incomplete | 197 | 6.6% |
| Accepted | 161 | 5.4% |
| Analysis failed | 3 | 0.1% |

Acceptance varies strongly across proposal coordinates:

- lowest speed quintile: 0% accepted;
- smallest configured-impact quintile: 10.3% accepted;
- largest configured-impact quintile: 1.0% accepted;
- one incoming-angle quintile: 25.8% accepted;
- three central/positive-angle quintiles: 0% accepted.

Plots made only from accepted trajectories are therefore conditioned on a
strong selection event. Correlations in that subset are not estimates of
unconditional physical relationships and may contain collider bias.

### 5. Most accepted trajectories are not outbound solutions

Legacy encounter states are selected using the first and last stored nodes
whose planet-relative distance exceeds a threshold. That is not equivalent to
an inbound and outbound boundary crossing.

In the Ortiz-labeled run:

- 21/161 accepted trajectories reached the terminal escape event;
- 140/161 ended at the maximum integration time;
- only 93/161 had outbound barycentric radius at least as large as the
  reported inbound radius;
- the outbound-to-inbound radius ratio ranged from 0.056 to 3.0.

Endpoint speed and energy differences are therefore frequently comparisons
between unequal, non-asymptotic states.

### 6. Stellar collisions contaminate the accepted population

The stellar-clearance filter is disabled in the Ortiz-labeled configuration.
The star has radius approximately 2,825,760 km, while accepted trajectories
reach as close as 8,855 km.

Among accepted trajectories:

- 41/161, or 25.5%, cross the stellar photosphere;
- 3/17 selected candidates cross the stellar photosphere.

Removing stellar-interior paths reduces the reported maximum turning
quadratic from 27,268 to 13,995 km²/s². This confirms that a substantial part
of the extreme turning signal is produced by physically invalid
stellar-interior motion.

### 7. Sample maxima are unstable

The report emphasizes best candidates and maxima. Randomly partitioning the
161 accepted Ortiz trajectories into ten groups produced group maximum scalar
speed changes ranging from about 4.5 to 44.2 km/s.

This is ordinary extreme-value instability. A sample maximum should not be
described as a physical limit, ceiling, or converged optimum without a
dedicated extreme-value or bounded optimization analysis.

### 8. The stored seed is null

The Ortiz-labeled run has `pipeline.seed: null`. The frozen configuration does
not reproduce the proposal draw. `results.pkl` preserves the realized states,
so the artifact can be inspected, but the run cannot be independently
regenerated from configuration plus code revision.

The legacy artifact also lacks a schema-v4 manifest, Git commit, validation
status, and package-version record.

## Legacy Metrics

### Scalar speed change

`delta_v = |v_f|-|v_i|` is descriptive but is not itself the specific-energy
change unless potential terms are equal or negligible at both endpoints.

For the accepted Ortiz trajectories, the Spearman correlation between scalar
speed change and endpoint total specific-energy change is only about 0.29.
This metric is not a reliable substitute for energy gain in the stored runs.

### Turning quadratic

`0.5 * |v_f-v_i|²` is a turning diagnostic. It can be very large for a
constant-speed reversal with zero kinetic-energy gain.

Legacy reports and plots repeatedly label it “scattering energy” and compare
its maxima as if it were an energy-gain ceiling. Those comparisons do not
support claims about acceleration or extracted orbital energy.

### Heuristic planet attribution

The legacy `energy_from_planet_orbit` subtracts a monopole endpoint difference
from the three-body endpoint difference. It is not a work integral and is
sensitive to unmatched endpoints.

In the Ortiz accepted population, its Spearman correlation with total endpoint
specific-energy change is about 0.06 and is not statistically distinguishable
from zero in that sample.

The tier labels and mechanism claims derived from this quantity should remain
legacy diagnostics only.

## Plot Audit

### Phase energy and heading maps

**Verdict: not scientifically interpretable in their current form.**

The phase grid has 140 impact bins and 120 angle bins: 16,800 cells. Only
3,663 cells are populated. Each populated cell contains exactly 20 values from
the hidden velocity grid.

The plotting function retains the maximum energy in each impact-angle cell:

```text
if sample_e[k] > e_grid[i, j]:
    e_grid[i, j] = sample_e[k]
```

It therefore projects a three-dimensional deterministic sweep
\((v,b,\alpha)\) into two dimensions by taking a maximum over \(v\).

Measured within populated cells:

- median relative energy standard deviation: approximately 170%;
- 90th-percentile relative standard deviation: approximately 181%;
- median within-cell energy range: approximately 146 km²/s²;
- 90th-percentile within-cell range: approximately 2,560 km²/s².

The square blocks are not merely a rendering problem. They reveal that the
plot is showing per-cell extrema over a highly variable hidden dimension.

Recommended replacement:

- facet by fixed \(v_\infty\), or put \(v_\infty\) on one axis;
- display event probability or median gain, not maxima;
- include a count/support panel;
- require enough independent samples per cell for the intended uncertainty;
- show binomial or bootstrap intervals;
- use adaptive bins only after defining the estimand;
- reserve smoothing for visualization and never treat it as added evidence.

### Trajectory-track gradient

**Verdict: illustrative only.**

The background hexbin assigns each trajectory’s single energy value to every
stored point along that trajectory. Adaptive solver nodes are not equal-time
or equal-arclength samples, so dense regions can reflect solver behavior and
trajectory duration rather than physical probability density.

Overlaid trajectories are useful examples, but the colored background is not
a spatial probability or energy field.

### Velocity phase plots

**Verdict: one panel is descriptive; the radial-normal panel is currently
incorrect.**

The radial basis is computed from `star_position - planet_position`, not
`particle_position - planet_position`. The resulting radial and normal
components are relative to the star-planet axis rather than the particle’s
planet-relative radius vector.

Even after correction, these plots describe one selected trajectory and do
not provide population-level statistical evidence.

### Energy CDF

**Verdict: invalid comparison.**

The coarse Monte Carlo curve uses `0.5 * scalar_delta_v²`, while the selected
rerun curve and baseline lines use `0.5 * |vector_delta_v|²`. The figure labels
them as one scattering-energy distribution.

These are different quantities and cannot share one CDF or ceiling
comparison.

### Sampling histograms

**Verdict: useful quality-control figures.**

Showing all proposals and accepted proposals is valuable. Improve them by
plotting acceptance probability with binomial intervals rather than overlaying
raw histograms whose areas differ by a factor of about twenty.

### Rejection breakdown

**Verdict: useful quality control, not a physical result.**

It identifies where computational effort is lost and how strongly the
analysis conditions the sample. It should be reported by speed bin and seed in
v4.

### Star-proximity distributions

**Verdict: essential physical-validity diagnostic.**

This figure correctly exposes stellar-interior trajectories. In future runs,
physical collision outcomes should be handled by root events before any
scientific metric is interpreted.

### Deflection and pairwise-correlation plots

**Verdict: exploratory only.**

They use accepted trajectories only, contain no uncertainty bands or density
normalization, and do not control for speed, impact parameter, phase, or
selection. They can motivate hypotheses but not establish mechanisms.

### Pareto and candidate-ranking plots

**Verdict: optimization diagnostics only.**

They visualize a selected tail of the run and include legacy metrics that are
not valid gain or mechanism estimators. They should not be presented as
population statistics.

### Best-candidate and multi-trajectory figures

**Verdict: useful illustrations after physical validation.**

They should be labeled as examples and accompanied by proposal probability,
outcome class, collision clearance, convergence checks, and cross-integrator
agreement.

## Legacy Uncertainty and Robustness

The June 19 run enabled parameter-posterior reruns but disabled robustness.
The uncertainty workflow is not sufficient for the research estimand:

- it reruns only preselected candidates;
- it does not reselect the population for each observational draw;
- it assumes independent Gaussian parameter distributions;
- it changes the system while retaining finite particle states;
- `a_planet_AU` is drawn but is not passed to the circular initializer;
- confidence bands therefore describe conditional candidate sensitivity, not
  uncertainty in encounter width or event probability.

The legacy convergence workflow tracks sample maxima from one nested random
sequence. Maxima are intrinsically sample-size dependent, and a less-than-1%
change between two sample sizes is not evidence of convergence.

The sensitivity workflow also compares maxima with one seed. It is vulnerable
to Monte Carlo noise and changes in filtering.

## v4 Methodology

### What is now meaningful

The v4 core makes several necessary corrections:

- Keplerian eccentric binary states;
- proper asymptotic \(v_\infty\) and signed-impact proposals;
- uniform binary mean anomaly;
- center-of-mass specific-energy gain;
- root-event collision and escape handling;
- moving-potential work integrals;
- a declared planar-width estimand;
- Wilson intervals;
- multiple seeds;
- compact per-sample proposal, outcome, and solver records;
- version-aware exclusion of legacy runs.

For the stated planar estimand, uniform signed impact sampling and
\(W=2b_{\max}N_{\rm event}/N\) are statistically coherent.

### Remaining v4 gaps

#### 1. No substantive v4 data exist yet

The method has tests and a tiny smoke run, but no research pilot is present
under `results/`. The project currently has a method, not empirical v4 support.

#### 2. The analytic-deflection gate does not test deflection

The validation routine computes analytic and numerical deflection, but the
gate checks only energy and angular-momentum conservation.

In the existing smoke validation:

- analytic deflection: approximately 3.0991 rad;
- numerical finite-boundary velocity-angle difference: approximately
  3.0819 rad.

That difference is much larger than the configured \(10^{-6}\) tolerance, yet
the gate passes. Either the finite-boundary observable must be analytically
corrected or the test must compare asymptotic directions.

#### 3. Numerical failures are counted as non-events

The width denominator is currently `len(rows)`. Integration failures and time
limits therefore reduce the width estimate as if they were verified misses.

Unresolved numerical outcomes should instead:

- fail a campaign validation gate above a very small declared rate;
- be reported separately;
- be rerun with stricter settings or larger time limits;
- produce lower and upper identification bounds if any remain unresolved.

#### 4. The tail check is a point estimate

If no threshold events occur, the tail check automatically passes. With few
events, zero observed tail events does not demonstrate a negligible tail.

The tail gate should use an upper confidence bound on tail contribution and
should trigger an automated larger-`b_max` rerun when it fails or is
uninformative.

#### 5. “Effective sample size” is hard-coded

It is currently reported as the number of trials. That is correct for equal,
independent Bernoulli proposals within a bin, but the implementation should
either call it `N` or compute ESS from weights if importance sampling is later
introduced.

#### 6. Seeds are restarted for every speed bin

Each speed bin initializes `default_rng(seed)` again. This creates common
random numbers across speeds: matching sample indices share the same impact,
direction, and phase quantiles.

This can be beneficial for paired speed comparisons, but it must be declared
and analyzed as a paired design. If independent bins are intended, derive
child streams from speed, seed, model, and campaign identifiers.

#### 7. Seed variability is not emphasized

Per-seed summaries are saved, but the report focuses on pooled Wilson
intervals. Add:

- seed-level points;
- between-seed variance;
- a heterogeneity diagnostic;
- paired Quinn-Ortiz differences when common random numbers are used.

#### 8. Simultaneous inference is absent

Seven speeds and six nested thresholds produce 42 related width estimates per
model. Pointwise 95% intervals are acceptable for descriptive work but do not
provide simultaneous 95% coverage.

Use simultaneous bootstrap bands or explicitly state that intervals are
pointwise and avoid binary significance language.

#### 9. Rare-event precision may be inadequate

With 10,000 pooled proposals per speed:

- 10 events imply roughly 65% relative interval half-width;
- 25 events imply roughly 40%;
- 100 events imply roughly 20%.

The fixed 10,000-sample design is adequate only when event probabilities are
not too small. Pilot event counts should determine final sample sizes for each
threshold and speed.

#### 10. Observational uncertainty is metadata only

Quinn and Ortiz are correctly separated as discrete models, but parameter
uncertainties are not propagated in v4.

In addition, Quinn metadata should preserve the published asymmetric
uncertainties. The current preset symmetrizes several values, including the
planet mass, radius, and eccentricity.

#### 11. Validation remains incomplete

The implemented gates do not yet demonstrate:

- boundary-radius convergence;
- tolerance convergence over representative trajectory classes;
- `b_max` convergence;
- zero-softening convergence against diagnostic softened runs;
- independent-integrator agreement in the exact published observables;
- campaign-level failure-rate control.

The REBOUND check compares the norm of the complete final state, which can be
dominated by large position coordinates. It should compare periapsis,
outbound direction, and specific-energy gain directly.

#### 12. Version metadata disagree

The v4 module and manifest report 4.0.0 while `pyproject.toml` still reports
3.0.0. Reproducibility metadata should come from one package-version source.

## Recommended v4 Figures

The v4 reporting layer should prioritize figures aligned with the estimand:

1. **Width versus \(v_\infty\)**  
   One curve per threshold, Wilson intervals, individual seed points, and
   separate Quinn/Ortiz panels.

2. **Outcome fractions versus \(v_\infty\)**  
   Escaped, stellar collision, planetary collision, time limit, and numerical
   failure with intervals.

3. **Tail support plot**  
   Event probability versus \(|b|/b_{\max}\), including binomial intervals and
   the outer-tail boundary.

4. **Gain distributions**  
   ECDFs or quantile intervals of dimensionless gain by speed, explicitly
   conditional on escape.

5. **Seed stability**  
   Seed-level width estimates and pooled estimates.

6. **Quinn-Ortiz paired differences**  
   If common random numbers are retained, plot paired width differences and
   paired confidence intervals.

7. **Work-energy diagnostics**  
   Closure residual distribution and signed stellar versus planetary work.

8. **Convergence dashboard**  
   Width changes under sample size, boundary radius, tolerance, `b_max`, and
   integrator.

9. **Conditional phase maps**  
   Fixed-speed facets showing event probability or median gain over signed
   impact parameter and relative binary phase. Include count and uncertainty
   panels. Do not maximize over a hidden variable.

10. **Representative trajectories**  
    A small, preregistered set of median, high-quantile, collision, and
    mechanism-extreme examples. These remain illustrations, not estimators.

## Priority Remediation

### P0 — Required before scientific campaign claims

1. Fix the analytic two-body deflection gate.
2. Add campaign gates for time limits and integration failures.
3. Add boundary-radius, tolerance, and `b_max` convergence.
4. Replace the point tail check with a confidence-bound rule.
5. Correct Quinn asymmetric uncertainty metadata.
6. Unify package and manifest versions.
7. Generate and archive true schema-v4 pilot runs for both presets.

### P1 — Required before graduate-level statistical presentation

1. Perform pilot-based sample-size calculations.
2. Add seed-level variance and simultaneous uncertainty bands.
3. Declare common-random-number pairing or make speed bins independent.
4. Propagate observational uncertainty within each discrete preset.
5. Add v4 estimand-aligned figures and remove legacy “energy” terminology.
6. Report numerical-outcome identification bounds if unresolved cases remain.

### P2 — Improvements for publication and extension

1. Add adaptive or stratified sampling for rare thresholds with correct
   weights.
2. Add response-surface diagnostics only after the probability estimand is
   defined.
3. Add full 3D isotropic sampling for physical area cross-sections.
4. Add external population models only after 3D cross-sections are validated.

## Final Assessment

### Legacy v3

**Research support: weak.**  
It is valuable as exploratory software and as evidence that the numerical
system can generate diverse scattering trajectories. It does not provide
defensible estimates of Kepler-432 encounter gain, physical limits,
cross-sections, or rates.

### v4 implementation

**Research design: promising but not yet publication-ready.**  
The primary planar-width estimand, asymptotic proposal, COM energy metric, and
event handling are appropriate foundations. Once the remaining validation,
failure handling, convergence, tail inference, and reporting issues are fixed,
the methodology can support a strong graduate-level planar study.

### Current project state

The project should be described as:

> A validated computational framework under active methodological
> development, with legacy exploratory results and an implemented but not yet
> fully campaigned v4 planar-width model.

It should not yet be described as having statistically established a
Kepler-432 slingshot efficiency, maximum gain, or event cross-section.

## Observational Sources

- Quinn et al., *Kepler-432: a red giant interacting with one of its two long
  period giant planets*: https://arxiv.org/abs/1411.4666
- Ortiz et al., *Kepler-432 b: a massive warm Jupiter in a 52-day eccentric
  orbit transiting a giant star*: https://arxiv.org/abs/1410.3000
