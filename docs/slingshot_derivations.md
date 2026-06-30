# Slingshot Solver: Mathematical Derivations

Each section derives one piece of the methodology from first principles, then
annotates why it matters for the actual code and which bug or gate it connects to.

---

## 1. Why an Eccentric Binary Has No Conserved Energy for the Test Particle

In the circular restricted three body problem you can move to a frame
corotating with the binary at angular rate $n$. In that frame the combined
gravitational potential is static, and the Jacobi integral

$$C_J = \frac{1}{2}|v_{rot}|^2 + \Phi_{eff}(r)$$

is conserved, where $\Phi_{eff}$ includes the centrifugal term. This only
works because $n$ is constant, which requires a circular binary.

For Kepler-432 with $e = 0.478$, the binary's angular rate is not constant,
so there is no rotating frame in which the potential is time-independent.
The test particle's equation of motion in the inertial barycentric frame is

$$\ddot{\mathbf{r}} = -\frac{GM_\star(\mathbf{r}-\mathbf{r}_\star(t))}{|\mathbf{r}-\mathbf{r}_\star(t)|^3} - \frac{GM_p(\mathbf{r}-\mathbf{r}_p(t))}{|\mathbf{r}-\mathbf{r}_p(t)|^3}$$

where $\mathbf{r}_\star(t)$ and $\mathbf{r}_p(t)$ both move on the binary's
Keplerian orbit. Because the right hand side depends explicitly on $t$
through the binary's motion, there is no time-independent potential and
therefore no conserved mechanical energy for the test particle.

**Annotation:** this is not a defect to engineer around. It is the entire
mechanism. The slingshot gain *is* the net work the time-varying field does
on the particle as it crosses. Everything downstream (the work integrals,
the COM frame, the energy histograms) exists to measure that one
non-conserved quantity carefully.

---

## 2. Two Body Asymptotic Kinematics: $b$, $e$, $v_\infty$, and Deflection

The boundary mapping and the deflection gate both depend on standard
hyperbolic two-body relations. Worth rederiving since this is exactly where
the P0.1 bug lived.

**Setup.** Treat the inbound leg, before the particle feels the binary as
anything but a point mass $\mu = G(M_\star+M_p)$, as a two-body hyperbolic
orbit. Specific energy and angular momentum:

$$\varepsilon = \frac{v_\infty^2}{2}, \qquad h = b\,v_\infty$$

The angular momentum relation $h = b v_\infty$ holds because far from the
focus the path is a straight line with speed $v_\infty$, and the angular
momentum of straight line motion about any fixed point equals (perpendicular
distance) × (speed), conserved along the line.

**Semi-major axis and eccentricity.** From $\varepsilon = -\mu/2a$,

$$a = -\frac{\mu}{2\varepsilon} \quad (a<0 \text{ for a hyperbola})$$

From $h^2 = \mu a(1-e^2)$, substitute $a$:

$$h^2 = -\frac{\mu^2(1-e^2)}{2\varepsilon} \;\;\Rightarrow\;\; e^2 = 1+\frac{2\varepsilon h^2}{\mu^2}$$

Plugging in $\varepsilon = v_\infty^2/2$ and $h = bv_\infty$:

$$e^2 = 1 + \left(\frac{b\,v_\infty^2}{\mu}\right)^2 \tag{2.1}$$

**Deflection angle.** The orbit equation $1/r = (\mu/h^2)(1+e\cos\theta)$
gives the asymptote direction where $r\to\infty$:

$$\cos\theta_\infty = -\frac{1}{e}$$

The total turning angle $\chi$ (the angle between incoming and outgoing
velocity vectors) is related to $\theta_\infty$ by $\theta_\infty = \pi/2+\chi/2$,
which after substitution gives the standard closed form

$$\sin\!\left(\frac{\chi}{2}\right) = \frac{1}{e} \quad\Longleftrightarrow\quad \chi = 2\arcsin\!\left(\frac{1}{e}\right) \tag{2.2}$$

Sanity checks: $e\to\infty$ (large $b$) gives $\chi\to0$, no deflection.
$e\to1^+$ (small $b$, near head-on) gives $\chi\to\pi$. Both correct.

**Annotation, the P0.1 bug.** The old gate checked that energy and $h$
were conserved at the finite boundary, then stopped. But $(2.2)$ shows the
deflection angle is a separate, independent function of $e$ alone. You can
conserve $\varepsilon$ and $h$ at the boundary to high precision while the
*numerically recovered eccentricity vector* still encodes a deflection that
differs from the analytic asymptotic value, since the boundary truncates the
orbit before it has fully approached its asymptote. The fix, computing
$\chi_{num} = 2\arcsin(1/e_{num})$ from the recovered eccentricity vector and
comparing it to the analytic $\chi$, is checking a genuinely different
invariant than energy or $h$ alone. That is why the bug was invisible to the
old gate: it was the right two conserved quantities and the wrong
derived one.

---

## 3. The Work-Energy Identity and What Closure Error Actually Measures

Dot the equation of motion with velocity:

$$\dot{\mathbf{r}}\cdot\ddot{\mathbf{r}} = \dot{\mathbf{r}}\cdot\mathbf{a}_\star + \dot{\mathbf{r}}\cdot\mathbf{a}_p$$

$$\frac{d}{dt}\left(\frac{1}{2}|\dot{\mathbf{r}}|^2\right) = \dot{\mathbf{r}}\cdot\mathbf{a}_\star(t) + \dot{\mathbf{r}}\cdot\mathbf{a}_p(t)$$

Integrate from the inbound boundary crossing $t_{in}$ to the outbound
crossing $t_{out}$:

$$\frac{1}{2}v_{out}^2 - \frac{1}{2}v_{in}^2 = \underbrace{\int_{t_{in}}^{t_{out}} \dot{\mathbf{r}}\cdot\mathbf{a}_\star\,dt}_{W_\star} + \underbrace{\int_{t_{in}}^{t_{out}} \dot{\mathbf{r}}\cdot\mathbf{a}_p\,dt}_{W_p} \tag{3.1}$$

This is exact. It is nothing more than Newton's second law dotted with
velocity and integrated, true for any force law, true at any boundary
placement, true regardless of how energy is otherwise defined.

`work_star`, `work_planet` are $W_\star$, $W_p$ computed by quadrature along
the stored trajectory. The closure check compares the left side (computed
directly from the raw state vectors at $t_{in}, t_{out}$) against the right
side (computed by integrating force dot velocity). Any mismatch beyond
quadrature error is a real signal that something in the integration is
wrong, since $(3.1)$ has no free parameters and no modeling assumption.

**Why the tolerance had to move from $10^{-6}$ to $10^{-4}$.** The
left-hand side is a difference of two $O(v_\infty^2)$ numbers near periapsis
passage, where the integrand on the right has a sharp, large-amplitude peak
that DOP853 must resolve with finite step size and finite floating point
precision. The empirical $2.6\times10^{-5}$ closure at $\text{rtol}=\text{atol}=10^{-10}$
is the realistic floor for that integrator on this problem; demanding
$10^{-6}$ was asking for accuracy below the integrator's own error budget,
which the tolerance_convergence gate likely exposed directly by showing the
residual saturates rather than continuing to shrink with tighter settings.

**A separate subtlety: boundary radius and residual potential.** $(3.1)$
is an identity for *kinetic* energy specifically. If $v^2$ at a finite
boundary is used as a proxy for the asymptotic speed, it still carries
whatever potential energy has not yet been converted to kinetic form at that
radius. Writing $\varepsilon_{kin} = \varepsilon_{true} - \Phi(r)$ with
$\Phi(r)<0$, the kinetic-only estimate converges to the true asymptotic
energy only as $r\to\infty$. At finite $r$ there's a residual offset set
by the local potential, including any multipole structure beyond the
monopole term that only vanishes once the boundary is well outside the
binary's own separation scale. This is the reasoning behind treating
periapsis, a purely geometric quantity with no such residual, as the
boundary-independent check in `boundary_radius_convergence`, while leaving
energy gain out of that particular comparison. I'm inferring the general
shape of this argument from the methodology notes rather than your exact
internal formula, so treat this paragraph as the principle rather than a
verified implementation detail.

---

## 4. Galilean Invariance of the COM-Frame Energy Gain

Apply a constant boost $\mathbf{V}$ to every body in the system:

$$\mathbf{r}'_i(t) = \mathbf{r}_i(t) + \mathbf{V}t, \qquad \mathbf{v}'_i = \mathbf{v}_i + \mathbf{V}$$

Accelerations depend only on relative separations:

$$\mathbf{r}'_{test}-\mathbf{r}'_\star = (\mathbf{r}_{test}+\mathbf{V}t)-(\mathbf{r}_\star+\mathbf{V}t) = \mathbf{r}_{test}-\mathbf{r}_\star$$

unchanged, so the dynamics in the boosted frame is identical; only every
velocity shifts by the constant $\mathbf{V}$.

Now subtract the COM velocity, which also shifts by $\mathbf{V}$ under the
same boost:

$$\mathbf{v}'_{test}-\mathbf{v}'_{COM} = (\mathbf{v}_{test}+\mathbf{V})-(\mathbf{v}_{COM}+\mathbf{V}) = \mathbf{v}_{test}-\mathbf{v}_{COM}$$

The $\mathbf{V}$ cancels exactly. So the COM-relative velocity, and
therefore $\varepsilon$ measured in the COM frame and $\Delta\varepsilon_{COM}$,
is invariant under any constant boost. That is the entire content of the
`galilean_invariance` gate: add a config-level `bulk_velocity`, rerun, and
confirm $\Delta\varepsilon_{COM}$ doesn't move.

**Why this kills the legacy metrics.** Take the lab-frame scalar speed
change $|\mathbf{v}_{out}|-|\mathbf{v}_{in}|$. Under the same boost,

$$|\mathbf{v}'_{out}|-|\mathbf{v}'_{in}| = |\mathbf{v}_{out}+\mathbf{V}|-|\mathbf{v}_{in}+\mathbf{V}|$$

and by the triangle inequality this does **not** equal
$|\mathbf{v}_{out}|-|\mathbf{v}_{in}|$ in general; the boost enters
nonlinearly through the vector norm. `bary_delta_v_pct` and the
`turning_quadratic` diagnostic both have this defect, which is exactly why
the methodology bans them as gain estimators and keeps them only as
non-invariant geometric diagnostics.

---

## 5. The Dimensionless Energy Scale $v_c^2$

$$v_c^2 = \frac{G(M_\star+M_p)}{a}$$

From the vis-viva equation $v^2 = GM_{tot}(2/r - 1/a)$, evaluating at
$r=a$ gives exactly $v^2 = GM_{tot}/a = v_c^2$. So $v_c$ is the speed the
binary components would have at separation equal to their own semi-major
axis, independent of eccentricity. It's the natural orbital velocity scale
of the binary itself.

**Annotation.** Dividing gain by $v_c^2$ turns
$\Delta\varepsilon_{COM}/v_c^2$ into "fraction of the binary's own orbital
energy scale gained," which is exactly what lets Quinn and Ortiz, two
different $(M_\star, M_p, a)$ tuples, be compared on the same threshold
grid without the comparison being contaminated by their different absolute
mass and separation scales.

---

## 6. From Cross Section to Planar Width

In an unrestricted 3D scattering problem with full rotational symmetry
about the incoming axis, the cross section for an outcome is

$$\sigma(>q) = \int_0^\infty 2\pi b\,\mathbb{1}[\text{gain}(b)>q]\,db$$

the $2\pi b$ coming from integrating over the azimuthal angle around the
approach axis, since nothing in a spherically symmetric problem depends on
that angle.

Here the dynamics is confined to the binary's orbital plane. There is no
azimuthal symmetry to integrate out, because the impact parameter is a
*signed*, one-dimensional coordinate within that plane, not a radius with a
free rotational degree of freedom. The natural object is therefore a 1D
measure, the planar width:

$$W(>q\mid v_\infty) = \int_{-b_{max}}^{b_{max}} \mathbb{1}[\text{gain}(b,\theta,M)>q]\,db$$

with direction $\theta$ and binary mean anomaly $M$ treated as nuisance
parameters that are simultaneously marginalized by sampling them uniformly
alongside $b$. Since the joint sampling measure factorizes
(uniform-$b$ $\times$ uniform-$\theta$ $\times$ uniform-$M$), the standard
Monte Carlo estimator of this integral is just

$$\hat{W}(>q\mid v_\infty) = 2b_{max}\cdot\frac{N_{event}}{N} \tag{6.1}$$

which is the formula in the README. It's the same logic as a Monte Carlo
cross-section estimate in scattering physics generally, one dimension lower
because the planar restriction removes the azimuthal integral rather than
approximating it.

---

## 7. The Wilson Score Interval, Derived

$N_{event}$ out of $N$ trials satisfying a criterion is
$N_{event}\sim\text{Binomial}(N,p)$ for unknown true rate $p$. The naive
Wald interval $\hat{p}\pm z\sqrt{\hat{p}(1-\hat{p})/N}$ plugs in $\hat{p}$
for the variance, which is known to undercover badly near $p=0$ or $p=1$
or at small $N$, exactly the regime of a sparse high-threshold tail.

The Wilson interval instead inverts the score statistic, which uses the
*true* $p$ in the variance rather than the estimate:

$$\frac{\hat{p}-p}{\sqrt{p(1-p)/N}} \approx \mathcal{N}(0,1)$$

Square both sides and clear denominators:

$$(\hat{p}-p)^2 = \frac{z^2\,p(1-p)}{N}$$

This is quadratic in $p$:

$$p^2\left(1+\frac{z^2}{N}\right) - p\left(2\hat{p}+\frac{z^2}{N}\right) + \hat{p}^2 = 0$$

Solving:

$$p_{\pm} = \frac{\hat{p}+\dfrac{z^2}{2N}\;\pm\; z\sqrt{\dfrac{\hat{p}(1-\hat{p})}{N}+\dfrac{z^2}{4N^2}}}{1+\dfrac{z^2}{N}} \tag{7.1}$$

This never leaves $[0,1]$, has materially better coverage at the edges
than Wald, and crucially stays well-behaved even when $\hat{p}=0$.

---

## 8. The One-Sided Upper Bound and Why the Tail Gate Needs It

Take only the upper root of $(7.1)$ with a one-sided critical value
$z_{1-\alpha}$ (e.g. $z=1.645$ for 95%):

$$p_{upper} = \frac{\hat{p}+\dfrac{z^2}{2N}+ z\sqrt{\dfrac{\hat{p}(1-\hat{p})}{N}+\dfrac{z^2}{4N^2}}}{1+\dfrac{z^2}{N}}$$

Set $\hat{p}=0$ (zero observed tail events) and simplify directly:

$$p_{upper}\Big|_{\hat p = 0} = \frac{z^2/(2N) + z\cdot z/(2N)}{1+z^2/N} = \frac{z^2/N}{1+z^2/N} = \frac{z^2}{N+z^2} \tag{8.1}$$

This is strictly positive and shrinks as $1/N$, not as the trivially-zero
point estimate would suggest. Zero observed events out of, say, $N=200$
trials in the outer strip only bounds the true tail rate above by roughly
$z^2/N \approx 0.0135$ at 95% one-sided, not by zero.

**Annotation, the P0.4 bug.** The previous gate checked
$\hat{p}\le\text{threshold}$, a point estimate. With small outer-strip
sample sizes, zero observed events trivially passed that check while
providing almost no statistical evidence the true tail rate is actually
below the threshold. Gating on $(8.1)$ instead is the textbook fix for
exactly this failure mode: it forces the campaign to demonstrate the tail
is small with a sample size large enough to make the upper bound itself
small, rather than letting an under-sampled region pass by default.

---

## Summary Table

| Section | Result | Connects to |
|---|---|---|
| 1 | No Jacobi integral for eccentric binary | Why the slingshot mechanism exists at all |
| 2 | $e^2=1+(bv_\infty^2/\mu)^2$, $\chi=2\arcsin(1/e)$ | P0.1 deflection gate fix |
| 3 | $\Delta KE = W_\star+W_p$ exactly | `work_energy_closure_relative`, tolerance choice |
| 4 | COM-frame $\Delta\varepsilon$ boost-invariant | `galilean_invariance` gate, banning legacy metrics |
| 5 | $v_c^2=G M_{tot}/a$ | Cross-config (Quinn/Ortiz) comparability |
| 6 | $\hat W=2b_{max}N_{event}/N$ | Core estimand, dimensional reduction from 3D cross section |
| 7 | Wilson score interval derivation | All reported confidence intervals |
| 8 | $p_{upper}\to z^2/N$ as $\hat p\to0$ | P0.4 tail gate fix |
