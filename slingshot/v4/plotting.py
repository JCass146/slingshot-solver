"""Diagnostic figures for planar-width campaigns.

All plots read from the CSV/JSON artifacts written by campaign.run_campaign().
Nothing here depends on results.pkl or v3 data structures.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Optional

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    _MPL = True
except ImportError:
    _MPL = False

AU_KM = 1.495978707e8


def _require_mpl():
    if not _MPL:
        raise ImportError("matplotlib is required for plotting")


def _load(run_dir: Path):
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    samples = list(csv.DictReader((run_dir / "samples.csv").open(encoding="utf-8")))
    summary = list(csv.DictReader((run_dir / "width_summary.csv").open(encoding="utf-8")))
    # Also parse the stored config.yaml for system/orbit/sampling parameters
    try:
        import yaml
        cfg_raw = yaml.safe_load((run_dir / "config.yaml").read_text(encoding="utf-8"))
    except Exception:
        cfg_raw = {}
    return manifest, samples, summary, cfg_raw


def _savefig(fig, path: Path, name: str):
    fig.savefig(path / name, dpi=120, bbox_inches="tight")
    plt.close(fig)


def _meta(cfg_raw: dict) -> dict:
    """Extract commonly-needed physical parameters from stored config.yaml."""
    sys = cfg_raw.get("system", {})
    smp = cfg_raw.get("asymptotic_sampling", {})
    return {
        "name": sys.get("name", "Unknown system"),
        "b_max_au": float(smp.get("b_max_au", 1.0)),
        "boundary_radius_au": float(smp.get("boundary_radius_au", 5.0)),
        "star_radius_rsun": float(sys.get("star_radius_rsun", 1.0)),
        "planet_radius_rjup": float(sys.get("planet_radius_rjup", 1.0)),
    }


# ──────────────────────────────────────────────────────────────────────────────
# 1. Width vs v∞  (primary estimand)
# ──────────────────────────────────────────────────────────────────────────────

def plot_width_vs_vinf(run_dir: Path, thresholds: Optional[list] = None) -> Path:
    """Planar width W(v∞) with Wilson CI bands, one curve per energy threshold.

    Also shows per-seed points to visualise between-seed variance.
    """
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, _, summary, cfg_raw = _load(run_dir)
    b_max_au = _meta(cfg_raw)["b_max_au"]

    combined = [r for r in summary if r["scope"] == "combined"
                and r["statistic"] == "energy_threshold"]
    seed_rows = [r for r in summary if r["scope"] == "seed"
                 and r["statistic"] == "energy_threshold"]
    all_thresh = sorted({float(r["threshold"]) for r in combined})
    if thresholds is not None:
        all_thresh = [t for t in all_thresh if t in thresholds]

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(all_thresh)))

    for color, thresh in zip(colors, all_thresh):
        rows = [r for r in combined if float(r["threshold"]) == thresh]
        rows.sort(key=lambda r: float(r["v_inf_kms"]))
        vs = [float(r["v_inf_kms"]) for r in rows]
        ws = [float(r["width_km"]) / AU_KM for r in rows]
        lo = [float(r["width_low_km"]) / AU_KM for r in rows]
        hi = [float(r["width_high_km"]) / AU_KM for r in rows]
        ax.fill_between(vs, lo, hi, alpha=0.18, color=color)
        ax.plot(vs, ws, "-o", color=color, lw=1.8, ms=5,
                label=f"Δε/vc²>{thresh:.3g}")
        # Per-seed points for thresh=0 only (avoids clutter)
        if thresh == 0.0:
            srows = [r for r in seed_rows if float(r["threshold"]) == 0.0]
            for r in srows:
                ax.scatter(float(r["v_inf_kms"]),
                           float(r["width_km"]) / AU_KM,
                           color=color, alpha=0.4, s=14, zorder=3)

    ax.axhline(2.0 * b_max_au, color="gray", ls="--", lw=1,
               label=f"2b_max={2*b_max_au:.1f} AU (ceiling)")
    ax.set_xlabel("v∞ (km/s)", fontsize=11)
    ax.set_ylabel("Effective planar width W (AU)", fontsize=11)
    ax.set_title(f"{_meta(cfg_raw)['name']}\n"
                 "Planar-width estimates — 95% Wilson CI + per-seed points (Δε>0)")
    ax.legend(fontsize=8, ncol=2)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    out = run_dir / "width_vs_vinf.png"
    _savefig(fig, run_dir, "width_vs_vinf.png")
    return out


# ──────────────────────────────────────────────────────────────────────────────
# 2. Outcome fractions vs v∞
# ──────────────────────────────────────────────────────────────────────────────

def plot_outcome_fractions(run_dir: Path) -> Path:
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)

    vinfs = sorted({float(r["v_inf_kms"]) for r in samples})
    outcomes = ["escaped", "star_collision", "planet_collision", "time_limit", "integration_failed"]
    colors = {"escaped": "#2ca02c", "star_collision": "#d62728",
              "planet_collision": "#ff7f0e", "time_limit": "#9467bd",
              "integration_failed": "#8c8c8c"}

    fracs = {o: [] for o in outcomes}
    for v in vinfs:
        rows = [r for r in samples if float(r["v_inf_kms"]) == v]
        n = len(rows)
        for o in outcomes:
            fracs[o].append(sum(1 for r in rows if r["outcome"] == o) / n)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bottom = np.zeros(len(vinfs))
    for o in outcomes:
        vals = np.array(fracs[o])
        if vals.max() < 1e-6:
            continue
        ax.bar(range(len(vinfs)), vals, bottom=bottom, color=colors[o],
               label=o.replace("_", " "), alpha=0.85)
        bottom += vals

    ax.set_xticks(range(len(vinfs)))
    ax.set_xticklabels([f"{v:.0f}" for v in vinfs])
    ax.set_xlabel("v∞ (km/s)", fontsize=11)
    ax.set_ylabel("Fraction of samples", fontsize=11)
    ax.set_title(f"{_meta(cfg_raw)['name']}\nOutcome fractions vs v∞")
    ax.legend(fontsize=9, loc="upper right")
    ax.set_ylim(0, 1.05)
    fig.tight_layout()
    _savefig(fig, run_dir, "outcome_fractions.png")
    return run_dir / "outcome_fractions.png"


# ──────────────────────────────────────────────────────────────────────────────
# 3. Tail-support plot (event probability vs |b|/b_max)
# ──────────────────────────────────────────────────────────────────────────────

def plot_tail_support(run_dir: Path, v_inf_kms: float = None) -> Path:
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)
    b_max_km = _meta(cfg_raw)["b_max_au"] * AU_KM

    if v_inf_kms is None:
        vinfs = sorted({float(r["v_inf_kms"]) for r in samples})
        v_inf_kms = vinfs[len(vinfs) // 2]  # middle bin

    rows = [r for r in samples if float(r["v_inf_kms"]) == v_inf_kms]
    b_vals = np.array([float(r["impact_parameter_km"]) for r in rows])
    escaped = np.array([r["outcome"] == "escaped" for r in rows])

    bins = np.linspace(0, 1.0, 21)
    b_frac = np.abs(b_vals) / b_max_km
    event_rate = []
    bin_centers = []
    ci_lo, ci_hi = [], []
    from statistics import NormalDist
    z = NormalDist().inv_cdf(0.975)
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (b_frac >= lo) & (b_frac < hi)
        n = mask.sum()
        k = escaped[mask].sum()
        if n == 0:
            continue
        p = k / n
        denom = 1 + z**2 / n
        center = (p + z**2 / (2 * n)) / denom
        hw = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
        event_rate.append(center)
        ci_lo.append(max(0.0, center - hw))
        ci_hi.append(min(1.0, center + hw))
        bin_centers.append(0.5 * (lo + hi))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.fill_between(bin_centers, ci_lo, ci_hi, alpha=0.25, color="#1f77b4")
    ax.plot(bin_centers, event_rate, "-o", color="#1f77b4", ms=4)
    ax.axvline(0.9, color="red", ls="--", lw=1.2, label="Outer 10% strip")
    ax.set_xlabel("|b| / b_max", fontsize=11)
    ax.set_ylabel("Escape event probability", fontsize=11)
    ax.set_title(f"{_meta(cfg_raw)['name']}\n"
                 f"Tail-support: event rate vs |b|/b_max  (v∞={v_inf_kms:.0f} km/s)")
    ax.legend(fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0)
    fig.tight_layout()
    _savefig(fig, run_dir, "tail_support.png")
    return run_dir / "tail_support.png"


# ──────────────────────────────────────────────────────────────────────────────
# 4. Energy-gain ECDF by speed bin
# ──────────────────────────────────────────────────────────────────────────────

def plot_gain_ecdf(run_dir: Path) -> Path:
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)

    vinfs = sorted({float(r["v_inf_kms"]) for r in samples})
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(vinfs)))

    fig, ax = plt.subplots(figsize=(8, 5))
    for color, v in zip(colors, vinfs):
        gains = np.array([float(r["energy_gain_dimensionless"])
                          for r in samples
                          if float(r["v_inf_kms"]) == v
                          and r["outcome"] == "escaped"
                          and r["energy_gain_dimensionless"] not in ("nan", "")])
        if gains.size < 2:
            continue
        gains.sort()
        ecdf = np.arange(1, len(gains) + 1) / len(gains)
        ax.step(gains, ecdf, color=color, lw=1.5, label=f"{v:.0f} km/s")

    ax.axvline(0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("Δε / vc²  (dimensionless COM energy gain)", fontsize=11)
    ax.set_ylabel("ECDF (conditional on escape)", fontsize=11)
    ax.set_title(f"{_meta(cfg_raw)['name']}\n"
                 "COM energy-gain ECDF by v∞ — escaped trajectories only")
    ax.legend(fontsize=8, ncol=2)
    ax.set_xlim(left=min(-0.05, ax.get_xlim()[0]))
    ax.set_ylim(0, 1.02)
    fig.tight_layout()
    _savefig(fig, run_dir, "gain_ecdf.png")
    return run_dir / "gain_ecdf.png"


# ──────────────────────────────────────────────────────────────────────────────
# 5. Seed stability
# ──────────────────────────────────────────────────────────────────────────────

def plot_seed_stability(run_dir: Path, threshold: float = 0.0) -> Path:
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, _, summary, cfg_raw = _load(run_dir)

    seed_rows = [r for r in summary
                 if r["scope"] == "seed"
                 and r["statistic"] == "energy_threshold"
                 and abs(float(r["threshold"]) - threshold) < 1e-9]
    combined = [r for r in summary
                if r["scope"] == "combined"
                and r["statistic"] == "energy_threshold"
                and abs(float(r["threshold"]) - threshold) < 1e-9]

    vinfs = sorted({float(r["v_inf_kms"]) for r in seed_rows})
    seeds = sorted({r["seed"] for r in seed_rows})

    fig, ax = plt.subplots(figsize=(8, 4.5))
    colors = plt.cm.tab10(np.linspace(0, 0.9, len(seeds)))

    for color, seed in zip(colors, seeds):
        rows = [r for r in seed_rows if r["seed"] == seed]
        rows.sort(key=lambda r: float(r["v_inf_kms"]))
        vs = [float(r["v_inf_kms"]) for r in rows]
        ws = [float(r["width_km"]) / AU_KM for r in rows]
        ax.plot(vs, ws, "o--", color=color, lw=1.1, ms=5, alpha=0.7, label=f"seed {seed}")

    # Pooled Wilson CI
    combined.sort(key=lambda r: float(r["v_inf_kms"]))
    vs_c = [float(r["v_inf_kms"]) for r in combined]
    ws_c = [float(r["width_km"]) / AU_KM for r in combined]
    lo_c = [float(r["width_low_km"]) / AU_KM for r in combined]
    hi_c = [float(r["width_high_km"]) / AU_KM for r in combined]
    ax.fill_between(vs_c, lo_c, hi_c, alpha=0.15, color="black")
    ax.plot(vs_c, ws_c, "k-", lw=2.2, label="Pooled (Wilson CI)", zorder=5)

    ax.set_xlabel("v∞ (km/s)", fontsize=11)
    ax.set_ylabel("Planar width W (AU)", fontsize=11)
    ax.set_title(f"{_meta(cfg_raw)['name']}\n"
                 f"Seed stability — Δε/vc²>{threshold:.3g}")
    ax.legend(fontsize=8, ncol=3)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    _savefig(fig, run_dir, "seed_stability.png")
    return run_dir / "seed_stability.png"


# ──────────────────────────────────────────────────────────────────────────────
# 6. Sampling parameter distributions
# ──────────────────────────────────────────────────────────────────────────────

def plot_sampling_distributions(run_dir: Path) -> Path:
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)
    b_max_au = _meta(cfg_raw)["b_max_au"]

    b_vals = np.array([float(r["impact_parameter_km"]) / AU_KM for r in samples])
    dir_vals = np.degrees([float(r["incoming_direction_rad"]) for r in samples])
    ma_vals = np.degrees([float(r["binary_mean_anomaly_rad"]) for r in samples])
    v_vals = np.array([float(r["v_inf_kms"]) for r in samples])
    escaped = np.array([r["outcome"] == "escaped" for r in samples])

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))

    # Impact parameter
    ax = axes[0, 0]
    ax.hist(b_vals, bins=60, color="#aec7e8", label="All proposals", density=True)
    ax.hist(b_vals[escaped], bins=60, color="#1f77b4", alpha=0.6,
            label="Escaped", density=True)
    ax.set_xlabel("Impact parameter b (AU)")
    ax.set_ylabel("Density")
    ax.set_title("Signed impact parameter")
    ax.legend(fontsize=8)

    # Incoming direction
    ax = axes[0, 1]
    ax.hist(dir_vals, bins=60, color="#c5b0d5", label="All", density=True)
    ax.hist(dir_vals[escaped], bins=60, color="#9467bd", alpha=0.6,
            label="Escaped", density=True)
    ax.set_xlabel("Incoming direction (°)")
    ax.set_title("Incoming direction (uniform in [0,360°])")
    ax.legend(fontsize=8)

    # Binary mean anomaly
    ax = axes[1, 0]
    ax.hist(ma_vals, bins=60, color="#c7e9c0", label="All", density=True)
    ax.hist(ma_vals[escaped], bins=60, color="#2ca02c", alpha=0.6,
            label="Escaped", density=True)
    ax.set_xlabel("Binary mean anomaly (°)")
    ax.set_title("Binary orbital phase")
    ax.legend(fontsize=8)

    # v∞ bins
    ax = axes[1, 1]
    vinfs_all = sorted({v for v in v_vals})
    counts_all = [np.sum(v_vals == v) for v in vinfs_all]
    counts_esc = [np.sum(v_vals[escaped] == v) for v in vinfs_all]
    x = np.arange(len(vinfs_all))
    ax.bar(x - 0.2, counts_all, 0.35, color="#aec7e8", label="All")
    ax.bar(x + 0.2, counts_esc, 0.35, color="#1f77b4", alpha=0.8, label="Escaped")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:.0f}" for v in vinfs_all])
    ax.set_xlabel("v∞ (km/s)")
    ax.set_ylabel("Count")
    ax.set_title("Sample counts per speed bin")
    ax.legend(fontsize=8)

    name = _meta(cfg_raw)["name"]
    fig.suptitle(f"{name} — Proposal & acceptance distributions", fontsize=12)
    fig.tight_layout()
    _savefig(fig, run_dir, "sampling_distributions.png")
    return run_dir / "sampling_distributions.png"


# ──────────────────────────────────────────────────────────────────────────────
# 7. Work-energy closure diagnostics
# ──────────────────────────────────────────────────────────────────────────────

def plot_work_energy_diagnostics(run_dir: Path) -> Path:
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)

    closure = np.array([float(r["work_energy_closure_relative"])
                        for r in samples
                        if r.get("work_energy_closure_relative", "nan") not in ("nan", "")])
    w_star = np.array([float(r["work_star"])
                       for r in samples if r.get("work_star", "nan") not in ("nan", "")
                       and r["outcome"] == "escaped"])
    w_planet = np.array([float(r["work_planet"])
                         for r in samples if r.get("work_planet", "nan") not in ("nan", "")
                         and r["outcome"] == "escaped"])

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    # Closure residual distribution
    ax = axes[0]
    finite = closure[np.isfinite(closure)]
    ax.hist(np.log10(finite + 1e-18), bins=60, color="#9ecae1")
    tol_line = manifest["validation"]["quick"].get("gates", [])
    ax.set_xlabel("log₁₀(work-energy closure relative error)")
    ax.set_ylabel("Count")
    ax.set_title("Work-energy closure residuals")
    ax.axvline(np.log10(1e-4), color="red", ls="--", lw=1.2, label="tol=1e-4")
    ax.legend(fontsize=8)

    # Signed stellar vs planetary work (escaped only)
    ax = axes[1]
    if w_star.size > 0:
        ax.scatter(w_star / (w_star + w_planet + 1e-30),
                   w_planet / (w_star + w_planet + 1e-30),
                   s=3, alpha=0.3, color="#ff7f0e")
    ax.set_xlabel("Star work fraction")
    ax.set_ylabel("Planet work fraction")
    ax.set_title("Signed work fractions (escaped)")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)

    # Total work vs energy gain (escaped only)
    ax = axes[2]
    gains = np.array([float(r["energy_gain_dimensionless"])
                      for r in samples
                      if r["outcome"] == "escaped"
                      and r.get("energy_gain_dimensionless", "nan") not in ("nan", "")])
    ws = np.array([float(r["work_planet"])
                   for r in samples
                   if r["outcome"] == "escaped"
                   and r.get("work_planet", "nan") not in ("nan", "")])
    if gains.size > 0 and ws.size == gains.size:
        ax.scatter(gains, ws, s=3, alpha=0.25, color="#2ca02c")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("Δε/vc² (dimensionless energy gain)")
    ax.set_ylabel("Planet work (km²/s²)")
    ax.set_title("Planet work vs energy gain (escaped)")

    name = _meta(cfg_raw)["name"]
    fig.suptitle(f"{name} — Work-energy diagnostics", fontsize=12)
    fig.tight_layout()
    _savefig(fig, run_dir, "work_energy_diagnostics.png")
    return run_dir / "work_energy_diagnostics.png"


# ──────────────────────────────────────────────────────────────────────────────
# 8. Periapsis distributions
# ──────────────────────────────────────────────────────────────────────────────

def plot_periapsis_distributions(run_dir: Path) -> Path:
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    name = _meta(cfg_raw)["name"]
    star_r_km = (_meta(cfg_raw)["star_radius_rsun"]
                 * 695700.0)
    planet_r_km = (_meta(cfg_raw)["planet_radius_rjup"]
                   * 71492.0)

    for ax, body, field, body_r, color in [
        (axes[0], "Planet", "periapsis_planet_km", planet_r_km, "#ff7f0e"),
        (axes[1], "Star",   "periapsis_star_km",   star_r_km,   "#d62728"),
    ]:
        vals = np.array([float(r[field]) for r in samples
                         if r.get(field, "nan") not in ("nan", "")])
        escaped_mask = np.array([r["outcome"] == "escaped" for r in samples
                                 if r.get(field, "nan") not in ("nan", "")])
        ax.hist(vals / body_r, bins=80, color="#c7c7c7",
                label="All", density=True)
        ax.hist(vals[escaped_mask] / body_r, bins=80, color=color,
                alpha=0.65, label="Escaped", density=True)
        ax.axvline(1.0, color="black", ls="--", lw=1.2, label=f"Surface (1 R)")
        ax.set_xlabel(f"r_min / R_{body.lower()}")
        ax.set_ylabel("Density")
        ax.set_title(f"{body} closest approach")
        ax.legend(fontsize=8)
        ax.set_xscale("log")

    fig.suptitle(f"{name} — Periapsis distributions", fontsize=12)
    fig.tight_layout()
    _savefig(fig, run_dir, "periapsis_distributions.png")
    return run_dir / "periapsis_distributions.png"


# ──────────────────────────────────────────────────────────────────────────────
# 9. Impact-parameter phase map (conditional mean gain)
# ──────────────────────────────────────────────────────────────────────────────

def plot_phase_map(run_dir: Path, v_inf_kms: float = None) -> Path:
    """Conditional mean Δε/vc² over (b, binary_mean_anomaly) for one speed bin.

    Shows count panel alongside to make support visible.
    """
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)
    b_max_au = _meta(cfg_raw)["b_max_au"]

    if v_inf_kms is None:
        vinfs = sorted({float(r["v_inf_kms"]) for r in samples})
        # Pick the middle bin
        v_inf_kms = vinfs[len(vinfs) // 2]

    rows = [r for r in samples
            if float(r["v_inf_kms"]) == v_inf_kms
            and r["outcome"] == "escaped"
            and r.get("energy_gain_dimensionless", "nan") not in ("nan", "")]

    if not rows:
        return run_dir / "phase_map.png"

    b_vals = np.array([float(r["impact_parameter_km"]) / AU_KM for r in rows])
    ma_vals = np.degrees([float(r["binary_mean_anomaly_rad"]) for r in rows])
    gains = np.array([float(r["energy_gain_dimensionless"]) for r in rows])

    nb, nm = 30, 24
    b_edges = np.linspace(-b_max_au, b_max_au, nb + 1)
    m_edges = np.linspace(0, 360, nm + 1)

    gain_grid = np.full((nm, nb), np.nan)
    count_grid = np.zeros((nm, nb))
    for b, m, g in zip(b_vals, ma_vals, gains):
        ib = np.searchsorted(b_edges, b) - 1
        im = np.searchsorted(m_edges, m) - 1
        if 0 <= ib < nb and 0 <= im < nm:
            count_grid[im, ib] += 1
            if np.isnan(gain_grid[im, ib]):
                gain_grid[im, ib] = g
            else:
                gain_grid[im, ib] += g

    mask = count_grid > 0
    gain_grid[mask] /= count_grid[mask]

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    extent = [-b_max_au, b_max_au, 0, 360]

    im0 = axes[0].imshow(gain_grid, origin="lower", aspect="auto",
                          extent=extent, cmap="RdBu_r",
                          vmin=-np.nanpercentile(np.abs(gain_grid[mask]), 95),
                          vmax=np.nanpercentile(np.abs(gain_grid[mask]), 95))
    plt.colorbar(im0, ax=axes[0], label="Mean Δε/vc²")
    axes[0].set_xlabel("Impact parameter b (AU)")
    axes[0].set_ylabel("Binary mean anomaly (°)")
    axes[0].set_title(f"Conditional mean energy gain\n(v∞={v_inf_kms:.0f} km/s, escaped only)")

    im1 = axes[1].imshow(count_grid, origin="lower", aspect="auto",
                          extent=extent, cmap="Blues")
    plt.colorbar(im1, ax=axes[1], label="Escaped sample count")
    axes[1].set_xlabel("Impact parameter b (AU)")
    axes[1].set_ylabel("Binary mean anomaly (°)")
    axes[1].set_title("Sample support (count)")

    name = _meta(cfg_raw)["name"]
    fig.suptitle(f"{name} — Phase map: b × binary phase", fontsize=11)
    fig.tight_layout()
    _savefig(fig, run_dir, "phase_map.png")
    return run_dir / "phase_map.png"


# ──────────────────────────────────────────────────────────────────────────────
# 10. Deflection distribution
# ──────────────────────────────────────────────────────────────────────────────

def plot_deflection_distribution(run_dir: Path) -> Path:
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)

    vinfs = sorted({float(r["v_inf_kms"]) for r in samples})
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(vinfs)))

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for color, v in zip(colors, vinfs):
        v_str = f"{v:.10g}"
        deflections = [float(r["deflection_rad"])
                       for r in samples
                       if f"{float(r['v_inf_kms']):.10g}" == v_str
                       and r["outcome"] == "escaped"
                       and r.get("deflection_rad", "nan") not in ("nan", "")]
        if not deflections:
            continue
        ax.hist(np.degrees(deflections), bins=60, density=True, histtype="step",
                color=color, lw=1.5, label=f"{v:.0f} km/s")

    ax.set_xlabel("Deflection angle (°)", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_title(f"{_meta(cfg_raw)['name']}\n"
                 "Deflection distribution by v∞ (escaped trajectories)")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    _savefig(fig, run_dir, "deflection_distribution.png")
    return run_dir / "deflection_distribution.png"


# ──────────────────────────────────────────────────────────────────────────────
# 11. Collision-width overlay (compare escape vs collision)
# ──────────────────────────────────────────────────────────────────────────────

def plot_collision_vs_escape_width(run_dir: Path) -> Path:
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, _, summary, cfg_raw = _load(run_dir)

    esc = sorted([r for r in summary if r["scope"] == "combined"
                  and r["statistic"] == "energy_threshold"
                  and float(r["threshold"]) == 0.0],
                 key=lambda r: float(r["v_inf_kms"]))
    col = sorted([r for r in summary if r["scope"] == "combined"
                  and r["statistic"] == "collision"],
                 key=lambda r: float(r["v_inf_kms"]))

    vs_e = [float(r["v_inf_kms"]) for r in esc]
    ws_e = [float(r["width_km"]) / AU_KM for r in esc]
    lo_e = [float(r["width_low_km"]) / AU_KM for r in esc]
    hi_e = [float(r["width_high_km"]) / AU_KM for r in esc]

    vs_c = [float(r["v_inf_kms"]) for r in col]
    ws_c = [float(r["width_km"]) / AU_KM for r in col]
    lo_c = [float(r["width_low_km"]) / AU_KM for r in col]
    hi_c = [float(r["width_high_km"]) / AU_KM for r in col]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.fill_between(vs_e, lo_e, hi_e, alpha=0.18, color="#2ca02c")
    ax.plot(vs_e, ws_e, "-o", color="#2ca02c", lw=2, ms=5, label="Escape (Δε>0)")
    ax.fill_between(vs_c, lo_c, hi_c, alpha=0.18, color="#d62728")
    ax.plot(vs_c, ws_c, "-s", color="#d62728", lw=2, ms=5, label="Any collision")
    ax.set_xlabel("v∞ (km/s)", fontsize=11)
    ax.set_ylabel("Effective planar width W (AU)", fontsize=11)
    ax.set_title(f"{_meta(cfg_raw)['name']}\n"
                 "Escape vs collision widths — 95% Wilson CI")
    ax.legend(fontsize=10)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    _savefig(fig, run_dir, "collision_vs_escape.png")
    return run_dir / "collision_vs_escape.png"


# ──────────────────────────────────────────────────────────────────────────────
# 12. Parameter correlations  (periapsis / deflection / star proximity vs gain)
# ──────────────────────────────────────────────────────────────────────────────

def plot_parameter_correlations(run_dir: Path) -> Path:
    """Scatter plots of key trajectory observables vs COM energy gain.

    Each panel shows escaped trajectories only, coloured by v∞.
    This replaces the v3 parameter_correlations_* figures using current scientific metrics.
    """
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)

    escaped = [r for r in samples
               if r["outcome"] == "escaped"
               and r.get("energy_gain_dimensionless", "nan") not in ("nan", "")]

    vinfs = sorted({float(r["v_inf_kms"]) for r in escaped})
    cmap = plt.cm.plasma
    norm = plt.Normalize(min(vinfs), max(vinfs))

    gain = np.array([float(r["energy_gain_dimensionless"]) for r in escaped])
    peri_p = np.array([float(r["periapsis_planet_km"]) / AU_KM for r in escaped])
    peri_s = np.array([float(r["periapsis_star_km"]) / AU_KM for r in escaped])
    defl = np.degrees([float(r["deflection_rad"]) for r in escaped])
    v_col = np.array([float(r["v_inf_kms"]) for r in escaped])

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    sc_kwargs = dict(s=4, alpha=0.35, cmap=cmap, norm=norm)

    axes[0].scatter(peri_p, gain, c=v_col, **sc_kwargs)
    axes[0].set_xlabel("Periapsis (planet) / AU", fontsize=10)
    axes[0].set_ylabel("Δε / vc²", fontsize=10)
    axes[0].set_title("Energy gain vs planet periapsis")
    axes[0].axhline(0, color="gray", lw=0.6)
    axes[0].set_xscale("log")

    axes[1].scatter(np.abs(defl), gain, c=v_col, **sc_kwargs)
    axes[1].set_xlabel("|Deflection| (°)", fontsize=10)
    axes[1].set_ylabel("Δε / vc²", fontsize=10)
    axes[1].set_title("Energy gain vs deflection angle")
    axes[1].axhline(0, color="gray", lw=0.6)

    sc = axes[2].scatter(peri_s, np.abs(defl), c=v_col, **sc_kwargs)
    axes[2].set_xlabel("Periapsis (star) / AU", fontsize=10)
    axes[2].set_ylabel("|Deflection| (°)", fontsize=10)
    axes[2].set_title("Deflection vs star proximity")
    axes[2].set_xscale("log")

    plt.colorbar(sc, ax=axes, label="v∞ (km/s)", shrink=0.8)
    fig.suptitle(f"{_meta(cfg_raw)['name']} — Parameter correlations (escaped)", fontsize=11)
    fig.subplots_adjust(right=0.88)
    _savefig(fig, run_dir, "parameter_correlations.png")
    return run_dir / "parameter_correlations.png"


# ──────────────────────────────────────────────────────────────────────────────
# 13. Velocity phase space  (initial vs final COM speed, coloured by gain)
# ──────────────────────────────────────────────────────────────────────────────

def plot_velocity_phase_space(run_dir: Path) -> Path:
    """Initial vs final COM-frame speed, coloured by dimensionless energy gain.

    Replaces the v3 velocity_phase_space_vx_vy figure using correct COM-frame
    speeds (boost-invariant).  The radial-normal panel is intentionally omitted
    because the v3 version had an incorrect radial basis (star-planet axis rather
    than particle-planet axis).
    """
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)

    escaped = [r for r in samples
               if r["outcome"] == "escaped"
               and r.get("initial_speed_com", "nan") not in ("nan", "")
               and r.get("final_speed_com", "nan") not in ("nan", "")
               and r.get("energy_gain_dimensionless", "nan") not in ("nan", "")]

    if not escaped:
        return run_dir / "velocity_phase_space.png"

    v_in = np.array([float(r["initial_speed_com"]) for r in escaped])
    v_out = np.array([float(r["final_speed_com"]) for r in escaped])
    gain = np.array([float(r["energy_gain_dimensionless"]) for r in escaped])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: v_in vs v_out coloured by gain
    vmax = np.percentile(np.abs(gain), 99)
    sc = axes[0].scatter(v_in, v_out, c=gain, s=3, alpha=0.35,
                         cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    plt.colorbar(sc, ax=axes[0], label="Δε / vc²")
    lim = max(v_in.max(), v_out.max()) * 1.05
    axes[0].plot([0, lim], [0, lim], "k--", lw=0.8, label="v_out = v_in")
    axes[0].set_xlabel("Initial COM speed (km/s)", fontsize=10)
    axes[0].set_ylabel("Final COM speed (km/s)", fontsize=10)
    axes[0].set_title("COM-frame speed: initial vs final")
    axes[0].legend(fontsize=8)

    # Panel 2: Δv_com (= v_out - v_in) distribution by v∞
    vinfs = sorted({float(r["v_inf_kms"]) for r in escaped})
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(vinfs)))
    for color, v in zip(colors, vinfs):
        v_str = f"{v:.10g}"
        dv = [float(r["delta_speed_com"])
              for r in escaped
              if f"{float(r['v_inf_kms']):.10g}" == v_str]
        if dv:
            axes[1].hist(dv, bins=50, density=True, histtype="step",
                         color=color, lw=1.5, label=f"{v:.0f} km/s")
    axes[1].axvline(0, color="gray", lw=0.8)
    axes[1].set_xlabel("Δ|v|_COM (km/s)  [final − initial COM speed]", fontsize=10)
    axes[1].set_ylabel("Density", fontsize=10)
    axes[1].set_title("COM speed-change distribution by v∞")
    axes[1].legend(fontsize=7, ncol=2)

    fig.suptitle(f"{_meta(cfg_raw)['name']} — COM-frame velocity phase space", fontsize=11)
    fig.tight_layout()
    _savefig(fig, run_dir, "velocity_phase_space.png")
    return run_dir / "velocity_phase_space.png"


# ──────────────────────────────────────────────────────────────────────────────
# 14. Top-N candidate ranking  (replaces candidate_ranking_* figures)
# ──────────────────────────────────────────────────────────────────────────────

def _float(row: dict, key: str, default: float = np.nan) -> float:
    value = row.get(key, default)
    if value in ("", None, "nan", "NaN"):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _candidate_top_n(cfg_raw: dict, default: int = 30) -> int:
    try:
        return int(cfg_raw.get("candidate_diagnostics", {}).get("top_n", default))
    except (TypeError, ValueError):
        return default


def _candidate_track_top_n(cfg_raw: dict, default: int = 10) -> int:
    try:
        return int(cfg_raw.get("candidate_diagnostics", {}).get("trajectory_top_n", default))
    except (TypeError, ValueError):
        return default


def _ensure_top_candidates(run_dir: Path, top_n: int = 30) -> list[dict]:
    from .candidates import load_top_candidates, write_top_candidates_for_run

    run_dir = Path(run_dir)
    rows = load_top_candidates(run_dir)
    if rows:
        return rows[:top_n]
    return write_top_candidates_for_run(run_dir, top_n=top_n)


def _integrate_candidate_row(run_dir: Path, row: dict):
    from .config import load_config
    from .dynamics import init_binary_barycentric, integrate_encounter
    from .sampling import state_at_inbound_boundary
    from .validation import physical_values
    from ..constants import G_KM

    config = load_config(run_dir / "config.yaml")
    values = physical_values(config)
    total_mu = G_KM * (values["star_mass_kg"] + values["planet_mass_kg"])
    pos, vel = state_at_inbound_boundary(
        _float(row, "v_inf_kms"),
        _float(row, "impact_parameter_km"),
        _float(row, "incoming_direction_rad"),
        values["boundary_radius_km"],
        total_mu,
    )
    binary = init_binary_barycentric(
        values["semi_major_axis_km"],
        config.orbit.eccentricity,
        _float(row, "binary_mean_anomaly_rad"),
        config.orbit.argument_periapsis_rad,
        values["star_mass_kg"],
        values["planet_mass_kg"],
        config.orbit.prograde,
        config.system.bulk_velocity_x_kms,
        config.system.bulk_velocity_y_kms,
    )
    bulk = np.array([config.system.bulk_velocity_x_kms, config.system.bulk_velocity_y_kms])
    initial_state = np.concatenate([binary, pos, vel + bulk, np.zeros(2)])
    integration = integrate_encounter(
        initial_state,
        values["star_mass_kg"],
        values["planet_mass_kg"],
        values["star_radius_km"],
        values["planet_radius_km"],
        values["boundary_radius_km"],
        config.numerical.max_time_sec,
        method=config.numerical.method,
        rtol=config.numerical.rtol,
        atol=config.numerical.atol,
        softening_km=config.numerical.softening_km,
        max_step_sec=config.numerical.max_step_sec,
    )
    return config, values, integration


def _unavailable_figure(run_dir: Path, filename: str, title: str, message: str) -> Path:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.text(0.5, 0.55, title, ha="center", va="center", fontsize=13, weight="bold")
    ax.text(0.5, 0.42, message, ha="center", va="center", fontsize=10)
    ax.set_axis_off()
    fig.tight_layout()
    _savefig(fig, run_dir, filename)
    return run_dir / filename


def _candidate_ranking_data(run_dir: Path, cfg_raw: dict, top_n: Optional[int] = None) -> dict:
    top_n = _candidate_top_n(cfg_raw) if top_n is None else int(top_n)
    top = _ensure_top_candidates(run_dir, top_n=top_n)
    if not top:
        return {"top": []}

    v_vals = np.array([_float(r, "v_inf_kms") for r in top])
    vmin = np.nanmin(v_vals)
    vmax = np.nanmax(v_vals)
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin, vmax = 0.0, max(1.0, vmin + 1.0 if np.isfinite(vmin) else 1.0)
    closure_raw = np.array([_float(r, "work_energy_closure_relative") for r in top])
    return {
        "top": top,
        "ranks": np.array([int(float(r.get("rank", i + 1))) for i, r in enumerate(top)]),
        "gain": np.array([_float(r, "energy_gain_dimensionless") for r in top]),
        "delta_e": np.array([_float(r, "delta_specific_energy_com") for r in top]),
        "turning": np.array([_float(r, "turning_quadratic") for r in top]),
        "deflection": np.abs(np.array([_float(r, "deflection_deg") for r in top])),
        "periapsis_planet_au": np.array([_float(r, "periapsis_planet_km") / AU_KM for r in top]),
        "closure": np.where(np.isfinite(closure_raw), np.maximum(closure_raw, 1e-16), np.nan),
        "v_vals": v_vals,
        "norm": plt.Normalize(vmin, vmax),
    }


_CANDIDATE_PANEL_FILENAMES = [
    "candidate_ranking_gain.png",
    "candidate_ranking_delta_energy.png",
    "candidate_ranking_turning.png",
    "candidate_ranking_deflection.png",
    "candidate_ranking_periapsis.png",
    "candidate_ranking_closure.png",
]


def _candidate_panel_specs(data: dict) -> list[dict]:
    return [
        {
            "filename": "candidate_ranking_gain.png",
            "title": "COM energy gain",
            "ylabel": "Delta epsilon / vc^2",
            "values": data["gain"],
            "zero_line": True,
        },
        {
            "filename": "candidate_ranking_delta_energy.png",
            "title": "True specific-energy change",
            "ylabel": "Delta epsilon (km^2/s^2)",
            "values": data["delta_e"],
            "zero_line": True,
        },
        {
            "filename": "candidate_ranking_turning.png",
            "title": "Turning diagnostic only",
            "ylabel": "0.5 |Delta V|^2",
            "values": data["turning"],
        },
        {
            "filename": "candidate_ranking_deflection.png",
            "title": "Deflection",
            "ylabel": "|Deflection| (deg)",
            "values": data["deflection"],
        },
        {
            "filename": "candidate_ranking_periapsis.png",
            "title": "Encounter depth",
            "ylabel": "Planet periapsis (AU)",
            "values": data["periapsis_planet_au"],
            "yscale": "log",
        },
        {
            "filename": "candidate_ranking_closure.png",
            "title": "Work-energy closure",
            "ylabel": "Closure relative error",
            "values": data["closure"],
            "yscale": "log",
        },
    ]


def _plot_candidate_panel(run_dir: Path, data: dict, spec: dict, cfg_raw: dict) -> Path:
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    sc = ax.scatter(
        data["ranks"],
        spec["values"],
        c=data["v_vals"],
        cmap=plt.cm.plasma,
        norm=data["norm"],
        s=42,
        zorder=3,
    )
    if spec.get("zero_line"):
        ax.axhline(0, color="gray", lw=0.7)
    if spec.get("yscale"):
        ax.set_yscale(spec["yscale"])
    ax.set_xlabel("Rank by energy gain")
    ax.set_ylabel(spec["ylabel"])
    ax.set_title(f"{_meta(cfg_raw)['name']} - {spec['title']}")
    ax.grid(True, alpha=0.2)
    plt.colorbar(sc, ax=ax, label="v_inf (km/s)")
    fig.tight_layout()
    _savefig(fig, run_dir, spec["filename"])
    return run_dir / spec["filename"]


def plot_candidate_ranking(run_dir: Path, top_n: Optional[int] = None) -> Path:
    """Top-candidate diagnostics using current v4 energy metrics."""
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, _, _, cfg_raw = _load(run_dir)
    data = _candidate_ranking_data(run_dir, cfg_raw, top_n)
    top = data["top"]
    if not top:
        return _unavailable_figure(
            run_dir,
            "candidate_ranking.png",
            "No eligible top candidates",
            "No escaped solver-successful samples with finite COM energy gain were available.",
        )

    cmap = plt.cm.plasma
    sc_kw = dict(c=data["v_vals"], cmap=cmap, norm=data["norm"], s=34, zorder=3)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()
    sc = axes[0].scatter(data["ranks"], data["gain"], **sc_kw)
    axes[0].set_ylabel("Delta epsilon / vc^2")
    axes[0].set_title("COM energy gain")
    axes[0].axhline(0, color="gray", lw=0.6)

    axes[1].scatter(data["ranks"], data["delta_e"], **sc_kw)
    axes[1].set_ylabel("Delta epsilon (km^2/s^2)")
    axes[1].set_title("True specific-energy change")
    axes[1].axhline(0, color="gray", lw=0.6)

    axes[2].scatter(data["ranks"], data["turning"], **sc_kw)
    axes[2].set_ylabel("0.5 |Delta V|^2")
    axes[2].set_title("Turning diagnostic only")

    axes[3].scatter(data["ranks"], data["deflection"], **sc_kw)
    axes[3].set_ylabel("|Deflection| (deg)")
    axes[3].set_title("Deflection")

    axes[4].scatter(data["ranks"], data["periapsis_planet_au"], **sc_kw)
    axes[4].set_yscale("log")
    axes[4].set_ylabel("Planet periapsis (AU)")
    axes[4].set_title("Encounter depth")

    axes[5].scatter(data["ranks"], data["closure"], **sc_kw)
    axes[5].set_yscale("log")
    axes[5].set_ylabel("Closure relative error")
    axes[5].set_title("Work-energy closure")

    for ax in axes:
        ax.set_xlabel("Rank by energy gain")
        ax.grid(True, alpha=0.2)

    plt.colorbar(sc, ax=axes, label="v_inf (km/s)", shrink=0.85)
    fig.suptitle(
        f"{_meta(cfg_raw)['name']} - Top-{len(top)} exploratory candidates",
        fontsize=12,
    )
    fig.subplots_adjust(right=0.9, hspace=0.35, wspace=0.35)
    _savefig(fig, run_dir, "candidate_ranking.png")
    return run_dir / "candidate_ranking.png"


def plot_candidate_ranking_panels(run_dir: Path, top_n: Optional[int] = None) -> list[Path]:
    """Write each top-candidate ranking diagnostic as its own PNG."""
    _require_mpl()
    run_dir = Path(run_dir)
    _, _, _, cfg_raw = _load(run_dir)
    data = _candidate_ranking_data(run_dir, cfg_raw, top_n)
    if not data["top"]:
        return [
            _unavailable_figure(
                run_dir,
                filename,
                "No eligible top candidates",
                "No escaped solver-successful samples with finite COM energy gain were available.",
            )
            for filename in _CANDIDATE_PANEL_FILENAMES
        ]
    return [
        _plot_candidate_panel(run_dir, data, spec, cfg_raw)
        for spec in _candidate_panel_specs(data)
    ]


def _barycentric_tracks(sol, values: dict) -> dict:
    total_mass = values["star_mass_kg"] + values["planet_mass_kg"]
    com_x = (
        values["star_mass_kg"] * sol.y[0, :]
        + values["planet_mass_kg"] * sol.y[4, :]
    ) / total_mass
    com_y = (
        values["star_mass_kg"] * sol.y[1, :]
        + values["planet_mass_kg"] * sol.y[5, :]
    ) / total_mass
    return {
        "star_x": (sol.y[0, :] - com_x) / AU_KM,
        "star_y": (sol.y[1, :] - com_y) / AU_KM,
        "planet_x": (sol.y[4, :] - com_x) / AU_KM,
        "planet_y": (sol.y[5, :] - com_y) / AU_KM,
        "test_x": (sol.y[8, :] - com_x) / AU_KM,
        "test_y": (sol.y[9, :] - com_y) / AU_KM,
    }


def _set_equal_limits(ax, xs: list[np.ndarray], ys: list[np.ndarray], min_span: float = 1e-3) -> None:
    x_arrays = [arr[np.isfinite(arr)] for arr in xs if np.size(arr) and np.any(np.isfinite(arr))]
    y_arrays = [arr[np.isfinite(arr)] for arr in ys if np.size(arr) and np.any(np.isfinite(arr))]
    if not x_arrays or not y_arrays:
        return
    x = np.concatenate(x_arrays)
    y = np.concatenate(y_arrays)
    x_mid = 0.5 * (np.nanmin(x) + np.nanmax(x))
    y_mid = 0.5 * (np.nanmin(y) + np.nanmax(y))
    span = max(np.nanmax(x) - np.nanmin(x), np.nanmax(y) - np.nanmin(y), min_span)
    half = 0.55 * span
    ax.set_xlim(x_mid - half, x_mid + half)
    ax.set_ylim(y_mid - half, y_mid + half)


def plot_best_candidate(run_dir: Path) -> Path:
    """Plot the strongest observed candidate in the binary barycentric frame."""
    _require_mpl()
    run_dir = Path(run_dir)
    _, _, _, cfg_raw = _load(run_dir)
    candidates = _ensure_top_candidates(run_dir, top_n=1)
    if not candidates:
        return _unavailable_figure(
            run_dir,
            "best_candidate.png",
            "No best candidate available",
            "No eligible escaped finite-gain candidate was available for plotting.",
        )
    row = candidates[0]
    try:
        _, values, integration = _integrate_candidate_row(run_dir, row)
    except Exception as exc:
        return _unavailable_figure(
            run_dir,
            "best_candidate.png",
            "Best candidate plot unavailable",
            f"Re-integration failed: {exc}",
        )

    tracks = _barycentric_tracks(integration.solution, values)

    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    ax.plot(tracks["test_x"], tracks["test_y"], color="#1f77b4", lw=1.5, label="Test particle")
    ax.scatter(tracks["test_x"][0], tracks["test_y"][0], marker="o", s=36, color="#2ca02c", label="Start")
    ax.scatter(tracks["test_x"][-1], tracks["test_y"][-1], marker="x", s=48, color="#d62728", label="End")
    ax.plot(tracks["star_x"], tracks["star_y"], color="#e15759", lw=1.0, alpha=0.75, label="Star")
    ax.plot(tracks["planet_x"], tracks["planet_y"], color="#f28e2b", lw=1.0, alpha=0.85, label="Planet")
    ax.scatter([0], [0], s=52, color="#222222", marker="+", zorder=5, label="Binary barycenter")
    ax.set_aspect("equal")
    ax.set_xlabel("Barycentric x (AU)")
    ax.set_ylabel("Barycentric y (AU)")
    ax.grid(True, alpha=0.25)
    _set_equal_limits(
        ax,
        [tracks["test_x"], tracks["star_x"], tracks["planet_x"]],
        [tracks["test_y"], tracks["star_y"], tracks["planet_y"]],
    )

    details = (
        f"rank={row.get('rank', '1')}  v_inf={_float(row, 'v_inf_kms'):.0f} km/s\n"
        f"gain={_float(row, 'energy_gain_dimensionless'):.4g}  "
        f"Delta eps={_float(row, 'delta_specific_energy_com'):.4g} km^2/s^2\n"
        f"deflection={_float(row, 'deflection_deg'):.3g} deg  "
        f"periapsis={_float(row, 'periapsis_planet_km')/AU_KM:.4g} AU\n"
        "binary barycentric frame; finite-sample example, not a converged optimum"
    )
    ax.text(
        0.02,
        0.98,
        details,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.82),
    )
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title(f"{_meta(cfg_raw)['name']} - Strongest observed candidate")
    fig.tight_layout()
    _savefig(fig, run_dir, "best_candidate.png")
    return run_dir / "best_candidate.png"


def plot_trajectory_tracks(run_dir: Path, top_n: Optional[int] = None) -> Path:
    """Re-integrate ranked top candidates and plot trajectories in barycentric frame."""
    _require_mpl()
    run_dir = Path(run_dir)
    _, _, _, cfg_raw = _load(run_dir)
    top_n = _candidate_track_top_n(cfg_raw) if top_n is None else int(top_n)
    candidates = _ensure_top_candidates(run_dir, top_n=top_n)
    if not candidates:
        return _unavailable_figure(
            run_dir,
            "trajectory_tracks.png",
            "No trajectory tracks available",
            "No eligible top candidates were available for re-integration.",
        )

    gain_vals = np.array([_float(r, "energy_gain_dimensionless") for r in candidates])
    gain_min, gain_max = np.nanmin(gain_vals), np.nanmax(gain_vals)
    if gain_min == gain_max:
        gain_max = gain_min + 1.0
    fig, ax = plt.subplots(figsize=(9, 9))
    cmap = plt.cm.RdYlGn
    norm = plt.Normalize(gain_min, gain_max)
    context_tracks = None
    limit_xs = []
    limit_ys = []

    for row in candidates:
        try:
            _, values, integration = _integrate_candidate_row(run_dir, row)
            tracks = _barycentric_tracks(integration.solution, values)
            if context_tracks is None:
                context_tracks = tracks
            color = cmap(norm(_float(row, "energy_gain_dimensionless")))
            ax.plot(tracks["test_x"], tracks["test_y"], lw=0.9, alpha=0.75, color=color)
            limit_xs.append(tracks["test_x"])
            limit_ys.append(tracks["test_y"])
        except Exception:
            continue

    if context_tracks is not None:
        ax.plot(context_tracks["star_x"], context_tracks["star_y"], color="#e15759", lw=1.0, alpha=0.8, label="Star path (rank 1)")
        ax.plot(context_tracks["planet_x"], context_tracks["planet_y"], color="#f28e2b", lw=1.0, alpha=0.9, label="Planet path (rank 1)")
        limit_xs.extend([context_tracks["star_x"], context_tracks["planet_x"]])
        limit_ys.extend([context_tracks["star_y"], context_tracks["planet_y"]])
    ax.scatter([0], [0], s=60, color="#222222", marker="+", zorder=5, label="Binary barycenter")
    _set_equal_limits(ax, limit_xs, limit_ys)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label="Delta epsilon / vc^2")
    ax.set_xlabel("Barycentric x (AU)")
    ax.set_ylabel("Barycentric y (AU)")
    ax.set_title(f"{_meta(cfg_raw)['name']} - Top-{len(candidates)} candidate tracks")
    ax.legend(fontsize=9)
    ax.set_aspect("equal")
    fig.tight_layout()
    _savefig(fig, run_dir, "trajectory_tracks.png")
    return run_dir / "trajectory_tracks.png"

def _pareto_mask(objectives: np.ndarray) -> np.ndarray:
    """Return boolean mask of Pareto-non-dominated points.

    Parameters
    ----------
    objectives:
        Array of shape (N, M) where each row is one solution and each column
        is one objective to MAXIMISE.
    """
    n = objectives.shape[0]
    dominated = np.zeros(n, dtype=bool)
    for i in range(n):
        if dominated[i]:
            continue
        # point i is dominated if some j dominates it on all objectives
        dominated[i] = np.any(
            np.all(objectives >= objectives[i], axis=1)
            & np.any(objectives > objectives[i], axis=1)
        )
    return ~dominated


def plot_pareto_front(run_dir: Path) -> Path:
    """Two-panel Pareto front using current scientific metrics.

    Panel 1: energy gain (↑) vs planet periapsis (↓ = closer = higher risk).
             The Pareto front marks trajectories where you cannot improve gain
             without accepting a closer planetary approach.

    Panel 2: energy gain (↑) vs |deflection| (↑).
             Shows whether high energy gain and high deflection are jointly
             achievable, or whether they trade off against each other.

    Both panels show all escaped trajectories as background scatter and the
    Pareto-optimal subset as highlighted points, coloured by v∞.

    Note: the v3 Pareto used ``energy_from_planet_orbit`` (Spearman ρ ≈ 0.06
    with true energy change) and ``bary_delta_v_pct`` (scalar speed at finite
    unmatched endpoints).  Those metrics are intentionally replaced here.
    """
    _require_mpl()
    run_dir = Path(run_dir)
    manifest, samples, _, cfg_raw = _load(run_dir)

    escaped = [r for r in samples
               if r["outcome"] == "escaped"
               and r.get("energy_gain_dimensionless", "nan") not in ("nan", "")
               and r.get("periapsis_planet_km", "nan") not in ("nan", "")
               and r.get("deflection_rad", "nan") not in ("nan", "")]
    if not escaped:
        return run_dir / "pareto_front.png"

    candidate_corrections = {}
    for row in _ensure_top_candidates(run_dir, top_n=_candidate_top_n(cfg_raw)):
        key = (
            str(row.get("seed", "")),
            str(row.get("sample_index", "")),
            str(row.get("v_inf_kms", "")),
        )
        candidate_corrections[key] = _float(row, "periapsis_planet_km")

    def _corrected_planet_periapsis(row: dict) -> float:
        key = (
            str(row.get("seed", "")),
            str(row.get("sample_index", "")),
            str(row.get("v_inf_kms", "")),
        )
        return candidate_corrections.get(key, float(row["periapsis_planet_km"]))

    gain = np.array([float(r["energy_gain_dimensionless"]) for r in escaped])
    peri = np.array([_corrected_planet_periapsis(r) / AU_KM for r in escaped])
    defl = np.abs(np.degrees([float(r["deflection_rad"]) for r in escaped]))
    v_col = np.array([float(r["v_inf_kms"]) for r in escaped])

    vinfs = sorted({float(r["v_inf_kms"]) for r in escaped})
    cmap = plt.cm.plasma
    norm = plt.Normalize(min(vinfs), max(vinfs))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ── Panel 1: gain (max) vs periapsis (max negative ↔ min positive) ──
    ax = axes[0]
    # Pareto objectives: maximise gain, maximise (−periapsis) ≡ minimise periapsis
    obj1 = np.column_stack([gain, -peri])
    pareto1 = _pareto_mask(obj1)

    ax.scatter(peri[~pareto1], gain[~pareto1], c=v_col[~pareto1],
               cmap=cmap, norm=norm, s=4, alpha=0.2, label="All escaped")
    sc = ax.scatter(peri[pareto1], gain[pareto1], c=v_col[pareto1],
                    cmap=cmap, norm=norm, s=30, alpha=0.9, zorder=5,
                    edgecolors="k", linewidths=0.4, label=f"Pareto front ({pareto1.sum()})")
    # Connect Pareto points with a step line (sort by periapsis)
    order = np.argsort(peri[pareto1])
    ax.step(peri[pareto1][order], gain[pareto1][order],
            color="black", lw=1.0, alpha=0.5, where="post")
    ax.axhline(0, color="gray", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("Planet periapsis (AU)  [log scale, closer →  left]", fontsize=10)
    ax.set_ylabel("Δε / vc²  (dimensionless energy gain)", fontsize=10)
    ax.set_title("Pareto: gain vs planet periapsis\n(top-left = more gain at closer approach)")
    ax.legend(fontsize=8)

    # ── Panel 2: gain (max) vs |deflection| (max) ──
    ax = axes[1]
    obj2 = np.column_stack([gain, defl])
    pareto2 = _pareto_mask(obj2)

    ax.scatter(defl[~pareto2], gain[~pareto2], c=v_col[~pareto2],
               cmap=cmap, norm=norm, s=4, alpha=0.2)
    sc = ax.scatter(defl[pareto2], gain[pareto2], c=v_col[pareto2],
                    cmap=cmap, norm=norm, s=30, alpha=0.9, zorder=5,
                    edgecolors="k", linewidths=0.4,
                    label=f"Pareto front ({pareto2.sum()})")
    order2 = np.argsort(defl[pareto2])
    ax.step(defl[pareto2][order2], gain[pareto2][order2],
            color="black", lw=1.0, alpha=0.5, where="post")
    ax.axhline(0, color="gray", lw=0.6)
    ax.set_xlabel("|Deflection angle|  (°)", fontsize=10)
    ax.set_ylabel("Δε / vc²  (dimensionless energy gain)", fontsize=10)
    ax.set_title("Pareto: gain vs deflection\n(top-right = more gain and more deflection)")
    ax.legend(fontsize=8)

    plt.colorbar(sc, ax=axes, label="v∞ (km/s)", shrink=0.8)
    fig.suptitle(f"{_meta(cfg_raw)['name']} — Pareto fronts (escaped trajectories)", fontsize=11)
    fig.subplots_adjust(right=0.88)
    _savefig(fig, run_dir, "pareto_front.png")
    return run_dir / "pareto_front.png"


# ──────────────────────────────────────────────────────────────────────────────
# Master entry point
# ──────────────────────────────────────────────────────────────────────────────

def generate_all_plots(run_dir: str | Path, verbose: bool = True) -> list[Path]:
    """Generate all diagnostic figures for a completed campaign run.

    Parameters
    ----------
    run_dir:
        Path to a run directory containing manifest.json, samples.csv,
        and width_summary.csv.
    verbose:
        Print progress messages.

    Returns
    -------
    List of paths to the generated PNG files.
    """
    _require_mpl()
    run_dir = Path(run_dir)
    if not (run_dir / "manifest.json").exists():
        raise FileNotFoundError(f"No manifest.json in {run_dir}")

    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 4:
        raise ValueError("generate_all_plots requires a current schema run directory")

    candidate_enabled = True
    candidate_top_n = 30
    candidate_track_top_n = 10
    try:
        from .candidates import write_top_candidates_for_run
        from .config import load_config

        config = load_config(run_dir / "config.yaml")
        candidate_enabled = bool(config.candidate_diagnostics.enabled)
        candidate_top_n = config.candidate_diagnostics.top_n
        candidate_track_top_n = config.candidate_diagnostics.trajectory_top_n
        if candidate_enabled:
            write_top_candidates_for_run(run_dir, top_n=candidate_top_n)
    except Exception as exc:
        if verbose:
            print(f"  Candidate table refresh skipped ({exc})")

    _, samples, _, cfg_raw = _load(run_dir)
    vinfs = sorted({float(r["v_inf_kms"]) for r in samples})
    mid_v = vinfs[len(vinfs) // 2]

    plots = [
        ("Width vs v-inf",            lambda: plot_width_vs_vinf(run_dir)),
        ("Outcome fractions",         lambda: plot_outcome_fractions(run_dir)),
        ("Tail support",              lambda: plot_tail_support(run_dir, mid_v)),
        ("Gain ECDF",                 lambda: plot_gain_ecdf(run_dir)),
        ("Seed stability",            lambda: plot_seed_stability(run_dir)),
        ("Sampling distributions",    lambda: plot_sampling_distributions(run_dir)),
        ("Work-energy diagnostics",   lambda: plot_work_energy_diagnostics(run_dir)),
        ("Periapsis distributions",   lambda: plot_periapsis_distributions(run_dir)),
        ("Phase map",                 lambda: plot_phase_map(run_dir, mid_v)),
        ("Deflection distribution",   lambda: plot_deflection_distribution(run_dir)),
        ("Collision vs escape widths",lambda: plot_collision_vs_escape_width(run_dir)),
        ("Parameter correlations",    lambda: plot_parameter_correlations(run_dir)),
        ("Velocity phase space",      lambda: plot_velocity_phase_space(run_dir)),
    ]

    if candidate_enabled:
        plots.extend([
            ("Best candidate",        lambda: plot_best_candidate(run_dir)),
            ("Candidate ranking",     lambda: plot_candidate_ranking(run_dir, candidate_top_n)),
            ("Candidate panels",      lambda: plot_candidate_ranking_panels(run_dir, candidate_top_n)),
        ])

    plots.append(("Pareto front", lambda: plot_pareto_front(run_dir)))

    if candidate_enabled:
        plots.append((
            "Trajectory tracks",
            lambda: plot_trajectory_tracks(run_dir, candidate_track_top_n),
        ))

    generated = []
    for name, fn in plots:
        if verbose:
            print(f"  {name}...", end=" ", flush=True)
        try:
            out = fn()
            if isinstance(out, (list, tuple)):
                generated.extend(out)
            else:
                generated.append(out)
            if verbose:
                print("done")
        except Exception as exc:
            if verbose:
                print(f"SKIPPED ({exc})")
    return generated
