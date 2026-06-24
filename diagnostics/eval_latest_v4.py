"""Print a structured evaluation of the most recent schema-v4 run."""
import csv
import json
import pathlib

AU_KM = 1.495978707e8
results = pathlib.Path("results")
dirs = sorted(results.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)

latest = None
for d in dirs:
    m = d / "manifest.json"
    if m.exists():
        manifest = json.loads(m.read_text(encoding="utf-8"))
        if manifest.get("schema_version") == 4:
            latest = d
            break

if latest is None:
    print("No v4 run found.")
    raise SystemExit(1)

print("=" * 72)
print(f"Run directory : {latest.name}")
print(f"Started       : {manifest['started_utc']}")
print(f"Duration      : {manifest['duration_sec']/3600:.2f} h")
print(f"Samples       : {manifest['sample_count']}  "
      f"({manifest['samples_per_bin']} / bin × {len(manifest['seeds'])} seeds × "
      f"{len(manifest['v_inf_kms'])} speeds)")
print(f"Package ver   : {manifest['package_version']}")
print(f"Git commit    : {manifest.get('git_commit','unknown')[:12]}")
print(f"Validation    : {manifest['validation_status'].upper()}")

val = manifest["validation"]
print()
print("─── Campaign gates ───")
gate_keys = [
    ("work_energy_passed",       "Work-energy closure"),
    ("work_energy_max_relative", "  max relative error"),
    ("tail_checks_passed",       "Tail CI checks"),
    ("time_limit_passed",        "Time-limit fraction gate"),
    ("time_limit_fraction",      "  fraction"),
    ("time_limit_count",         "  count"),
    ("numerical_failure_passed", "Numerical failure gate"),
    ("numerical_failure_fraction","  fraction"),
]
for k, label in gate_keys:
    v = val.get(k, "n/a")
    if isinstance(v, float):
        print(f"  {label:<34} {v:.3e}")
    else:
        print(f"  {label:<34} {v}")

print()
print("─── Quick gates ───")
for gate in val["quick"]["gates"]:
    status = "PASS" if gate["passed"] else "FAIL"
    req = "" if gate.get("required", True) else " [diagnostic]"
    print(f"  {gate['name']:<40} {status}{req}")

# Width summary
rows = list(csv.DictReader((latest / "width_summary.csv").open(encoding="utf-8")))
combined = [r for r in rows
            if r.get("scope") == "combined" and r.get("statistic") == "energy_threshold"]

print()
print("─── Planar-width estimates (combined, threshold Δε/vc²=0.0) ───")
print(f"  {'v∞ (km/s)':>10}  {'events':>7}  {'trials':>7}  "
      f"{'W (AU)':>9}  {'CI low':>9}  {'CI high':>9}  "
      f"{'tail_UB':>8}  {'pass':>5}")
for r in combined:
    if float(r["threshold"]) != 0.0:
        continue
    w  = float(r["width_km"]) / AU_KM
    wl = float(r["width_low_km"]) / AU_KM
    wh = float(r["width_high_km"]) / AU_KM
    tub = float(r.get("tail_fraction_upper_bound", "nan"))
    print(f"  {float(r['v_inf_kms']):>10.0f}  {int(r['events']):>7}  "
          f"{int(r['trials']):>7}  {w:>9.5f}  {wl:>9.5f}  {wh:>9.5f}  "
          f"{tub:>8.4f}  {r['tail_check_passed']:>5}")

print()
print("─── Planar-width estimates (combined, all thresholds, v∞=60 km/s) ───")
print(f"  {'thresh':>7}  {'events':>7}  {'W (AU)':>9}  {'CI low':>9}  {'CI high':>9}")
for r in combined:
    if float(r["v_inf_kms"]) != 60.0:
        continue
    w  = float(r["width_km"]) / AU_KM
    wl = float(r["width_low_km"]) / AU_KM
    wh = float(r["width_high_km"]) / AU_KM
    print(f"  {float(r['threshold']):>7.3f}  {int(r['events']):>7}  "
          f"{w:>9.5f}  {wl:>9.5f}  {wh:>9.5f}")

# Seed variance
svar = [r for r in rows if r.get("scope") == "seed_variance"
        and float(r.get("threshold", -1)) == 0.0]
if svar:
    print()
    print("─── Seed-level variance (threshold=0.0) ───")
    print(f"  {'v∞':>5}  {'seeds':>5}  {'mean W (AU)':>11}  "
          f"{'std W (AU)':>10}  {'heterogeneity':>13}")
    for r in svar:
        mw = float(r["seed_mean_width_km"]) / AU_KM
        sw = float(r["seed_std_width_km"]) / AU_KM
        h  = float(r["seed_heterogeneity"])
        print(f"  {float(r['v_inf_kms']):>5.0f}  {int(r['n_seeds']):>5}  "
              f"{mw:>11.5f}  {sw:>10.5f}  {h:>13.4f}")

# Outcome breakdown
sample_rows = list(csv.DictReader((latest / "samples.csv").open(encoding="utf-8")))
outcomes = {}
for s in sample_rows:
    outcomes[s["outcome"]] = outcomes.get(s["outcome"], 0) + 1
print()
print("─── Outcome breakdown (all samples) ───")
for outcome, count in sorted(outcomes.items(), key=lambda x: -x[1]):
    pct = 100.0 * count / len(sample_rows)
    print(f"  {outcome:<25} {count:>6}  ({pct:5.2f}%)")

# Collision widths
col_rows = [r for r in combined[:0]]  # placeholder
col_rows = [r for r in rows
            if r.get("scope") == "combined" and r.get("statistic") == "collision"]
if col_rows:
    print()
    print("─── Collision widths (combined) ───")
    print(f"  {'v∞':>5}  {'events':>7}  {'W (AU)':>9}  {'CI low':>9}  {'CI high':>9}")
    for r in col_rows:
        w  = float(r["width_km"]) / AU_KM
        wl = float(r["width_low_km"]) / AU_KM
        wh = float(r["width_high_km"]) / AU_KM
        print(f"  {float(r['v_inf_kms']):>5.0f}  {int(r['events']):>7}  "
              f"{w:>9.5f}  {wl:>9.5f}  {wh:>9.5f}")
