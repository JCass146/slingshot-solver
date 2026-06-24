"""One-shot script to fix stale manifest accesses in slingshot/v4/plotting.py."""
import pathlib

src = pathlib.Path("slingshot/v4/plotting.py").read_text(encoding="utf-8")

replacements = [
    # _load unpackings
    ("manifest, _, summary = _load(run_dir)",
     "manifest, _, summary, cfg_raw = _load(run_dir)"),
    ("manifest, samples, _ = _load(run_dir)",
     "manifest, samples, _, cfg_raw = _load(run_dir)"),
    ("manifest, samples, summary = _load(run_dir)",
     "manifest, samples, summary, cfg_raw = _load(run_dir)"),
    # system name
    ("manifest['observational_metadata']['system']['name']",
     "_meta(cfg_raw)['name']"),
    ('manifest["observational_metadata"]["system"]["name"]',
     '_meta(cfg_raw)["name"]'),
    # b_max_au
    ("manifest['observational_metadata']['asymptotic_sampling']['b_max_au']",
     "_meta(cfg_raw)['b_max_au']"),
    ('manifest["observational_metadata"]["asymptotic_sampling"]["b_max_au"]',
     '_meta(cfg_raw)["b_max_au"]'),
    # star_radius_rsun
    ('manifest["observational_metadata"]["system"]["star_radius_rsun"]',
     '_meta(cfg_raw)["star_radius_rsun"]'),
    # planet_radius_rjup
    ('manifest["observational_metadata"]["system"]["planet_radius_rjup"]',
     '_meta(cfg_raw)["planet_radius_rjup"]'),
]

for old, new in replacements:
    count = src.count(old)
    src = src.replace(old, new)
    print(f"  {count}x: {old[:60]!r}")

pathlib.Path("slingshot/v4/plotting.py").write_text(src, encoding="utf-8")
remaining = [ln for ln in src.splitlines() if "observational_metadata" in ln]
if remaining:
    print("REMAINING STALE ACCESSES:")
    for ln in remaining:
        print(" ", ln)
else:
    print("All stale accesses replaced.")
print(f"Total lines: {src.count(chr(10))}")
