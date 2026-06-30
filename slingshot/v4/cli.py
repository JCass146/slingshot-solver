"""Command-line interface for planar-width campaigns."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .campaign import run_campaign
from .config import load_config
from .statistics import summarize_tail_gate_status
from .validation import run_quick_validation


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Slingshot Solver effective-planar-width campaign"
    )
    subparsers = parser.add_subparsers(dest="command")
    run_parser = subparsers.add_parser("run", help="Run a campaign")
    run_parser.add_argument("config")
    run_parser.add_argument("--output-dir", "-o")
    run_parser.add_argument("--samples-per-bin", type=int)
    run_parser.add_argument("--seeds", help="Comma-separated integer seeds")
    run_parser.add_argument("--quiet", action="store_true")

    validate_parser = subparsers.add_parser(
        "validate", help="Run deterministic validation gates"
    )
    validate_parser.add_argument("config")

    plot_parser = subparsers.add_parser(
        "plot", help="Generate diagnostic figures for an existing run directory"
    )
    plot_parser.add_argument("run_dir", help="Path to a completed run directory")
    plot_parser.add_argument("--quiet", action="store_true")
    return parser



def _read_csv_rows(path: Path) -> list[dict]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _refresh_report_for_run(run_dir: str | Path, generated: list[Path]) -> None:
    """Refresh manifest artifact metadata and regenerate REPORT.md."""
    from .report import generate_report

    run_path = Path(run_dir)
    config = load_config(run_path / "config.yaml")
    manifest_path = run_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    summary_rows = _read_csv_rows(run_path / "width_summary.csv")
    candidate_rows = _read_csv_rows(run_path / "top_candidates.csv")
    tail_status = summarize_tail_gate_status(summary_rows)
    validation = manifest.setdefault("validation", {})
    validation.update(tail_status)
    validation["passed"] = (
        validation.get("quick", {}).get("passed", False)
        and validation.get("work_energy_passed", False)
        and validation.get("tail_checks_passed", False)
        and validation.get("time_limit_passed", False)
        and validation.get("numerical_failure_passed", False)
    )
    manifest["validation_status"] = "passed" if validation["passed"] else "failed"

    manifest["candidate_diagnostics"] = config.candidate_diagnostics.model_dump(mode="json")
    manifest["candidate_count"] = len(candidate_rows)
    manifest["best_observed_gain"] = (
        float(candidate_rows[0]["energy_gain_dimensionless"])
        if candidate_rows else None
    )

    artifacts = manifest.setdefault("artifacts", [])
    for artifact in ["top_candidates.csv", "REPORT.md"] + [Path(path).name for path in generated]:
        if (run_path / artifact).exists() and artifact not in artifacts:
            artifacts.append(artifact)

    with manifest_path.open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2)
    generate_report(run_path, config, summary_rows, manifest)

def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "validate":
        config = load_config(args.config)
        validation = run_quick_validation(config)
        for gate in validation["gates"]:
            print(f"{gate['name']}: {'PASS' if gate['passed'] else 'FAIL'}")
        raise SystemExit(0 if validation["passed"] else 1)
    if args.command == "run":
        seeds = (
            [int(value) for value in args.seeds.split(",")]
            if args.seeds
            else None
        )
        result = run_campaign(
            config_path=args.config,
            output_dir=args.output_dir,
            samples_per_bin=args.samples_per_bin,
            seeds=seeds,
            verbose=not args.quiet,
        )
        if not args.quiet:
            print(f"Results: {Path(result['output_dir'])}")
        return
    if args.command == "plot":
        from .plotting import generate_all_plots
        generated = generate_all_plots(args.run_dir, verbose=not args.quiet)
        _refresh_report_for_run(args.run_dir, generated)
        if not args.quiet:
            print(f"Generated {len(generated)} figures and refreshed REPORT.md in {args.run_dir}")
        return
    parser.print_help()
    raise SystemExit(2)
