"""Command-line interface for planar-width campaigns."""

from __future__ import annotations

import argparse
from pathlib import Path

from .campaign import run_campaign
from .config import load_config
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
        if not args.quiet:
            print(f"Generated {len(generated)} figures in {args.run_dir}")
        return
    parser.print_help()
    raise SystemExit(2)
