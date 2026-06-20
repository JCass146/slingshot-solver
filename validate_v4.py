#!/usr/bin/env python
"""Run publication-grade v4 validation gates."""

import argparse
import json

from slingshot.v4.config import load_config
from slingshot.v4.gates import run_publication_validation


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result = run_publication_validation(load_config(args.config))
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        for gate in result["gates"]:
            status = "SKIP" if gate.get("skipped") else ("PASS" if gate["passed"] else "FAIL")
            print(f"{gate['name']}: {status}")
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
