"""Validation-gate tests for the v4 physical kernel."""

from pathlib import Path


ROOT = Path(__file__).parent.parent


def test_quick_validation_passes_for_kepler_presets():
    from slingshot.v4.config import load_config
    from slingshot.v4.validation import run_quick_validation

    for name in ("v4_kepler432_quinn.yaml", "v4_kepler432_ortiz.yaml"):
        validation = run_quick_validation(load_config(ROOT / "configs" / name))
        assert validation["passed"]
        assert all(
            gate["passed"]
            for gate in validation["gates"]
            if gate["name"] != "rebound_available"
        )
