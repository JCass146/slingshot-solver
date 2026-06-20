"""Publication-gate regression tests."""

from pathlib import Path


ROOT = Path(__file__).parent.parent


def test_publication_validation_required_gates_pass():
    from slingshot.v4.config import load_config
    from slingshot.v4.gates import run_publication_validation

    config = load_config(ROOT / "configs" / "v4_kepler432_quinn.yaml")
    result = run_publication_validation(config)
    failures = [
        gate
        for gate in result["gates"]
        if gate.get("required", True)
        and not gate.get("skipped", False)
        and not gate["passed"]
    ]
    assert not failures
    assert result["passed"]
