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


def test_two_body_deflection_uses_asymptotic_comparison():
    """Verify the deflection gate compares asymptotic angles, not finite-boundary angles."""
    from slingshot.constants import G_KM, M_SUN, M_JUP, AU_KM
    from slingshot.v4.validation import numerical_two_body_deflection

    # Sun + Jupiter mass central body, 1 AU boundary
    mu = G_KM * (1.0 * M_SUN + 1.0 * M_JUP)
    result = numerical_two_body_deflection(
        v_inf_kms=30.0,
        impact_parameter_km=0.1 * AU_KM,
        boundary_radius_km=1.0 * AU_KM,
        mu_km3_s2=mu,
    )
    # Must have the asymptotic deflection field, NOT the old finite-boundary field
    assert "numerical_asymptotic_deflection_rad" in result
    assert "numerical_boundary_deflection_rad" not in result
    # Asymptotic deflection must agree with analytic to within configured tolerance
    tol = 1e-6
    assert result["deflection_relative_error"] <= tol, (
        f"Asymptotic deflection relative error {result['deflection_relative_error']:.2e} "
        f"exceeds tolerance {tol}"
    )
