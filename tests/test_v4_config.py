"""Configuration and observational-preset tests for the research core."""

from pathlib import Path

import pytest


ROOT = Path(__file__).parent.parent


def test_config_roundtrip(tmp_path):
    from slingshot.v4.config import load_config, save_config

    config = load_config(ROOT / "configs" / "kepler432_quinn.yaml")
    output = tmp_path / "roundtrip.yaml"
    save_config(config, output)
    loaded = load_config(output)
    assert loaded == config
    assert "python/tuple" not in output.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    ("filename", "mass", "semi_major_axis", "eccentricity"),
    [
        ("kepler432_quinn.yaml", 1.32, 0.301, 0.5134),
        ("kepler432_ortiz.yaml", 1.35, 0.303, 0.478),
    ],
)
def test_observational_presets(filename, mass, semi_major_axis, eccentricity):
    from slingshot.v4.config import load_config

    config = load_config(ROOT / "configs" / filename)
    assert config.schema_version == 4
    assert config.system.star_mass_msun == pytest.approx(mass)
    assert config.orbit.semi_major_axis_au == pytest.approx(semi_major_axis)
    assert config.orbit.eccentricity == pytest.approx(eccentricity)
    assert config.metadata.citation.startswith("https://")
    assert config.metadata.legacy_science_model is False


def test_sampling_boundary_must_exceed_impact_domain():
    from slingshot.v4.config import AsymptoticSamplingConfig

    with pytest.raises(ValueError):
        AsymptoticSamplingConfig(b_max_au=2.0, boundary_radius_au=1.0)
