"""Tiny end-to-end artifact test for the v4 campaign runner."""

import json


def test_tiny_campaign_writes_versioned_artifacts(tmp_path):
    from slingshot.v4.campaign import run_campaign
    from slingshot.v4.config import (
        AsymptoticSamplingConfig,
        MetadataConfig,
        NumericalConfig,
        OrbitConfig,
        SystemConfig,
        V4Config,
        save_config,
    )

    config = V4Config(
        system=SystemConfig(
            name="tiny",
            star_mass_msun=1.0,
            star_radius_rsun=1.0,
            planet_mass_mjup=1.0,
            planet_radius_rjup=1.0,
        ),
        orbit=OrbitConfig(
            model="circular",
            semi_major_axis_au=0.1,
            eccentricity=0.0,
        ),
        asymptotic_sampling=AsymptoticSamplingConfig(
            v_inf_kms=[120.0],
            b_max_au=0.05,
            boundary_radius_au=0.3,
            samples_per_bin=2,
            seeds=[7],
        ),
        numerical=NumericalConfig(max_time_sec=2.0e6),
        metadata=MetadataConfig(
            case_name="tiny",
            description="test",
            parameter_source="synthetic",
            citation="internal",
        ),
    )
    config_path = tmp_path / "config.yaml"
    output_path = tmp_path / "run"
    save_config(config, config_path)
    result = run_campaign(config_path, output_path, verbose=False)
    assert len(result["samples"]) == 2
    for filename in (
        "config.yaml",
        "samples.csv",
        "width_summary.csv",
        "manifest.json",
        "REPORT.md",
    ):
        assert (output_path / filename).exists()
    manifest = json.loads((output_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == 4
    assert manifest["legacy_science_model"] is False
    assert manifest["seeds"] == [7]
    assert manifest["v_inf_kms"] == [120.0]
    assert manifest["samples_per_bin"] == 2
    assert manifest["integrator"] == "DOP853"
    assert manifest["validation_status"] in {"passed", "failed"}
    assert manifest["observational_metadata"]["case_name"] == "tiny"
