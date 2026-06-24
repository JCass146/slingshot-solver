"""Configuration models for the defensible planar v4 research core."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class SystemConfig(BaseModel):
    """Physical star and planet properties."""

    model_config = ConfigDict(extra="forbid")

    name: str
    star_mass_msun: float = Field(gt=0.0)
    star_radius_rsun: float = Field(gt=0.0)
    planet_mass_mjup: float = Field(gt=0.0)
    planet_radius_rjup: float = Field(gt=0.0)
    bulk_velocity_x_kms: float = 0.0
    bulk_velocity_y_kms: float = 0.0


class OrbitConfig(BaseModel):
    """Relative Keplerian orbit of the star-planet binary."""

    model_config = ConfigDict(extra="forbid")

    model: Literal["circular", "keplerian"] = "keplerian"
    semi_major_axis_au: float = Field(gt=0.0)
    eccentricity: float = Field(default=0.0, ge=0.0, lt=1.0)
    mean_anomaly_rad: float = 0.0
    argument_periapsis_rad: float = 0.0
    prograde: bool = True

    @model_validator(mode="after")
    def normalize_circular(self):
        if self.model == "circular":
            self.eccentricity = 0.0
        return self


class AsymptoticSamplingConfig(BaseModel):
    """Planar inbound-flux proposal."""

    model_config = ConfigDict(extra="forbid")

    v_inf_kms: List[float] = Field(
        default_factory=lambda: [10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 120.0]
    )
    b_max_au: float = Field(default=1.0, gt=0.0)
    boundary_radius_au: float = Field(default=5.0, gt=0.0)
    samples_per_bin: int = Field(default=2000, ge=1)
    seeds: List[int] = Field(default_factory=lambda: [42, 43, 44, 45, 46])
    sample_incoming_direction: bool = True
    sample_binary_mean_anomaly: bool = True

    @field_validator("v_inf_kms")
    @classmethod
    def validate_speeds(cls, values: List[float]):
        if not values or any(value <= 0.0 for value in values):
            raise ValueError("v_inf_kms must contain positive speeds")
        return values

    @field_validator("seeds")
    @classmethod
    def validate_seeds(cls, values: List[int]):
        if not values or any(seed < 0 for seed in values):
            raise ValueError("seeds must contain non-negative integers")
        return values

    @model_validator(mode="after")
    def validate_boundary(self):
        if self.boundary_radius_au <= self.b_max_au:
            raise ValueError("boundary_radius_au must exceed b_max_au")
        return self


class PlanarWidthConfig(BaseModel):
    """Statistical estimands and quality gates for effective planar width."""

    model_config = ConfigDict(extra="forbid")

    dimensionless_energy_thresholds: List[float] = Field(
        default_factory=lambda: [0.0, 0.01, 0.03, 0.1, 0.3, 1.0]
    )
    confidence_level: float = Field(default=0.95, gt=0.0, lt=1.0)
    tail_fraction: float = Field(default=0.10, gt=0.0, lt=0.5)
    max_tail_event_fraction: float = Field(default=0.01, ge=0.0, le=1.0)

    @field_validator("dimensionless_energy_thresholds")
    @classmethod
    def validate_thresholds(cls, values: List[float]):
        if not values or any(value < 0.0 for value in values):
            raise ValueError("dimensionless energy thresholds must be non-negative")
        return sorted(set(float(value) for value in values))


class NumericalConfig(BaseModel):
    """Campaign integrator settings."""

    model_config = ConfigDict(extra="forbid")

    method: Literal["DOP853", "Radau"] = "DOP853"
    rtol: float = Field(default=1e-10, gt=0.0)
    atol: float = Field(default=1e-10, gt=0.0)
    max_time_sec: float = Field(default=3e8, gt=0.0)
    softening_km: float = Field(default=0.0, ge=0.0)
    max_step_sec: Optional[float] = Field(default=None, gt=0.0)


class ValidationConfig(BaseModel):
    """Publication gates for the numerical core."""

    model_config = ConfigDict(extra="forbid")

    binary_elements_relative_tolerance: float = 1e-10
    two_body_deflection_relative_tolerance: float = 1e-6
    jacobi_relative_tolerance: float = 1e-8
    work_energy_relative_tolerance: float = 1e-6
    integrator_agreement_relative_tolerance: float = 1e-5
    run_rebound_if_available: bool = True
    # Campaign-level failure gates (P0.2)
    max_time_limit_fraction: float = Field(default=0.05, ge=0.0, le=1.0)
    max_numerical_failure_fraction: float = Field(default=0.01, ge=0.0, le=1.0)


class MetadataConfig(BaseModel):
    """Observational provenance attached to every run."""

    model_config = ConfigDict(extra="allow")

    case_name: str
    description: str
    parameter_source: str
    citation: str
    uncertainties: Dict[str, Any] = Field(default_factory=dict)
    legacy_science_model: bool = False


class V4Config(BaseModel):
    """Complete schema-version-4 configuration."""

    model_config = ConfigDict(extra="forbid")

    schema_version: Literal[4] = 4
    system: SystemConfig
    orbit: OrbitConfig
    asymptotic_sampling: AsymptoticSamplingConfig = Field(
        default_factory=AsymptoticSamplingConfig
    )
    planar_width: PlanarWidthConfig = Field(default_factory=PlanarWidthConfig)
    numerical: NumericalConfig = Field(default_factory=NumericalConfig)
    validation: ValidationConfig = Field(default_factory=ValidationConfig)
    metadata: MetadataConfig


def load_config(path: str | Path) -> V4Config:
    """Load a strict v4 YAML or JSON configuration."""
    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as stream:
        if config_path.suffix.lower() in {".yaml", ".yml"}:
            data = yaml.safe_load(stream)
        elif config_path.suffix.lower() == ".json":
            data = json.load(stream)
        else:
            raise ValueError(f"Unsupported config format: {config_path.suffix}")
    return V4Config.model_validate(data or {})


def save_config(config: V4Config, path: str | Path) -> None:
    """Safely serialize a v4 configuration without Python-specific YAML tags."""
    config_path = Path(path)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    data = config.model_dump(mode="json")
    with config_path.open("w", encoding="utf-8") as stream:
        if config_path.suffix.lower() in {".yaml", ".yml"}:
            yaml.safe_dump(data, stream, sort_keys=False)
        elif config_path.suffix.lower() == ".json":
            json.dump(data, stream, indent=2)
        else:
            raise ValueError(f"Unsupported config format: {config_path.suffix}")
