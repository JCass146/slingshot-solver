"""Physically defined planar asymptotic-flux initial conditions."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class AsymptoticSample:
    """One proposal defined at infinity and mapped to an inbound boundary."""

    v_inf_kms: float
    impact_parameter_km: float
    incoming_direction_rad: float
    binary_mean_anomaly_rad: float
    position_km: np.ndarray
    velocity_kms: np.ndarray
    specific_energy_km2_s2: float
    angular_momentum_km2_s: float


def state_at_inbound_boundary(
    v_inf_kms: float,
    impact_parameter_km: float,
    incoming_direction_rad: float,
    boundary_radius_km: float,
    total_mu_km3_s2: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Map asymptotic speed and signed impact parameter to a finite boundary."""
    if v_inf_kms <= 0.0:
        raise ValueError("v_inf_kms must be positive")
    if boundary_radius_km <= 0.0 or total_mu_km3_s2 <= 0.0:
        raise ValueError("boundary radius and gravitational parameter must be positive")

    impact_abs = abs(float(impact_parameter_km))
    if impact_abs == 0.0:
        incoming_unit = np.array(
            [np.cos(incoming_direction_rad), np.sin(incoming_direction_rad)]
        )
        position = -boundary_radius_km * incoming_unit
        speed = np.sqrt(v_inf_kms**2 + 2.0 * total_mu_km3_s2 / boundary_radius_km)
        return position, speed * incoming_unit

    angular_momentum_abs = impact_abs * v_inf_kms
    eccentricity = np.sqrt(
        1.0 + impact_abs**2 * v_inf_kms**4 / total_mu_km3_s2**2
    )
    semi_latus_rectum = angular_momentum_abs**2 / total_mu_km3_s2
    cosine_true_anomaly = (
        semi_latus_rectum / boundary_radius_km - 1.0
    ) / eccentricity
    if cosine_true_anomaly < -1.0 or cosine_true_anomaly > 1.0:
        raise ValueError(
            "Impact parameter does not intersect the configured inbound boundary"
        )

    true_anomaly = -np.arccos(np.clip(cosine_true_anomaly, -1.0, 1.0))
    speed_scale = np.sqrt(total_mu_km3_s2 / semi_latus_rectum)
    position_perifocal = boundary_radius_km * np.array(
        [np.cos(true_anomaly), np.sin(true_anomaly)]
    )
    velocity_perifocal = speed_scale * np.array(
        [-np.sin(true_anomaly), eccentricity + np.cos(true_anomaly)]
    )
    asymptote_true_anomaly = -np.arccos(-1.0 / eccentricity)
    asymptote_velocity = speed_scale * np.array(
        [
            -np.sin(asymptote_true_anomaly),
            eccentricity + np.cos(asymptote_true_anomaly),
        ]
    )
    if impact_parameter_km < 0.0:
        position_perifocal[1] *= -1.0
        velocity_perifocal[1] *= -1.0
        asymptote_velocity[1] *= -1.0

    asymptote_angle = np.arctan2(asymptote_velocity[1], asymptote_velocity[0])
    rotation_angle = incoming_direction_rad - asymptote_angle
    cosine = np.cos(rotation_angle)
    sine = np.sin(rotation_angle)
    rotation = np.array([[cosine, -sine], [sine, cosine]])
    position = rotation @ position_perifocal
    velocity = rotation @ velocity_perifocal
    if float(np.dot(position, velocity)) >= 0.0:
        raise RuntimeError("Constructed state is not inbound")
    return position, velocity


def draw_samples(
    rng: np.random.Generator,
    count: int,
    v_inf_kms: float,
    b_max_km: float,
    boundary_radius_km: float,
    total_mu_km3_s2: float,
    sample_incoming_direction: bool = True,
    sample_binary_mean_anomaly: bool = True,
    fixed_direction_rad: float = 0.0,
    fixed_mean_anomaly_rad: float = 0.0,
) -> list[AsymptoticSample]:
    """Draw uniform signed-impact proposals for one fixed v-infinity bin."""
    samples = []
    for _ in range(count):
        impact_parameter = float(rng.uniform(-b_max_km, b_max_km))
        direction = (
            float(rng.uniform(0.0, 2.0 * np.pi))
            if sample_incoming_direction
            else fixed_direction_rad
        )
        mean_anomaly = (
            float(rng.uniform(0.0, 2.0 * np.pi))
            if sample_binary_mean_anomaly
            else fixed_mean_anomaly_rad
        )
        position, velocity = state_at_inbound_boundary(
            v_inf_kms=v_inf_kms,
            impact_parameter_km=impact_parameter,
            incoming_direction_rad=direction,
            boundary_radius_km=boundary_radius_km,
            total_mu_km3_s2=total_mu_km3_s2,
        )
        samples.append(
            AsymptoticSample(
                v_inf_kms=v_inf_kms,
                impact_parameter_km=impact_parameter,
                incoming_direction_rad=direction,
                binary_mean_anomaly_rad=mean_anomaly,
                position_km=position,
                velocity_kms=velocity,
                specific_energy_km2_s2=0.5 * v_inf_kms**2,
                angular_momentum_km2_s=impact_parameter * v_inf_kms,
            )
        )
    return samples
