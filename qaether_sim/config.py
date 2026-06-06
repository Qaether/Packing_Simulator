from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import List, Optional


@dataclass
class ExperimentConfig:
    ell_q: float = 1.0
    radius: float = 0.5
    k_core: float = 1.0
    damping: float = 1.0
    dt: float = 0.01
    relax_dt: Optional[float] = None
    relax_steps: int = 80
    relax_min_steps: int = 20
    relax_convergence_window: int = 10
    relax_msd_tolerance: float = 1.0e-14
    relax_force_tolerance: float = 1.0e-6
    dynamics_steps: int = 120
    snapshot_stride: int = 10
    epsilon_contact: float = 0.06
    epsilon_contact_on: float = 0.01
    epsilon_contact_off: float = 0.02
    epsilon_cycle: float = 0.10
    epsilon_planar: float = 0.12
    epsilon_volume: float = 1.0e-4
    epsilon_center: float = 0.15
    epsilon_perp: float = 0.25
    phi_targets: List[float] = field(
        default_factory=lambda: [0.20, 0.30, 0.40, 0.50, 0.58, 0.64]
    )
    omega_q: float = 1.0
    lambda_phase: float = 0.0

    def to_dict(self) -> dict:
        return asdict(self)
