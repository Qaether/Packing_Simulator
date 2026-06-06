from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import h5py
import numpy as np


@dataclass
class QaetherState:
    # Computational embedding coordinates for graph construction, not a
    # fundamental background-space position of a particle.
    positions: np.ndarray
    # Effective relaxation variable used by numerical protocols, not a literal
    # microscopic particle velocity.
    velocities: np.ndarray
    box: np.ndarray
    radius: float = 0.5
    ell_q: float = 1.0
    theta: Optional[np.ndarray] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int(self.positions.shape[0])

    @property
    def volume(self) -> float:
        return float(np.prod(self.box))

    @property
    def phi(self) -> float:
        sphere_volume = 4.0 * np.pi * self.radius**3 / 3.0
        return self.n * sphere_volume / self.volume

    def copy(self) -> "QaetherState":
        return QaetherState(
            positions=self.positions.copy(),
            velocities=self.velocities.copy(),
            box=self.box.copy(),
            radius=self.radius,
            ell_q=self.ell_q,
            theta=None if self.theta is None else self.theta.copy(),
            metadata=dict(self.metadata),
        )

    def wrap(self) -> None:
        self.positions = np.mod(self.positions, self.box)

    def save_h5(self, path: str) -> None:
        with h5py.File(path, "w") as h5:
            h5.create_dataset("positions", data=self.positions)
            h5.create_dataset("velocities", data=self.velocities)
            h5.create_dataset("box", data=self.box)
            if self.theta is not None:
                h5.create_dataset("theta", data=self.theta)
            h5.attrs["radius"] = self.radius
            h5.attrs["ell_q"] = self.ell_q
            h5.attrs["metadata_json"] = json.dumps(self.metadata, sort_keys=True)

    @classmethod
    def load_h5(cls, path: str) -> "QaetherState":
        with h5py.File(path, "r") as h5:
            theta = h5["theta"][()] if "theta" in h5 else None
            metadata = json.loads(h5.attrs.get("metadata_json", "{}"))
            return cls(
                positions=h5["positions"][()],
                velocities=h5["velocities"][()],
                box=h5["box"][()],
                radius=float(h5.attrs["radius"]),
                ell_q=float(h5.attrs["ell_q"]),
                theta=theta,
                metadata=metadata,
            )
