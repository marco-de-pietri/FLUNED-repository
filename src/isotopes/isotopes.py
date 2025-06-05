from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union
import json
import numpy as np

EPS = 1e-5

@dataclass(frozen=True)
class Isotope:
    name: str
    decay_constant: float
    e_lines: list[float]
    p_lines: list[float]
    emitting_particle: str
    data_source: str
    # derived fields
    branching_ratio: float = field(init=False)
    p_lines_norm: np.ndarray = field(init=False)
    e_bins: np.ndarray = field(init=False)
    p_bins: np.ndarray = field(init=False)

    def __post_init__(self):
        if len(self.e_lines) != len(self.p_lines):
            raise ValueError(f"{self.name}: e_lines and p_lines length mismatch")
        object.__setattr__(self, "branching_ratio", sum(self.p_lines))
        if self.branching_ratio == 0:
            raise ValueError(f"{self.name}: branching_ratio is zero")
        p_norm = np.asarray(self.p_lines) / self.branching_ratio
        e_bins, p_bins = _bins_from_lines(self.e_lines, self.p_lines)
        object.__setattr__(self, "p_lines_norm", p_norm)
        object.__setattr__(self, "e_bins", e_bins)
        object.__setattr__(self, "p_bins", p_bins)

def _bins_from_lines(e_lines: list[float],
                     p_lines: list[float],
                     eps: float = EPS) -> tuple[np.ndarray, np.ndarray]:
    e_bins = [0.0]
    for e in e_lines:
        delta = max(abs(e) * eps, 1e-6)
        e_bins.extend([e - delta, e + delta])
    p_bins = [x for p in p_lines for x in (0.0, p)] + [0.0]
    return np.asarray(e_bins), np.asarray(p_bins)


def load_isotopes(path: Union[str, Path] = "isotopes.json") -> dict[str, Isotope]:
    p = Path(path)
    if not p.is_absolute():
        p = Path(__file__).parent / p
    raw = json.loads(p.read_text(encoding="utf-8"))
    return {k: Isotope(**v) for k, v in raw.items()}


