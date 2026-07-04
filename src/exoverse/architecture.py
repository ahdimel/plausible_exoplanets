"""System-architecture knobs: intra-system radius correlation (sigma_r) and
mutual-inclination dispersion (sigma_i), the two axes of the Kepler-dichotomy
decomposition study (docs/phase0_audit.md; math note in docs/sigma_r_note.md).

Radius correlation is a Gaussian copula, NOT a literal lognormal hierarchy:
each system draws a latent scale z_sys ~ N(0,1); each small planet draws
z = (z_sys + sigma_r * eps) / sqrt(1 + sigma_r^2),  eps ~ N(0,1),
so z is marginally N(0,1) for every sigma_r, and the radius is the marginal
distribution's quantile at Phi(z). The one-point (marginal) radius
distribution is therefore preserved *exactly* as sigma_r varies — the
renormalization the study design requires — while sibling radii share the
intra-system correlation rho = 1 / (1 + sigma_r^2):

    sigma_r -> 0    : rho -> 1, perfect peas-in-a-pod (identical quantiles)
    sigma_r -> inf  : rho -> 0, independent draws (the baseline null)

sigma_r is dimensionless (latent scatter in units of the marginal's width).
Giant-branch planets keep independent draws: the peas-in-a-pod signal is a
small-planet phenomenon (Weiss+ 2018), and real systems hosting giants are
not size-uniform either.

Mutual inclinations stay Rayleigh about the system plane (system.py applies
them at a random nodal azimuth) with sigma_i now free; an optional
two-component mixture (f_hot drawing sigma_i_hot) implements the classic
dichotomous population, and isotropic=True gives the exact independent
isotropic limit (per-planet cos i uniform, no shared plane).

DEFAULTS ARE BIT-FOR-BIT NEUTRAL: with Architecture() (or arch=None) every
rng draw matches the post-geometry-fix baseline exactly — enforced by
tests/test_architecture.py::test_default_arch_bit_for_bit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

DEG = "deg"


@dataclass(frozen=True)
class Architecture:
    sigma_r: float | None = None    # None = independent radii (baseline null)
    sigma_i: float = 1.5            # Rayleigh scale of mutual inclinations [deg]
    f_hot: float = 0.0              # fraction of systems drawing sigma_i_hot
    sigma_i_hot: float = 30.0       # hot-component Rayleigh scale [deg]
    isotropic: bool = False         # independent isotropic inclinations

    @property
    def is_default(self) -> bool:
        return self == DEFAULT_ARCH

    def meta_items(self) -> dict[str, str]:
        return {"arch_sigma_r": "" if self.sigma_r is None else repr(self.sigma_r),
                "arch_sigma_i": repr(self.sigma_i),
                "arch_f_hot": repr(self.f_hot),
                "arch_sigma_i_hot": repr(self.sigma_i_hot),
                "arch_isotropic": "1" if self.isotropic else "0"}

    @classmethod
    def from_meta(cls, get) -> "Architecture":
        """Rebuild from a get(key, default) accessor (WorldDB.get_meta)."""
        sr = get("arch_sigma_r", "")
        return cls(sigma_r=float(sr) if sr else None,
                   sigma_i=float(get("arch_sigma_i", "1.5")),
                   f_hot=float(get("arch_f_hot", "0.0")),
                   sigma_i_hot=float(get("arch_sigma_i_hot", "30.0")),
                   isotropic=get("arch_isotropic", "0") == "1")


DEFAULT_ARCH = Architecture()


def normal_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def correlated_quantile(rng: np.random.Generator, z_sys: float,
                        sigma_r: float) -> float:
    """Quantile u in (0, 1) for one planet, correlated with its siblings
    through the shared latent z_sys. Marginally uniform for any sigma_r."""
    eps = rng.standard_normal()
    z = (z_sys + sigma_r * eps) / math.sqrt(1.0 + sigma_r * sigma_r)
    return normal_cdf(z)
