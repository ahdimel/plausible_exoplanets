"""exoverse: procedurally generated, physics-validated exoplanetary systems."""

__version__ = "0.1.0"

from .stars import Star, generate_star
from .planets import Planet, generate_planet
from .system import StellarSystem, generate_system
from .transits import compute_geometry, model_light_curve, transit_flux
from .observatories import observe
from .database import WorldDB
from .generate import generate_population

__all__ = [
    "Star", "generate_star", "Planet", "generate_planet",
    "StellarSystem", "generate_system", "compute_geometry",
    "model_light_curve", "transit_flux", "observe", "WorldDB",
    "generate_population",
]
