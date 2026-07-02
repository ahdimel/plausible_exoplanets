"""Physical and astronomical constants (SI unless noted).

CODATA 2018 / IAU 2015 nominal values.
"""

# Universal
G = 6.67430e-11            # m^3 kg^-1 s^-2
SIGMA_SB = 5.670374419e-8  # W m^-2 K^-4
C_LIGHT = 2.99792458e8     # m s^-1

# Sun (IAU nominal)
M_SUN = 1.98892e30         # kg
R_SUN = 6.957e8            # m
L_SUN = 3.828e26           # W
TEFF_SUN = 5772.0          # K
MBOL_SUN = 4.74            # IAU absolute bolometric magnitude
RHO_SUN = 1408.0           # kg m^-3 (mean density)

# Earth
M_EARTH = 5.9722e24        # kg
R_EARTH = 6.3781e6         # m (equatorial, IAU nominal)

# Jupiter
M_JUP = 1.89813e27         # kg
R_JUP = 7.1492e7           # m (equatorial, IAU nominal)

# Distances / time
AU = 1.495978707e11        # m
PC = 3.0856775814913673e16 # m
DAY = 86400.0              # s
YEAR = 365.25 * DAY        # s (Julian year)
GYR = 1e9 * YEAR

# Derived conveniences
M_JUP_IN_EARTH = M_JUP / M_EARTH      # ~317.8
DEUTERIUM_BURNING_MASS = 13.0 * M_JUP  # planet / brown-dwarf boundary (conventional)
