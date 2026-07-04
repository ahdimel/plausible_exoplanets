# Appendix note: the σ_R radius-correlation knob and marginal preservation

*Phase 1 of the dichotomy-decomposition study (docs/phase0_audit.md).
Implementation: `architecture.py`, `planets.py::small_radius_ppf`.*

## Requirement

The study sweeps intra-system radius correlation while holding the
population's one-point (marginal) radius distribution fixed. If the sweep
changed the marginal, the (σ_R, σ_i) surface would confound correlation
strength with the radius function itself.

## Why a literal lognormal hierarchy fails here

The brief's naive construction — draw a system scale from the population
radius function, scatter each planet lognormally with width σ_R about it —
does *not* preserve the marginal: the marginal of that hierarchy is the
population distribution *convolved* with the lognormal scatter, so the
one-point distribution broadens as σ_R grows, and de-convolving the base
distribution to compensate has no exact solution for our marginal, which is
a period-dependent two-lognormal mixture with clip atoms at 0.4 and 4.0 R⊕
(`sample_small_radius`): deconvolution against a mixture-with-atoms is
ill-posed.

## The Gaussian-copula construction (exact)

Work in latent-Gaussian space and let the quantile function carry all the
structure of the marginal:

1. Per system: `z_sys ~ N(0, 1)`.
2. Per (small) planet: `ε ~ N(0, 1)` and

       z = (z_sys + σ_R · ε) / sqrt(1 + σ_R²)

3. `u = Φ(z)`; the radius is `R = F⁻¹_P(u)` where `F_P` is the *exact* CDF
   of `sample_small_radius` at the planet's period P
   (`small_radius_cdf/ppf`; clip mass maps to the boundary atoms, matching
   `np.clip` semantics).

**Marginal preservation.** z_sys and ε are independent standard normals, so
z is N(0, 1) *exactly* for every σ_R; hence u is Uniform(0, 1) and R has
the marginal F_P exactly. No renormalization step is needed — the
construction is marginal-preserving by identity, not by correction.

**Correlation structure.** Two siblings share z_sys, giving latent
correlation

       ρ = corr(z_i, z_j) = 1 / (1 + σ_R²)

- σ_R → 0: ρ → 1. Siblings get identical quantiles; radii differ only
  through the period-dependent valley shift (factor (P/10)^(−0.05), a few
  per cent across a typical system) — the perfect peas-in-a-pod limit.
- σ_R → ∞: ρ → 0. Independent draws; **bit-identical code path to the
  baseline is recovered by σ_R = None** (the knob off), which skips the
  latent draw entirely.

**Interpretation of σ_R.** σ_R is dimensionless: the intra-system scatter
in units of the *population's* latent width. For a single-lognormal
marginal with log-width s, the observable sibling dispersion is

       std(log R_i − log R_j) = s · sqrt(2) · σ_R / sqrt(1 + σ_R²)

i.e., for small σ_R it reads as a lognormal scatter of width s·σ_R about a
system scale — the brief's intended parameter — while automatically
saturating at the independent-draw dispersion for large σ_R. Our marginal
is a bimodal mixture, so the mapping to |Δlog R| is evaluated numerically
(validation plots, `analysis/phase1_validation_plots.py`); the monotonicity
and the two limits are unit-tested (`test_sigma_r_uniformity_monotonic`,
`test_marginal_radius_invariant_under_sigma_r`).

## Scope decisions (documented, revisit in Phase 4 if needed)

- **Only small-branch planets (R < 4 R⊕ population) are correlated.** The
  giant coin flip and giant radii stay independent: peas-in-a-pod is an
  observed small-planet phenomenon (Weiss+ 2018), and real giant-hosting
  systems are not size-uniform. Detected giant–small pairs therefore dilute
  the |Δlog R| statistic identically in data and simulation.
- **The radius–period valley coupling is preserved** (quantiles are mapped
  through the period-dependent CDF), so σ_R does not distort the
  radius–period plane.
- The class coin (giant vs small) is uncorrelated within a system; if
  Phase 4 shows the |Δlog R| statistic is sensitive to giant clustering,
  correlating the coin through the same latent z_sys is a one-line
  extension.
