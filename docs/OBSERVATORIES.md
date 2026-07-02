# Observatory models: assumptions, sources, and volatility

Documentation for each modeled observatory, with emphasis on the two
**future** facilities whose specifications WILL change as the missions
evolve. Last researched: **2026-07-02**. Anything marked (volatile) should
be re-verified before being relied on.

## Operating / archival

| Facility | Key assumptions | Basis |
|---|---|---|
| TESS | 1-hr precision ~150 ppm @ T=8 rising to ~6000 @ T=14, 60 ppm floor; 27.4-d sector; SNR>=7.1 + >=2 transits | Stassun+ 2018 noise curves; SPOC pipeline thresholds |
| Kepler | 6.5-h CDPP 30 ppm @ Kp=12, photon-scaled, 20 ppm floor; 4-yr baseline, 92% duty | Kepler mission performance; retired 2018 (results are "what Kepler would have seen") |
| JWST NIRISS SOSS | ~20 ppm/hr @ J=8, 10 ppm floor, saturates J<~6.5; single targeted transit | ERS/GO program performance (e.g. WASP-39b ERS papers) |
| JWST NIRSpec Prism | ~12 ppm/hr @ J=11, 10 ppm floor, saturates J<~10.5 | Prism bright limit + published transit precisions |
| Ground 1-m | 2 mmag/point + 1 mmag correlated floor; 90 nights, 33% duty | Typical small-telescope survey performance |

## Nancy Grace Roman Space Telescope — Galactic Bulge Time Domain Survey

**Status (2026-07)**: construction complete, launch expected late 2026,
science operations from 2027. (volatile: schedule)

What we implemented, and where each number comes from:

- **Telescope/survey**: 2.4 m aperture; Wide Field Instrument; GBTDS observes
  bulge fields in F146 (0.93–2.0 µm) every ~12 minutes, with ~6 seasons of
  ~72 days across the ~5-yr primary mission. Sensitive for stars F146 ≲ 21.
- **Yields**: predicted ~60,000–200,000 transiting planets (mostly giants,
  a < 0.3 AU) plus ~7,000–12,000 small (<4 R⊕) planets — Wilson et al. 2023
  ([arXiv:2305.16204](https://arxiv.org/abs/2305.16204)); ~1300 small planets
  around mid-M dwarfs — [arXiv:2303.09959](https://arxiv.org/abs/2303.09959).
- **Our noise model**: photon-scaled around ~700 ppm/hr at F146=16, 200 ppm
  bright-end floor, usable 8 < F146 < 21, F146 ≈ J as a color proxy.
  These precision anchors are OUR calibration to reproduce the published
  yield regime, not official mission numbers. (volatile)
- **Caveat baked into every result**: the GBTDS points at ~2 deg² of Galactic
  bulge. Our generated systems are local. Roman rows therefore answer the
  counterfactual "if this system sat in a Roman bulge field" — flagged in
  the observation note.

## Habitable Worlds Observatory (HWO)

**Status (2026-07)**: pre-Phase-A concept maturation (technology programs,
trade-space studies); no confirmed architecture. Everything below is
requirement-era and WILL change. Target launch: **2040s**.

- **Concept**: ~6–8 m class segmented UV/optical/NIR telescope; prime goal
  to directly image and characterize ≥25 potentially habitable exo-Earths
  (Astro2020 decadal recommendation).
- **Coronagraph**: raw contrast ~1e-10 goal at the relevant working angles;
  post-processed sensitivity requirement often quoted to flux ratios
  ~1e-11. Sources: HWO concept/technology papers
  ([arXiv:2601.11803](https://arxiv.org/abs/2601.11803),
  [arXiv:2510.02547](https://arxiv.org/abs/2510.02547)),
  [Wikipedia/HWO](https://en.wikipedia.org/wiki/Habitable_Worlds_Observatory).
- **Our model** (all volatile): detection when
  - host V < 8 and d < 30 pc (nearby-bright target-list proxy),
  - quadrature separation a/d within [60, 500] mas (IWA ≈ 3λ/D at V for
    ~6.5 m),
  - reflected-light contrast C = A_g·Φ(90°)·(Rp/a)² > 3e-11, with A_g = 0.3
    and Lambertian Φ(90°) = 1/π. (Earth analog at 10 pc: C ≈ 1.7e-10 at
    100 mas — comfortably detected, as intended.)
- **Not modeled**: spectroscopy/biosignature yield, exozodi confusion,
  observation scheduling, UV transit spectroscopy mode.

## Sources

- [Wilson et al. 2023, Roman GBTDS transit yields (arXiv:2305.16204)](https://arxiv.org/abs/2305.16204)
- [Roman GBTDS ultra-cool dwarf transit yields (arXiv:2303.09959)](https://arxiv.org/abs/2303.09959)
- [NASA JPL: Roman construction complete (2026)](https://www.jpl.nasa.gov/news/nasa-completes-nancy-grace-roman-space-telescope-construction/)
- [HWO concept & technology maturation (arXiv:2601.11803)](https://arxiv.org/abs/2601.11803)
- [HWO coronagraph ConOps (arXiv:2510.02547)](https://arxiv.org/abs/2510.02547)
- [Habitable Worlds Observatory — Wikipedia](https://en.wikipedia.org/wiki/Habitable_Worlds_Observatory)
