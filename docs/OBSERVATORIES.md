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
  - host V < 11 and d < 30 pc (evaluation envelope; the faint end is then
    governed by the photon-limited floor below, not a hard cut — V ≥ 11
    hosts are excluded because exposure times become impractical),
  - quadrature separation a/d within [60, 500] mas (IWA ≈ 3λ/D at V for
    ~6.5 m),
  - reflected-light contrast C = A_g·Φ(90°)·(Rp/a)² above the
    brightness-dependent post-processed floor
    **3e-11 · 10^(0.2·max(V−7, 0))** — systematics-limited 3e-11 on bright
    hosts, photon-limited degradation fainter than V ≈ 7.
  - A_g drawn **per planet** from class- and temperature-dependent
    distributions (`atmospheres.sample_geometric_albedo`), stored in the DB
    (`atmospheres.geometric_albedo`); Lambertian Φ(90°) = 1/π. The old fixed
    per-class values remain only as fallbacks for planets without a draw.
    Distribution anchors (all V band):
    - **H/He envelopes** follow the Sudarsky cloud sequence bracketed by
      measurements: cold ammonia-cloud giants ~0.5 (Jupiter 0.52, Saturn
      0.47, Neptune 0.44); a bright water-cloud regime near Teq ~ 250 K
      (mean 0.65); clear alkali-absorbing atmospheres 700–1500 K are DARK —
      the measured Kepler/CHEOPS/TESS hot-Jupiter population sits at
      A_g ≈ 0.03–0.11 (TrES-2b < 0.04); ultra-hot silicate-cloud objects
      recover ~0.3 (Kepler-7b 0.35). Lognormal ×~1.4 scatter, clip
      [0.02, 0.85].
    - **Secondary atmospheres**: dark branch N(0.24, 0.08) centered on
      Earth's re-measured visual A_g of 0.242 (2026 phase-curve analysis;
      classic yield studies assume a flat 0.2), plus a 22% Venus-like
      branch N(0.65, 0.07) — a global cloud deck nearly triples the
      reflected signal (Venus A_g ≈ 0.7).
    - **Steam worlds**: no measurements exist; 0.30 lognormal, [0.08, 0.6].
    - **Airless**: regolith-dark N(0.12, 0.05) (Moon 0.07, Mercury 0.14,
      Mars 0.17), with a 15% icy-bright branch (mean 0.55; Europa 0.67,
      Enceladus > 1) for Teq < 200 K.
    (Earth analog at 10 pc with the branch means: C ≈ 1.4e-10 at 100 mas
    for A_g = 0.24 — still comfortably detected; a Venus-branch twin
    reaches ≈ 3.7e-10.)
- **Headline metric**: `stats` reports `hwo_exo_earth_candidates` — HZ
  planets of 0.8–1.4 R⊕ that HWO can image — mirroring the Astro2020
  "≥25 exo-Earths" goal. The neighborhood population is generated at the
  model-expected count of real systems within 30 pc (~7,100; see
  `stars.expected_systems_within`), so its count reads as an absolute yield
  estimate; `hwo_exo_earth_candidates_scaled` normalizes any other
  population to reality (the 300 pc transit population contains almost no
  HWO-reachable systems by construction).
- **Not modeled**: spectroscopy/biosignature yield, exozodi confusion,
  orbital-phase scheduling (quadrature assumed), UV transit spectroscopy
  mode, multi-visit orbit determination.

## Sources

- [Wilson et al. 2023, Roman GBTDS transit yields (arXiv:2305.16204)](https://arxiv.org/abs/2305.16204)
- [Roman GBTDS ultra-cool dwarf transit yields (arXiv:2303.09959)](https://arxiv.org/abs/2303.09959)
- [NASA JPL: Roman construction complete (2026)](https://www.jpl.nasa.gov/news/nasa-completes-nancy-grace-roman-space-telescope-construction/)
- [HWO concept & technology maturation (arXiv:2601.11803)](https://arxiv.org/abs/2601.11803)
- [HWO coronagraph ConOps (arXiv:2510.02547)](https://arxiv.org/abs/2510.02547)
- [Habitable Worlds Observatory — Wikipedia](https://en.wikipedia.org/wiki/Habitable_Worlds_Observatory)
