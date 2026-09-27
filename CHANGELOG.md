# Second-draft change log — 27 September 2026

- Corrected radiator area by eight orders of magnitude and distinguished total
  emitting area, one-face disk radius, and spherical radius.
- Recomputed all benchmark fluxes, exact two-band slopes, and conditional horizons
  from a shared physical model; retained the original correct isotropic flux normalization.
- Replaced the Moran-process description with a synchronous four-state cellular
  automaton; made competing updates explicit and reproducible with per-run seeds.
- Added the actual discrete inward loss probability, colonization competition,
  and remnant/expansionist occupancy to the balance relation. Replaced unsupported
  R0 values and percolation thresholds with a clearly limited reproduction bound.
- Ran two parameter settings, ten seeds each, 20,000 steps per seed. Reported
  late-time means and between-seed standard deviations instead of universal abundance claims.
- Defined irreversible erasure separately from gross computation, separated bath
  and radiator temperatures, and added cooling-area and irradiation limitations.
- Corrected Wien/Rayleigh–Jeans terminology, Snu/Slambda peak conventions, and
  FIR/submillimeter coverage. Removed claims that spectral slope proves artificiality.
- Updated GERONIMO to PRIMAger and PRIMA's dated development status; replaced
  undocumented sensitivity curves and integration-time promises with explicit
  required noise and SNR. Kept bandpass and confusion limitations visible.
- Removed the unsupported FIRAS exclusion; corrected the assumed-budget limit to
  1.90e45 bits/s and its separation from the benchmark to 5.28 dex.
- Regenerated three companion figures, supplied machine-readable results, added
  15 regression/identity tests and an original-code bug reproduction, and made
  `main.tex` independent of external image and bibliography files.
- Built-in manuscript compilation remains unverified because the compiler failed
  with a platform-directory error. No manuscript PDF is represented as compiled.
