# Provenance and reproducibility limits

Reference for assessing and adapting [Particle_Filter](../README.md). The historical source baseline is [588308f](https://github.com/rhymesg/Particle_Filter/tree/588308fadbe7de7463b51cdae339b29e619da739). Current changes correct APF importance weights and the reported sample mean; the research limitations below remain.

## Differences from the paper

The [published method](https://doi.org/10.1109/TAES.2017.2741878) and the supplied implementation are related but not identical.

| Paper | Supplied implementation | Implication |
|---|---|---|
| Algorithm 2 draws past states conditioned on the last accepted and current states, Eqs. (8)–(9) | [OOSM.m](../OOSM.m) uses stored prior clouds; `last.particle` is recorded in the entry point but not passed to the helper | Conditional sampling is absent |
| Algorithm 3 accepts an OOSM candidate only when its covariance determinant improves | `dcov_oosm` is calculated but never compared with input `dcov_pl`; the final candidate is returned | The additional acceptance gate is absent |
| Algorithm 3's pseudocode uses a strict decrease for the current update | The entry point skips only a strict increase | Equal determinants are accepted by the code |
| Eq. (7) normalizes covariance by particle count | Calls to `cov` use default sample normalization | Covariance magnitudes differ |
| Section IV-C uses 1000 particles and 100 Monte Carlo runs | The settings block supplies a smaller demonstration configuration | Defaults do not reproduce the experiment |
| Section IV compares five methods, including RHKF | Only four PF variants are supplied | The full comparison cannot be run here |
| Table I specifies three noise cases | One active configuration uses `sig_z = 15 + 4.71` and its own process/bias settings | A case-by-case reproduction requires reconciling settings |

The MPF uses critical-bandwidth mode estimation and nearest-centre assignment, rather than the included mean-shift helper. The MPF resamples each mode to its previous particle count and then discards posterior mode mass when resetting global weights. A repair must preserve particle and mode weights through clustering, or reallocate particles by posterior mode mass; its current multi-step posterior is not reliable. The APF now includes its lookahead likelihood-ratio correction.

## Runtime and numerical limits

- [DEM_height.m](../DEM_height.m) assumes every particle stays inside the [terrain interpolation boundary](simulation.md#terrain-data); out-of-range indices cause errors.
- Likelihood densities can still underflow. Resampling and APF second-stage weighting reject invalid totals rather than silently returning a cloud; the other direct normalizations have no recovery policy.
- `OOSM.m` requires a nonempty queue; otherwise `particle_star` is undefined, and its `oosmSucceed` output is never updated even for a nonempty queue.
- [FindCriticalBW.m](../FindCriticalBW.m) can reference uninitialized `r`/`c` if its initial bandwidth search reaches the fallback before finding a mode transition.
- [dskensity2d.m](../dskensity2d.m) assumes independent coordinate marginals; [Significance.m](../Significance.m) divides by marginal variances without handling zero variance.
- The full entry point assumes all filter outputs and its original time grid are present; see [settings](simulation.md#settings-and-randomness).
- [HGMeanShiftCluster.m](../HGMeanShiftCluster.m) is unused by the main run; its Gaussian branch references `gaussfun`, which is not supplied.

## Provenance and attribution

- The canonical repository is [rhymesg/Particle_Filter](https://github.com/rhymesg/Particle_Filter); Git history and the original entry-point header identify Youngjoo Kim as the code author.
- The local journal PDF and [KAIST publication record](https://pure.kaist.ac.kr/en/publications/utilizing-out-of-sequence-measurement-for-ambiguous-update-in-par/) establish the [citation metadata](../README.md#citation); the paper DOI identifies the publication, not a software release.
- The repository's [MIT license](../LICENSE) is retained; `HGMeanShiftCluster.m` also carries Han Gong and Bart Finkston copyright notices, which remain intact, but no separate third-party license file is supplied.
- The paper describes SRTM terrain near 38° N, 128° E; the repository does not supply the exact source tile, preprocessing history, or a dataset-specific license for `DB_part.mat`.
- Local reference PDFs remain in ignored `ref/` and are not distributed with these documentation changes.

## Verification status

- Publication metadata and Algorithm 2–3 mappings were checked against the supplied journal PDF; source paths, helper calls, and run-output descriptions were inspected.
- `CITATION.cff` passes the CFF 1.2.0 schema; local documentation links and anchors resolve, and the terrain structure was inspected with SciPy.
- Terrain data are unchanged. [APF regression checks](../tests/integration/auxiliary/README.md) cover the changed weight algebra; these checks have not run natively.
- MATLAB and Octave are unavailable in the review environment; neither the simulation nor the synthetic helper command has been run.
- No validated numerical tolerances, deterministic full-run baseline, or reproduction of published figures is available.

## Suggested repository metadata

These drafts are for the GitHub About section; remote metadata has not been changed.

- Description: `MATLAB particle filters for terrain-referenced navigation: PF, APF, MPF and OOSM research code accompanying an IEEE TAES paper.`
- Topics: `matlab`, `particle-filter`, `sequential-monte-carlo`, `out-of-sequence-measurements`, `terrain-referenced-navigation`, `terrain-aided-navigation`, `state-estimation`, `aerospace`, `bayesian-filtering`.
