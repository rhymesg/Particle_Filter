# Implementation and provenance

Reference for adapting the terrain-navigation filters and connecting them to the [IEEE TAES paper](../README.md#citation).

## Differences from the paper

| Component | Source procedure |
|---|---|
| OOSM candidates | [OOSM.m](../OOSM.m) reuses stored prior clouds; paper Eqs. (8)–(9) instead condition past-state sampling on accepted and current states. |
| OOSM acceptance | The helper returns its final candidate; the paper adds a covariance-determinant acceptance test. The entry point accepts equal current-update determinants. |
| Covariance | MATLAB `cov` uses sample normalization; paper Eq. (7) normalizes by particle count. |
| Experiment | The settings block selects a demonstration configuration; use the paper's Section IV and Table I for its experiment settings. |

The APF applies the lookahead likelihood-ratio correction. The MPF uses critical-bandwidth mode estimation and nearest-center assignment, resamples each mode to its existing count, then resets global weights uniformly. Thus mode population counts determine the next step's prior mass; preserve that behavior when comparing historical outputs.

## Numerical contracts

Terrain particles must remain inside the [interpolation neighborhood](simulation.md#terrain-data). Likelihood normalization requires a positive finite total, and `OOSM` requires a nonempty queue. `dskensity2d` treats coordinate marginals independently; the significance calculation requires positive marginal variances.

`FindCriticalBW` returns centers from the final density grid when its search reaches the minimum bandwidth. A transition found on the final iteration retains the preceding bandwidth and centers.

## Provenance and attribution

Git history and the original entry-point header identify Youngjoo Kim as the author. Use the [paper citation](../README.md#citation) and preserve the [MIT license](../LICENSE) and existing third-party notices, including those in `HGMeanShiftCluster.m`.

The paper describes SRTM terrain near 38° N, 128° E; the [data guide](simulation.md#terrain-data) describes the supplied `DB_part.mat` structure. Confirm dataset and third-party reuse terms with the respective rights holders.

## Checks

[APF and bandwidth checks](../tests/integration/auxiliary/README.md) exercise importance weights, ancestor correspondence, and critical-bandwidth termination. The [simulation guide](simulation.md) specifies dependencies, entry points, and output variables.
