# Particle filtering for terrain-referenced navigation

Implementation reference for [Particle_Filter](../README.md), accompanying the [IEEE TAES paper](https://doi.org/10.1109/TAES.2017.2741878). Use the [simulation guide](simulation.md) to run the code and the [citation section](../README.md#citation) when adapting it.

## Model and common operations

The state is a two-component horizontal position; the code uses local distances in metres and terrain heights in metres. A particle array has shape `2 × N`, with one position per column.

The implemented process and measurement models correspond to paper Eqs. (10)–(11):

$$x_k^{(i)-} = x_{k-1}^{(i)} + v_{k-1}^{\mathrm{meas}}\Delta t + \sigma_{\mathrm{proc}}\epsilon_k^{(i)}\Delta t, \qquad z_k = h(x_k) + \eta_k.$$

Here `epsilon` is a two-component standard Gaussian draw; [DEM_height.m](../DEM_height.m) evaluates `h` by local bicubic interpolation. The filter likelihood in [likelihood.m](../likelihood.m) is scalar Gaussian:

$$L(z\mid x)=\frac{1}{\sqrt{2\pi}\sigma_{\mathrm{meas}}}\exp\left[-\frac{(z-h(x))^2}{2\sigma_{\mathrm{meas}}^2}\right].$$

- [main_OOSM.m](../main_OOSM.m), `%% PF`: propagate, multiply weights by `L`, normalize, resample, and average the resampled positions (Algorithm 1).
- [Resample.m](../Resample.m): independently sample indices from the cumulative weights and return uniformly weighted particles; this is multinomial resampling.
- Covariance calls use MATLAB's default sample covariance (`N-1` normalization); paper Eq. (7) uses `N`.

## Postponing ambiguous measurements

The paper addresses occasional local measurement ambiguity, rather than persistent global localization ambiguity. Algorithm 3 compares prior and posterior covariance determinants; this is the criterion implemented in `main_OOSM.m`, not a multimodality test.

1. Save predicted particles as `particle_mi`, compute current likelihoods, and save the unnormalized weights as `weight_un`.
2. Resample a candidate posterior and compute `dcov_pl = det(cov(particle_pl'))`.
3. If `dcov_pl > dcov_mi`, restore the prior particles and append the measurement and prior cloud to `skipped`.
4. Otherwise, pass any queued measurements to [OOSM.m](../OOSM.m), then clear the queue.
5. Average the resulting particle positions to obtain the current estimate.

For each stored measurement `z_a`, `OOSM.m` cumulatively multiplies current weights by likelihoods evaluated at the stored particles of time `a`, normalizes, and resamples the current prior cloud. It returns the final resampled cloud.

The helper reuses stored particles and returns the final candidate. Paper Algorithms 2–3 specify conditional sampling in Eqs. (8)–(9) and an additional covariance acceptance test; see the [procedure comparison](implementation-notes.md#differences-from-the-paper).

## Comparison filters

| Method | Implemented procedure | Paper connection |
|---|---|---|
| APF | Propagate with measured velocity; use the current likelihood to resample; add process noise; apply the lookahead likelihood-ratio correction and compute a weighted estimate | Section IV-B1 describes auxiliary look-ahead resampling; the second stage now divides by the selected ancestor's lookahead likelihood |
| MPF | Estimate modes; assign particles to nearest centres; update and resample each cluster; combine cluster estimates by total likelihood | Section IV-B2 describes mixture filtering; clustering here uses critical-bandwidth mode estimation |

The MPF path is [NumMode.m](../NumMode.m) → [FindCriticalBW.m](../FindCriticalBW.m) / [Significance.m](../Significance.m) → [dskensity2d.m](../dskensity2d.m) → [Cluster.m](../Cluster.m). `NumMode.m` cites Silverman's 1981 mode-estimation method, also reference [25] in the paper.

`dskensity2d.m` multiplies separate one-dimensional kernel density estimates, assuming independent coordinates. The included [HGMeanShiftCluster.m](../HGMeanShiftCluster.m) is not called by this path.

## Metrics and contracts

| Function | Input → output | Requirements |
|---|---|---|
| `DEM_height(pos, DEM)` | `2 × 1` position → scalar height | A valid five-by-five terrain patch; see [data contract](simulation.md#terrain-data) |
| `likelihood(Z_est, Z_mea, sig_meas)` | Scalar heights and standard deviation → density | Positive `sig_meas` |
| `Resample(particle, weight)` | `d × N` particles, `1 × N` weights → same-sized particles and uniform weights | Finite nonnegative weights with positive sum |
| `OOSM(...)` | Current prior, weights, queued records, terrain → `2 × N` cloud and unchanged `oosmSucceed` | Nonempty queue with matching particle-column order and count |
| `RMSE(x_err, length, numMonte)` | `2 × K × M` errors → three `K × 1` RMSE series | Coordinate and Euclidean-distance RMSE; Eq. (12) |
| `covAnal(C)` | `2 × 2` covariance → two scalar summaries | Both entries algebraically reduce to `sqrt(trace(C))`; Eq. (13) is the run-averaged measure |

These are assumptions of the existing code, not validated input checks. Degenerate weights, boundary positions, and mode-estimation failures are covered in [implementation details](implementation-notes.md#numerical-contracts).

`Resample` optionally returns ancestor indices as its third output. [auxiliary_weights.m](../auxiliary_weights.m) normalizes the ratios `p(z|x_new)/p(z|lookahead_ancestor)` in log space after the densities are evaluated. Invalid or completely underflowed weights fail explicitly. [Regression checks](../tests/integration/auxiliary/README.md) cover ancestor correspondence and the zero-process-noise case.
