# Particle_Filter

## Overview

MATLAB particle filtering examples for terrain-referenced navigation (TRN, also called terrain-aided navigation), comparing standard, auxiliary, mixture, and out-of-sequence measurement particle filters for nonlinear position estimation.

This is Youngjoo Kim's research code accompanying **“Utilizing Out-of-Sequence Measurement for Ambiguous Update in Particle Filtering,” published in IEEE Transactions on Aerospace and Electronic Systems (2018)**, an [established peer-reviewed journal covering aerospace systems, navigation, and target tracking](https://ieee-aess.org/publications/taes). See the [paper and citation](#citation) and [canonical repository](https://github.com/rhymesg/Particle_Filter).

Use the source and [method guide](docs/particle-filtering.md) to study ambiguity handling in terrain navigation and compare particle-filter update strategies.

## Method

When the terrain likelihood is ambiguous, retain the measurement and reconsider it after later observations provide more context. The supplied simulation compares this idea with standard, auxiliary, and mixture particle filtering.

### Algorithms and source

All four filters run in [main_OOSM.m](main_OOSM.m); the [algorithm reference](docs/particle-filtering.md) maps the paper to implementation details.

| Method | Paper location | Source entry points |
|---|---|---|
| Standard PF with sequential importance resampling (SIR) | Section II-B, Algorithm 1 | `%% PF` section; [Resample.m](Resample.m), [likelihood.m](likelihood.m) |
| Out-of-sequence measurement PF (OOSMPF) | Section III, Algorithms 2–3 | `%% OOSM` section; [OOSM.m](OOSM.m) |
| Auxiliary PF (APF) | Section IV-B1 | `%% APF` section |
| Mixture PF (MPF) | Section IV-B2; mode analysis in IV-D | `%% MPF` section; [NumMode.m](NumMode.m), [Cluster.m](Cluster.m) |
| Terrain model and error metrics | Section IV-A, IV-C, Eqs. (10)–(13) | [DEM_height.m](DEM_height.m), [RMSE.m](RMSE.m), [covAnal.m](covAnal.m) |

## Examples

Run from the repository root with MATLAB, Statistics and Machine Learning Toolbox (`ksdensity`), and Image Processing Toolbox (`imregionalmax`). The terrain file [DB_part.mat](DB_part.mat) is included; shell commands use `matlab -batch` (R2019a or later).

Run the four-filter comparison without opening figure windows:

```bash
matlab -batch "set(groot,'defaultFigureVisible','off'); main_OOSM"
```

For visible plots, select the repository root as MATLAB's current folder and enter `main_OOSM`. The script writes or overwrites `result.mat` and `result_mode.mat` in that folder.

See the [simulation guide](docs/simulation.md) for settings, data layout, outputs, randomness, and a small synthetic helper example. Keep all `RUN_*` switches enabled for the complete script; its final plots depend on all four filter outputs.

## Implementation scope

The standard and auxiliary filters provide estimation examples; the APF includes the likelihood-ratio weight correction. The OOSM helper differs from the published conditional update, and MPF mode-weight handling needs correction before quantitative reuse; see [implementation differences](docs/limitations.md#differences-from-the-paper). The paper's receding-horizon Kalman comparison is not included, and native MATLAB execution remains unverified.

### Checks

Run the focused [APF importance-weight checks](tests/integration/auxiliary/README.md), which use base MATLAB:

```bash
matlab -batch "addpath('tests/integration/auxiliary'); verify_auxiliary"
```

## Citation

If you use or adapt this work, please cite:

> Youngjoo Kim, Kyungwoo Hong, and Hyochoong Bang. “Utilizing Out-of-Sequence Measurement for Ambiguous Update in Particle Filtering.” *IEEE Transactions on Aerospace and Electronic Systems*, 54(1), 493–501, February 2018. [doi:10.1109/TAES.2017.2741878](https://doi.org/10.1109/TAES.2017.2741878).

[CITATION.cff](CITATION.cff) provides machine-readable software and publication metadata.

## License and provenance

The repository includes an [MIT license](LICENSE). See [provenance and limitations](docs/limitations.md#provenance-and-attribution) for the source revision and third-party notices.
