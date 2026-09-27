# Particle_Filter

## Overview

MATLAB particle filtering examples for terrain-referenced navigation (TRN, also called terrain-aided navigation), comparing standard, auxiliary, mixture, and out-of-sequence measurement particle filters for nonlinear position estimation.

This is Youngjoo Kim's research code accompanying **“Utilizing Out-of-Sequence Measurement for Ambiguous Update in Particle Filtering,” published in IEEE Transactions on Aerospace and Electronic Systems (2018)**, a peer-reviewed aerospace and electronic systems journal. See the [paper and citation](#citation) and [canonical repository](https://github.com/rhymesg/Particle_Filter).

The proposed method postpones locally ambiguous terrain measurements and reuses them in a later update.

Use this code as a source-level reference for the four-filter terrain simulation. Before adapting it, review the [differences from the published algorithm](docs/limitations.md#differences-from-the-paper) and [license and provenance](#license-and-provenance); exact reproduction of the paper's results and suitability for operational navigation have not been established.

## Installation

- Install MATLAB with Statistics and Machine Learning Toolbox ([`ksdensity`](https://www.mathworks.com/help/stats/ksdensity.html)) and Image Processing Toolbox ([`imregionalmax`](https://www.mathworks.com/help/images/ref/imregionalmax.html)) for the full comparison.
- The shell command uses [`matlab -batch`](https://www.mathworks.com/help/matlab/ref/matlabmacos.html), available from R2019a; put the MATLAB executable on PATH.
- The required terrain file, [DB_part.mat](DB_part.mat), is included; the ignored `ref/` folder is not required.
- No MATLAB release or Octave compatibility has been validated for this checkout.

Clone the repository:

```bash
git clone https://github.com/rhymesg/Particle_Filter.git
```

Enter the repository root:

```bash
cd Particle_Filter
```

## Usage

Run the four-filter comparison without opening figure windows:

```bash
matlab -batch "set(groot,'defaultFigureVisible','off'); main_OOSM"
```

For visible plots, select the repository root as MATLAB's current folder and enter `main_OOSM`. The script writes or overwrites `result.mat` and `result_mode.mat` in that folder.

See the [simulation guide](docs/simulation.md) for settings, data layout, outputs, randomness, and a small synthetic helper example. Keep all `RUN_*` switches enabled for the complete script; its final plots depend on all four filter outputs.

## Development

No automated MATLAB test suite is supplied. MATLAB execution and the published numerical results have not been verified in this documentation update; [verification and reuse limits](docs/limitations.md) distinguish source inspection from reproduction.

Report issues through the [issue tracker](https://github.com/rhymesg/Particle_Filter/issues), including the commit, MATLAB/toolbox versions, settings, and error or unexpected output.

## Algorithms and source

All four filters run in [main_OOSM.m](main_OOSM.m); the [algorithm reference](docs/particle-filtering.md) maps the paper to implementation details.

| Method | Paper location | Source entry points |
|---|---|---|
| Standard PF with sequential importance resampling (SIR) | Section II-B, Algorithm 1 | `%% PF` section; [Resample.m](Resample.m), [likelihood.m](likelihood.m) |
| Out-of-sequence measurement PF (OOSMPF) | Section III, Algorithms 2–3 | `%% OOSM` section; [OOSM.m](OOSM.m) |
| Auxiliary PF (APF) | Section IV-B1 | `%% APF` section |
| Mixture PF (MPF) | Section IV-B2; mode analysis in IV-D | `%% MPF` section; [NumMode.m](NumMode.m), [Cluster.m](Cluster.m) |
| Terrain model and error metrics | Section IV-A, IV-C, Eqs. (10)–(13) | [DEM_height.m](DEM_height.m), [RMSE.m](RMSE.m), [covAnal.m](covAnal.m) |

The paper's receding-horizon Kalman filter comparison is not included in this repository.

## Citation

If you use or adapt this work, please cite:

> Youngjoo Kim, Kyungwoo Hong, and Hyochoong Bang. “Utilizing Out-of-Sequence Measurement for Ambiguous Update in Particle Filtering.” *IEEE Transactions on Aerospace and Electronic Systems*, 54(1), 493–501, February 2018. [doi:10.1109/TAES.2017.2741878](https://doi.org/10.1109/TAES.2017.2741878).

For code reuse, please also reference [rhymesg/Particle_Filter](https://github.com/rhymesg/Particle_Filter) and the commit or release used. [CITATION.cff](CITATION.cff) provides machine-readable software and publication metadata; citation requests are separate from license obligations.

## License and provenance

The repository includes an [MIT license](LICENSE). [Provenance and limitations](docs/limitations.md#provenance-and-attribution) identify the source revision, preserved third-party notices, and unresolved terrain-data provenance.
