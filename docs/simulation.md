# Running the particle-filter simulation

Guide to the inputs and outputs of [main_OOSM.m](../main_OOSM.m). Follow the [installation and run commands](../README.md#installation) first; algorithm details are in the [technical reference](particle-filtering.md).

## Settings and randomness

- Edit simulation duration, timestep, particle count, and Monte Carlo count in `main_OOSM.m` under `%% setting`; that block is the canonical source of defaults.
- `x_init` and `x_init_err` set the true starting position and initial estimate offset; the sensor and filter blocks separately set simulated errors and assumed uncertainty.
- Keep `RUN_DATA` enabled for a fresh run; disabling it requires existing compatible arrays in the workspace.
- Keep the four filter switches enabled for the complete run; final plotting and reporting use their variables unconditionally.
- The plotting/reporting code assumes the original time grid, including indices `51:151` and the first 70 samples; changing duration or timestep also requires adapting those consumers.
- The script calls `rng('shuffle')` inside data-generation and filter loops; setting a seed before calling it does not make the full run deterministic.
- A reproducible adaptation must replace those reseeds with a documented seed policy, seed before all random draws, and record MATLAB/toolbox versions and settings; no seed policy is imposed by this documentation.

## Terrain data

[DB_part.mat](../DB_part.mat) supplies `DEM`; no preprocessing is performed by the entry point.

- `DEM.DB` is a finite `4801 × 1000` double terrain-height matrix; its first index is the local x coordinate and its second index is y.
- `DEM.resolution` is 30 metres, the grid spacing used to convert metre coordinates to grid indices.
- [DEM_height.m](../DEM_height.m) uses `floor(pos / DEM.resolution)` without adding a one-based offset, then interpolates a five-by-five neighbourhood.
- Each floored index must be at least 3 and at most the corresponding matrix dimension minus 2; there is no clipping or boundary recovery.
- [plot_terrain.m](../plot_terrain.m) additionally requires `x_true` from the simulation workspace and changes two local plot samples; its contour values are not untouched source data.

## Expected outputs

These outputs are identified from source inspection; they are not measured results from a verified run.

| Output | Contents |
|---|---|
| Console | Start/completion and Monte Carlo progress messages for OOSM, PF, MPF, APF; final `ePF`, `eAPF`, `eOOSM`, `eMPF` summaries |
| `result.mat` | Workspace saved before final plotting/reporting, including estimated positions, errors, particles, RMSE and covariance summaries |
| `result_mode.mat` | MPF mode significance (`signi`) and PF covariance-increase counts (`covIncrease_PF`) |
| Figures | MPF mode significance overlaid with PF ambiguity markers, four-method distance RMSE, and covariance summary versus time |

`x_est_*` and `x_err_*` have shape `2 × K × M`; `d_err_RMS_*` has shape `K × 1`. The suffix identifies the method, `K` is the number of simulated samples, and `M` is `numMonte`.

The mode-significance plot combines different filters; it does not establish a within-filter relationship between modes and covariance increases. With multiple Monte Carlo runs, `signi` retains the final MPF run while `covIncrease_PF` accumulates PF counts.

The `e*` summaries sum 101 samples and divide by 100; treat them as the script's reported statistic, not a conventional sample mean. Published plots and Table II are not expected numerical outputs of the default run.

## Small synthetic helper example

From the repository root, run this base-MATLAB example without the terrain file or additional toolboxes:

```bash
matlab -batch "rng(1); p=likelihood(100,100,2); assert(abs(p-1/sqrt(8*pi))<1e-12); [x,w]=Resample([1 2 3;4 5 6],[0 1 0]); assert(isequal(x,repmat([2;5],1,3))); assert(max(abs(w-1/3))<1e-12); disp('Particle filter helper checks passed')"
```

This demonstrates the Gaussian likelihood at zero residual and resampling with all mass on one particle. Expected output is `Particle filter helper checks passed`; this command has not been executed here and does not validate the full filter or reproduce the paper.
