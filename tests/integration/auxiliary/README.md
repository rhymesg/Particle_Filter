# Numerical regression checks

A seeded zero-process-noise population checks that auxiliary resampling does not count the observation twice. A nonidentity ancestor permutation checks the denominator correspondence.

From the repository root, with base MATLAB:

```bash
matlab -batch "addpath('tests/integration/auxiliary'); verify_auxiliary"
```

Check bandwidth fallback centers and a transition on the final iteration:

```bash
matlab -batch "addpath('tests/integration/auxiliary'); verify_bandwidth"
```

The bandwidth check substitutes deterministic density and regional-maxima functions in a temporary folder, covering search termination with base MATLAB. Both checks restore the caller's environment and require no external data or plotting.
