# Numerical regression checks

A seeded zero-process-noise population checks that auxiliary resampling does not count the observation twice. A nonidentity ancestor permutation checks the denominator correspondence.

From the repository root, with base MATLAB:

```bash
matlab -batch "addpath('tests/integration/auxiliary'); verify_auxiliary"
```

These checks require no external data or plotting and cover the numerical routines described above.
