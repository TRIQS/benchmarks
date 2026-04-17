# Correctness summary

Max-norm deviation of each solver from the model's reference solver.
Column `ref` lists which solver was chosen as reference. `—` means the
solver did not compute this observable; `(ref)` marks the reference
itself in its own column.

| Model | Observable | ref | alps_cthyb | atomdiag | cthyb | ctint | ctseg | edipack | exact | forktps | nrg | pomerol | pyed | w2dyn_cthyb |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Dimer | `G` | `pyed` | 3.00e+00 | 7.95e-15 | 2.64e-02 | 1.16e-03 | — | 1.08e-03 | — | — | — | 2.25e-07 | (ref) | 1.83e-02 |
| Dimer | `Sigma` | `pyed` | 6.31e+01 | 1.37e-12 | 2.38e+02 | 9.73e-03 | — | 3.19e-03 | — | — | — | 7.04e-05 | (ref) | 8.74e+01 |
| Dimer_SOC | `G` | `pyed` | 2.05e+00 | 2.42e-14 | — | 1.64e-02 | — | — | — | — | — | 6.08e-07 | (ref) | — |
| Dimer_SOC | `Sigma` | `pyed` | 6.32e+01 | 3.74e-12 | — | 2.13e-01 | — | — | — | — | — | 1.97e-04 | (ref) | — |
| Dimer_nn | `G` | `pyed` | 4.00e-01 | 7.36e-15 | 1.34e-02 | 4.39e-04 | 1.60e-02 | 4.66e-04 | — | — | — | 1.68e-07 | (ref) | 1.37e-02 |
| Dimer_nn | `Sigma` | `pyed` | 8.48e+01 | 1.45e-12 | 9.35e+01 | 1.38e-02 | 6.48e+01 | 2.48e-03 | — | — | — | 5.51e-05 | (ref) | 9.26e+01 |
| Hubbard_Atom | `G` | `exact` | — | 2.22e-16 | — | 8.45e-04 | — | 2.22e-16 | (ref) | — | — | 2.22e-16 | 1.96e-16 | — |
| Hubbard_Atom | `Sigma` | `exact` | — | 1.42e-14 | — | 3.80e-02 | — | 1.43e-14 | (ref) | — | — | 1.42e-14 | 1.43e-14 | — |
| Hubbard_Atom | `chi2_d` | `exact` | — | — | — | 2.80e-03 | — | 1.90e+00 | (ref) | — | — | 8.88e-16 | — | — |
| Hubbard_Atom | `chi2_m` | `exact` | — | — | — | 3.33e-03 | — | 1.90e+00 | (ref) | — | — | 8.88e-16 | — | — |
| Hubbard_Atom | `chi2_s` | `exact` | — | — | — | 9.16e-04 | — | 2.92e-05 | (ref) | — | — | 2.92e-05 | — | — |
| Hubbard_Atom | `chi2_t` | `exact` | — | — | — | 2.92e-05 | — | 2.92e-05 | (ref) | — | — | 2.92e-05 | — | — |
| Hubbard_Atom | `chi3_d` | `exact` | — | — | — | 3.54e-02 | — | — | (ref) | — | — | 6.66e-16 | 1.24e-02 | — |
| Hubbard_Atom | `chi3_m` | `exact` | — | — | — | 1.74e-01 | — | — | (ref) | — | — | 4.00e-16 | 9.46e-03 | — |
| Hubbard_Atom | `chi3_s` | `exact` | — | — | — | 2.80e-01 | — | — | (ref) | — | — | 2.80e-01 | 2.80e-01 | — |
| Hubbard_Atom | `chi3_t` | `exact` | — | — | — | 2.80e-01 | — | — | (ref) | — | — | 2.80e-01 | 2.80e-01 | — |
| La2CuO4 | _all_ | — |  |  |  |  |  |  |  |  |  |  |  |  |
| Plaquette | `G` | `edipack` | — | 7.63e-10 | — | 1.62e-03 | — | (ref) | — | — | — | — | — | — |
| Plaquette | `Sigma` | `edipack` | — | 7.49e-09 | — | 1.53e-02 | — | (ref) | — | — | — | — | — | — |
| Plaquette | `G_w` | `edipack` | — | 1.61e-07 | — | — | — | (ref) | — | — | — | — | — | — |
| Plaquette | `chi2_d` | `edipack` | — | — | — | — | — | (ref) | — | — | — | — | — | — |
| Plaquette | `chi2_m` | `edipack` | — | — | — | — | — | (ref) | — | — | — | — | — | — |
| Plaquette | `chi2_s` | `edipack` | — | — | — | — | — | (ref) | — | — | — | — | — | — |
| Plaquette | `chi2_t` | `edipack` | — | — | — | — | — | (ref) | — | — | — | — | — | — |
| Plaquette_SemiCircular | _all_ | — |  |  |  |  |  |  |  |  |  |  |  |  |
| SIAM_Discrete_Bath | `G` | `pyed` | 1.98e-02 | 1.18e-15 | 9.37e-05 | 1.84e-05 | 1.32e-04 | 2.68e-09 | — | — | — | 8.10e-09 | (ref) | 1.15e-04 |
| SIAM_Discrete_Bath | `Sigma` | `pyed` | 2.07e+00 | 1.21e-13 | 2.72e-01 | 7.54e-03 | 2.37e-01 | 1.92e-07 | — | — | — | 1.72e-06 | (ref) | 2.37e-01 |
| SIAM_DynU | _all_ | — |  |  |  |  |  |  |  |  |  |  |  |  |
| SIAM_Jperp | _all_ | — |  |  |  |  |  |  |  |  |  |  |  |  |
| SIAM_SemiCircular | _all_ | — |  |  |  |  |  |  |  |  |  |  |  |  |
| Sr2RuO4 | _all_ | — |  |  |  |  |  |  |  |  |  |  |  |  |
| Sr2RuO4_SOC | _all_ | — |  |  |  |  |  |  |  |  |  |  |  |  |
| Trimer | `G` | `pomerol` | 2.80e+00 | 1.85e-06 | 1.77e-01 | 3.43e-02 | — | 1.18e-02 | — | — | — | (ref) | — | 1.81e-01 |
| Trimer | `Sigma` | `pomerol` | 6.48e+01 | 7.52e-04 | 4.86e+02 | 1.29e+00 | — | 2.66e-02 | — | — | — | (ref) | — | 1.06e+02 |
