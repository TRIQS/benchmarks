# Solver coverage

Wall time (seconds) per solver/model. `✓` means the solver ran but
did not report `run_time` in its `Solver_Info` group; `—` means the
solver has no result for this model.

| Model | alps_cthyb | atomdiag | cthyb | ctint | ctseg | edipack | exact | forktps | nrg | pomerol | pyed | w2dyn_cthyb |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Dimer | 63 | ✓ | 322 | 8 | — | 0 | — | — | — | ✓ | ✓ | 151 |
| Dimer_SOC | 63 | ✓ | — | 4 | — | — | — | — | — | ✓ | ✓ | — |
| Dimer_nn | 63 | ✓ | 175 | 337 | 236 | 0 | — | — | — | ✓ | ✓ | 137 |
| Hubbard_Atom | — | ✓ | — | 171 | — | 0 | ✓ | — | — | ✓ | ✓ | — |
| La2CuO4 | 64 | — | 630 | 354 | 668 | — | — | 90 | — | — | — | 57 |
| Plaquette | — | ✓ | — | 129 | — | 1 | — | — | — | — | — | — |
| Plaquette_SemiCircular | 67 | — | 997 | 66 | — | — | — | — | — | — | — | — |
| SIAM_Discrete_Bath | 63 | ✓ | 229 | 161 | 330 | 0 | — | — | — | ✓ | ✓ | 67 |
| SIAM_DynU | — | — | — | 8 | 190 | — | — | — | — | — | — | — |
| SIAM_Jperp | — | — | — | 63 | 234 | — | — | — | — | — | — | — |
| SIAM_SemiCircular | 63 | — | 28 | 46 | 73 | — | — | — | 12 | — | — | 33 |
| Sr2RuO4 | 66 | — | 1162 | — | — | — | — | — | — | — | — | 39 |
| Sr2RuO4_SOC | 65 | — | — | — | — | — | — | — | — | — | — | — |
| Trimer | 63 | ✓ | 1303 | 120 | — | 12 | — | — | — | ✓ | — | 473 |
