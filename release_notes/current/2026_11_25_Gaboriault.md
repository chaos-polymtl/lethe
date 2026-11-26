## [Master] - 2026/11/25

### Changed

- MINOR The application tests cls_phase_change_darcy_with_density and time_dependent_boundaries_heat_transfer failed intermittently because their non-linear solvers stopped at a tolerance of 1e-8, which left results that differed slightly from the expected output from one run to the next. Their non-linear tolerances and linear-solver minimum residuals are now 1e-11, so the results no longer change between runs and the tests pass consistently. [#XXXX](https://github.com/chaos-polymtl/lethe/pull/XXXX)
