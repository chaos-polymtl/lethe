# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2026/03/06

### Added

- MINOR This PR adds the Fichera oven validation test for the time-harmonic Maxwell solver. The test keep track of the DPG error norm at each h-refinement step and also output a .dat that contains the number of cell, interior dofs (even if its not the number of dofs in the linear system) the maximum error and the L2 error.[#1931](https://github.com/chaos-polymtl/lethe/pull/1931)

### Added

- MINOR This PR adds a new mortar manager class that allows a linear (straight) interface, which is useful for debugging. At the moment, the implementation supports only two-dimensional cases, and it assumes that the mortar interface is parallel to the y axis.[#1926](https://github.com/chaos-polymtl/lethe/pull/1926)

### Fixed

- MINOR The p-refinement in the multigrid preconditioner was not working with the mortar feature. When calling the function `compute_n_subdivisions_and_radius` to rotate the mapping at each level, the level triangulation was being passed as argument. This was not compatible with the p-refinement structure, since no new triangulation is created for levels with the same refinement but different p order. This PR fixes this by computing the `interface_radius` outside the levels loop, which is suitable because such variable is the same for all levels. [#1927](https://github.com/chaos-polymtl/lethe/pull/1927)

- MAJOR This refactoring changes the way the h-adaptivity is decided for the .prm by changing the `kelly` `type` to `adaptive` and adding a new parameter in the mesh adaptation subsection : `error estimator`. This new parameter follow the same logic as the variable parameter and therefore enable the use of different types of error estimators for different physics. This refactor was necessary for the Time-harmonic Maxwell solver that has a DPG built-in error indicator. This new error indicator is also introduced in the PR along with h-adaptivity for the time-harmonic auxiliary physics. Note that I also needed to change the type of element from FE_NedelecSZ to FE_Nedelec, because there is a bug with the SZ version when using multiple core in deal.ii library currently.[#1922](https://github.com/chaos-polymtl/lethe/pull/1922)

- MAJOR The CFD-DEM solvers (`CFDDEMSolver` and `CFDDEMMatrixFree`) have been refactored to centralize DEM sub-iteration management within the `SubSimulationControlDEM` class. The `coupling_frequency` and `dem_time_step`  member variables have been removed and replaced by a shared `SubSimulationControlDEM` object from which the state of the DEM steps is queried. This also fixes a bug in which particle location within the cells would not be forced at the last DEM time-step if the DEM time step was based on a fraction of the Rayleigh time step. [#1929](https://github.com/chaos-polymtl/lethe/pull/1929)

## [Master] - 2026/03/03

### Added

- MAJOR The double contact PR between solid surface and particle still had a problem. The contact_info was being passed by copy in the contact_record which caused the tangential displacement and the rolling resistance spring torque to be reset instead of properly updated. Another bug was that contact_info was not cleared when a particle was no longer in contact with a given triangle. As a result, during a subsequent contact with the same triangle, the particle could start with a non-zero tangential displacement and rolling resistance spring torque. [#1928](https://github.com/chaos-polymtl/lethe/pull/1928s)

## [Master] - 2026/02/26

### Added

- MINOR This PR adds a new parameter to the multigrid setup, named `set mg p min coarsening degree`. Combined with the p-multigrid coarsening type `decrease by one`, it allows the user to limit the lowest polynomial degree to a value greater than 1 (which was previously the default).[#1925](https://github.com/chaos-polymtl/lethe/pull/1925)

### Fixed

- MINOR The matrix-free CFD-DEM solver (lethe-fluid-particles-matrix-free) requires the hesisan of the velocity field in both the matrix and the right-hand side. Without the hessian, the solver does not converge adequately. This PR adds a check that ensures that hessians are enabled when using this solver, otherwise the hessian of the velocity is used uninitialized and this leads to NaNs. [#1924](https://github.com/chaos-polymtl/lethe/pull/1924)

### Fixed

- MAJOR The dynamic time-stepping functionality of the unresolved CFD-DEM did not work adequately because the set dem iteration control = fraction of rayleigh time could not be parsed adequately. This was due to a small spelling mistake. This PR fixes this and the mode works adequately. [#1923](https://github.com/chaos-polymtl/lethe/pull/1923)

## [Master] - 2026/02/24

### Added

- MINOR This PR adds three new simulation parameters: `set output qcriterion`, `set output vorticity`, and `set output velocity gradient`. They allow the user to enable/disable the output of the respective fields within the `pvd`/`vtu` files. As default, all of them are set to `true`. [#1919](https://github.com/chaos-polymtl/lethe/pull/1919)

## [Master] - 2026/02/23

### Changed
- MAJOR This PR refactor the grid.cc file so it now have a special enum argument for custom Lethe mesh. The type associated to custom grid is now ``lethe`` and the previous option are now handled with the ``grid_type`` following a similar convention as deal.II. The PR also add two new grid that will be used in future time-harmonic benchmark examples (i.e the fichera_oven and the uniform_channel_with meshed_cylinder). [#1914](https://github.com/chaos-polymtl/lethe/pull/1914)


## [Sample] - YYYY/MM/DD

### Added

- MAJOR/MINOR/PATCH Description (#PR).

### Changed

- MAJOR/MINOR/PATCH Description (#PR).

### Deprecated

- MAJOR/MINOR/PATCH Description (#PR).

### Removed

- MAJOR/MINOR/PATCH Description (#PR).

### Fixed

- MAJOR/MINOR/PATCH Description (#PR).

