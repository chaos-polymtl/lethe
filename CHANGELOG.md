# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2026/05/02

### Fixed

- MINOR Fixed several Clang warnings encountered when compiling Lethe with Apple Clang. The unused-parameter warning in `LetheGridTools::find_cells_around_edge` is silenced with `[[maybe_unused]]` (the parameter is used when compiled with OpenCascade). The unused-private-field warnings in `SubSimulationControlDEM` are resolved by reading the constructor's stored members (via `this->`) instead of the shadowing parameters. The obsolete `#if __GNUC__ > 7` guards around `std::filesystem` usage in `create_output_folder`, `delete_output_folder`, and `delete_vtu_and_pvd_files` are removed; they prevented these functions from working on Apple Clang (which defines `__GNUC__ == 4`) and are no longer needed since Lethe requires C++20. The implicit `double`-to-`int` narrowing in `TimeHarmonicMaxwell::should_solve_auxiliary_physics` is made explicit with `static_cast<int>`. [#1981](https://github.com/chaos-polymtl/lethe/pull/1981)

## [Master] - 2026/05/01

### Fixed

- MINOR This PR removes the obligation to have a constant mortar cell height in the direction parallel to the rotation axis. [#1978](https://github.com/chaos-polymtl/lethe/pull/1978)

- MAJOR This PR adds a missing factor in the uniform and Gaussian laser models to account for the reduce heat flux when there is an angle of incidence between the laser and the CLS interface. The factor corresponds to the scalar product between the normal of the interface and the laser direction unit vector. [#1977](https://github.com/chaos-polymtl/lethe/pull/1977)

## [Master] - 2026/04/30

### Added

- MAJOR This PR adds the postprocessing of the melt volume for cases with phase change. Setting `calculate geometric melt volume` to true enables the computation of the geometric melt volume. The region where the temperature of the fluid with phase change is above the `melting temperature` describes the melt subdomain. The implementation is compatible with CLS multiphase flows as long as only 1 of the 2 fluids has phase change. `monitored fluid with phase change` set the fluid monitored. `geometric melt volume name` sets the name of the file where computed melt volume values are saved.

  This PR also changes the computation of the liquid fraction with `calculate liquid fraction` to `calculate algebraic melt volume`. The liquid fraction in the domain is no longer outputted. Instead, a liquid fraction weighted melt volume is outputted in a file named with `algebraic melt volume name`. In CLS simulations, the fluid of interest must be specified with `monitored fluid with phase change`. [#1972](https://github.com/chaos-polymtl/lethe/pull/1972)

## [Master] - 2026/04/28

### Added

- MINOR Added application test (MMS case) for 3D mortar implementation. [#1976](https://github.com/chaos-polymtl/lethe/pull/1976)

### Fixed

- MINOR This PR removes the obligation to start the domain at z=0 when using the mortar method in 3D cases. The domain can now be shifted along the rotation axis without affecting the results. [#1975](https://github.com/chaos-polymtl/lethe/pull/1975)

## [Master] - 2026/04/22

### Added

- MINOR Added checkpoint functionality for the restart of the time-harmonic Maxwell solver. [#1974](https://github.com/chaos-polymtl/lethe/pull/1974)

- MINOR This PR adds a new mesh to simulate a channel with an object in it. Here the object is a rectangular prism with a square base. The name of this new mesh is "uniform_channel_with_meshed_square_prism". This will be used for benchmark and test the multiphysic microwave-heating capabilities of Lethe. The square prism can rotate and the size of its side can be adjusted in addition to the channel dimension. [#1968](https://github.com/chaos-polymtl/lethe/pull/1968)

## [Master] - 2026/04/21

### Fixed

- MINOR This PR fixes the computation of cell radial distance for 3D mortar cases, and improves the identification of the rotor/stator domains by using a material_id() flag (as opposed to comparing the radial distance between the cell center and the computed interface radius). The linear velocity computation is also generalized for 3D cases with rotation axis other than z. [#1967](https://github.com/chaos-polymtl/lethe/pull/1967)

### Added

- MINOR Added unit tests for the custom cylinder grid generator (`CylinderGrid`) and the periodic hills grid generator (`PeriodicHillsGrid`). The cylinder test exercises the four cylinder types (`cylinder_classic`, `cylinder_balanced`, `cylinder_squared`, `cylinder_regularized`) in 3D, while the periodic hills test exercises a default and a stretched configuration in both 2D and 3D (the stretched case is uniformly refined to make the stretching transformation visible in the statistics). Both tests report the number of active cells, vertices, mesh volume and boundary-id face counts, and write a VTK file per case for visual inspection. [#1973](https://github.com/chaos-polymtl/lethe/pull/1973)

## [Master] - 2026/04/20

### Added

- MINOR Added a unit test for the Birmingham fluidized bed grid generator (`GridBirminghamFluidizedBed`). The test exercises the 3D mesh with and without the chimney and with a non-zero inlet offset, reports the number of active cells, vertices, mesh volume and boundary-id face counts, and writes the generated mesh both in UCD format to the deallog and in VTK format for visual inspection. [#1971](https://github.com/chaos-polymtl/lethe/pull/1971)

## [Master] - 2026/04/16

### Changed

- MINOR The mortar TGV tests were not robust enough, and some were failing in debug only but not in release mode. This PR decreases the tolerance of both linear and non-linear solvers in the 3 mortar TGV tests. [#1970](https://github.com/chaos-polymtl/lethe/pull/1970)

## [Master] - 2026/04/10

### Fixed

- MINOR The dynamic_with_sparse_contacts load balance method could crash if load balancing was triggered at the first iteration of a restart (yes this is oddly specific). This PR fixes this by refreshing the ASC structures when the DEM is restarted. [#1964](https://github.com/chaos-polymtl/lethe/pull/1964)

## [Master] - 2026/04/09

### Added

- MAJOR This PR couples the time-harmonic electromagnetic solver with heat transfer to perform microwave-heating multiphysics simulations. In addition, it enables more complex time integration coupling strategies between the fluid dynamics and auxiliary physics solvers. At the moment, only the time-harmonic electromagnetic solver does not use the default behaviour of solving the physics at every time step. This physics has three time coupling strategies implemented: none, iteration, and time. The choice of strategy instructs the multiphysics interface on when to compute the electromagnetic field during time integration. [#1961](https://github.com/chaos-polymtl/lethe/pull/1961)

## [Master] - 2026/04/08

### Fixed

- MINOR Fixed the VANS 2D particle-location curl handling so it matches deal.II's curl result type and remains compatible with the supported deal.II versions. [#1963](https://github.com/chaos-polymtl/lethe/pull/1963)

### Added

- MINOR Sharp-IB now mirrors DEM-style particle output and post-processing more closely: IB particle PVD output follows the fluid output cadence, including step 0, and DEM force-chain and lagrangian post-processing can now be rebuilt from Sharp-IB particles. Application tests and documentation were added for the new output paths. [#1958](https://github.com/chaos-polymtl/lethe/pull/1958)

## [Master] - 2026/03/31

### Changed

- MAJOR Following PRs [#1937](https://github.com/chaos-polymtl/lethe/pull/1937), [#1938](https://github.com/chaos-polymtl/lethe/pull/1938), [#1944](https://github.com/chaos-polymtl/lethe/pull/1944) and [#1950](https://github.com/chaos-polymtl/lethe/pull/1950) this PR renames all occurrences of "interface regularization", "algebraic interface reinitialization", and "phase fraction" with "interface reinitialization", "pde-based interface reinitialization", and "phase indicator" in the code.[#1954](https://github.com/chaos-polymtl/lethe/pull/1954).

### Added

- MINOR This PR adds the possibility to scale the solution of the time-harmonic Maxwell solver (which is always solved in a dimensionless form) with three different method `electric_field`, `magnetic_field` or `power` to recover the dimensionality according to user input. Prior to this PR the outputted solution was always dimensionless. This behavior can still be obtain by using the scaling method `none`.[#1942](https://github.com/chaos-polymtl/lethe/pull/1942)

- MINOR This PR adds the option to assign a different weight for the cells located at the mortar interface, which is taken into consideration for the workload partitioning. [#1953](https://github.com/chaos-polymtl/lethe/pull/1953)

- MINOR This PR adds the capability to output the ratio of the minimum cell size to the maximum particle diameter in DEM and CFD-DEM simulations. It also updates the output related to diamond-shaped cells and the expand particle–wall contact parameter. Tests and documentation have been updated accordingly. [#1955](https://github.com/chaos-polymtl/lethe/pull/1955)

## [Master] - 2026/03/30

### Added

- MINOR Added optional distance-based coarsening for Sharp-IB local mesh refinement through the parameters `enable distance based coarsening` and `coarsen mesh outside radius factor`, together with documentation and an application test. [#1956](https://github.com/chaos-polymtl/lethe/pull/1956)

## [Master] - 2026/03/26

### Added

- MAJOR This PR adds a 2 dimensional example using the mortar method, named "Paddle Mixer Using the Mortar Method". [#1939](https://github.com/chaos-polymtl/lethe/pull/1939)

### Fixed

- MINOR The two prototype "matrix_based_advection_diffusion" and "matrix_based_non_linear_poisson" would not compile with the current deal.II master version.  Essentialy the changes to the Trilinos sparsity pattern prevented compilation. This has been fixed by adopting the process used in step-50. This has no consequence on the main Lethe software stack. [#1951](https://github.com/chaos-polymtl/lethe/pull/1951).

## [Master] - 2026/03/24

### Changed

- MAJOR Following PRs [#1937](https://github.com/chaos-polymtl/lethe/pull/1937), [#1938](https://github.com/chaos-polymtl/lethe/pull/1938), and [#1944](https://github.com/chaos-polymtl/lethe/pull/1944), this PR renames all occurrences of "VOF" with "CLS" in the code. [#1950](https://github.com/chaos-polymtl/lethe/pull/1950).

## [Master] - 2026/03/19

### Removed

- MINOR The square gas–solid fluidized bed example is removed, as PR [#1844] introduced a more comprehensive cylindrical fluidized bed simulation. This closes issue [#1860]. [#1948](https://github.com/chaos-polymtl/lethe/pull/1948).

- MAJOR The way we parametrized the display precision on the terminal in Lethe was incredibly messy. The log_precision parameter in the SimulationControl is supposed to control it, but we also have another display_precision in the non-linear solver which can overwrite it. This led to a lot of confusion. The simple thing should be that one parameter controls the log precision and that's it. This PR achieves this. It deletes the display_precision parameter from the non-linear solver parameters. Now, the usage of the log_precision is made uniform and the precision is set once at the start of the simulation. The documentation has been adapted as well.  [#1947](https://github.com/chaos-polymtl/lethe/pull/1947).

## [Master] - 2026/03/18

### Changed

- MAJOR Following PRs [#1937](https://github.com/chaos-polymtl/lethe/pull/1937) and [#1938](https://github.com/chaos-polymtl/lethe/pull/1938), this PR renames all occurrences of "VOF", "interface regularization", "algebraic interface reinitialization", and "phase fraction" with "CLS", "interface reinitialization", "pde-based interface reinitialization", and "phase indicator" respectively in the example section of the documentation. [#1944](https://github.com/chaos-polymtl/lethe/pull/1944).

### Fixed

- MAJOR The lethe-fluid-particles-matrix-free solver would not restart adequately. This was mainly caused by three factors. The first is that the vector used to project the forces was never reset to zero, and so the projection of the forces and the other variables became dependent on the history of this projection since this would affect the GMRES solver behavior. The second was that the void fraction time-history would get partially erased if a BDF2 time-stepping scheme was used. The last one was that some of the multigrid operators being used to transfer the particle-fluid force/drag information were not the correct ones. All these three things have been fixed. They don't changes tests (except a small one for some floating point error), but now the code restarts to machine precision. [#1943](https://github.com/chaos-polymtl/lethe/pull/1943).

## [Master] - 2026/03/17

### Changed

- MAJOR Following PR [#1937](https://github.com/chaos-polymtl/lethe/pull/1937), this PR renames all occurrences of "VOF", "interface regularization", "algebraic interface reinitialization", and "phase fraction" with "CLS", "interface reinitialization", "pde-based interface reinitialization", and "phase indicator" respectively in the parameters and theory sections of the documentation. [#1938](https://github.com/chaos-polymtl/7lethe/pull/1938).

## [Master] - 2026/03/15

### Added

- MINOR This PR adds to the utilities a bash script (`copy_files_for_postprocessing.sh`) that helps with copying only simulation result files of interest in a specific folder. This script can be used for identifying and copying files of importance before transferring to another machine or simply to clean up irrelevant files. This can be especially useful for parametric studies with more than necessary outputs for postprocessing. [#1941](https://github.com/chaos-polymtl/lethe/pull/1941)

## [Master] - 2026/03/12

### Changed

- MAJOR This PR renames all occurrences of "VOF", "interface regularization", "algebraic interface reinitialization", and "phase fraction" with "CLS", "interface reinitialization", "pde-based interface reinitialization", and "phase indicator" respectively in:
  - parameter names, values and descriptions;
  - parameter subsection names;
  - console output strings (including warning and exception messages) that are not related to in-code variables or classes names, and
  - .prm, .tpl, .output files in the code (application tests and examples).
  
  All modified parameter names have an alias with their previous names.
  The objective of this PR is to update all in-code user-seen aspects when running the code that refer to VOF with CLS. All related documentation will be updated in a subsequent PR. In the same mindset, all the changes in the following:
  - Filenames;
  - Classes, Structs, variables, and functions names;
  - VTK output fields, and
  - post-processing data
    
  will also be done in subsequent PRs. [#1937](https://github.com/chaos-polymtl/lethe/pull/1937)

## [Master] - 2026/03/11

### Added

- MINOR Added the `birmingham_fluidized_bed` Lethe grid type. This custom mesh generates the Birmingham fluidized bed geometry, composed of a bottom cylinder, a truncated cone, and a top cylinder merged along the x-axis with appropriate manifold and boundary IDs.

## [Master] - 2026/03/10

### Added

- MINOR This PR adds a new parameter to the .prm file: `comment message`. It belongs to the head level of parameters just like `dimension` and `print parameters`. It allows users to comment a message on the console at the very beginning of the simulation. When combined with `print parameters`, it helps to dissociate console outputs from the parameters file (if it were to ever be lost). Indeed, it allows user to describe more explicitly what they are simulating/studying/analyzing. This can be especially useful when looking back at older simulations after a while. This PR also adds Sphinx documentation for the head level parameters. [#1932](https://github.com/chaos-polymtl/lethe/pull/1932)

- MAJOR Added support to pressure boundary conditions in the matrix-free solvers including two application tests. Improved the documentation on the pressure boundary conditions. [#1933](https://github.com/chaos-polymtl/lethe/pull/1933)

### Fixed

- MINOR To avoid division by zero in navier_stokes_vof_assembler, a tolerance of DBL_MIN was used. However, since DBL_MIN is approx. 1e-300, it didn't really avoid a division by zero. It led to issues for fine meshes. This PR fixes this issue by changing DBL_MIN for 1e-14. [#1934](https://github.com/chaos-polymtl/lethe/pull/1934)

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
