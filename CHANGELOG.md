
# Change Log
All notable changes to the Lethe project will be documented in this file.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2024-03-13

### Changed

- MINOR In a similar way to "VOF Barycenter", the output of the "VOF Mass Conservation" is now displayed at every time iteration on the terminal when `calculate mass conservation` in the `post-processing` is enabled (default) and the `verbosity` is set to `verbose`. [#1062](https://github.com/lethe-cfd/lethe/pull/1062)

- MINOR A new boolean argument in `make_table` functions has been added, namely `display_scientific_notation` to decide if the display of table values is displayed at a fixed decimal or in scientific notation. [#1062](https://github.com/lethe-cfd/lethe/pull/1062)

### Added

- MINOR Assertion on fields required for physical property models were added to the code. This is done to ensure that if a physical property model necessitating a certain field is employed and that field is not defined, the exception message thrown is comprehensible. [#1065](https://github.com/lethe-cfd/lethe/pull/1065)

- MINOR A new `make_table_scalars_vectors` function has been added in `utilities.h`(`.cc`). [#1062](https://github.com/lethe-cfd/lethe/pull/1062)

## [Master] - 2024-03-11

## Changed

- MINOR The eigenvalue estimation for the multigrid preconditioners in the lethe-fluid-matrix-free solver is now done internally by the PreconditionRelaxation class instead of using a PreconditionChebyshev and an estimate omega function [#1064](https://github.com/lethe-cfd/lethe/pull/1064).

## [Master] - 2024-03-11

### Added

- MINOR Temperature-dependent `stasis constraint` is now featured in the Melting Cavity heat transfer example. [#1061](https://github.com/lethe-cfd/lethe/pull/1061)

## [Master] - 2024-03-05

### Added

- MINOR The new temperature-dependent solid domain constraints for two-phase Volume of Fluid (VOF) simulations has been implemented (`constrain_stasis_with_temperature_vof()`). To select cells on which the constraints are applied, we check if the filtered phase fraction at quadrature points in the cell are within a range of accepted values. This range of accepted values is defined with the `phase fraction tolerance` parameter. This parameter defines an absolute tolerance on filtered phase fraction. Related documentation was updated and a new application test was also added. [#1048](https://github.com/lethe-cfd/lethe/pull/1048)

- MINOR Added the capability to locally refine the mesh in the vicinity of the boundary conditions as specified by a list of boundary conditions. [#1056](https://github.com/lethe-cfd/lethe/pull/1056)

### Changed

- MINOR Previously implemented temperature-dependant solid domain constraints for one-fluid simulations (`constrain_solid_domain()`) was renamed to `constrain_stasis_with_temperature()` and its content changed a bit to avoid copied lines in `constrain_stasis_with_temperature_vof()`. [#1048](https://github.com/lethe-cfd/lethe/pull/1048)

- MINOR The parameter subsection `constrain solid domain` was renamed to `constrain stasis`. Related documentation was also updated. [#1048](https://github.com/lethe-cfd/lethe/pull/1048)

### Added - 2024-03-04

- MINOR The "file" insertion method has been added to the DEM and CFD-DEM solvers.[#1054](https://github.com/lethe-cfd/lethe/pull/1054)

### Fixed

- MINOR In LPBF simulations, recoil pressure formulation $`p_\mathrm{rec} = 0.55p_\mathrm{sat} + [1/\rho]m_\mathrm{dot}^2`$ was accounting 2 times for the term $`[1/\rho]m_\mathrm{dot}^2`$ because the latter is included in the term $`0.55p_\mathrm{sat}`$. The formulation is now corrected and reads $`p_\mathrm{rec} = 0.55p_\mathrm{sat}`$.

## [Master] - 2024-03-02

### Fixed

- MINOR The previously implemented weak compressibility ([#815](https://github.com/lethe-cfd/lethe/pull/815)) was not working with the heat transfer (HT) physic. The missing pressure field was added to scratch data. Two application tests were added to test the compressibility with HT for single-fluid and VOF simulations. [#1051](https://github.com/lethe-cfd/lethe/pull/1051)

## [Master] - 2024-02-27

### Added

- MINOR Temperature-dependant solid domain constraints were added for one-fluid simulations. In order to improve calculation time and linear system conditioning, homogeneous constraints are imposed on velocity DOFs in a defined solid domain through a temperature range. Additionally, pressure DOFs in the same solid domain that are not connected to fluid cells are also imposed with homogeneous constraints. The temperature range is specified through a new parameter subsection, namely `constrain solid domain`. The use of a `dynamic_zero_constraints` AffineConstraints object, rather than reconstruction of the `zero_constraints` one at every time step was determined to be optimal after some profiling tests. Documentation and a new application test have been also added. [#1038](https://github.com/lethe-cfd/lethe/pull/1038)

### Changed

- MINOR The `Variable` enum class was moved from `parameters.h` to `multiphysics.h`. Using variables from the `Variable` enum class, the constrained and constraining fields could be parametrized in the future according to user needs. [#1038](https://github.com/lethe-cfd/lethe/pull/1038)

### Fixed

- MINOR In `establish_solid_domain()` the pressure was not constrained in solid cells that were not connected to a fluid due to a small logic mistake. This was corrected. [#1038](https://github.com/lethe-cfd/lethe/pull/1038)

## [Master] - 2024-02-26

### Fixed

- MAJOR The first level of the fine search candidates for particle-particle, particle-wall and particle-floating walls was never cleared even when the value of the unordered_map became an empty unordered_map. In simulations with load balancing or particles that were moving significantly, this could potentially lead to a scenario where the size of the first level of the fine search candidate became equal to the number of total particles in the simulation, potentially leading to a crash.

## [Master] - 2024-02-09

### Fixed

- MAJOR Restarts using lethe-fluid-particles with "void fraction time derivative = true" would be incoherent because the void fraction was reinitialized when restarting, leading to a wrong time derivative of the void fraction. This has been patched.

## [Master] - 2024-01-31

### Changed

- MAJOR Robin boundary condition was renamed to incorporate the imposed heat flux (Neumann). The corresponding type parameter was ``convection-radiation`` and the current is ``convection-radiation-flux``. All application tests were updated accordingly.

## [Master] - 2024-01-28

### Fixed

- MAJOR All application_tests that use files have now mpirun=1 in their output to ensure that they are run from a mpirun=1 folder. Since https://github.com/dealii/dealii/pull/16551 deal.II 9.6 uses a serial folder for serial tests. However, deal.II 9.5 does not. To maintain portability between the two versions, we manually force a folder called mpirun=1 when running tests and adapt the copy of the files to this folder. This is a temporary fix that can be reverted once we drop support for deal.II 9.5

## [Master] - 2024-01-24

### Added

- MINOR Darcy-like penalization is now extended to two-fluid VOF simulations with the `PhaseChangeDarcyVOFAssembler` assembler class. Related documentation has been updated and a new application test was added. [#990](https://github.com/lethe-cfd/lethe/pull/990)

- MINOR Extrapolation in time of temperature and temperature gradient values has been added. This is used for source terms in the Navier-Stokes equations requiring them (Marangoni and evaporation recoil force terms). [#994](https://github.com/lethe-cfd/lethe/pull/994)

### Changed

- MINOR Velocity extrapolation has been moved to `include/solvers/vof_scratch_data.h`. Solution extrapolations are now located in `ScratchData` to avoid constraint on assemblers order. [#994](https://github.com/lethe-cfd/lethe/pull/994)

## [Master] - 2024-01-24

### Fixed

- MINOR In the previous implementation of the heat flux postprocessor ([#953](https://github.com/lethe-cfd/lethe/pull/953)), there was a confusion related to the `material_id` of a cell. As a cell's `material_id` is defined in the mesh, fluids moving from a cell to another did not have the right `material_id` when outputted. Furthermore, some deal.II meshes have an intrinsic way of numbering `material_ids` (e.g. colorized subdivided_hyper_rectangle) that were leading to inadequate comparison of the `material_ids`. As of now, in the heat flux postprocessor, the `material_id` is only used to indentify cell material when solids are present. When simulating with fluids and solids, cells of the mesh in which the fluids lie should have a `material_id` of `0`. For _fluids-only_ simulations, their heat fluxes are outputted on the entire domain, as the current `DataPostprocessor` can only hold 1 `dof_handler`. During postprocessing, the user could clip the domain of interest of each fluid using the `phase` or `filtered_phase` fields outputted. [#988](https://github.com/lethe-cfd/lethe/pull/988)

## [Master] - 2023-12-27

### Changed

- MINOR  Phase change thermal expansion now considers the solid thermal expansion within the mushy zone (between Tsolidus and Tliquidus) instead of using a blending rule. This leads to more stable results in the melting cavity benchmark

## [Master] - 2023-12-22

### Fixed

- MINOR  PVD handles for solid DEM objects were written with the wrong text format ("_" separators instead of ".")

### [Master] - Changed - 2023-12-17

- MINOR The "insertion random number range" and "insertion random number seed" parameters got renamed to "insertion maximum offset" and "insertion prn seed" respectively. The old names didn't make sense, as they're not random because they're defined in the parameter file. [#970](https://github.com/lethe-cfd/lethe/pull/896)

## [Master] - 2023-12-11

### Fixed

- MINOR The number of remaining particle to insert of each type is being checkpointed adequately. This means that no modification are required to the "number of particle" parameter after restarting a simulation. [#964](https://github.com/lethe-cfd/lethe/pull/964)

### Fixed

- MINOR Solid objects can now be restarted adequately in DEM. They will resume at the position they had at the end of the simulation [#959](https://github.com/lethe-cfd/lethe/pull/959)

## [Master] - 2023-11-27

### Added

- MINOR A few tables were omitted from being "checkpointed" causing a lost of valuable information when simulations were restarted (Issue: [#916](https://github.com/lethe-cfd/lethe/issues/916)). The missing tables have been added to the `write_checkpoint()` and `read_checkpoint()` of each physic. [#938](https://github.com/lethe-cfd/lethe/pull/938)

## [Master] - 2023-11-27

### Removed

- MINOR The average diameter of the uniform size distribution with the DEM module was specified with an "average diameter" parameter. It is now specified directly from the diameter parameter. This is correctly documented. [#940](https://github.com/lethe-cfd/lethe/pull/940).

### Fixed

- MINOR The DEM time step verification was outputting the most permissive time step (the biggest) and not the most restrictive (the smallest). This bugfix doesn't affect the uniform particle size simulation. [#939](https://github.com/lethe-cfd/lethe/pull/939).


## [Master] - 2023-11-23

### Fixed

- MINOR The plane insertion for the DEM was only supporting the uniform diameter distribution. Now it supports all types of distribution.

## [Master] - 2023-11-16

### Changed

- MINOR The maximum number of boundary conditions for all physics was fixed to 14 since the boundary conditions had to be declared before being parsed. A new mechanism is now in place which parses the "number" parameter for each physics and keeps the maximal value. Then, this maximal value is used to pre-declare the boundary conditions. This enables much more robust sanity checking of the input parameter. The major drawback of this (and this is a major one) is that if we ever have another parameter with the name "number" then this parameter would also be parsed and used to establish the maximum number of boundary conditions. In this case, the best approach would be to replace "number" with "number of boundary conditions" in the parameter file. I (B-saurus-rex) did not want to do this at the time of this change to not have a massive PR which breaks every parameter files.

- MAJOR The "number" parameter within "subsection lagrangian physical properties" and "particle type n" was changed to "number of particles" to prevent confusions with the "number" used for boundary conditions. The "number" for boundary conditions will be changed to "number of boundary conditions" in the near future.

## [Master] - 2023-11-12

### Deprecated

- MINOR The uniform insertion method had been removed. The non-uniform insertion method has been renamed to volume method to remain coherent with the plane method. If you want to use an insertion method equivalent to the uniform insertion method, use the volume method with an "insertion random number range " equal to zero. [#926](https://github.com/lethe-cfd/lethe/pull/926)


## [Master] - 2023-10-01

### Fixed

- MINOR The calculation of the velocity on rotating walls was calculated inadequately [#896](https://github.com/lethe-cfd/lethe/pull/896): The velocity of rotating boundary conditions (e.g. boundaries of the mesh) was calculated inadequately in the case where the normal vector of the wall was aligned with the rotation axis. The whole calculation procedure was slightly messed up and only worked for cylinders. This has been fixed and made general, paving the way for a full refactor of the calculation of the particle-wall contact force.

## [Master] - 2023-09-30

### Fixed

- MINOR Affine constraints used to bound the void fraction and used as boundary conditions within the heat transfer and the block Navier-stokes solver were corrected [#885](https://github.com/lethe-cfd/lethe/pull/895): There was an error in the application of the affine constraints used to clip the void fraction in the lethe-fluid-vans and lethe-fluid-particles solvers. This led to assertions being thrown in debug. This has been corrected by reiniting the constraints with the appropriate size. For the heat transfer and the block Navier-stokes, the issue was that the constraints were never reinited with the correct size, so they contained ghosted elements. This was caught by a new assert introduced in 2023-09 within deal.II master.

## [Master] - 2023-09-20

### Deprecated

- MINOR The calculation of the source term was enabled using a parameter called "enable". This parameter was used in some physics, not in others and was poorly implemented. We deprecate the usage of this parameter and always enable source term, considering the fact that the default value of the source term is zero anyway. This prevents the false perception that source terms could be enabled or disabled, while the behavior was inconsistent across physics.


## [Master] - 2023-09-19

### Changed

- MAJOR All the applications were renamed [#882](https://github.com/lethe-cfd/lethe/pull/882): `gls_navier_stokes` is now `lethe-fluid`, `gd_navier_stokes` is now `lethe-fluid-block`, `nitsche_navier_stokes` is now `lethe-fluid-nitsche`, `gls_sharp_navier_stokes` is now `lethe-fluid-sharp`, `gls_vans` is now `lethe-fluid-vans`, `mf_navier_stokes` is now `lethe-fluid-matrix-free`, `dem` is now `lethe-particles`, `cfd_dem_coupling` is now `lethe-fluid-particles`, `rpt3d` is now `lethe-rpt-3d`, `rpt_cell_reconstruction_3d` is now `lethe-rpt-cell-reconstruction-3d`, `rpt_fem_reconstruction_3d` is now `lethe-rpt-fem-reconstruction-3d`, and `rpt_l2_projection_3d` is now `lethe-rpt-l2-projection-3d`.


## [Master] - 2023-10-30

- MINOR The rotational vector for the rotational boundary condition in the lethe-particles solver is now define with one line in the parameters file. [#920](https://github.com/lethe-cfd/lethe/pull/920)



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
