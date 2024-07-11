
# Change Log
All notable changes to the Lethe project will be documented in this file.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2024-07-09

### Added

- MINOR The multigrid output now also prints the workload imbalance and vertical communication efficiency of the multigrid hierarchy being used. In addition, the mulrigrid timers were corrected. [#1194](https://github.com/chaos-polymtl/lethe/pull/1194)

## [Master] - 2024-07-09

### Fixed

- MINOR The lethe-fluid-nitsche solver could not checkpoint once it had been restarted. This has been fixed by ensuring that the particles are properly checkpointed and read from the restart files. [#1192](https://github.com/chaos-polymtl/lethe/pull/1192)

## [Master] - 2024-07-08

### Changed

- MINOR  The timer output for the lethe-fluid and lethe-fluid-matrix-free applications contains now more detailed elements. Including post-processing capabilities that were not timed, the time spent in reading mesh and manifolds, and setting the initial conditions. [#1187](https://github.com/chaos-polymtl/lethe/pull/1187) 

## [Master] - 2024-07-05

### Changed

- MINOR  The NS-VOF assembler was not optimized as the NS assembler. The if conditions on the dof components were added to not compute multiplications that are a priori 0 in the matrix and rhs assembly. [#1188](https://github.com/chaos-polymtl/lethe/pull/1188)

## [Master] - 2024-07-04

### Changed

- MINOR  The `entry_string_to_tensor` functions were duplicated; they were in `utilities.h` and `parameters.h`. In one case, it was templated with `spacedim` and the other not. This PR merges both functions in `utilities.h` as `value_string_to_tensor`. Since these are small functions, they are now inline functions. [#1186](https://github.com/chaos-polymtl/lethe/pull/1186)

## [Master] - 2024-07-03

### Fixed

- MINOR The term associated with the viscosity jump in the strong residual for the NS-VOF assembler added in [#1049] (https://github.com/chaos-polymtl/lethe/pull/1149) led to an ill-posed formulation in solid-liquid phase change cases where the solid is represented by a highly viscous fluid. Hence, this term was removed from the current formulation. [#1180](https://github.com/chaos-polymtl/lethe/pull/1180)

### Changed

- MINOR LPBF benchmark example prm files are up to date with the current working version. [#1180](https://github.com/chaos-polymtl/lethe/pull/1180)

### Added

- MINOR There was no appropriate application test for lpbf/phase change with vof cases that could have allowed us to identify this ill-posed formulation prior to the merge of [#1049] (https://github.com/chaos-polymtl/lethe/pull/1149). A new application test is now included to avoid future mistake in lpbf/phase change with vof cases. [#1180](https://github.com/chaos-polymtl/lethe/pull/1180)

## [Master] - 2024-07-02

### Fixed

- MINOR Box refinement was applied when restarting simulations, changing the restarting triangulation and leading to undesirable behavior. A boolean stating if the simulation restarts is now checked before applying the box refinement. [#1184](https://github.com/chaos-polymtl/lethe/pull/1184) 

## [Master] - 2024-06-16

### Changed

- MINOR Some application-test restart files have been updated using p4est 2.3.6. Some test results have changed for lethe-fluid-particles and lethe-fluid-vans, since the DEM solver has slightly changed since the previous restart files generation, and it is now impossible to regenerate the exact same initial condition. [#1181](https://github.com/chaos-polymtl/lethe/pull/1181)


## [Master] - 2024-06-16

### Changed

- MINOR A new subsection ``solid surfaces`` needs to be used when defining solid objects in the parameter file. For more information, see the solid object documentation. [#1169](https://github.com/chaos-polymtl/lethe/pull/1169)
- MINOR The center of rotation of a solid object is no longer being defined using a subsection and three parameters. It is now defined with one parameter and a list of doubles. For more information, see the solid object documentation. [#1169](https://github.com/chaos-polymtl/lethe/pull/1169)

## [Master] - 2024-06-13

### Added

- MINOR The DEM solver supports deprecated parameters when 3 individual component parameters are changed to a list of values parameter. [#1171](https://github.com/chaos-polymtl/lethe/pull/1171)

### Changed

- MINOR The parameters <gx>, <gy> and <gz> of the DEM solver are changed to the parameter <g>. [#1171](https://github.com/chaos-polymtl/lethe/pull/1171)

## [Master] - 2024-06-05

### Added

- MINOR A cylindrical manifold type has been added in the manifold parameter file option. Input arguments arg1, arg2, arg3, etc. used to describe the manifold geometry were modified to point coordinates and direction vector. A manifold section has been added to the documentation [#1167](https://github.com/chaos-polymtl/lethe/pull/1167)

## [Master] - 2024-05-27

### Changed

- MINOR The "insertion file name" parameter has been renamed to "list of input files". It now supports a list of input files that will be used in order when inserting particles. [#1135](https://github.com/chaos-polymtl/lethe/pull/1135)

## [Master] - 2024-05-26

### Fixed

- MINOR Simulation restart would not work adequately when a boundary condition was time-dependent and this would lead unstable simulation restarts. [#1158](https://github.com/chaos-polymtl/lethe/pull/1158) 

## [Master] - 2024-05-23

### Added

- MINOR dependency on muParser. This requires to load an additional module "muparser/2.3.2" when compiling deal.II and Lethe on clusters. [#1143](https://github.com/chaos-polymtl/lethe/pull/1143) 

## [Master] - 2024-05-20

### Removed

- MINOR The initial conditions executable which could be used to test some initial conditions has now been removed completely. This executable was never used and was unstable to the multiphysics interface. [#1148](https://github.com/chaos-polymtl/lethe/pull/1148)

## [Master] - 2024-05-13

### Changed

- MAJOR The source term specifications for the Navier-Stokes equations were using the label "xyz" which made no sense. This has been changed to "fluid dynamics" [#1130](https://github.com/chaos-polymtl/lethe/pull/1130)

## [Master] - 2024-05-09

### Added

- MAJOR Added advection of particles in the Adaptive Sparse Contacts (ASC) feature in DEM and CFD-DEM (Disabling Contacts is renamed Adaptive Sparse Contacts). [#1113](https://github.com/chaos-polymtl/lethe/pull/1113)

## [Master] - 2024-05-09

### Added

- MINOR Added a condition to the generation of .vtu and .pvd files in addition to the output generation frequency. If specified in the parameter file, only the results within the specified simulation time interval will be generated. [#1120] (https://github.com/chaos-polymtl/lethe/pull/1120) 


## [Master] - 2024-05-02

### Fixed

- MAJOR Periodic boundary conditions for the DEM were not working if they weren't the last boundary conditions being declare in the parameter file. Now, every boundary condition work in which even order they are being declared. [#1110](https://github.com/chaos-polymtl/lethe/pull/1110)

## [Master] - 2024-04-30

### Fixed

- MINOR The lethe-fluid-nitsche solver was unable to restart when the immersed triangulation was made of simplices. This has been fixed, however, mesh adaptation crashes when it is done after the restart process. I (problembär) have an idea why (the previous particle_handler of the previous checkpoint is still registered somehow in the triangulation), but I will need more time to come up with an adequate solution. [#1106](https://github.com/chaos-polymtl/lethe/pull/1106)

## [Master] - 2024-04-24

### Added

- MAJOR A new parameter to reuse the preconditioner between Newton iterations was added. Therefore, the logic of the Newton solver was also modified as it required a different call for the assembly of the matrix and the set up of the preconditioner. These calls were separated for all relevant physics solvers. [#1102](https://github.com/chaos-polymtl/lethe/pull/1102)

## [Master] - 2024-04-24

### Changed

- MAJOR For the multigrid preconditioners of the matrix-free application a new class was added, where the constructor sets the operators, constraints and transfers, and the initialize function sets the smoother, coarse-grid solver and final multigrid object. This allows to reuse the operators, constraints and transfers whenever possible reducing computational cost. [#1102](https://github.com/chaos-polymtl/lethe/pull/1102)

## [Master] - 2024-04-18

### Fixed

- MINOR Checkpoints were not established correctly with the lethe-fluid-particles solver. The void fraction was delayed in time and the restart file would not yield the same result because of this, resulting in crashes when the time derivative of the void fraction was used in the calculations. [#1096](https://github.com/chaos-polymtl/lethe/pull/1096)

## [Master] - 2024-04-17

### Added

- MINOR Added Cahn-Hilliard equations energy computation and output. The documentation was updated accordingly. [#1095](https://github.com/chaos-polymtl/lethe/pull/1095)

## [Master] - 2024-04-16

### Fixed

- MINOR Simulations with adaptive time-stepping would checkpoint correctly, but would not restart with the correct time-step due to the way the time-step was read from the checkpointing files. This is now fixed and tested. [#1093](https://github.com/chaos-polymtl/lethe/pull/1093)

## [Master] - 2024-04-15

### Fixed

- MINOR The previous flow rate post processing tool would add a line to the .dat file for each defined boundary (in the boundary condition subsection of the .prm file) for every time iteration. This is not the expected output as one would expect a single line to be added for each time step with the flow rates of every boundary. This is now fixed. [#1092](https://github.com/chaos-polymtl/lethe/pull/1092)

## [Master] - 2024-04-07

### Added

- MINOR The lethe-fluid-matrix-free solver would not work adequately on mesh that were not refined homogenously if all the boundary conditions of the domain were periodic or dirichlet boundary condition. This is fixed by adding a new parameter to fix the pressure constant in the coarse grid solver. This prevents the coarse grid linear solver from failing. [#1089](https://github.com/chaos-polymtl/lethe/pull/1089)

## [Master] - 2024-04-05

### Changed

- MINOR The output at every iteration in the terminal of the non-linear solver now displays the norms of the correction for both the velocity and pressure contributions for the Navier-Stokes solver when `verbosity` in the `non-linear solver` subsection is set to `verbose`. [#1076](https://github.com/chaos-polymtl/lethe/pull/1076)

## [Master] - 2024-04-04

### Added

- MINOR Added the turbulent Taylor-Couette incompressible flow example. [#1083](https://github.com/chaos-polymtl/lethe/pull/1083)

## [Master] - 2024-03-28

### Changed

- MINOR The ``velocity x/y/z`` parameters have been replaced by one parameter named ``initial velocity``. The ``omega x/y/z`` parameters have been replaced by one parameter named ``initial angular velocity``. [#1078](https://github.com/chaos-polymtl/lethe/pull/1078)
 [#1078](https://github.com/chaos-polymtl/lethe/pull/1078)


## [Master] - 2024-03-20

### Changed

- MINOR The ``insertion box minimum/maximum x/y/z`` parameters have been replaced by one parameter named ``insertion box points coordinates``. The ``insertion first/second/third direction `` parameters have been replaced by one parameter named ``insertion direction sequence``. [#1074](https://github.com/chaos-polymtl/lethe/pull/1074)

## [Master] - 2024-03-19

### Added

- MINOR The mass conservation feature of the VOF model now also monitors the momentum of both phases. This should be useful for some cases and improve the sharpening mechanism. [#1073](https://github.com/chaos-polymtl/lethe/pull/1073)

## [Master] - 2024-03-16

### Fixed

- MINOR The recently implemented VOF stasis constraint ([#1048](https://github.com/chaos-polymtl/lethe/pull/1048)) was not compatible with adaptive mesh refinements and checkpointing system. When refining the mesh or when restarting a simulation, the filtered phase fraction was reinitialized, but not filled with the appropriate values. This issue is now solved by adding `apply_phase_filter()` to `post_mesh_adaptation()` and `read_checkpoint()` VOF functions. [#1072](https://github.com/chaos-polymtl/lethe/pull/1072)

## [Master] - 2024-03-13

### Changed

- MINOR In a similar way to "VOF Barycenter", the output of the "VOF Mass Conservation" is now displayed at every time iteration on the terminal when `calculate mass conservation` in the `post-processing` is enabled (default) and the `verbosity` is set to `verbose`. [#1062](https://github.com/chaos-polymtl/lethe/pull/1062)

- MINOR A new boolean argument in `make_table` functions has been added, namely `display_scientific_notation` to decide if the display of table values is displayed at a fixed decimal or in scientific notation. [#1062](https://github.com/chaos-polymtl/lethe/pull/1062)

### Added

- MINOR Assertion on fields required for physical property models were added to the code. This is done to ensure that if a physical property model necessitating a certain field is employed and that field is not defined, the exception message thrown is comprehensible. [#1065](https://github.com/chaos-polymtl/lethe/pull/1065)

- MINOR A new `make_table_scalars_vectors` function has been added in `utilities.h`(`.cc`). [#1062](https://github.com/chaos-polymtl/lethe/pull/1062)

## [Master] - 2024-03-11

### Changed

- MINOR The eigenvalue estimation for the multigrid preconditioners in the lethe-fluid-matrix-free solver is now done internally by the PreconditionRelaxation class instead of using a PreconditionChebyshev and an estimate omega function [#1064](https://github.com/chaos-polymtl/lethe/pull/1064).

## [Master] - 2024-03-11

### Added

- MINOR Temperature-dependent `stasis constraint` is now featured in the Melting Cavity heat transfer example. [#1061](https://github.com/chaos-polymtl/lethe/pull/1061)

## [Master] - 2024-03-05

### Added

- MINOR The new temperature-dependent solid domain constraints for two-phase Volume of Fluid (VOF) simulations has been implemented (`constrain_stasis_with_temperature_vof()`). To select cells on which the constraints are applied, we check if the filtered phase fraction at quadrature points in the cell are within a range of accepted values. This range of accepted values is defined with the `phase fraction tolerance` parameter. This parameter defines an absolute tolerance on filtered phase fraction. Related documentation was updated and a new application test was also added. [#1048](https://github.com/chaos-polymtl/lethe/pull/1048)

- MINOR Added the capability to locally refine the mesh in the vicinity of the boundary conditions as specified by a list of boundary conditions. [#1056](https://github.com/chaos-polymtl/lethe/pull/1056)

### Changed

- MINOR Previously implemented temperature-dependant solid domain constraints for one-fluid simulations (`constrain_solid_domain()`) was renamed to `constrain_stasis_with_temperature()` and its content changed a bit to avoid copied lines in `constrain_stasis_with_temperature_vof()`. [#1048](https://github.com/chaos-polymtl/lethe/pull/1048)

- MINOR The parameter subsection `constrain solid domain` was renamed to `constrain stasis`. Related documentation was also updated. [#1048](https://github.com/chaos-polymtl/lethe/pull/1048)

### Added - 2024-03-04

- MINOR The "file" insertion method has been added to the DEM and CFD-DEM solvers.[#1054](https://github.com/chaos-polymtl/lethe/pull/1054)

### Fixed

- MINOR In LPBF simulations, recoil pressure formulation $`p_\mathrm{rec} = 0.55p_\mathrm{sat} + [1/\rho]m_\mathrm{dot}^2`$ was accounting 2 times for the term $`[1/\rho]m_\mathrm{dot}^2`$ because the latter is included in the term $`0.55p_\mathrm{sat}`$. The formulation is now corrected and reads $`p_\mathrm{rec} = 0.55p_\mathrm{sat}`$.

## [Master] - 2024-03-02

### Fixed

- MINOR The previously implemented weak compressibility ([#815](https://github.com/chaos-polymtl/lethe/pull/815)) was not working with the heat transfer (HT) physic. The missing pressure field was added to scratch data. Two application tests were added to test the compressibility with HT for single-fluid and VOF simulations. [#1051](https://github.com/chaos-polymtl/lethe/pull/1051)

## [Master] - 2024-02-27

### Added

- MINOR Temperature-dependant solid domain constraints were added for one-fluid simulations. In order to improve calculation time and linear system conditioning, homogeneous constraints are imposed on velocity DOFs in a defined solid domain through a temperature range. Additionally, pressure DOFs in the same solid domain that are not connected to fluid cells are also imposed with homogeneous constraints. The temperature range is specified through a new parameter subsection, namely `constrain solid domain`. The use of a `dynamic_zero_constraints` AffineConstraints object, rather than reconstruction of the `zero_constraints` one at every time step was determined to be optimal after some profiling tests. Documentation and a new application test have been also added. [#1038](https://github.com/chaos-polymtl/lethe/pull/1038)

### Changed

- MINOR The `Variable` enum class was moved from `parameters.h` to `multiphysics.h`. Using variables from the `Variable` enum class, the constrained and constraining fields could be parametrized in the future according to user needs. [#1038](https://github.com/chaos-polymtl/lethe/pull/1038)

### Fixed

- MINOR In `establish_solid_domain()` the pressure was not constrained in solid cells that were not connected to a fluid due to a small logic mistake. This was corrected. [#1038](https://github.com/chaos-polymtl/lethe/pull/1038)

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

- MINOR Darcy-like penalization is now extended to two-fluid VOF simulations with the `PhaseChangeDarcyVOFAssembler` assembler class. Related documentation has been updated and a new application test was added. [#990](https://github.com/chaos-polymtl/lethe/pull/990)

- MINOR Extrapolation in time of temperature and temperature gradient values has been added. This is used for source terms in the Navier-Stokes equations requiring them (Marangoni and evaporation recoil force terms). [#994](https://github.com/chaos-polymtl/lethe/pull/994)

### Changed

- MINOR Velocity extrapolation has been moved to `include/solvers/vof_scratch_data.h`. Solution extrapolations are now located in `ScratchData` to avoid constraint on assemblers order. [#994](https://github.com/chaos-polymtl/lethe/pull/994)

## [Master] - 2024-01-24

### Fixed

- MINOR In the previous implementation of the heat flux postprocessor ([#953](https://github.com/chaos-polymtl/lethe/pull/953)), there was a confusion related to the `material_id` of a cell. As a cell's `material_id` is defined in the mesh, fluids moving from a cell to another did not have the right `material_id` when outputted. Furthermore, some deal.II meshes have an intrinsic way of numbering `material_ids` (e.g. colorized subdivided_hyper_rectangle) that were leading to inadequate comparison of the `material_ids`. As of now, in the heat flux postprocessor, the `material_id` is only used to indentify cell material when solids are present. When simulating with fluids and solids, cells of the mesh in which the fluids lie should have a `material_id` of `0`. For _fluids-only_ simulations, their heat fluxes are outputted on the entire domain, as the current `DataPostprocessor` can only hold 1 `dof_handler`. During postprocessing, the user could clip the domain of interest of each fluid using the `phase` or `filtered_phase` fields outputted. [#988](https://github.com/chaos-polymtl/lethe/pull/988)

## [Master] - 2023-12-27

### Changed

- MINOR  Phase change thermal expansion now considers the solid thermal expansion within the mushy zone (between Tsolidus and Tliquidus) instead of using a blending rule. This leads to more stable results in the melting cavity benchmark

## [Master] - 2023-12-22

### Fixed

- MINOR  PVD handles for solid DEM objects were written with the wrong text format ("_" separators instead of ".")

### [Master] - Changed - 2023-12-17

- MINOR The "insertion random number range" and "insertion random number seed" parameters got renamed to "insertion maximum offset" and "insertion prn seed" respectively. The old names didn't make sense, as they're not random because they're defined in the parameter file. [#970](https://github.com/chaos-polymtl/lethe/pull/896)

## [Master] - 2023-12-11

### Fixed

- MINOR The number of remaining particle to insert of each type is being checkpointed adequately. This means that no modification are required to the "number of particle" parameter after restarting a simulation. [#964](https://github.com/chaos-polymtl/lethe/pull/964)

### Fixed

- MINOR Solid objects can now be restarted adequately in DEM. They will resume at the position they had at the end of the simulation [#959](https://github.com/chaos-polymtl/lethe/pull/959)

## [Master] - 2023-11-27

### Added

- MINOR A few tables were omitted from being "checkpointed" causing a lost of valuable information when simulations were restarted (Issue: [#916](https://github.com/chaos-polymtl/lethe/issues/916)). The missing tables have been added to the `write_checkpoint()` and `read_checkpoint()` of each physic. [#938](https://github.com/chaos-polymtl/lethe/pull/938)

## [Master] - 2023-11-27

### Removed

- MINOR The average diameter of the uniform size distribution with the DEM module was specified with an "average diameter" parameter. It is now specified directly from the diameter parameter. This is correctly documented. [#940](https://github.com/chaos-polymtl/lethe/pull/940).

### Fixed

- MINOR The DEM time step verification was outputting the most permissive time step (the biggest) and not the most restrictive (the smallest). This bugfix doesn't affect the uniform particle size simulation. [#939](https://github.com/chaos-polymtl/lethe/pull/939).


## [Master] - 2023-11-23

### Fixed

- MINOR The plane insertion for the DEM was only supporting the uniform diameter distribution. Now it supports all types of distribution.

## [Master] - 2023-11-16

### Changed

- MINOR The maximum number of boundary conditions for all physics was fixed to 14 since the boundary conditions had to be declared before being parsed. A new mechanism is now in place which parses the "number" parameter for each physics and keeps the maximal value. Then, this maximal value is used to pre-declare the boundary conditions. This enables much more robust sanity checking of the input parameter. The major drawback of this (and this is a major one) is that if we ever have another parameter with the name "number" then this parameter would also be parsed and used to establish the maximum number of boundary conditions. In this case, the best approach would be to replace "number" with "number of boundary conditions" in the parameter file. I (B-saurus-rex) did not want to do this at the time of this change to not have a massive PR which breaks every parameter files.

- MAJOR The "number" parameter within "subsection lagrangian physical properties" and "particle type n" was changed to "number of particles" to prevent confusions with the "number" used for boundary conditions. The "number" for boundary conditions will be changed to "number of boundary conditions" in the near future.

## [Master] - 2023-11-12

### Deprecated

- MINOR The uniform insertion method had been removed. The non-uniform insertion method has been renamed to volume method to remain coherent with the plane method. If you want to use an insertion method equivalent to the uniform insertion method, use the volume method with an "insertion random number range " equal to zero. [#926](https://github.com/chaos-polymtl/lethe/pull/926)


## [Master] - 2023-10-01

### Fixed

- MINOR The calculation of the velocity on rotating walls was calculated inadequately [#896](https://github.com/chaos-polymtl/lethe/pull/896): The velocity of rotating boundary conditions (e.g. boundaries of the mesh) was calculated inadequately in the case where the normal vector of the wall was aligned with the rotation axis. The whole calculation procedure was slightly messed up and only worked for cylinders. This has been fixed and made general, paving the way for a full refactor of the calculation of the particle-wall contact force.

## [Master] - 2023-09-30

### Fixed

- MINOR Affine constraints used to bound the void fraction and used as boundary conditions within the heat transfer and the block Navier-stokes solver were corrected [#885](https://github.com/chaos-polymtl/lethe/pull/895): There was an error in the application of the affine constraints used to clip the void fraction in the lethe-fluid-vans and lethe-fluid-particles solvers. This led to assertions being thrown in debug. This has been corrected by reiniting the constraints with the appropriate size. For the heat transfer and the block Navier-stokes, the issue was that the constraints were never reinited with the correct size, so they contained ghosted elements. This was caught by a new assert introduced in 2023-09 within deal.II master.

## [Master] - 2023-09-20

### Deprecated

- MINOR The calculation of the source term was enabled using a parameter called "enable". This parameter was used in some physics, not in others and was poorly implemented. We deprecate the usage of this parameter and always enable source term, considering the fact that the default value of the source term is zero anyway. This prevents the false perception that source terms could be enabled or disabled, while the behavior was inconsistent across physics.


## [Master] - 2023-09-19

### Changed

- MAJOR All the applications were renamed [#882](https://github.com/chaos-polymtl/lethe/pull/882): `gls_navier_stokes` is now `lethe-fluid`, `gd_navier_stokes` is now `lethe-fluid-block`, `nitsche_navier_stokes` is now `lethe-fluid-nitsche`, `gls_sharp_navier_stokes` is now `lethe-fluid-sharp`, `gls_vans` is now `lethe-fluid-vans`, `mf_navier_stokes` is now `lethe-fluid-matrix-free`, `dem` is now `lethe-particles`, `cfd_dem_coupling` is now `lethe-fluid-particles`, `rpt3d` is now `lethe-rpt-3d`, `rpt_cell_reconstruction_3d` is now `lethe-rpt-cell-reconstruction-3d`, `rpt_fem_reconstruction_3d` is now `lethe-rpt-fem-reconstruction-3d`, and `rpt_l2_projection_3d` is now `lethe-rpt-l2-projection-3d`.


## [Master] - 2023-10-30

- MINOR The rotational vector for the rotational boundary condition in the lethe-particles solver is now define with one line in the parameters file. [#920](https://github.com/chaos-polymtl/lethe/pull/920)



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
