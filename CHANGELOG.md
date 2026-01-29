
# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2026/01/27

### Changed

- MINOR The analytical Jacobian computation in the SignedDistanceSolver involves symmetric matrices and the implementation was not taking advantage of that. Additionally, intermediary matrices where reallocated each time the function was called. This PR improves the performance of compute_analytical_jacobian() by using the symmetry property in matrix-matrix products. It also includes the intermediary matrices as attributes of the SignedDistanceSolver class to avoid multiple re-instanciation. The gain in terms of computational time is about 20% in 3D. [#1893](https://github.com/chaos-polymtl/lethe/pull/1893)

### Fixed

- MINOR Added a missing term in the volume-filtered momentum conservation equation stress tensor. The term $\varepsilon_{\mathrm{f}} \nu \partial_j \partial_i u_j$ was omitted in the previous implementation. The outputs of the `lethe-fluid-particles` and `lethe-fluid-vans` application tests changed only slightly, except for MMS test case 2, which was updated so that the source term accounts for this correction. The outputs of the affected tests were updated accordingly. [#1890](https://github.com/chaos-polymtl/lethe/pull/1890)

## [Master] - 2026/01/26

### Added

- MINOR This PR adds the option to set a simulation initial time-step that respects the capillary time-step constraint through enabling (setting to `true`) the parameter `respect capillary time-step constraint` ($\Delta t_\sigma$) in the `simulation control` subsection. This constraint is especially useful when simulating capillary flows. With the parameter `capillary time-step ratio` ($N_\text{CTR}$), users can specify a multiple of the default capillary time-step constraint as defined by [1]. By default, the $N_\text{CTR}$ is set to 1. A $N_\text{CTR}$ value < 1 can improve stability for highly dynamic cases, and using a $N_\text{CTR}$ value > 1 can lead to a significant decrease in computation time in some cases. 

  Additionally, when `respect capillary time-step constraint` is  paired with `adapt` (for adaptive time-stepping), the time-step will be dynamically updated to respect the capillary time-step ($N_\text{CTR} \times \Delta t_\sigma $) at every time-step. A warning message is printed on the console when `adapt` is disabled to inform the user that the capillary time-step constraint might not be respected if either of the mesh cell size, the densities or the surface tension coefficient evolves dynamically.

  Finally, a new overload of the `announce_string` function has also been added to allow multiple line banners. [#1886](https://github.com/chaos-polymtl/lethe/pull/1886)

  [1] F. Denner, F. Evrard, and B. van Wachem “Breaching the capillary
  time-step constraint using a coupled VOF method with implicit surface tension,” J. Comput. Phys., vol. 459, p. 111128, Jun. 2022, doi: 10.1016/j.jcp.2022.111128

- MINOR This PR is a follow-up for the time harmonic electromagnetic auxiliary physics rudimentary implementation #1852. It links its physical material properties #1862 to the assembly, creates the waveguide port boundary conditions, update previous application tests related to the physics and add a new one for a dissipative medium. It also removes .prm in the examples/multiphysics/waveguide as this example will come in a follow up PR. [#1882](https://github.com/chaos-polymtl/lethe/pull/1882) 

## [Master] - 2026-01-23

### Changed

- MAJOR Some initial tests created for the mortar feature became obsolete as the code developed. Hence, this PR cleans up the mortar tests and moves the more representative ones (that are not used anymore) to the prototypes folder. [#1887](https://github.com/chaos-polymtl/lethe/pull/1887)

## [Master] - 2026-01-22

### Fixed

- MINOR The hopper example post-processing script was not working anymore since the change in the definition of the floating wall. The script would try to parse "nx" from the .prm file, but this variable did not exist anymore. This PR fixes this issue by removing this parameter parsing. Furthermore, we have identified an issue with the automatic validation script which can keep on reusing old pdf of the results if the generation of the new ones did not work. This is also fixed by automatically removing all PDF that were already present in the case folder. [#1889](https://github.com/chaos-polymtl/lethe/pull/1889)

## [Master] - 2026-01-21

### Changed

- MINOR The cylindrical gas–solid fluidized bed example, introduced in PR [#1844], has been added to the Lucille validation test suite. The selected validation case uses the matrix-free solver with a semi-implicit drag scheme and a QCM filter. The corresponding post-processing Python script has been modified to allow users to select the simulation configuration to be post-processed. [#1880](https://github.com/chaos-polymtl/lethe/pull/1880)

### Fixed

- MINOR Calls to convert_vector_dealii_to_trilinos() were failling due to new restriction on ghosted Trilinos vectors in dealii. This PR fixes this using a temporary locally owned vector. [#1888](https://github.com/chaos-polymtl/lethe/pull/1888)

## [Master] - 2026-01-20

### Changed

- MINOR This PR improves the performances of the signed distance solver. Previously, the bottleneck was the calls to FEPointEvaluation.reinit(cell, points_vector) in compute_second_neighbors(). Hence, the cost of these calls has been reduce in two ways: 1)  the numerical jacobian of the minimization problem was replace by the analytic jacaobian to reduce the number of points to reinit for and 2) the MappingFE(Fe_Q) was replaced by a MappingQ(fe.degree) to take the fast path in FEPointEvaluation.reinit() (those two mapping are equivalent, but MappingFe was not recognized for the fast path). For the 3D simulation of the rising bubble benchmark, we have a speed up of approx. 7x for the signed distance computations. Hence, the geometric reinitialization is not the bottleneck of the simulation anymore. [#1885](https://github.com/chaos-polymtl/lethe/pull/1885)

## [Master] - 2026-01-19

### Changed

- MINOR The latest Lethe release is compatible with deal.ii 9.7, and we don't have CI checks with 9.6 version (or below) anymore. Hence, this PR removes deprecated checks for deal.II versions prior to 9.7. [#1884](https://github.com/chaos-polymtl/lethe/pull/1884)

## [Master] - 2026-01-16

### Changed

- MINOR In the cylindrical packed bed example, the pressure drop across the bed is compared with the Ergun correlation. This was done in PR [#1840], where a bed porosity of 0.45 was used in the correlation. In this PR, the bed porosity is calculated using the bed height that is extracted from the simulation instead, leading to a better estimate for the pressure drop. [#1881] (https://github.com/chaos-polymtl/lethe/pull/1881)

- MINOR This PR adds the option to declare a particle size distribution in the parameter file as volume-based for the ``normal`` and ``lognormal`` distribution. The ``custom`` distribution is the only distribution that is not supporting the ``number-based`` distribution. The will get fixed in a follow up PR. [#1875](https://github.com/chaos-polymtl/lethe/pull/1875)

## [Master] - 2026/01/15

### Changed

- MAJOR The algebraic reinitialization PDE had some issues converging towards steady-state. To solve them, this PR: 
  - changes the weak formulation of the problem solved (the curvature field is no more required);
  - changes how the artificial time-step is computed;
  - changes how the steady-state criterion is evaluated;
  - removes the filtered VOF solution field for the computation of the projected phase gradient used to compute the unit interface normal vector and the curvature;
  - renames the `reinitialization CFL` parameter to `artificial time-step factor` (an alias was made for `reinitialization CFL`, to avoid parameter file version issues), and;
  - changes default values of algebraic reinitialization parameters.

  This PR also updates the rising bubble example. Other VOF examples with surface tension force (static bubble, capillary wave, and Rayleigh-Plateau instability) will have to be checked and updated in a future PR. [#1879](https://github.com/chaos-polymtl/lethe/pull/1879)

## [Master] - 2026-01-13

### Changed

- MINOR This PR replaces the duplicate ``set_insertion_type`` function definitions in ``dem.cc``, ``cfd_dem_coupling.cc``, and ``cfd_dem_coupling_matrix_free.cc`` with a single definition in ``set_insertion_method.cc``. [#1877](https://github.com/chaos-polymtl/lethe/pull/1877)

## [Master] - 2026/01/12

## Fixed

- Minor The lethe-particles/distribution_normal.prm application test was using a log-normal distribution. The parameter was changed back to ``normal``. As a result, the output of the test did change. [#1876](https://github.com/chaos-polymtl/lethe/pull/1876) 

- Minor The convection term in the VANS equations for models A and B is changed from ``local_matrix_ij += ((phi_u_j * void_fraction * velocity_gradient * phi_u_i) + (grad_phi_u_j * void_fraction * velocity * phi_u_i));`` to ``local_matrix_ij += ((velocity_gradient * phi_u_j * void_fraction * phi_u_i) + (grad_phi_u_j * void_fraction * velocity * phi_u_i))`` in the matrix assembly, to match the correct algebraic formula. The output of the ``lethe-fluid-particles`` and ``lethe-fluid-vans`` application tests changed only slightly. Therefore, their outputs were updated in this PR as well. [#1874](https://github.com/chaos-polymtl/lethe/pull/1874) 


## [Master] - 2026/01/07

## Added

- MINOR The particle projector can be used to project the particle-fluid force onto the CFD mesh. This requires the assembly of a matrix and a right-hand side. On smaller meshes, this assembly is quite cheap, but when using the matrix-free methods, it becomes a significant cost for large parallel simulations. This PR optimizes this projection step by ensuring that the matrix and the preconditioner are only recalculated and reinitialized when it is necessary. In the case when the projected field does not have a Neumann boundary condition, then the matrix is only assembled once, whenever the degrees of freedom are set up, instead of at every iteration. This greatly diminishes the cost of the explicit, semi-implicit and implicit coupling. There are still optimizations remaining (for example, to carry out all the projections in a single step instead of one by one), but these will be addressed in a follow-up PR. This already reduces the cost by 20-25% of the projection. Furthermore, the projection was not adequately timed, this has been fixed in this PR. [#1873](https://github.com/chaos-polymtl/lethe/pull/1873)

## Fixed

- MINOR  The use of an automatic pointer was creating a compilation error on modern C++ compiler. The communicator is explicitely declared instead of using an automatic pointer. This does not affect anything, but makes compilation more stable accross platforms. [#1872](https://github.com/chaos-polymtl/lethe/pull/1872)

## [Master] - 2026/01/06

## Fixed

- MINOR The test read_mortar_data_02 constantly fails in the CI release, even if it passes on debug mode and also in local office computers. Since read_mortar_data_01 and read_mortar_data_03 are very similar to such test (only with different grids), this PR solves this issue by removing the (unnecessary) read_mortar_data_02 test. [#1871](https://github.com/chaos-polymtl/lethe/pull/1871)

## [Master] - 2025/12/19

### Added

- MAJOR This PR adds a rudimentary version of the new time harmonic electromagnetic auxiliary physics solver. It puts in place the time_harmonic_maxwell.h and .cc files with the minimum required functions to solve the DPG linear system with two .prm files as input. The electromagnetic multiphysics, its most basic boundary conditions, the initial condition, the analytical solution, and the post processing of the interior fields function skeletons are set in place. The other functionality implementation will come in future PRs. The solver is tested with two new application tests based on the .prm. [#1852](https://github.com/chaos-polymtl/lethe/pull/1852)

## [Master] - 2025/12/23

### Added

- MINOR This PR adds a gas-solid cylindrical fluidized bed example to the Lethe example suite. Different filtering and drag coupling approaches are compared using a robust benchmark problem: the dependence of the pressure drop on the superficial gas inlet velocity. VANS Model A is used with both QCM and cell-based filters using the matrix-based solver. The same model is also tested with the QCM filter using the matrix-free solver. For each configuration, three drag coupling strategies are considered: explicit, semi-implicit, and implicit. [#1844](https://github.com/chaos-polymtl/lethe/pull/1844)

- MAJOR The unresolved CFD-DEM matrix-free implementation could not allow dynamic load balancing during the simulation since this feature had not been implemented yet. This adds the load balance implementation as well as make sure that it works well in both releaseand debug mode. An application test which is the sedimentation of a particle in parallel with dynamic load balancing is added to ensure that the test is stable in parallel since it was quite difficult to figure out all of the ghost logic of the vectors considering that both Trilinos and deal.II vectors are used in the ParticleProjector class. [#1869](https://github.com/chaos-polymtl/lethe/pull/1869)

## [Master] - 2025/12/19

### Fixed

- MINOR This is the follow-up to PR 1855, where we fixed the infinite while loop related to solid surface. This PR adds an application that was falling previously to this fix. This PR also changes how the main loop on every solid object is built during the force calculation which results in a significant speed-up of the code. [#1865](https://github.com/chaos-polymtl/lethe/pull/1865)

## [Master] - 2025/12/18**

### Fixed

- MINOR This PR fixes boundary list that were added in #1697. Boundary lists were not being parsed properly for physics other than fluid dynamic. The parameter parsing now expected a list of integers instead of an integer. [#1866](https://github.com/chaos-polymtl/lethe/pull/1866)

## [Master] - 2025/12/17

### Changed

- MAJOR Floating walls declaration in a parameter file got changed so that it requires fewer parameters for the ``point on wall`` and the ``normal vector``. A bug was also find during this change: the ``normal vector`` wasn't normalized after being parsed. As a results, if the user did not declare a unit vector in the parameter file, the normal overlap between particles and this ``floating wall`` is off by some factor. This problem got fixed, but resulted in the change of an application-test output (which was declaring the normal vector as a non-unit vector). [#1850](https://github.com/chaos-polymtl/lethe/pull/1850)

- MINOR The verification of the number of cells at the stator and rotor interfaces was being done inside the loops for dealii and gmsh grids. This PR adjusts this by doing the verification only once, after the initial mesh refinement is performed. [#1854](https://github.com/chaos-polymtl/lethe/pull/1854)

### Fixed

- MINOR This is a follow-up of #1863. Although #1863 fixed the restart and prevented the simulations from crashing, the time history of the void fraction was not store appropriately. This stemmed from a confusion since in the matrix-free solver, it is the deal.II distributed vectors for the void fraction which are checkpointed and not the Trilinos ones. This PR reads the correct deal.II vector when reading a checkpoint, but also ensures that the values in the Trilinos vectors matches that of the deal.II vectors. This allows reproducing the norm of the residuals to machine accuracy when restarting. [#1864](https://github.com/chaos-polymtl/lethe/pull/1864)

- MINOR The Matrix-free CFD-DEM solver could not restart adequately. This is because the solution vector for the fluid was not sized accordingly. The deal.II vectors require that the locally relevant dofs be provided to the vector before reading a checkpoint, which is not the case for the Trilinos vectors. This would prevent restarts in parallel. [#1863](https://github.com/chaos-polymtl/lethe/pull/1863)

## [Master] - 2025/12/16

### Added

- MINOR This PR add the material properties that will be used by the time-harmonic electromagnetic solver which include the electric conductivity, the complex electric permittivity and the complex magnetic permeability. At the moment, those properties only support a "constant" field. A unit test for the electromagnetic constant properties was added, along with the update of the already existing "physical_properties_manager" tests. [#1862](https://github.com/chaos-polymtl/lethe/pull/1862)

## [Master] - 2025/12/15

### Fixed

- MINOR The tolerance adopted in the radius computation at the rotor-stator mortar interface was a hard-coded value, which was not ideal. This PR fixes this by introducing a radius tolerance parameter. [#1853](https://github.com/chaos-polymtl/lethe/pull/1853)

## [Master] - 2025/12/14

### Changed

- MAJOR A new message is added at the start of every DEM and unresolved CDF-DEM simulation which informs the user about the kind of distribution being used for every particle type. When using a normal or lognormal distribution, an extra message is written about the diameter cutoff values relative to the entire distribution. This change is MAJOR since every application test using the `lethe-particles`, `lethe-fluid-particles` and `lethe-fluid-particles-matrix-free` solver had to be updated. Also in this PR, the documentation relative to the `update-golden.tl` script got updated. [#1849](https://github.com/chaos-polymtl/lethe/pull/1849)

### Added

- MINOR This PR adds the possibility to use the semi-implicit and implicit drag couplings with the CFD-DEM matrix-based solver while also using the QCM filter to project the particle forces onto the fluid. At this stage, it has been tested with the single particle sedimentation case. [#1845](https://github.com/chaos-polymtl/lethe/pull/1845)

## [Master] - 2025/12/13

### Fixed

- MAJOR The new implementation of the solid-particle contact detection could lead to an infinite loop for one of the cases that happens only on deformed mesh. The issue was that two continue statement could be called without having the iterator increased. This lead to the aforementioned bug. Running any example with a solid object (bunny drill or granular mixer) would directly lead to this infinite loop being encountered no matter the number of cores. The solution is simple (increment the iterator), but it also demonstrates that we need more robust testing for this part of the code. [#1855](https://github.com/chaos-polymtl/lethe/pull/1855)

## [Master] - 2025/12/11

### Fixed

- MINOR The documentation for the simulation control and the analytical solution was improved. [#1843](https://github.com/chaos-polymtl/lethe/pull/1843)

### Removed

- MAJOR Officially deprecate the RPT applications and remove them. The applications are now available at https://github.com/chaos-polymtl/lethe-rpt. This significantly decreases the size of the repository and the length of the CI procedure. [#1851](https://github.com/chaos-polymtl/lethe/pull/1851) 

## [Master] - 2025/12/10

### Added

- MINOR This PR adds the possibility to change the diffusion factor of the VOF DCDD stabilization in the prm file. [#1847](https://github.com/chaos-polymtl/lethe/pull/1847)

### Fixed

- MINOR The extract-slice-from-vtu.py tool was not working properly with --np > 1 because an argument (file path) was missing in the call of extract_slice. This PR fixes this bug and add 2 functionalities. 1) When the flag group files in the .prm (parallel output of Lethe) was set to a value larger than 1, the extract-slice-from-vtu.py would not slice all .vtus associated with a .pvtu. Since pyvista can read .pvtus correctly, the list of files to read is now based on .pvtus instead of .vtus. 2) There was no sliced .pvd file generated, so the slices could not be used easily in other post-processing tools using .pvds. The .pvd file of the slices is now generated. [#1848] (https://github.com/chaos-polymtl/lethe/pull/1848)

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

