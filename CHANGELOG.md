# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2026/06/01

### Fixed

- MINOR This PR refactor the verbosity of the THM linear solver and add a memory consumption diagnostic for the case where the verbosity of the linera solver is set to extra verbose. Additionally, when running the diagnostics, it was noticed that the memory consumption of the dynamic sparsity pattern (dsp) was large. From discussion with @blaisb we established that recasting the dsp to a sparsity pattern will help the memory consumption. However, the `TrilinosWrappers::SparsityPattern` do not have the `memory_consumption` function and it is therefore not clear if the cost is really reduced as the proxy used to estimate it (`sparsity_pattern.n_nonzero_elements() * sizeof(TrilinosWrappers::types::int_type) * bytes_to_gb;`) actually show an increase in memory usage in comparison to the value obtained from `dsp.memory_consumption()`. Finally, the sparsity pattern and the sparse matrix are now only set up when the function `should_solve_auxiliary_physics()` returns true. Again this is to free memory when the THM solver is not required to be called. [#1998](https://github.com/chaos-polymtl/lethe/pull/1998)

## [Master] - 2026/06/01

### Added

- MINOR A new feature to the `grid_uniform_channel_with_meshed_cylinder` and `grid_uniform_channel_with_meshed_square_prism` is added so it is now possible to mesh or not the obstacle inside the channel. The default behavior of the meshes is also change so the obstacle is not meshed.[#1997](https://github.com/chaos-polymtl/lethe/pull/1997)

## [Master] - 2026/05/29

### Added

- MINOR The Method of Manufactured Solution example now features results with P3-P2 and P3-P2 elements. The right order of convergence is recovered. [#2003](https://github.com/chaos-polymtl/lethe/pull/2003)

## [Master] - 2026/05/27

### Changed

- MINOR The project version in `CMakeLists.txt` has been updated with the latest released version. [#2001](https://github.com/chaos-polymtl/lethe/pull/2001)

### Deprecated

- MAJOR The FEM subsection parameters `velocity order`, `pressure order`, `temperature order`, `tracer order`, `cls order` (and its alias `VOF order`), `void fraction order`, `phase cahn hilliard order`, `potential cahn hilliard order`, `electromagnetics trial order`, and `electromagnetics test order` have been renamed to use `degree` instead of `order` (e.g., `velocity degree`, `pressure degree`). The parameter name `degree` is more adequate, since it refers to the degree of the underlying shape functions that are used. The old names are kept as deprecated aliases and will continue to work but will trigger a deprecation warning. All examples, application tests, and documentation have been updated to use the new names. [#1994](https://github.com/chaos-polymtl/lethe/pull/1994)

## [Master] - 2026/05/18

### Fixed

- MINOR Fixed an intermittent (hard to reproduce) failure of the `particle_projector_03` MPI unit test that occurred with the Ninja build system. The `exchange_ghost_particles` call was moved out of `generate_particle_grid` so it no longer overlaps with the destructor of the local `particle_triangulation`. A missing `MPI_Barrier` was also added to the per-rank output loop to match the pattern established in `particle_projector_02`. [#1992](https://github.com/chaos-polymtl/lethe/pull/1992)

## [Master] - 2026/05/14

### Added

- MAJOR The DEM and CFD-DEM solver could only manage one set of periodic boundary condition. This PR extends the periodicity to be able to manage multiple periodic boundary conditions such that up to three pair of periodic boundary faces can be declared in 3D. [#1983](https://github.com/chaos-polymtl/lethe/pull/1983)

## [Master] - 2026/05/12

### Fixed

- MAJOR Fixed the usage of matrix free solver with solid material id in the domain for both the GCMG and LSMG multigrid preconditioners. This has been tested with the `grid_uniform_channel_with_meshed_cylinder` and required the mesh to have an additional option to disable the transfinite manifold as it does not work for multi-grid preconditioners with the matrix-free solvers. Finally, a test for the mesh `grid_uniform_channel_with_meshed_cylinder` has also been created to keep track of its changes in the future. [#1984](https://github.com/chaos-polymtl/lethe/pull/1984)

### Added

- MAJOR This PR enables the LSMG (local-smoothing multigrid) preconditioner for the VANS matrix-free solver (`lethe-fluid-vans-matrix-free`). When `set preconditioner = lsmg` is selected for fluid dynamics, the auxiliary particle-projected fields (void fraction, particle-fluid force, drag, particle velocity, momentum transfer coefficient) now have multigrid level DoFs distributed via the `ParticleProjector`, and one `MGTransferMatrixFree` per field is built to interpolate the fine vectors to every triangulation level. The VANS operator is extended to handle level cell iterators on coarse multigrid levels (gathering values via the mg dof indices instead of the active-cell `get_function_values` path). Only pure h-coarsening is supported on the LSMG path; mortar coupling and the FE_Q_iso_Q1 coarse-grid option are blocked with explicit assertions. Simplex meshes are now also explicitly rejected by the VANS matrix-free solver at construction to ensure that the code is a bit more defensive. [#1991](https://github.com/chaos-polymtl/lethe/pull/1991)

### Fixed

- MINOR The GCMG (global-coarsening multigrid) preconditioner of the VANS matrix-free solver was building the particle-fluid force transfer twice and never building the particle-fluid drag transfer (copy-paste bug). This left the drag contribution missing on coarse multigrid levels. The correct transfer is now built for each auxiliary field. [#1991](https://github.com/chaos-polymtl/lethe/pull/1991)

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
