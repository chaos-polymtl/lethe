# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

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
