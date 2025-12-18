
# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

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

