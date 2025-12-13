
# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2025/12/13

### Fixed

- MAJOR The new implementation of the solid-particle contact detection could lead to an infinite loop for one of the cases that happens only on deformed mesh. The issue was that two continue statement could be called without having the iterator increased. This lead to the aforementioned bug. Running any example with a solid object (bunny drill or granular mixer) would directly lead to this infinite loop being encountered no matter the number of cores. The solution is simple (increment the iterator), but it also demonstrates that we need more robust testing for this part of the code. [#1855](https://github.com/chaos-polymtl/lethe/pull/1855)

## [Master] - 2025/12/11

### Fixed

- MINOR The documentation for the simulation control and the analytical solution was improved. [#1843](https://github.com/chaos-polymtl/lethe/pull/1843)

## [Master] - 2025/12/10

### Added

- MINOR This PR adds the possibility to change the diffusion factor of the VOF DCDD stabilization in the prm file. [#1847](https://github.com/chaos-polymtl/lethe/pull/1847)

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

## [Master] - 2025/12/10

### Fixed

- MINOR The extract-slice-from-vtu.py tool was not working properly with --np > 1 because an argument (file path) was missing in the call of extract_slice. This PR fixes this bug and add 2 functionalities. 1) When the flag group files in the .prm (parallel output of Lethe) was set to a value larger than 1, the extract-slice-from-vtu.py would not slice all .vtus associated with a .pvtu. Since pyvista can read .pvtus correctly, the list of files to read is now based on .pvtus instead of .vtus. 2) There was no sliced .pvd file generated, so the slices could not be used easily in other post-processing tools using .pvds. The .pvd file of the slices is now generated. [#1848] (https://github.com/chaos-polymtl/lethe/pull/1848)
