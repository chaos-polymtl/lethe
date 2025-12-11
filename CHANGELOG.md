
# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2025/12/11

### Changed

- MAJOR Floating walls declaration in a parameter file got changed so that it requires fewer parameters for the ``point on wall`` and the ``normal vector``. A bug was also find during this change: the ``normal vector`` wasn't normalized after being parsed. As a results, if the user did not declare a unit vector in the parameter file, the normal overlap between particles and this ``floating wall`` is off by some factor. This problem got fixed, but resulted in the change of an application-test output (which was declaring the normal vector as a non-unit vector). [#1850](https://github.com/chaos-polymtl/lethe/pull/1850)

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
