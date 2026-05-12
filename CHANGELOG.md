# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## [Master] - 2026/05/11

### Fixed

- MAJOR Fixed the usage of matrix free solver with solid material id in the domain for both the GCMG and LSMG multigrid preconditioners. This has been tested with the `grid_uniform_channel_with_meshed_cylinder` and required the mesh to have an additional option to disable the transfinite manifold has it does not work for multi-grid which have also been added to this PR. Finally, a test for the mesh `grid_uniform_channel_with_meshed_cylinder` has also been created to keep track of its changes in the future. [#1984](https://github.com/chaos-polymtl/lethe/pull/1984)

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
