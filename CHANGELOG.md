# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

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

