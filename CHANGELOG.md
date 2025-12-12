
# Change Log
All notable changes to the Lethe project will be documented in this file.
The changelog for the previous releases of Lethe are located in the release_notes folder.
The format is based on [Keep a Changelog](http://keepachangelog.com/).

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
