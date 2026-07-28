## [Master] - 2026/07/06

### Changed

- MINOR This PR fixes issue [#2057](https://github.com/chaos-polymtl/lethe/pull/2057) which consists in the cleanup of the insertion_volume.cc file. Tests output using the volume insertion had to change since the second random number used for the insertion offset is no longer the same. This changes the insertion location of the particle, hence changes the test outputs. [#2058](https://github.com/chaos-polymtl/lethe/pull/2058)
