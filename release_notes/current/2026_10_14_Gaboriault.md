## [Master] - 2026/10/14

### Fixed

- MINOR The CFD-DEM checkpoint files of the particle handler and of the insertion object were written by every rank into the same file. They are now written by rank 0 only, as in the DEM solver. This could be the cause of the problem we were experiencing with the file system on the clusters.  [#2163](https://github.com/chaos-polymtl/lethe/pull/2163)
