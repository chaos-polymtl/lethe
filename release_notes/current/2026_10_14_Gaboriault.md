## [Master] - 2026/10/14

### Fixed

- MINOR The CFD-DEM checkpoint files of the particle handler and of the insertion object were written by every rank into the same file. They are now written by rank 0 only, as in the DEM solver. This should help  [#XXXX](https://github.com/chaos-polymtl/lethe/pull/XXXX)
