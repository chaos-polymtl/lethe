## [Master] - 2026/11/23

### Changed

- MINOR The R-tree of the locally owned background cells used to map the solid surfaces onto the background mesh in the DEM solver is now built once and rebuilt only when the background mesh changes (load balancing or grid motion), instead of at every mapping and for every solid. The mapping is unchanged. [#2232](https://github.com/chaos-polymtl/lethe/pull/2232)
