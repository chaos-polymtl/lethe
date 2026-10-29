## [Master] - 2026/10/29

### Changed

- MINOR The particle-point and particle-line contacts of the DEM use contiguous containers. The contact pairs are stored in `ankerl::unordered_dense` maps instead of `std::unordered_map`, and the contact candidates, the boundary cells with a point and the boundary cells with a line are stored in vectors instead of maps keyed by a counter or by the cell id string. The simulation results are unchanged. [#XXXX](https://github.com/chaos-polymtl/lethe/pull/XXXX)
