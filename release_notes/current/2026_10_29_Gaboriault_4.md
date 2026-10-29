## [Master] - 2026/10/29

### Changed

- MINOR The searches for boundary cells and for cell neighbors of the DEM, which run at the start of a simulation and after every load balancing, no longer scale quadratically with the mesh size. The search for boundary faces no longer compares every face to all the faces already found, the search for the boundary neighbors of cells with a boundary line no longer scans all the boundary cells for every neighbor, and the cell neighbor searches flag cells in vectors indexed by the active cell index instead of storing them in `std::set`. The boundary lines are gathered per cell in a map keyed by the active cell index instead of the cell id string. The simulation results are unchanged. [#XXXX](https://github.com/chaos-polymtl/lethe/pull/XXXX)
