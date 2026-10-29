## [Master] - 2026/10/29

### Changed

- MINOR The DEM contact search no longer iterates tree containers or copies them. The boundary cells and faces used by the particle-wall broad search are stored in a vector sorted like the map they came from, the floating-wall broad search reads the cells of each floating wall by reference instead of copying them, and the displacement of particles across periodic boundaries iterates a vector of the pairs of periodic cells. The ongoing particle-wall collisions of the collision statistics are stored in an `ankerl::unordered_dense` map instead of a `std::map`. The simulation results are unchanged. [#XXXX](https://github.com/chaos-polymtl/lethe/pull/XXXX)
