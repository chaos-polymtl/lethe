## [Master] - 2026/10/18

### Fixed

- MINOR The particle-triangle contact detection of solid surfaces is faster. The check meant to skip the projection of a particle on a triangle when the particle is far from the plane of the triangle was always false, so every particle-triangle contact candidate was projected on its triangle at every time step. The check now compares the distance between the particle center and the plane of the triangle to the largest distance at which a contact force is calculated. The contact candidates of each triangle are also stored in a contiguous container, and the broad search no longer copies the neighbor list of each cell. The simulation results are unchanged. [#XXXX](https://github.com/chaos-polymtl/lethe/pull/XXXX)
