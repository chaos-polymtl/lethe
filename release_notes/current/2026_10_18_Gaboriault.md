## [Master] - 2026/10/18

### Fixed

- MINOR The bunny-drill example was not working since we check the minimum cell size to the maximum particle diameter ratio. The cells at the center of the cylinder was too small relative to the particles. The mesh was adjusted to pass this new check. Since the mesh is now bigger, the number of particles inserted was increased and the simulation time was also increased so that all the particles have the time to be inserted in the loading simulation. The example documentation was changed accordingly. [#2171](https://github.com/chaos-polymtl/lethe/pull/2171) 
