## [Master] - 2026/10/18

### Fixed

- MINOR The bunny-drill example was not working since we check the minimum cell size to the maximum particle diameter ratio. The mesh was adjusted to pass this check. Since the mesh is not bigger, the number of particles inserted was increased and the simulation time was also increase so that all the particles have the time to be inserted. The documentation was changed accordingly. [#2171](https://github.com/chaos-polymtl/lethe/pull/2171) 
