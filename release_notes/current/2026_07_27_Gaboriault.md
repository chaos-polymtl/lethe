## [Master] - 2026/07/27

### Fixed

- MINOR The ``Minimum cell size to maximum particle diameter ratio`` was not properly displayed at the beginning of a restart DEM simulation. The ``report_cell_size_to_particle_diameter_ratio`` was called before setting the ``triangulation`` for a restart in ``read_checkpoint``, thus the minimum cell size used in the report function was a default value instead of the real one. [#2070](https://github.com/chaos-polymtl/lethe/pull/2070)
