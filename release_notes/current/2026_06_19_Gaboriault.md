## [Master] - 2026/06/19

### Fixed

- MINOR The ``particle_sizes`` vector was kept as an attribute of the distribution class and wasn't resized after a particle insertion. This would result in a big ram usage. Now, this vector was removed, since the class does not need to keep track of the diameter values generated and the ``assign_properties`` function now passes a vector by reference. [#2041](https://github.com/chaos-polymtl/lethe/pull/2041)