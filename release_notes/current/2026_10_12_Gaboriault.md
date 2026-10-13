## [Master] - 2026/10/12

### Fixed

- MAJOR The packed insertion method now evaluates all contacts with the positions of the particles at the beginning of each pseudo time-step and then displaces each particle by the sum of the displacements of its contacts. Previously, particles were moved while the contacts were evaluated, which made the displacements depend on the evaluation order and on the partitioning of the domain, and made the convergence slower as the number of MPI processes increased. Particle-point and particle-line contacts are ignored by the packed insertion method. [#PR](link)
