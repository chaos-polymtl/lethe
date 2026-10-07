## [Master] - 2026/10/07

### Fixed

- MAJOR The periodic particle-particle contact force, torque and heat transfer were applied twice to the same particle in simulations combining two or more periodic directions with two or more MPI processes. With several periodic directions, a corner or edge cell can touch two principal periodic boundaries at once, so both cells of a periodic pair were mapped as main cells in `find_cell_periodic_neighbors` and each one discovered the other. The same contact was then stored both as a local-ghost and as a ghost-local contact, and both apply their contribution to the local particle. Of the two cells of such a pair, only the one with the smallest `dealii::CellId` now records it; since a `CellId` is identical on every process, both processes elect the same cell without any communication. Simulations with a single periodic direction or on a single process are unaffected. [#PR](link)
