## [Master] - 2026/09/16

### Fixed

- MAJOR The periodic particle-particle contact force, torque and heat
  transfer were applied twice to the same particle in simulations that
  combine two or more periodic directions with two or more MPI processes.
  `find_cell_periodic_neighbors` maps periodic cell pairs starting from the
  principal side of each periodic boundary, and feeds three containers
  accordingly: local-local, local-ghost and ghost-local. With a single
  periodic direction, a cell's periodic partner always lies on the matching
  (non-principal) boundary and is never itself mapped as a main cell, so
  each periodic pair is discovered from one side only. With two or more
  periodic directions, a corner or edge cell can touch two principal
  periodic boundaries at once, so both cells of a periodic pair can be
  mapped as main cells, and each one discovers the other when it sweeps its
  periodic vertices. The same physical contact was then stored both as a
  local-ghost and as a ghost-local contact, and since both of these apply
  their contribution to the local particle of the pair, that particle
  received the contact force, the contact torque and, in thermal
  simulations, the heat transfer twice per timestep. The check meant to
  prevent this could not: the local-owned branch performed no cross-branch
  check at all when its periodic neighbor was a ghost cell, and the ghost
  branch tested its locally owned neighbor against a set that, by
  construction, only ever held ghost cells and could therefore never
  contain it. Of the two cells of such a pair, only the one with the
  smallest `dealii::CellId` now records it. `CellId` identifies a cell
  identically on every process, ghost cells included, so the two processes
  sharing a cross-process periodic pair elect the same cell without
  exchanging anything, and each process applies its contribution exactly
  once to its own particle. Simulations with a single periodic direction,
  and simulations on a single process, are unaffected and their results are
  unchanged. [#PR](link)
