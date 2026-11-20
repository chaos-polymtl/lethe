## [Master] - 2026/11/18

### Added

- MAJOR This PR adds the transport of the moments of the quadrature method of moments (QMOM) auxiliary physics in physical space, with their source terms. All the QMOM components share a single scalar discontinuous Galerkin element and DoFHandler, with one vector per component. They are advected by the fluid velocity with an upwind discontinuous Galerkin discretization and a BDF time integration, with a single matrix and preconditioner (ILU or AMG) shared by all the components. The moments are constant in phase space, and the scalars, the Strang splitting and the phase-space mechanisms are rejected with an error until they are implemented. [#2216](https://github.com/chaos-polymtl/lethe/pull/2216)

### Changed

- MAJOR The `none` type of the `boundary conditions qmom` subsection is renamed `noflux`, which states that no particle crosses these boundaries (walls). It remains the default type. [#2216](https://github.com/chaos-polymtl/lethe/pull/2216)

### Removed

- MAJOR The `qmom` subsection of the `non-linear solver` subsection is removed, since the QMOM physics is solved with the linear solver strategy. [#2216](https://github.com/chaos-polymtl/lethe/pull/2216)
