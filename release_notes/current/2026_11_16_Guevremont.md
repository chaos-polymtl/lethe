## [Master] - 2026/11/16

### Added

- MAJOR This PR adds a Dirichlet boundary condition for the tracer on the surface of the particles of the sharp-interface immersed boundary method (`lethe-fluid-sharp`), for example to model a surface reaction that is limited by the diffusion of the tracer. It is set in the new subsection `tracer boundary condition` of the subsection `particles` and applies to every particle. As for the velocity, the equation of the degrees of freedom shared by a cut cell and a fluid cell is replaced by a polynomial stencil that extrapolates the solution of a neighbouring fluid cell to the surface of the particle. The tracer is no longer solved inside the particles: there, it is a harmonic extension of the solution in the fluid. The stencils and the extension are implemented in the new class `IBScalarConstraints`, which applies to any scalar field discretized with continuous elements. [#PR](https://github.com/chaos-polymtl/lethe/pull/PR)
