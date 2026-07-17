## [Master] - 2026/06/30

### Added

- MINOR This PR implements Neumann traction boundary condition assembler for both matrix-based and matrix-free solvers to solve for the incompressible Navier-Stokes equations. The feature has been tested with the `method of manufactured solutions - 2d problem` adapted to have Neumann boundary condition on the top edge of the domain. Finally, a test is added for both Matrix-based and Matrix-free Navier-Stokes solvers solving the adapted mms_2d_problem with the `mms_2d_fe2_neumann_traction_navierstokes.prm` to keep track of its changes in the future. [#1965](https://github.com/chaos-polymtl/lethe/pull/1965)