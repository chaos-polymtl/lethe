## [Master] - 2026/11/18

### Added

- MAJOR The constants of the shift contact models of the packed insertion method are now parameters of the `insertion info` subsection: `packed particle-particle overlap factor` (default 0.501), `packed particle-wall overlap factor` (default 1.02, used for the walls, the boundary points and lines, and the solid surfaces) and `packed clearance ratio` (default 0.01). The defaults are the previous hard-coded values, so the results are unchanged. A new application test, `insertion_packed_solid_surface`, checks the packed insertion of particles around a solid surface. [#PR](link)

### Fixed

- MAJOR The solid objects no longer move during the packed insertion. The pseudo time-steps of the packed insertion moved the solid objects that have a translational or angular velocity, so the particles were arranged around a solid object that was no longer at its initial position, and the solid object could push particles out of the domain. [#PR](link)
