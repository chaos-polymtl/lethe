## [Master] - 2026/11/18

### Added

- MAJOR The constants of the shift contact models of the packed insertion method are now parameters of the `insertion info` subsection: `packed particle-particle overlap factor` (default 0.501), `packed particle-wall overlap factor` (default 1.02) and `packed clearance ratio` (default 0.01). The defaults are the previous hard-coded values. A new application test, `insertion_packed_solid_surface`, checks the packed insertion of particles between a solid surface and a wall. [#PR](link)

### Fixed

- MAJOR The solid objects no longer move during the packed insertion. The pseudo time-steps of the packed insertion moved the solid objects that have a translational or angular velocity, so the particles were arranged around a solid object that was no longer at its initial position, and the solid object could push particles out of the domain. [#PR](link)
- MAJOR The particles can no longer cross the solid surfaces and the walls during the packed insertion. The displacements of the particle-particle contacts of a particle could push it through a solid surface, where it remained on the wrong side, or through a wall, where it was deleted. The displacement of the particle-particle contacts are now limited to the particle radius at each pseudo time-step, before the displacements of the particle-wall and particle-solid contacts are added. This changes the final particle positions obtained with the packed insertion. [#PR](link)
