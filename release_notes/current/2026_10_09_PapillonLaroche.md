## [Master] - 2026/10/09

### Fixed

- MINOR This PR fixes a bug in the timer output of the auxiliary physics. Previously, the timer outputs were called in each physic postprocess() routine, which is only called according to the calculation frequency parameter in the postprocess subsection. Hence, if different than 1, the timer tables were not printed. This PR removes the timer outputs from axiliary physic postprocess() calls and moves them to the new function output_per_iteration_timer() in the physics. The new functions are called in navier_stokes_base finish_time_step() via the multiphysics_interface new output_per_iteration_timer() function. [#2154](https://github.com/chaos-polymtl/lethe/pull/2154)
