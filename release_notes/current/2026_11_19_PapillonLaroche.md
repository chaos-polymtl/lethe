## [Master] - 2026/10/09

### Changed

- MINOR This PR improve the signed distance unit test (interface_tools_signed_distance.cc) by adding a mesh convergence analysis, distributed triangulation, and adaptive mesh refinement. The 3D test does not reach the asymptotic range to keep the simulation time low, but with two additional refinements the convergence orders are correct, and the 2D test is correct. [#2225](https://github.com/chaos-polymtl/lethe/pull/2225)
