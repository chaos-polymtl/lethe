## [Master] - 2026/08/18

### Added

- MINOR Add the computations of the geometric barycenter position and velocity for CLS simulations. The functions are in the InterfaceTool namespace, and can be re-used for other simulation types (e.g., compute the barycenter of a region enclosed by a given level of the temperature field). Additionally, the number of significant digits in tables for the mass conservation and barycenter in the file cls.cc are set to the `log_precision` of the simulation control subsection introduced in PR [#1947](https://github.com/chaos-polymtl/lethe/pull/1947). [#2095](https://github.com/chaos-polymtl/lethe/pull/2095)