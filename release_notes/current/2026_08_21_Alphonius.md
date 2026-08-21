## [Master] - 2026/08/21

### Added

- MINOR The `gaussian_heat_flux_cls_interface` and `uniform_heat_flux_cls_interface` laser heat sources have an angle of incidence factor. The factor can now be disabled with the new parameter `enable angle of incidence dependence`. By default, its value is set to `true` to not alter current application tests. A regression application test, `heat_transfer_cls_lpbf_benchmark_rotated_laser_desabled_angle_of_incidence_factor` has been added to test the parameter. [#2098](https://github.com/chaos-polymtl/lethe/pull/2098)