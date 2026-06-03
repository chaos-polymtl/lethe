## [Master] - 2026/06/03

### Added

- MAJOR Adds the new parameter `calculation frequency` to the `post-processing` subsection. It indicates the frequency at which the post-processed quantities are computed. The `cls-velocity-extrapolation` application-test was updated to include this feature. This new parameter requires the `output frequency` to be a multiple of the `calculation frequency`. [#2010](https://github.com/chaos-polymtl/lethe/pull/2010)