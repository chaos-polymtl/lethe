## [Master] - 2026/11/24

### Changed

- MINOR The unit test serial_solid_intersection_01 maps a coarser solid surface in a coarser background mesh, near the top of the ball, so that its expected output lists 102 cell pairs instead of 24,275. It still checks which background cells each triangle of the solid surface is mapped to. [#PR](link)
