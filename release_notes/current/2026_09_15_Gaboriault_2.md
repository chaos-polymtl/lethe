## [Master] - 2026/09/15

### Changed

- MINOR The periodic particle-particle fine search no longer scans a combinatorial list of periodic translations (up to 26 entries in 3D) to find the nearest periodic image of a particle. Since periodic directions are axis-aligned and independent, the nearest image is now computed directly per direction with the minimum image convention. This is an internal implementation change with no impact on simulation results. [#PR](link)
