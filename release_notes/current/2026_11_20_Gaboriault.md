## [Master] - 2026/11/20

### Changed

- MINOR The mapping of the solid surfaces in the background triangulation is faster, and its result is unchanged. It computed the distance between the center of every locally owned background cell and every cell of the solid surface. An R-tree of the locally owned background cells now selects, for each solid cell, the background cells whose center can be closer to it than their diameter, and the distance is only computed for them. [#PR](link)
