## [Master] - 2026/11/20

### Changed

- MINOR The broad search of the particle-solid surface contacts is faster, and its result is unchanged. It visited the neighbor cells of every background cell to which a solid cell is mapped, so a cell that neighbors several of these background cells was visited several times for the same solid cell. Each neighbor cell is now visited once per solid cell. [#PR](link)
