## [Master] - 2026/10/14

### Added

- MAJOR The DEM `file` insertion method now reads the particle output of a previous simulation directly when a file of the `list of input files` ends with `.pvtu` or `.vtu`, without extracting the particles in a text file first. The `.vtu` files referred to by a `.pvtu` file are split among the processes, so each process only reads its share of the particles. These compressed files are much smaller than the equivalent text files, which could reach several GB and were read by every process. [#PR](link)
