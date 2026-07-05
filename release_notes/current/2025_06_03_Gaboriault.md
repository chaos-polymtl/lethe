## [Master] - 2026/05/29

### Added

- MAJOR The ``insertion acceptance function`` feature has been added. When using the ``volume`` insertion, this feature allows to insert particle from an insertion box of an arbitrary shape by defining a function. Insertion points from the original insertion box are kept if the function returns a value greater then ``0.`` relative to its coordinates before applying the insertion offset. [#2045](https://github.com/chaos-polymtl/lethe/pull/2045)
