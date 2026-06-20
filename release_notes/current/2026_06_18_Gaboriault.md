## [Master] - 2026/06/17

### Removed

- MINOR Adds `-Wno-unknown-warning-option` to the CMake arguments in `run_clang_tidy.sh` to suppress spurious warnings caused by GCC-specific flags in deal.II that Clang does not recognize. [#2036](https://github.com/chaos-polymtl/lethe/pull/2036)