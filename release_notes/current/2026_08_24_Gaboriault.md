## [Master] - 2026/08/24

### Fixed

- MINOR The `save()` function used to checkpoint the simulation control was
  using the default stream precision of 6 significant digits. This precision was
  insufficient for small time steps, which could result in a loss of precision
  when saving the simulation time and time steps. The precision is now set to 17
  significant digits.