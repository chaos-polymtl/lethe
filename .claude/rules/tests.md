---
paths:
  - "tests/**/*"
  - "applications_tests/**/*"
---

# Tests in Lethe

A new feature needs unit and/or application tests. A bug fix needs a test that
fails without the fix. Prefer unit tests when the feature can be tested in
isolation.

## Unit tests (`tests/<module>/`)

- Start from `tests/TEMPLATE`: include `<../tests/tests.h>`, call `initlog()` in
  `main`, and write the results to `deallog`. `<name>.output` holds the
  expected `deallog` content. A test that needs MPI creates a
  `Utilities::MPI::MPI_InitFinalize` and uses `<name>.mpirun=N.output`.
- `TEST_LIBRARIES` in `tests/<module>/CMakeLists.txt` follows the library
  graph (`tests/dem` links `lethe-core lethe-dem`). A unit test cannot use a
  library that its module does not depend on.
- Data files read by a test are copied in the folder's `CMakeLists.txt` with
  `file(COPY ... DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/<test>.${_build_type}")`
  (see `tests/core/CMakeLists.txt`).

## Application tests (`applications_tests/<app>/`)

- `<name>.prm` + `<name>.output`, the expected standard output of the
  application. `.prm` files start with the SPDX lines (`#` comments) and are
  formatted with `prmindent`.
- `<name>.mpirun=N.output` runs on N ranks. Stay at N ≤ 3 (existing tests use
  1 to 3).
- `<name>.with_qmomkit=on.output` runs only when Lethe is built with QMOMKit.
- Restart tests: the `.prm` that writes the checkpoint goes in
  `applications_tests/<app>/generators/`. The checkpoint files are committed in
  `<test>_restart_files/` and copied to
  `${CMAKE_CURRENT_BINARY_DIR}/<test>.${_build_type}/mpirun=N/` by the folder's
  `CMakeLists.txt`. Regenerate them when the checkpoint format changes.

## For every test

- Deterministic across MPI ranks and platforms: no wall-clock timings, no
  unordered output, no random seed without a fixed value.
- Fast: a small mesh and few time steps. A test that needs a longer timeout in
  Debug gets `set_tests_properties(<app>/<test>.mpirun=N.debug PROPERTIES
  TIMEOUT <s>)` (no `.mpirun=N` for a serial test) in the folder's
  `CMakeLists.txt`. Justify any new or longer timeout in the PR: it can hide a
  performance regression.
- Golden files: `contrib/utilities/update-golden.tl <build-dir>` (needs TXR; `-n`
  for a dry run) copies the build outputs over the `.output` files. Use it only
  after checking that the new output is correct, and explain in the PR why it
  changed. A change in user-facing output can touch many `.output` files.
