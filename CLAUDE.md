# CLAUDE.md

Guidance for Claude Code (and other AI assistants) working in Lethe. More
detailed rules in `.claude/rules/` load when a matching file is read or edited:
`cpp.md` (naming, clang-tidy, `-Werror`), `tests.md`, `parameters-and-docs.md`.

Lethe is an open-source CFD / DEM / CFD-DEM solver built on **deal.II** (≥ 9.7,
with MPI, p4est, Trilinos and magic_enum). Assume an MPI-parallel, template-heavy,
high-performance C++ codebase. Docs: https://chaos-polymtl.github.io/lethe/

## Layout

- `include/`, `source/` — headers (`.h`) and sources (`.cc`), each split into
  `core/`, `dem/`, `solvers/`, `fem-dem/` (one library each).
- `applications/` — one executable per folder; minimal main files (< 150 lines)
  that instantiate the templated solvers for a given `dim`. Driven by `.prm`.
- `applications_tests/`, `tests/` — functional tests (per application) and unit
  tests (per module).
- `doc/` (Sphinx site), `examples/` (documented cases), `prototypes/`
  (throwaway, unstable), `contrib/utilities/` (maintenance scripts),
  `release_notes/`.

## Library dependencies — never violate

```
core ──► dem ──────┐
  └────► solvers ──┴──► fem-dem
```

- `core` includes no other Lethe library. `dem` and `solvers` include only
  `core` and **never each other**. `fem-dem` may include all three.
- Include Lethe headers as `<module/file.h>`, never `"file.h"`: `checkdeps`
  only parses the angle-bracket form. `qmomkit/` headers are allowed only in
  `solvers` and `fem-dem`, under `#ifdef LETHE_WITH_QMOMKIT`.
- `target_link_libraries` mirrors the graph (`lethe-dem`, `lethe-solvers` →
  `lethe-core`; `lethe-fem-dem` → all three). Never add a link to make a build pass.
- When a lower library needs something from a higher one, move the shared code
  down to the lowest library that needs it (often `core`), or invert the
  dependency with an interface or callback declared in the lower library. A
  forward declaration does not make an upward dependency acceptable.
- Run `contrib/utilities/checkdeps` after changing includes — **CI does not run
  it**. It must report no new line; the OpenCASCADE `*.hxx` includes in
  `core/shape.*` and `<ankerl/unordered_dense.h>` in `dem/data_containers.h`
  are known exceptions.

## Before editing

- Read the relevant header, source, `CMakeLists.txt`, existing tests, and
  parameter documentation before changing code. Before writing or reviewing
  C++, read `doc/source/contributing/coding-standards.rst`.
- Prefer matching a nearby solver/model/test pattern over inventing a new
  structure.
- Check whether the change affects `.prm` parameters, restart files, output
  text, MPI behavior, or release notes.

Parameters are declared in `include/core/parameters*.h` and
`source/core/parameters*.cc`, and documented in `doc/source/parameters/`.

## Build, format and check

**Do not build unless asked** — a full build is long, and the developer
compiles themselves. Prefer reasoning from the source. When asked:

```bash
cmake ../ -DCMAKE_BUILD_TYPE=Debug -G Ninja && ninja -j<N>
ctest --output-on-failure -j<N>        # -R <regex> for a subset
```

**Never hand-format.** Formatting needs clang-format **16.0 exactly** (other
versions are rejected): install it once with
`contrib/utilities/download_clang_format`. Then, from the repository root,
`contrib/utilities/indent` formats the C++ files changed since master
(`indent-all`: every file), `prmindent -i <file.prm>` formats parameter files,
and `pre-commit.sh` does both and bumps copyright years.

Cheap checks that are usually acceptable: `checkdeps`,
`check_parameter_files_indentation.sh`, and a targeted `ctest -R <regex>` when
a build already exists. Expensive, only when asked: `run_clang_tidy.sh
<src-dir>` from an empty build directory (clang 20), and
`check_parameter_files_in_examples.sh` (needs an installed `lethe-parameter-check`).

## What CI runs

`Check_indent`, `Check_prm_indent` and the Sphinx build (`-W`: any warning is
an error) run first; the jobs below start only if they pass, and are skipped
when only `doc/` and `release_notes/` change.

- Debug build + full `ctest`.
- Release build + `ctest` against deal.II master **and v9.7.1**: guard newer
  deal.II APIs with `#if DEAL_II_VERSION_GTE(9, 8, 0)`.
- GCC Debug with `-Werror -pedantic-errors`, **including `prototypes/`**.
- `-DLETHE_USE_LDV=ON`: use `GlobalVectorType` / `GlobalBlockVectorType`
  (`include/core/vector.h`), never a hard-coded Trilinos vector.
- QMOMKit build with `-Werror` and the `qmomkit*` tests.
- clang-tidy 20 on the four libraries, every warning an error.
- `lethe-parameter-check` on every example `.prm`; any deprecated parameter or
  value fails it.

## Coding standards (summary)

`doc/source/contributing/coding-standards.rst` is authoritative and wins over
this summary; naming details are in `.claude/rules/cpp.md`.

- **Names** are complete `snake_case` words (`kinematic_viscosity`, not `nu`);
  a variable holding a parameter mirrors its `.prm` name. Types are
  `UpperCamelCase`. Options are an `enum class` with `snake_case` values, never
  ints, bools or strings. Constants are named; no magic numbers.
- **Validate early**: `Patterns::` ranges in `declare_parameters`, consistency
  checks with `AssertThrow` in `parse_parameters`. Handle every `enum` value; an
  unsupported combination throws instead of silently doing nothing.
- **`Assert` vs `AssertThrow`**: `Assert`, `AssertDimension`, `AssertIndexRange`
  and `AssertIsFinite` vanish in Release; use them in loops and hot paths
  (cells, quadrature points, DoFs, particles, matrix-free operators).
  Elsewhere (setup, parsing, I/O, once per time step or solve) use
  `AssertThrow`, internal checks included. User-triggerable errors are always
  `AssertThrow`, checked once outside the loops.
- **An `AssertThrow` must fail on every rank at once**: reduce a condition on
  local data (e.g. `Utilities::MPI::logical_or`) before checking it, or the
  other ranks hang in the next collective call.
- **Error messages** state what is wrong, the offending value and how to fix it.
  Shared exceptions live in `include/core/exceptions.h`; reuse them.
- **Smart pointers**: no `new`/`delete`; `std::make_unique` by default,
  `std::shared_ptr` only for truly shared ownership. Take a
  `const std::shared_ptr<T> &` only to store it, otherwise `const T &`.
- **`const`** every variable that does not change and initialize it where it is
  declared. Pass deal.II objects and `Parameters::*` structs by `const &`, and
  prefer a whole parameters struct over many scalars.
- **Files**: new `.cc`, `.h` and `.prm` files start with
  `// SPDX-FileCopyrightText: Copyright (c) <year> The Lethe Authors` and
  `// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later`
  (`#` in `.prm`). Bump the end year of modified files (`update_copyright.sh`
  does C++ only). Header guards are `lethe_<file_name>_h`.
- **Headers**: include a header for every symbol used and none you don't use.
  No new file-scope `using namespace` in headers.
- **Classes** are implemented in the `.cc` (trivial accessors excepted), with
  explicit instantiations at its end (typically `dim = 2` and `3`).

## CMakeLists.txt

- List sources and headers explicitly; no `file(GLOB)` (except `prototypes/`).
- Adding/removing a file means updating the relevant `CMakeLists.txt`; lists are
  sorted, sources before headers.
- `target_link_libraries` lists only **direct** dependencies, never transitive
  ones (`lethe-fluid` links only `lethe-solvers`).

## What a change must include

- **Doxygen**: `@brief`, a `@param` per argument and a `@tparam` per template
  parameter, in the header, for every function.
- **Tests**: a new feature needs unit and/or application tests; a bug fix needs
  a test reproducing the bug and a PR description explaining the root cause.
- **Docs**: new/changed input-file parameters update `doc/source/parameters/`
  and every affected `.prm`; new features get a doc page and often an example.
- **Release note**: `release_notes/current/YYYY_MM_DD_LASTNAME.md` following
  `release_notes/template.md`, tagged `MAJOR` (new feature, or any parameter
  name change) or `MINOR` (fix with no impact on results). Never edit the
  versioned `release_notes/v*.md` files.

## deal.II and numerical correctness

- **MPI**: after writing into a distributed vector or matrix, call `compress()`
  with the correct `VectorOperation` (`insert` vs `add`).
- **Ghosted vs. non-ghosted vectors**: reading off-processor entries requires a
  ghosted vector; writing requires a non-ghosted (writable) one.
- **Ownership**: loops over cells/DoFs respect `cell->is_locally_owned()` and
  skip artificial cells.
- **Parallel output** goes through `pcout` / `ConditionalOStream`, never raw
  `std::cout`/`std::cerr`, which every rank would emit.
- **Floating point**: never compare with `==`; use a named tolerance. Guard or
  regularize divisions by quantities that can vanish (cell measure, density,
  time step, gradient norm).
- **Performance**: avoid copies of large deal.II objects and repeated
  reinit/allocation in cell loops (reuse `FEValues`, scratch/copy data,
  WorkStream).
- **Deprecated API**: use the current deal.II equivalent; `-Werror` fails on it.

## Pull requests

One idea per PR (one fix, one feature, or one restructuring) plus its test.
Title < 60 characters, starting with an action verb. Rebase on `master`
(`git push --force-with-lease`). Commits are squashed, so the PR description is
the permanent record: use the matching `.github/PULL_REQUEST_TEMPLATE/`
template and fill in its **Use of LLMs** section with what Claude generated.
Reviews: 2 reviewers for small changes, 3 for large or architectural ones.
Replying to reviewers is the developer's responsibility
(`doc/source/contributing/pull-requests.rst`): draft replies, never post them.

When reviewing code, prioritize bugs, MPI/distributed correctness, numerical
robustness, performance regressions, missing tests/docs/release notes, and
dependency violations. Do not nitpick formatting handled by Lethe's tools.

## Working conventions

- Don't reformat code by hand or touch regions unrelated to the change.
- Don't touch build artifacts (`build/`, `CMakeFiles/`, `CMakeCache.txt`, …)
  or simulation outputs unless explicitly asked.
- When adding a file, do all of it: source, header, CMakeLists entry, Doxygen,
  test, docs, release note.
- When unsure about the physics or numerics, say so rather than guessing — a
  plausible-looking but wrong discretization is worse than no change.
