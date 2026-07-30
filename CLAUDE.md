# CLAUDE.md

Guidance for Claude Code (and other AI assistants) working in the Lethe repository.

Lethe is an open-source CFD / DEM / CFD-DEM solver built on **deal.II** (≥ 9.7,
with MPI, p4est and Trilinos). Assume an MPI-parallel, template-heavy,
high-performance C++ codebase. Docs: https://chaos-polymtl.github.io/lethe/

## Layout

- `include/`, `source/` — headers (`.h`) and sources (`.cc`), each split into
  `core/`, `dem/`, `solvers/`, `fem-dem/`.
- `applications/` — one executable per folder; minimal main files (< 100 lines)
  that instantiate the templated solvers for a given `dim`. Driven by `.prm`.
- `applications_tests/`, `tests/` — functional tests (per application) and unit
  tests (per module).
- `doc/` (Sphinx site), `examples/` (documented cases), `prototypes/`
  (throwaway, unstable), `contrib/utilities/` (maintenance scripts),
  `release_notes/`.

**Library dependency hierarchy — never violate.** core depends on nothing; DEM
and solvers depend only on core; FEM-DEM depends on core, DEM and solvers.
`contrib/utilities/checkdeps` must exit 0.

## Before editing

- Read the relevant header, source, `CMakeLists.txt`, existing tests, and
  parameter documentation before changing code.
- Prefer matching a nearby solver/model/test pattern over inventing a new
  structure.
- Check whether the change affects `.prm` parameters, restart files, output
  text, MPI behavior, or release notes.
- Do not modify generated build files or simulation outputs unless explicitly
  asked.

Useful places to search:

- Parameter declarations/parsing: `include/core/parameters*.h`,
  `source/core/parameters*.cc`, and solver-specific parameter files.
- Application tests: `applications_tests/<app>/`.
- Unit tests: `tests/<module>/`.
- Parameter docs: `doc/source/parameters/`.
- Examples: `examples/`.

## Build and test

**Do not build unless asked** — a full build is long, and the developer
compiles themselves. Prefer reasoning from the source. When asked:

```bash
cmake ../ -DCMAKE_BUILD_TYPE=Debug -G Ninja && ninja -j<N>
ctest --output-on-failure -j<N>        # -R <regex> for a subset
```

Tests diff stdout against golden `.output` files; regenerate them (only after
verifying the new output is correct) with `contrib/utilities/update-golden.tl <build-dir>`.
Be careful changing user-facing output: many tests compare stdout against
`.output` files. If output changes, explain why the new output is intended
before regenerating goldens.

Cheap checks that are usually acceptable:

- `contrib/utilities/checkdeps`.
- `contrib/utilities/prmindent -i <file.prm>`.
- `contrib/utilities/check_parameter_files_indentation.sh`.
- Targeted `ctest -R <regex> --output-on-failure` only when a build already
  exists or the user asks.

CI on every PR: Debug build + full ctest, warnings-as-errors GCC build,
clang-tidy, clang-format, `prmindent`, and example-parameter-file checks.

## Style

Never hand-format. Run `contrib/utilities/indent-all` (sources),
`contrib/utilities/prmindent -i file.prm` (parameter files), or
`contrib/utilities/pre-commit.sh` (both).

- SPDX header on every file; header guards named `lethe_<file_name>_h`.
- Don't add new file-scope `using namespace` in headers (existing
  `using namespace dealii;` in old headers stays).
- `const` correctness everywhere: pass deal.II objects (`Vector`, `Tensor`,
  `Point`, matrices, `FEValues` data) and `Parameters::*` structs by `const &`.
  Prefer passing a whole parameters struct over threading many scalars.
- No raw `new`/`delete` — use smart pointers.
- Implement classes in the `.cc`, not inline in the header (trivial accessors
  excepted), and explicitly instantiate every used instantiation at the end of
  the `.cc` (typically `dim = 2` and `3`).
- Include a header for **every** symbol used, and none you don't use — no
  reliance on transitive includes.

## CMakeLists.txt

- List sources and headers explicitly; no `file(GLOB)` (except `prototypes/`).
- Adding/removing a file means updating the relevant `CMakeLists.txt`; lists are
  sorted, sources before headers.
- `TARGET_LINK_LIBRARIES` lists only **direct** dependencies, never transitive
  ones (`lethe-fluid` links only `lethe-solvers`).

## What a change must include

- **Doxygen**: `@brief` and a `@param` per argument, in the header, for every
  function.
- **Tests**: a new feature needs unit and/or application tests; a bug fix needs
  a test reproducing the bug, and a PR description explaining the root cause.
  Unit tests are a `.cc` + `.output` + `.debug`/`.release` marker directory in
  `tests/<module>/`. Application tests are a `.prm` + `.output` sharing a prefix
  in `applications_tests/<app>/`; `foo.mpirun=2.output` runs on 2 ranks. Stay at
  ≤ 3 ranks and keep tests fast. Tests should be deterministic across MPI ranks
  and platforms. Avoid tests that depend on wall-clock timing, unordered output,
  random seeds without fixed values, or excessive mesh/problem sizes.
- **Docs**: new/changed input-file parameters update
  `doc/source/parameters/`; new features get a doc page and often an example.
- **Release note**: `release_notes/current/YYYY_MM_DD_LASTNAME.md` following
  `release_notes/template.md`, tagged `MAJOR` (new feature, or any parameter
  name change) or `MINOR` (fix with no impact on results). Never edit
  `CHANGELOG.md` — it is assembled from these at release time. Use one of:
  Added, Changed, Deprecated, Removed, Fixed, Disabled, Documentation.

When adding, renaming, or changing an input-file parameter:

- Update the declaring/parsing code in the relevant `Parameters::*` class.
- Update `doc/source/parameters/`.
- Update affected examples/tests `.prm`.
- Run `prmindent`.
- Add a `MAJOR` release note if the parameter name or behavior changes.

## deal.II and numerical correctness

The mistakes that matter most here:

- **MPI / distributed correctness**: after writing into a distributed vector or
  matrix, call `compress()` with the correct `VectorOperation` (`insert` vs
  `add`).
- **Ghosted vs. non-ghosted vectors**: reading off-processor entries requires a
  ghosted vector; writing requires a non-ghosted (writable) one.
- **Locally owned vs. locally relevant**: loops over cells/DoFs must respect
  ownership — `cell->is_locally_owned()` — and skip artificial cells.
- **Parallel output**: user-facing prints go through a rank-restricted stream
  (`pcout` / `ConditionalOStream`), never raw `std::cout`/`std::cerr`, which
  every rank would emit. This also keeps test output deterministic.
- **`Assert` vs `AssertThrow`**: `Assert` is compiled out in Release. Use it
  only for internal invariants that cannot fail in correct code. Validation of
  user input or of runtime conditions must use `AssertThrow`.
- **Floating point**: never compare with `==`; use a tolerance.
- **Numerical robustness**: guard or regularize divisions by quantities that can
  vanish (cell measure, density, time step, gradient norm). Avoid magic numbers
  and hard-coded tolerances; give them names.
- **Performance**: avoid copies of large deal.II objects, avoid repeated
  reinit/allocation inside cell loops (reuse `FEValues`, scratch/copy data,
  WorkStream), and hoist allocations out of tight loops.
- **Deprecated API**: prefer the current deal.II equivalent over deprecated
  functions; the warnings-as-errors CI build will catch them.

## Pull requests

One idea per PR (one fix, one feature, or one restructuring) plus its test.
Title < 60 characters, starting with an action verb. Rebase on `master`
(`git push --force-with-lease`). Commits are squashed, so the PR description is
the permanent record — use the appropriate PR template from
`.github/PULL_REQUEST_TEMPLATE/` (new feature, bug fix, documentation,
refactoring, or example). Reviews: 2 reviewers for small changes, 3 for large
features or architectural ones.

When reviewing code, prioritize bugs, MPI/distributed correctness, numerical
robustness, performance regressions, missing tests/docs/release notes, and
dependency violations. Do not nitpick formatting handled by Lethe's tools.

## Working conventions

- Don't reformat code by hand or touch regions unrelated to the change.
- Don't touch build artifacts (`build/`, `CMakeFiles/`, `CMakeCache.txt`, …),
  which may litter the working tree.
- When adding a file, do all of it: source, header, CMakeLists entry, Doxygen,
  test, docs, release note.
- When unsure about the physics or numerics, say so rather than guessing — a
  plausible-looking but wrong discretization is worse than no change.
