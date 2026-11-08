---
paths:
  - "include/**/*.h"
  - "source/**/*.cc"
  - "applications/**/*.cc"
  - "prototypes/**/*.{cc,h}"
  - "tests/**/*.cc"
---

# Writing C++ in Lethe

These rules come from `doc/source/contributing/coding-standards.rst` (authoritative),
`.clang-tidy`, `.clang-format` and the CI builds. When in doubt, imitate the
surrounding code.

## Naming

- Variables, functions and files: `snake_case`, complete words that describe the
  physical or numerical meaning (`pressure_scaling_tolerance`, not `tol`). Short
  names only for deal.II-style local indices (`q`, `i`, `j`).
- Counts start with `n_` (`n_q_points`, `n_dofs`).
- A variable storing a parameter mirrors its `.prm` name: `set number of fluids`
  → `number_of_fluids`.
- Booleans read as a statement: `is_`, `has_`, `enable_` variables;
  `density_is_constant()`, `is_non_newtonian()` functions.
- Functions start with a verb: `setup_dofs`, `assemble_system_matrix`, `get_`/`set_`.
- Members have no prefix or suffix (no `m_`, no trailing `_`). In class templates,
  access base-class members through `this->`.
- Constants are `constexpr`/`const`, `snake_case`, declared near their use.
  `ALL_CAPS` is reserved for macros.
- Classes, structs, aliases and namespaces: `UpperCamelCase`, acronyms kept
  (`CFDDEMSolver`, `GLSNavierStokesAssemblerCore`). Abstract bases end with
  `Base`; a derived class names what it specializes (`CahnHilliardFilterBase` →
  `CahnHilliardTanhFilter`).
- A file is named after its main class: `FluidDynamicsMatrixBased` lives in
  `fluid_dynamics_matrix_based.{h,cc}`, guarded by
  `lethe_fluid_dynamics_matrix_based_h`.
- Options are an `enum class` with `snake_case` values
  (`TimeSteppingMethod::bdf2`). Declare enum parameters with `declare_enum_entry`
  and parse them with `string_to_enum` (`include/core/utilities.h`): the `.prm`
  strings are the enumerator names.
- Template parameters: `dim`, `spacedim`, and `UpperCamelCase` types
  (`VectorType`, `Number`). `DeclException` names are `UpperCamelCase`, after the
  error they report (`FileDoesNotExist`).

## clang-tidy 20 (CI, every warning is an error)

It runs on `lethe-core`, `lethe-dem`, `lethe-solvers`, `lethe-fem-dem` and their
headers, not on applications or tests. The enabled checks and the reasons for
each exclusion are in `.clang-tidy`. Common triggers:

- `performance-unnecessary-value-param`: pass non-trivial parameters by
  `const &`. `modernize-pass-by-value` is disabled: do not switch to by-value +
  `std::move`.
- `performance-unnecessary-copy-initialization`, `performance-for-range-copy`:
  `const auto &x = object.get_x();`, `for (const auto &cell : ...)`.
- `performance-inefficient-vector-operation`: `reserve()` before a loop of
  `push_back`/`emplace_back` of known length.
- `performance-move-const-arg`: no `std::move` of a `const` object or of a
  trivially copyable value.
- `modernize-*`: `nullptr`; `using` aliases, never `typedef`; `= default` /
  `= delete`; `std::make_unique` / `std::make_shared`; range-based `for` when the
  index is not needed; `emplace_back`; `<cmath>`, not `<math.h>`; `true`/`false`
  for booleans; no `(void)` parameter lists.
- Signed/unsigned comparisons (also GCC `-Wsign-compare`): use matching types,
  or `std::cmp_less` and friends.
- `bugprone-integer-division`: write `1. / 2.`, not `1 / 2`, when the result is
  used as a floating-point value.
- `bugprone-throw-keyword-missing`: a constructed exception is thrown or passed
  to `AssertThrow`, never left as a bare statement.
- `cppcoreguidelines-pro-type-static-cast-downcast`: downcast polymorphic types
  with `dynamic_cast`, not `static_cast`.
- `readability-qualified-auto`: `auto *` / `const auto *` when `auto` deduces a
  pointer.
- `mpi-*`: the buffer type matches the `MPI_Datatype`. Prefer the
  `Utilities::MPI::` wrappers (`sum`, `max`, `logical_or`, ...) to raw MPI calls.

## GCC `-Werror` builds

The warnings build is Debug with `-Wall -Wextra -Wpedantic -pedantic-errors
-Werror`, **`prototypes/` included**; the QMOMKit build also uses `-Werror`.

- Unused parameter: comment out its name (`const unsigned int /*component*/`)
  or mark it `[[maybe_unused]]`, as the surrounding code does.
- Variable only read inside an `Assert` (removed in Release): `[[maybe_unused]]`.
- Unused variables and functions, signed/unsigned comparisons, a missing
  `return`, a stray `;` after a function body and deprecated deal.II calls are
  all errors.

## Code that compiles in every CI configuration

- Generic code uses `GlobalVectorType` / `GlobalBlockVectorType`
  (`include/core/vector.h`), never `TrilinosWrappers::MPI::Vector` or
  `LinearAlgebra::distributed::Vector` directly: CI builds with
  `-DLETHE_USE_LDV=ON`.
- A deal.II API newer than 9.7.1 goes under `#if DEAL_II_VERSION_GTE(9, 8, 0)`
  with a fallback (see `include/solvers/fluid_dynamics_matrix_free.h`).
- QMOMKit code goes under `#ifdef LETHE_WITH_QMOMKIT`, only in `solvers` and
  `fem-dem`.

## Doxygen and includes

- Document every function in the header, in the form used in
  `include/core/utilities.h`. clang-format does not reflow lines that start with
  these tags (`CommentPragmas` in `.clang-format`).

  ```cpp
  /**
   * @brief Return the enumerator whose identifier is @p name.
   *
   * @tparam EnumType Enumeration type.
   *
   * @param[in] name Identifier of the enumerator.
   *
   * @return The enumerator of @p EnumType named @p name.
   */
  ```

- clang-format regroups and sorts includes (`IncludeBlocks: Regroup`): Lethe
  modules first, then deal.II, Boost, other third-party, and the standard
  library. Add an include anywhere in the block and run
  `contrib/utilities/indent`.
