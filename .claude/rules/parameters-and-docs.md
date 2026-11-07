---
paths:
  - "doc/**/*"
  - "examples/**/*"
  - "**/*.prm"
  - "include/core/parameters*.h"
  - "source/core/parameters*.cc"
  - "release_notes/**/*"
---

# Parameters, examples, documentation and release notes

## Adding, renaming or changing an input-file parameter

1. Declare it in the relevant `Parameters::*` class (`include/core/parameters*.h`
   and `source/core/parameters*.cc`, or the solver-specific parameter file) with
   the most restrictive pattern: `Patterns::Double(0.)` for a quantity that
   cannot be negative, `Patterns::Double(0., 1.)` for one in [0, 1]. An option
   is an `enum class` declared with `declare_enum_entry` and parsed with
   `string_to_enum` (`include/core/utilities.h`).
2. Check its consistency with the other entries with `AssertThrow` in
   `parse_parameters`, so that the error is reported before the mesh is read.
3. Name the C++ member after the entry: `number of fluids` → `number_of_fluids`.
4. Document it in `doc/source/parameters/`: the `text` code block with the
   default values, then the explanation of each entry. Keep the defaults
   identical to `declare_entry`.
5. Update every `.prm` file that uses it, in `applications_tests/` and
   `examples/`, and the parameter blocks quoted in `doc/source/examples/`,
   which are copied rather than included.
6. Keeping old names working: `prm.declare_alias(new_name, old_name, ...)` for
   an entry (see `source/core/parameters_lagrangian.cc`), and a
   `DeprecatedEnumNames` list passed to `declare_enum_entry` and
   `string_to_enum` for enum strings, which prints a deprecation warning.
   Example files must use the new names: `lethe-parameter-check` fails CI on
   any deprecation message.
7. Run `contrib/utilities/prmindent -i` on every `.prm` file you touched.
8. Add a `MAJOR` release note for a new parameter, a renamed one, or a change
   of behavior.

## Parameter files

- New `.prm` files start with the two SPDX lines as `#` comments. Bump the end
  year of the ones you modify by hand (`update_copyright.sh` only handles C++).
- CI runs `check_parameter_files_indentation.sh`, which fails if `prmindent`
  changes any `.prm` file.
- CI runs `check_parameter_files_in_examples.sh` on every example `.prm` with
  `lethe-parameter-check`. The solver family comes from the folder:
  - `incompressible-flow`, `multiphysics` and `sharp-immersed-boundary` →
    `lethe-fluid`;
  - `dem` and `dem-mp` → `lethe-particles`;
  - `unresolved-cfd-dem*` → `lethe-fluid-particles`, or `lethe-particles` if the
    file name contains `particles`.

## Sphinx documentation (`doc/`)

- CI builds with `SPHINXOPTS = -W --keep-going`: every warning fails the build
  and blocks every other CI job. Typical causes are a broken `:doc:`/`:ref:`
  target, a page missing from a `toctree`, and bad RST indentation. Build
  locally with `cd doc && make html` (requirements in `doc/requirements.txt`).
- Titles and subtitles use Chicago-style capitalization. The other format rules
  are in `doc/source/contributing/documentation.rst`.
- A new example has its files in `examples/<category>/<name>/` and its page in
  `doc/source/examples/<category>/<name>/`, listed in the category's `toctree`.

## Release notes (`release_notes/current/`)

- One file per change, named `YYYY_MM_DD_LASTNAME.md` after the date and the
  author's last name. Add `_2`, `_3`, … when that name is already taken,
  including by another open PR.
- Content follows `release_notes/template.md`: a `## [Master] - YYYY/MM/DD`
  title, one heading (Added, Changed, Deprecated, Removed, Fixed, Disabled or
  Documentation) and entries of the form `- MAJOR|MINOR <what this PR does>. [#PR](link)`.
- `MAJOR`: a new significant feature, or any change to a parameter name.
  `MINOR`: a fix or correction with no impact on the results.
- Never edit the versioned `release_notes/v*.md` files; they are assembled
  from `current/` at release time.
