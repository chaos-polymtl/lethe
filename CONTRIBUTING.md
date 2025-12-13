# Lethe

Lethe is an open source software, and we try to follow the deal.II mentality of having it
developed in an accessible environment. Consequently, all scientific
developments made within Lethe are available to everyone even as they are
currently being written. Anybody that wishes to use and contribute to a feature
is welcomed to do so.

This document aims at establishing guidelines for contributions within Lethe.

Lethe has its own documentation (user guide and developer guide) which are available on: https://chaos-polymtl.github.io/lethe/index.html


# How can I contribute?

Contributions through code or documentation should be done through pull
request at the official GitHub repository of Lethe : https://github.com/chaos-polymtl/lethe.
We recommend that users either request access to the Lethe repository or
create their own fork of Lethe on their own GitHub account. This can then be
used to open pull requests

# What should be the content of a pull request (PR)?

A pull request should contain the following elements:

- A brief title (less than 60 characters) that describes the goal of the pull
request.
- A detailed description of the content added by the pull request:
  * If the PR adds a new feature:\
    The feature should be documented, tested (with unit tests and/or application tests) 
    and if the feature adds or alters parameters in the input file, the appropriate section 
    of the documentation should be updated accordingly.
  * If the pull request corrects a bug:\
    the source of the bug should be explained and the way it was identified should be briefly described. 
    A unit test or application test that reproduces the bug should be added.

See [this page](https://chaos-polymtl.github.io/lethe/documentation/contributing.html#pull-requests) for more information about the pull request process.


# Good practices

- Before making a pull request, the branch should be rebased on master to ensure
a linear history and to make merging as easy as possible. It is the responsibility
of the branch owner to ensure that the pull request goes as smoothly as possible.
Pushing to GitHub after a rebase requires a force push. You can accomplish
this by using `git push --force-with-lease`
- All functions should be documented in the `.h`. All functions should contain
a `@brief` that describes the function and a `@param` for each argument.
- Classes should be implemented in their `.cc`. All possibilities of the classes
should be instantiated at the end of the `.cc` files.
- Lethe follows strict indentation guidelines. Automatic indentation can be
carried out by launching the `/contrib/utilities/indent-all.sh` script as
long as a valid version of clang-format is installed. An installation script for
clang is provided in Lethe and within the deal.II repository. Lethe follows the
same indentation guidelines as deal.II to ensure continuity and compatibility.

# Specific guidelines

## Include appropriate headers for all symbols in use

Always include appropriate headers for all symbols used in source files
and headers (and conversely, don't include headers that don't declare
any symbols in use).
Although the build may still succeed otherwise, for example if a
required header is included by another included header, it is at risk of
failing whenever headers change, for instance when updating
dependencies.

Unfortunately, keeping headers up to date is an onerous task — because
keeping track of all symbols used in a source file and the headers they
are declared in is nontrivial — and thus some files do not follow this
guideline.
([Include What You Use][] should aid in this endeavor, but it currently
fails to build Lethe.)

[Include What You Use]: https://github.com/include-what-you-use/include-what-you-use

## Adding and removing files

All source files and headers need to be explicitly listed in the
`CMakeLists.txt` files' `ADD_EXECUTABLE` and `ADD_LIBRARY` commands.
(In other words, the `CMakeLists.txt` files eschew CMake's `file(GLOB)`
command.)
If any source files or headers are omitted, the build may fail, or even
worse, succeed.

When adding and removing files, update the file lists in the relevant
`CMakeLists.txt` files, and make sure the lists are sorted, with sources
appearing before headers.
See [the `CMakeLists.txt` of `lethe-core`](source/core/CMakeLists.txt)
for an example.

## Specify only direct dependencies in linked libraries

Always list the direct dependencies of libraries and executables (as
determined by the headers included in the source files and headers of
those libraries and executables) in the `CMakeLists.txt` files'
`TARGET_LINK_LIBRARIES` calls, and omit any transitive dependencies.

For example, the `lethe-particles` application
[links to `lethe-core` and `lethe-dem`](applications/lethe-particles/CMakeLists.txt)
because it
[includes core and DEM headers](applications/lethe-particles/dem.cc),
but the `lethe-fluid` application
[links only to `lethe-solvers`](applications/lethe-fluid/CMakeLists.txt)
because it
[includes only solvers headers](applications/lethe-fluid/fluid_dynamics_matrix_based.cc),
and `lethe-core` is a transitive dependency of `lethe-solvers`.

## Dependencies between Lethe's libraries

Lethe consists of four libraries: core, DEM, solvers and FEM-DEM.
They all depend on deal.II (and Boost), but they also depend on each
other, as follows:

- Core depends on none of the others.
- DEM depends only on core.
- Solvers depends only on core.
- FEM-DEM depends on core, DEM, and solvers.

So for example the DEM library must not include headers from the solvers
library.

Don't introduce additional dependencies between Lethe's libraries.
The `contrib/utilities/checkdeps` utility detects such situations:
Run it from the Lethe's top-level directory and ensure its exit status
is 0 (otherwise it prints the offending include directives, which you
should fix).
