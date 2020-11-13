# Lethe

Lethe is an open source and we try to follow the deal.II mentality of having it
developed an open accessible environment. Consequently, all scientific
development made within Lethe are accessible to everyone even as they are
currently being written. Anybody that wishes to use and contribute to a feature
is welcomed to do so.

This document aims at establishing from guidelines for contributions within
Lethe.


# How can I contribute?

Contributions through code or documentation should be done through pull
request at the official Github repository of Lethe : https://github.com/lethe-cfd/lethe.
We recommend that users either request access to the Lethe repository or
create their own fork of Lethe on their own github account. This can then be
used to open pull requests

# What should be the content of a pull request (PR)?

A pull request should contain the following elements:

- A brief title (less than 60 characters) that describes the goal of the pull
requests
- A detailed description of the content added by the pull request:
i) If the PR adds a new feature. The feature should be documented, tests (with
  a unit tests and/or an application test) and if the features adds or alters
  parameters in the input file, the appropriate section of the wiki should
  be updated accordingly.
ii) If the pull requests corrects a bug, the source of the bug should be
explained and the way it was identified should be briefy described. A unit
test or application test that reproduces the bug should be added.

# Good practices

- Before making a pull request, the branch should be rebased on master to ensure
a linear history and to make merging as easy as possible. It is the responsability
of the branch owner to ensure that the pull request goes as smoothly as possible.
Pushing to github after a rebase requires a force push. You can accomplish
this by using `git push --force-with-lease`
- All functions should be documented in the `.h`. All functions should contain
a `@brief` that describes the function and a `@param` for each argument.
- Classes should be implemented in their `.cc`. All possibilities of the Classes
should be instantiated at the end of the `.cc` files.
- Lethe follows strict indentation guidelines. Automatic indentation can be
carried out by launching the `/contrib/utilities/indent-all.sh` script as
long as a valid version of clang-format is installed. An installation script for
clan is provided in Lethe and within the deal.II respository. Lethe follows the
same indentation guidelines as deal.II to ensure continuity and compatibility.
