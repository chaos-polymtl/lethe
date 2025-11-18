#!/bin/sh
##-----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
##-----------------------------------------------------------------------------
#
# This is a script that is used by the continuous integration servers
# to make sure that the currently checked out version of a git repository
# satisfies our our identation standars for parameter files (.prm).
# It can also be run locally to check for indentation errors before
# committing code.

echo "Running parameter files indentation test"

cd "$(git rev-parse --show-toplevel)"

# Run formatter
find . -name '*.prm' -type f \
  -exec "$PWD/contrib/utilities/prmindent" -i {} \;

# Check for changes
if ! git diff --quiet; then
  echo "The following files were modified by prmindent:"
  git diff --name-only
  echo ""
  echo "Please run prmindent locally and commit the results."
  exit 1
fi
