#!/bin/sh
# SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# To install:
#   
#   (cd .git/hooks && ln -sf ../../contrib/utilities/pre-commit.sh pre-commit)

# Indent all prm
find . -name '*.prm' -type f -execdir "$PWD/contrib/utilities/prmindent" -i '{}' \;

# Indent all .cc and .h files
contrib/utilities/download_clang_format >/dev/null 2>&1
contrib/utilities/indent-all

# Update copyright years in modified files
contrib/utilities/update_copyright.sh
