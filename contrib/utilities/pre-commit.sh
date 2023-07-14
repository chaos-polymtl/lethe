#!/bin/sh
# SPDX-FileCopyrightText: Copyright (C) 2023 The Lethe Authors
# SPDX-License-Identifier: LGPL-2.1-or-later

# To install:
#   
#   (cd .git/hooks && ln -sf ../../contrib/utilities/pre-commit.sh pre-commit)

# Indent all prm
find . -name '*.prm' -type f -execdir "$PWD/contrib/utilities/prmindent" -i '{}' \;

# Indent all .cc and .h files
contrib/utilities/download_clang_format >/dev/null 2>&1
contrib/utilities/indent-all
