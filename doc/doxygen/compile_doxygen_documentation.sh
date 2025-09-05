#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

################################################################################
# Bash script that generates the appropriate header file and compiles doxygen
# documentation
#
# * IMPORTANT *
#   - Make sure that you have the proper permissions to execute the file before
#     calling it. If not, you can use "chmod" to change the files permissions.
################################################################################

# Get template header, footer and stylesheet files
doxygen -w html header.html dummy.html dummy.css

# Remove unnecessary files
rm -f dummy*

# Add favicon to the header template file
sed -i '/<\/head>/i <link rel="icon" type="image/x-icon" href="lethe_favicon.ico"/>' header.html

# Compile doxygen documentation
doxygen doxygen-config-file