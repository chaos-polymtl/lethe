#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

################################################################################
# Bash script that generates the appropriate header file, software version and
# compiles doxygen documentation
#
# * IMPORTANT *
#   - Make sure that you have the proper permissions to execute the file before
#     calling it. If not, you can use "chmod" to change the files permissions.
################################################################################

set -euo pipefail

# ------------------------------------------
# Extract Lethe version from CMakeLists.txt
# ------------------------------------------

CMAKE_FILE="../../CMakeLists.txt"

# Check if file exists
if [[ ! -f "$CMAKE_FILE" ]]
then
  echo "ERROR: $CMAKE_FILE cannot be found in the current directory: $(pwd)" >&2
  exit 1
fi

# Get version number
PROJECT_VERSION="$(
  sed -nE 's/.*project[[:space:]]*\([^)]*VERSION[[:space:]]+([0-9]+(\.[0-9]+){1,2}).*/\1/ip' "$CMAKE_FILE" \
  | head -n 1
)"
echo "Detected Lethe version = $PROJECT_VERSION"

# ---------------------
# Generate header file
# ---------------------

# Get template header, footer and stylesheet files
doxygen -w html header.html dummy.html dummy.css

# Remove unnecessary files
rm -f dummy*

# Add favicon to the header template file
sed -i '/<\/head>/i <link rel="icon" type="image/x-icon" href="lethe_favicon.ico"/>' header.html

# Add toggle button for dark mode and initialize
sed -i '/<\/head>/i <script type="text/javascript" src="$relpath^doxygen-awesome-darkmode-toggle.js"></script>' header.html
sed -i '/<\/head>/i <script type="text/javascript"> DoxygenAwesomeDarkModeToggle.init() </script>' header.html

# Add GitHub link
sed -i '/<\/head>/i <script type="text/javascript" src="$relpath^doxygen-github-widget.js"></script>' header.html

# --------------------------------------------------------------------------
# Build doxygen documentation, overriding PROJECT_NUMBER from CMake version
# --------------------------------------------------------------------------

# Compile doxygen documentation
DOXFILE="doxygen-config-file"

doxygen - <<EOF
@INCLUDE = $DOXFILE
PROJECT_NUMBER = "version $PROJECT_VERSION"
EOF