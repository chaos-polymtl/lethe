#!/usr/bin/env bash
# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
#
# Update SPDX copyright years in files changed since the split from master/main.
# To run this script in debug mode, run DEBUG=1 contrib/utilities/update_copyright.sh from the 
# repository root.

set -euo pipefail

CURRENT_YEAR=$(date +%Y)
DEBUG=${DEBUG:-0}
updated_any=0

# --- Find a reasonable base commit (try origin/master, master, origin/main, main, else initial commit) ---
BASE_COMMIT=""
for cand in origin/master master origin/main main; do
  if git show-ref --verify --quiet "refs/remotes/${cand}" 2>/dev/null; then
    BASE_COMMIT=$(git merge-base HEAD "refs/remotes/${cand}" 2>/dev/null || true)
    [ -n "$BASE_COMMIT" ] && break
  elif git show-ref --verify --quiet "refs/heads/${cand##*/}" 2>/dev/null; then
    BASE_COMMIT=$(git merge-base HEAD "${cand##*/}" 2>/dev/null || true)
    [ -n "$BASE_COMMIT" ] && break
  fi
done

# fallback to initial commit if merge-base detection failed
if [ -z "${BASE_COMMIT:-}" ]; then
  BASE_COMMIT=$(git rev-list --max-parents=0 HEAD)
fi

[ "$DEBUG" -eq 1 ] && printf 'DEBUG: using BASE_COMMIT=%s\n' "$BASE_COMMIT"

# Quick check: are there any changed files since base? (non-NUL check is fine here)
if [ -z "$(git diff --name-only "$BASE_COMMIT" --)" ]; then
  [ "$DEBUG" -eq 1 ] && echo "DEBUG: no files changed since base commit."
  exit 0
fi

# --- Main loop: read NUL-separated file list directly from git (do NOT store into a variable) ---
while IFS= read -r -d '' file; do
  [ -n "$file" ] || continue

  # skip deleted / non-regular files
  if [ ! -f "$file" ]; then
    [ "$DEBUG" -eq 1 ] && printf 'DEBUG: skipping non-regular: %s\n' "$file"
    continue
  fi

  # skip binary files
  if ! grep -Iq . -- "$file"; then
    [ "$DEBUG" -eq 1 ] && printf 'DEBUG: skipping binary: %s\n' "$file"
    continue
  fi

  # only process files containing the SPDX marker
  if ! grep -q 'SPDX-FileCopyrightText' -- "$file"; then
    [ "$DEBUG" -eq 1 ] && printf 'DEBUG: no SPDX marker in: %s\n' "$file"
    continue
  fi

  # Update the year range using a robust perl regex
  env CURRENT_YEAR="$CURRENT_YEAR" \
  perl -0777 -pe '
    s{^( [ \t]* // [ \t]* SPDX-FileCopyrightText: [ \t]* Copyright [ \t]* \(c\) [ \t]* ([0-9]{4}) (?:-([0-9]{4}))? [ \t]+ The [ \t]+ Lethe [ \t]+ Authors \b ) }
     {
       my $full = $1;
       my ($start,$end) = ($2, $3 // $2);
       if ($end < $ENV{CURRENT_YEAR}) {
         "// SPDX-FileCopyrightText: Copyright (c) $start-$ENV{CURRENT_YEAR} The Lethe Authors"
       } else {
         $full
       }
     }exm;
  ' -i -- "$file"


# Note: the process substitution below feeds the NUL-separated list directly to the loop
done < <(git diff --name-only -z "$BASE_COMMIT" --)

if [ "$updated_any" -eq 1 ]; then
  echo "Copyright years updated."
fi

exit 0
