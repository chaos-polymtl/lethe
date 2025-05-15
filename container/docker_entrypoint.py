#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) 2022 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import os
import subprocess
import sys
from typing import List

LETHE_TITLE = f"""
██╗     ███████╗████████╗██╗  ██╗███████╗
██║     ██╔════╝╚══██╔══╝██║  ██║██╔════╝
██║     █████╗     ██║   ███████║█████╗
██║     ██╔══╝     ██║   ██╔══██║██╔══╝
███████╗███████╗   ██║   ██║  ██║███████╗
╚══════╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚══════╝
"""

HELP_MESSAGE = f"""
{LETHE_TITLE}

Usage: <program name> [<program arg>, ...]

Available programs:
===================
"""

LETHE_INSTALL_DIR = os.getenv("LETHE_INSTALL_DIR", "/opt/lethe/")
LETHE_EXECUTABLES = sorted(os.listdir(f"{LETHE_INSTALL_DIR}/bin/"))


def print_help() -> None:
    print(HELP_MESSAGE)

    for program in LETHE_EXECUTABLES:
        print(f"▸ {program}")


def main(args: List[str]) -> None:
    # Print help if program not specified
    if len(args) < 2:
        print_help()
        sys.exit(-1)

    # Check if program is valid
    program = args[1]
    if program not in LETHE_EXECUTABLES:
        print(f"[Error] Unknown executable: {program}", file=sys.stderr)
        print_help()
        sys.exit(-1)

    # Call selected program and forward args
    subprocess.call(args[1:])


if __name__ == '__main__':
    main(sys.argv)
