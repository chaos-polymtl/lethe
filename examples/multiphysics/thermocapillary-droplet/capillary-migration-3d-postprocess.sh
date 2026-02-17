# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#!/bin/bash
mkdir -p figures

python capillary-migration-3d.py  -p capillary-migration-3d-AMR-pro-f001-eps06/output/ -g capillary-migration-3d-AMR-geo-f001-eps06/output/ -a capillary-migration-3d-AMR-pde-f001-eps04/output/ -l /home/hepap/lethe/lethe/ -fs true

mv *.png figures/.
mv *.svg figures/.

