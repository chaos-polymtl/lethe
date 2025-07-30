# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Data
dt = np.array([0.2, 0.1, 0.05, 0.01])
err = np.array([1.629277e-3, 2.232885e-4, 2.940144e-5, 2.509748e-7])

# log-log scale
log_dt = np.log10(dt)
log_err = np.log10(err)

# Linear regression
slope, intercept, r_value, p_value, std_err = linregress(log_dt, log_err)
regression_line = slope * log_dt + intercept

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(dt, err, 'o', label='Data', markersize=8)
plt.plot(dt, 10**regression_line, '--', label=f"Linear Regression (slope = {slope:.2f})")

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\Delta t$', fontsize=14)
plt.ylabel(r'L2 Error at $t=1.0$ s', fontsize=14)
plt.title('L2 Error vs $\Delta t$ (SDIRK33)', fontsize=16)
plt.grid(True, which="both", ls="--", lw=0.5)
plt.legend()
plt.tight_layout()
plt.show()
