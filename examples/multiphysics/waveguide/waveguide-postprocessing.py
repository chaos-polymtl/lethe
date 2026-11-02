# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for the waveguide example.
This code reads the table of the L2 errors on the electric and magnetic fields
that Lethe prints at the end of the simulation, and compares their convergence
with the expected order O(h^{p+1}).
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

# Plot colors
colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the waveguide')
parser.add_argument("-l", "--log", type=str, help="File in which the output of the simulation was written", required=True)
parser.add_argument("-p", "--degree", type=int, help="Degree of the trial space of the electromagnetics", default=1)
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers = parser.parse_known_args()

# Read the last table of errors on the electromagnetic fields in the log. It
# has the following form, where every error is followed by its convergence
# rate:
# cells   error_E_real      error_E_imag      error_H_real      error_H_imag
#     4 8.631365e-02    - 9.611845e-02    - 9.351336e-02    - 8.562457e-02    -
#    32 2.115071e-02 2.03 2.437166e-02 1.98 2.252773e-02 2.05 2.192158e-02 1.97
with open(args.log) as file:
    lines = file.read().splitlines()

header_indices = [i for i, line in enumerate(lines)
                  if line.split()[:2] == ["cells", "error_E_real"]]
if (len(header_indices) == 0):
    raise RuntimeError(f"No table of errors on the electromagnetic fields was found in {args.log}")

header_index = header_indices[-1]
fields = lines[header_index].split()[1:]

cells = []
errors = []
for line in lines[header_index + 1:]:
    tokens = line.split()
    if (len(tokens) != 1 + 2 * len(fields) or not tokens[0].isdigit()):
        break
    cells.append(int(tokens[0]))
    errors.append([float(value) for value in tokens[1::2]])

cells = np.array(cells)
errors = np.array(errors)

# Element size relative to that of the coarsest mesh, in 3D
h = (cells / cells[0]) ** (-1. / 3.)

# Convergence rates between consecutive meshes
rates = np.log(errors[:-1] / errors[1:]) / np.log(h[:-1, None] / h[1:, None])

for i, field in enumerate(fields):
    print(f"{field}: rates = {rates[:, i]}")

# Plot the errors with the expected order of convergence
labels = [r"$\Re(\mathbf{E})$", r"$\Im(\mathbf{E})$", r"$\Re(\mathbf{H})$", r"$\Im(\mathbf{H})$"]
markers = ["o", "s", "^", "D"]
plt.rcParams['font.size'] = '14'
fig0 = plt.figure(figsize=(8, 6))
ax0 = fig0.add_subplot(111)
for i in range(len(fields)):
    ax0.loglog(h, errors[:, i], marker=markers[i], mfc='none', color=colors[i], label=labels[i])
order = args.degree + 1
ax0.loglog(h, errors[0, :].max() * h ** order, "--k", label=fr"$\mathcal{{O}}(h^{order})$")
ax0.set_xlabel(r"$h/h_0$")
ax0.set_ylabel(r"$\|e\|_{L^2}$")
ax0.legend(loc="lower right")
plt.tight_layout()

if (not args.validate):
    plt.savefig("waveguide-convergence.png", dpi=300)
    plt.show()
else:
    plt.savefig("waveguide-convergence.pdf")
    # Store the errors and the convergence rates so that they can be compared
    # between validation runs
    np.savetxt("solution-errors.dat",
               np.column_stack((cells, errors)),
               header="cells " + " ".join(fields))
    np.savetxt("solution-convergence-rates.dat",
               np.column_stack((cells[1:], rates)),
               header="cells " + " ".join(["rate_" + field.replace("error_", "") for field in fields]))
