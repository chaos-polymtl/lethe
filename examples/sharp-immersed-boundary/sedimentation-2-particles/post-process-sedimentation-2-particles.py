# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for the contact between two sedimenting spheres (drafting,
kissing and tumbling).
This code compares the particles with Fig. 12 of Barbeau et al. [1] and with
the results of Nie and Lin [2] and of Apte et al. [3] that are shown on it.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Plot colors
colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']

# Definition of the case (cm, s)
RADIUS = 1. / 12.
INITIAL_HEIGHTS = (3.5, 3.16)
TIME_END = 0.7

# The definitions below are the same for all discretizations and for the curves
# of the paper, so that their metrics can be compared.
# The spheres kiss when their vertical gap is below one cell of 1/128 cm. The
# vertical gap is used because the paper only gives the vertical velocities.
KISSING_GAP = 1. / 128.
# The curves stop when a sphere is this close to the bottom
BOTTOM_GAP = 1. / 256.
# End of the window of the errors, which precedes the kissing of all the curves
DRAFTING_END = 0.28
# End of the window in which the peak velocity of p0 (at kissing) is searched
PEAK_WINDOW_END = 0.4
# Time step of the sampling of the data compared between validation runs
SAMPLING_STEP = 0.01

DISCRETIZATIONS = {
    "q1q1-h128": ("out-q1q1-h128", "Q1Q1 1/128"),
    "q2q1-h128": ("out-q2q1-h128", "Q2Q1 1/128"),
}
LITERATURE = (("nie-lin-2020", "Nie & Lin (2020)", 'o'),
              ("apte-2009", "Apte et al. (2009)", 's'))

script_folder = Path(__file__).resolve().parent

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the sedimentation of two particles')
parser.add_argument("-f", "--folder", type=Path, help="Folder which contains the results of the simulation (default: the output path of the discretization)")
parser.add_argument("-d", "--discretization", choices=DISCRETIZATIONS, default="q1q1-h128", help="Discretization of the simulation and of the curves of the paper it is compared with")
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers = parser.parse_known_args()
default_folder, label = DISCRETIZATIONS[args.discretization]
folder = args.folder or Path(default_folder)


def read_reference(name):
    return np.loadtxt(script_folder / "reference" / f"{name}.dat", unpack=True)


def integrate(time, velocity, initial_height):
    increments = 0.5 * (velocity[1:] + velocity[:-1]) * np.diff(time)
    return initial_height + np.concatenate(([0.], np.cumsum(increments)))


def first_time(time, condition):
    indices = np.flatnonzero(condition)
    return time[indices[0]] if indices.size else np.nan


def get_events(time, velocities, heights):
    """Return the metrics that can be computed from vertical motions only."""
    return {
        "peak_velocity_p0": np.nanmin(velocities[0][time <= PEAK_WINDOW_END]),
        "kissing_time": first_time(time, heights[0] - heights[1] - 2 * RADIUS <= KISSING_GAP),
        "overtaking_time": first_time(time, heights[0] < heights[1]),
        "bottom_time": first_time(time, np.fmin(heights[0], heights[1]) - RADIUS <= BOTTOM_GAP),
    }


def rms_error(time, velocity, time_reference, velocity_reference):
    drafting = time_reference <= error_window_end
    error = np.interp(time_reference[drafting], time, velocity) - velocity_reference[drafting]
    return np.sqrt(np.mean(error**2))


# Velocity and position of the particles written by Lethe. The initial state,
# which Lethe does not write, is added.
particles = [pd.read_csv(folder / f"ib_force.{p:02d}.dat", sep=r"\s+").dropna() for p in range(2)]
# During a simulation, a file can be one time step ahead of the other
n_steps = min(len(particle) for particle in particles)
particles = [particle.iloc[:n_steps] for particle in particles]
t = np.insert(particles[0]["time"].to_numpy(), 0, 0.)


def get_column(p, name, initial_value):
    return np.insert(particles[p][name].to_numpy(), 0, initial_value)


v_z = [get_column(p, "v_z", 0.) for p in range(2)]
x = [get_column(0, "p_x", 0.5 + 1. / 300.), get_column(1, "p_x", 0.5)]
y = [get_column(p, "p_y", 0.5) for p in range(2)]
z = [get_column(p, "p_z", INITIAL_HEIGHTS[p]) for p in range(2)]
# During a simulation, the errors only cover the simulated time
error_window_end = min(DRAFTING_END, t[-1])
gap = np.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2) - 2 * RADIUS

# Curves of the paper. The two spheres do not always stop at the same time, so
# they are sampled on a common time, with NaN after their end.
paper = [read_reference(f"barbeau-2024-{args.discretization}-p{p}") for p in range(2)]
paper = [(np.insert(t_p, 0, 0.), np.insert(v_p, 0, 0.)) for t_p, v_p in paper]
t_paper = np.arange(0., max(t_p[-1] for t_p, _ in paper), 1e-3)
v_paper = [np.interp(t_paper, t_p, v_p, right=np.nan) for t_p, v_p in paper]
z_paper = [integrate(t_paper, v_paper[p], INITIAL_HEIGHTS[p]) for p in range(2)]

events = get_events(t, v_z, z)
events_paper = get_events(t_paper, v_paper, z_paper)
# The curves of the paper end when a sphere reaches the bottom
events_paper["bottom_time"] = max(t_p[-1] for t_p, _ in paper)

# Rows of the metrics: value of Lethe and value of the paper
metrics = {name: (events[name], events_paper[name]) for name in events}
metrics["minimum_gap"] = (gap.min(), np.nan)
drafting_times = np.arange(0., error_window_end, 1e-3)
for p in range(2):
    metrics[f"rms_error_paper_p{p}"] = (
        rms_error(t, v_z[p], drafting_times, np.interp(drafting_times, *paper[p])), 0.)
    for name, _, _ in LITERATURE:
        reference = read_reference(f"{name}-p{p}")
        metrics[f"rms_error_{name}_p{p}"] = (rms_error(t, v_z[p], *reference),
                                             rms_error(*paper[p], *reference))

print(f"Simulated up to t = {t[-1]:.3f} s, RMS errors up to t = {error_window_end:.3f} s")
print(f"{'metric':32s} {'Lethe':>10s} {'paper':>10s}")
for name, (value, value_paper) in metrics.items():
    print(f"{name:32s} {value:10.4f} {value_paper:10.4f}")

shown = t <= (events["bottom_time"] if np.isfinite(events["bottom_time"]) else t[-1])

plt.rcParams['font.size'] = '13'

# Vertical velocity, one panel per sphere
fig0, axes0 = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
for p, ax in enumerate(axes0):
    ax.plot(t[shown], v_z[p][shown], color=colors[p], lw=2, label=f"Lethe {label}")
    ax.plot(*paper[p], '--', color='k', lw=2, label=f"Barbeau et al. (2024) {label}")
    # The curves of the paper for the other discretizations
    for other, (_, other_label) in DISCRETIZATIONS.items():
        if other != args.discretization:
            ax.plot(*read_reference(f"barbeau-2024-{other}-p{p}"), ':', color='k', lw=2,
                    label=f"Barbeau et al. (2024) {other_label}")
    for name, author, marker in LITERATURE:
        ax.plot(*read_reference(f"{name}-p{p}"), marker, mfc='none', color='k', label=author)
    ax.set_title(f"p{p} ({'upper' if p == 0 else 'lower'} sphere)")
    ax.set_xlabel("Time [s]")
    ax.set_xlim(0, TIME_END)
    ax.legend(loc="upper center", fontsize=10)
axes0[0].set_ylabel(r"$v_z$ [cm/s]")
fig0.tight_layout()

# Gap and trajectories
fig1, (ax_gap, ax_x, ax_z) = plt.subplots(1, 3, figsize=(15, 6))
ax_gap.plot(t[shown], gap[shown] / (2 * RADIUS), color=colors[2], lw=2)
ax_gap.axhline(0, color='k', lw=1)
ax_gap.set_ylabel(r"Gap / $d_p$ [-]")
for p in range(2):
    ax_x.plot(t[shown], x[p][shown], color=colors[p], lw=2, label=f"p{p}, x")
    ax_x.plot(t[shown], y[p][shown], '--', color=colors[p], lw=2, label=f"p{p}, y")
    ax_z.plot(t[shown], z[p][shown], color=colors[p], lw=2, label=f"p{p}, Lethe")
    ax_z.plot(t_paper, z_paper[p], '--' if p == 0 else ':', color='k', lw=2,
              label=f"p{p}, paper (integrated)")
ax_x.set_ylabel("Horizontal position [cm]")
ax_z.set_ylabel("z [cm]")
for ax in (ax_gap, ax_x, ax_z):
    ax.set_xlabel("Time [s]")
    ax.set_xlim(0, TIME_END)
ax_x.legend(fontsize=10)
ax_z.legend(fontsize=10)
fig1.tight_layout()

if (not args.validate):
    fig0.savefig(f"velocity-comparison-{args.discretization}.png", dpi=300)
    fig1.savefig(f"trajectories-{args.discretization}.png", dpi=300)
    plt.show()
else:
    # The files of every discretization have their own suffix, so that the
    # validation can process several discretizations
    suffix = args.discretization
    fig0.savefig(f"two-particles-velocity-{suffix}.pdf")
    fig1.savefig(f"two-particles-trajectories-{suffix}.pdf")
    # The data is sampled at fixed times, so that the files always have the
    # same number of entries, which compare_validation_runs.py requires.
    # The times that a simulation did not reach are NaN.
    samples = np.linspace(0., TIME_END, round(TIME_END / SAMPLING_STEP) + 1)
    t_end = t.copy()
    # The last time of Lethe can be short of the end time by round-off
    if abs(t_end[-1] - TIME_END) < 1e-6:
        t_end[-1] = TIME_END
    for p in range(2):
        np.savetxt(f"solution-particle-{p}-{suffix}.dat",
                   np.column_stack([samples] + [np.interp(samples, t_end, value, right=np.nan)
                                                for value in (v_z[p], x[p], y[p], z[p])]),
                   fmt="%.6f", header="t v_z x y z")
    with open(f"solution-metrics-{suffix}.dat", "w") as output:
        output.write("# metric lethe paper\n")
        for name, (value, value_paper) in metrics.items():
            output.write(f"{name} {value:.5f} {value_paper:.5f}\n")

# [1] L. Barbeau, S. Golshan, J. Deng, S. Etienne, C. Beguin and B. Blais,
# "High-order moving immersed boundary and its application to a resolved
# CFD-DEM model", Computers & Fluids 268 (2024) 106094,
# doi:10.1016/j.compfluid.2023.106094.
# [2] D. Nie and J. Lin, "Simulation of sedimentation of two spheres with
# different densities in a square tube", J. Fluid Mech. 896 (2020) A12,
# doi:10.1017/jfm.2020.291.
# [3] S. V. Apte, M. Martin and N. A. Patankar, "A numerical method for fully
# resolved simulation (FRS) of rigid particle-flow interactions in complex
# flows", J. Comput. Phys. 228 (2009) 2712-2738.
