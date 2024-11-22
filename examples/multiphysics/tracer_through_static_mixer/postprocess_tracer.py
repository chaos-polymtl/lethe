# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# We read the output file
with open('./output/tracer_flow_rate.dat', 'r') as input_file:
    lines = input_file.readlines()
    tracer_time       = []
    tracer_flow_in    = []
    tracer_flow_out   = []
    for line in lines[1:]:
        newLine = line.strip(' ').split()
        tracer_time      .append(float(newLine[0]))
        tracer_flow_in   .append(float(newLine[1]))
        tracer_flow_out  .append(float(newLine[2]))

tracer_time     = np.array(tracer_time)
tracer_flow_in  = np.array(tracer_flow_in)
tracer_flow_out = np.array(tracer_flow_out)


# We plot the tracer flow rates
fig,  ax  = plt.subplots()
ax.plot(tracer_time,\
    np.abs(tracer_flow_in),\
    label="Inlet",\
    linestyle="-",\
    c = "blue")

ax.plot(tracer_time,\
    np.abs(tracer_flow_out),\
    label="Outlet",\
    linestyle="-",\
    c = "black")

ax.set_xlabel("Time [s]")
ax.set_ylabel("Tracer flow rate [g/s]")
ax.set_xticks([0, 100, 200, 300, 400, 500])
ax.set_yticks([0, 250, 500, 750, 1000, 1250])
ax.grid(False)
ax.legend(loc='upper right',
      fancybox=True, shadow=True, ncol=1)
file_output = "./tracer_flow_rates.png"
fig.savefig(file_output,dpi=300,bbox_inches='tight')
file_output = "./tracer_flow_rates.svg"
fig.savefig(file_output,dpi=300,bbox_inches='tight')

