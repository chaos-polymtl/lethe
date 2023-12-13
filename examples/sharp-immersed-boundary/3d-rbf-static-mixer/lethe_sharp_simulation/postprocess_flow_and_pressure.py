import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# We read the output files
with open('./output/flow_rate.dat', 'r') as input_file:
    lines = input_file.readlines()
    flow_in    = []
    flow_out   = []
    for line in lines[1:]:
        newLine = line.strip(' ').split()
        flow_in   .append(float(newLine[1]))
        flow_out  .append(float(newLine[2]))
        
with open('./output/pressure_drop.dat', 'r') as input_file:
    lines = input_file.readlines()
    time_steps    = []
    pressure_drop = []
    for line in lines[1:]:
        newLine = line.strip(' ').split()
        time_steps   .append(float(newLine[0]))
        pressure_drop.append(float(newLine[1]))



time_steps    = np.array(time_steps)
flow_in       = np.array(flow_in)
flow_out      = np.array(flow_out)
pressure_drop = np.array(pressure_drop)

mass_conservation = flow_in + flow_out


# We plot the mass conservation and pressure drop
fig,  ax  = plt.subplots()

ax.plot(time_steps,\
    np.abs(mass_conservation/np.max(mass_conservation)),\
    label="Mass conservation [%]",\
    linestyle="-",\
    marker="3",\
    c = "blue")

ax.plot(time_steps,\
    pressure_drop/1e4,\
    label="Pressure drop [kPa]",\
    linestyle="-",\
    marker="4",\
    c = "black")

ax.set_xlabel("Time (s)")
ax.set_ylabel("[-]")
ax.set_title("Mass Conservation and Pressure Drop in Static Mixer")
ax.set_yscale('log')
ax.set_xticks([0, 0.002, 0.004])
ax.grid(True, which="both", ls="-")
ax.legend(loc='upper left', bbox_to_anchor=(1,1),
      fancybox=True, shadow=True, ncol=1)
file_output = "./mass_and_pressure_drop.png"
fig.savefig(file_output,dpi=300,bbox_inches='tight')
file_output = "./mass_and_pressure_drop.svg"
fig.savefig(file_output,dpi=300,bbox_inches='tight')
