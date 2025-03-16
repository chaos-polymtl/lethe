import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm

matplotlib.rcParams.update({'font.size': 22})

df = pd.read_csv("angle_data.csv")
print(df)

experimental_angle = 26 # degree at 60 Hz on the machine 
df['Error'] = (df['Angle'] - experimental_angle)
df['Abs-error'] = np.abs((df['Angle'] - experimental_angle))


fig, ax = plt.subplots(1, 2, figsize=(22, 10))

ax[0].scatter(df["Sliding friction"],df["Error"], c=df["Rolling friction"],edgecolors= "black", cmap='viridis_r', s=80)
ax[0].set_xlabel("Sliding friction")

colormap = ax[1].scatter(df["Coefficient of restitution"],df["Error"], c=df["Rolling friction"],edgecolors= "black", cmap='viridis_r', s=80)
ax[1].set_xlabel("Coefficient of restitution")

cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
fig.colorbar(colormap, cax = cbar_ax, label='Rolling friction')  # Add a colorbar to show the scale
for i in ax:
    i.set_ylabel("Error = Simulation - experiments) [degree]")
plt.savefig("error_analysis.png",dpi=300)
plt.show()





fig2, ax2 = plt.subplots(1, 1, figsize=(19, 10))
colormap = ax2.scatter(df["Sliding friction"],df["Rolling friction"], c=df["Abs-error"],edgecolors= "black", cmap='viridis_r', s=140-5*df["Abs-error"])

cbar_ax = fig2.add_axes([0.92, 0.15, 0.02, 0.7])
fig2.colorbar(colormap, cax = cbar_ax, label='Absolute-Error')
for i, v in enumerate(df["Abs-error"]):
    ax2.annotate(str(v)[0:3], xy=(df["Sliding friction"][i],df["Rolling friction"][i]), xytext=(-7,7), textcoords='offset points')
ax2.set_xlabel("Sliding friction")
ax2.set_ylabel("Rolling friction")
plt.savefig("error_analysis_2.png",dpi=300)
plt.show()
