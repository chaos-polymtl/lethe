import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



# Read the CSV files into Pandas DataFrames
df_epsd = pd.read_csv('results_epsd_noloadb_out.pvd.csv')
df_constant = pd.read_csv('results_constant_noloadb_out.pvd.csv')
df_viscous = pd.read_csv('results_viscous_noloadb_out.pvd.csv')

# Convert the DataFrames to NumPy arrays
time = df_epsd['time'].to_numpy()
height_epsd = df_epsd['height'].to_numpy()
height_constant = df_constant['height'].to_numpy()
height_viscous = df_viscous['height'].to_numpy()


plt.plot(time,height_epsd,label='epsd_resistance')
plt.plot(time,height_constant,label='constant_resistance')
plt.plot(time,height_viscous,label='viscous_resistance')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Height of the pile (m)')
plt.yticks(np.arange(0.1, 0.32, 0.02))
plt.savefig("Evolution of the height of the pile")
plt.show()
