import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

rollingmethod = ["constant" , "viscous" , "epsd"]

for x in rollingmethod :
    # Read the CSV files into Pandas DataFrames
    df_sim   = pd.read_csv('height_' + rollingmethod + '.csv')
    df_paper = pd.read_csv('extraction_model_' + rollingmethod + '.csv')

    # Convert the DataFrames to NumPy arrays
    time_sim     = df_sim['time'].to_numpy()
    time_paper   = df_paper['x'].to_numpy()
    height_sim   = df_sim['height'].to_numpy()
    height_paper = df_sim['Curve1'].to_numpy()

    # Plot the heights of the simulation and of the paper
    plt.plot(time_sim,height_sim,label= "Lethe-DEM " + rollingmethod)
    plt.plot(time_paper,height_paper,label= "Ai2010 " + rollingmethod)


plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Height of the pile (m)')
plt.yticks(np.arange(0.1, 0.32, 0.02))
plt.xticks(np.arange(0, 60, 10))
plt.savefig("height-comparison")
plt.show()
