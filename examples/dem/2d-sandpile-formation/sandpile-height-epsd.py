import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

prn_seeds = ['20','73','103','113','1669','2971']

# Initialize the array of heights
df_sim   = pd.read_csv('data/prn20/height_epsd.csv')
time_sim     = df_sim['time'].to_numpy()
height_sim = np.zeros((len(prn_seeds),len(time_sim)))


plt.figure(figsize=(15, 5))

# Plot the heights for the different prn seeds
plt.subplot(1,2,1)
for i, p in enumerate(prn_seeds):
    # Read the CSV files into Pandas DataFrames
    df_sim   = pd.read_csv('data/prn' + p + '/height_epsd.csv')
    height_sim[i,:] += df_sim['height']

    # Plot the heights of the simulation
    plt.plot(time_sim,height_sim[i,:],label= "Lethe-DEM epsd with prn = " + p)

plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Height of the pile (m)')
plt.yticks(np.arange(0.1, 0.32, 0.02))
plt.xticks(np.arange(0, 60, 10))
plt.grid()
plt.title("Evolution of the height of the pile with model epsd and different prn seeds")



plt.subplot(1,2,2)

#Plot the results from J. Ai et al.
df_paper = pd.read_csv('data/reference/extraction_model_epsd.csv')
time_paper   = df_paper['x'].to_numpy()
height_paper = df_paper['Curve1'].to_numpy()
plt.plot(time_paper,height_paper,  '--', label= "Ai2010 epsd")

# Plot the mean height for the different prn seeds and the margin of error
mean_height = np.mean(height_sim, axis=0)
std_dev = np.std(height_sim, axis=0 , ddof=1)
Z = 1.96
margin_error = Z * (std_dev / np.sqrt(len(prn_seeds)))
plt.plot(time_sim, mean_height, label="Mean Height", color="blue")
plt.fill_between(time_sim, mean_height - margin_error, mean_height + margin_error, color='blue', alpha=0.3, label="Margin of Error (95%)")


plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Height of the pile (m)')
plt.yticks(np.arange(0.1, 0.32, 0.02))
plt.xticks(np.arange(0, 60, 10))
plt.grid()
plt.title("Evolution of the height of the pile with model epsd")



plt.savefig("different_prn_seeds")
plt.show()
