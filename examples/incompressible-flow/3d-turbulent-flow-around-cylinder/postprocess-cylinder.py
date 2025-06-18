import numpy as np
import matplotlib.pyplot as plt
import argparse
from cycler import cycler

from cylinder_functions import compute_cp_average, compute_strouhal_from_lift, compute_drag_coefficient

# Set plot parameters
plt.rcParams['lines.markersize'] = '11'
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.fancybox'] = False
plt.rcParams['font.size'] = 16
plt.rcParams['legend.handlelength'] = 2
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['figure.figsize'] = (10,8)
colors = ['#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02']
plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
# Common parameters
rho = 1.0
U_inf = 1.0
center = (0, 0)
theta_array = np.linspace(0, 361 * np.pi/360, 362)
L = 1.5*np.pi
z_array = np.linspace(0,L,200)
n_final_timesteps = 2000
D = 1.0
frac = 50

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 3d turbulent flow around a cylinder')
parser.add_argument("-f","--folder", type=str, nargs='+', help="Paths the output files of the simulation", required=True)
parser.add_argument("-l","--label", type=str, nargs='+', help="Label for the plot", required=True)
args, leftovers=parser.parse_known_args()

folders = args.folder
labels = args.label
print(f"Folders: {folders}")
print(f"Labels: {labels}")
# Check if the number of folders and labels match
if len(folders) != len(labels):
    raise ValueError("The number of folders and labels must match.")

# Populate the simulations dictionary
simulations = {}
for folder, label in zip(folders, labels):
    simulations[folder] = label
    
# Loop through all simulations
for folder, label in simulations.items():
    theta_array = np.linspace(0, 181 * np.pi/180, 182)
    # Strouhal number
    force_file = f"{folder}/force.02.dat"
    st_peak = compute_strouhal_from_lift(force_file, D, U_inf, method="peak", min_peak_distance=10, purcentage=frac)
    st_fft = compute_strouhal_from_lift(force_file, D, U_inf, method="fft", min_peak_distance=10, purcentage=frac)
    print(f"[{label}] Strouhal number (peak): St = {st_peak:.4f}")
    print(f"[{label}] Strouhal number (fft):  St = {st_fft:.4f}")

    # Drag coefficient
    cd = compute_drag_coefficient(force_file, D, L, rho, U_inf, n_final_timesteps)
    print(f"[{label}] Mean drag coefficient Cd = {cd:.4f}")

# Experimental data
deg_exp = [9.158187511596124, 18.482886695571118, 28.806661605164106, 38.96392484237822, 48.78815335214513, 58.44587588087789, 68.76965648181557, 78.09435566579059, 88.75114822879661, 99.57446953821567, 108.89916872219067, 119.38945530416254, 129.2136951966189, 140.03699374065909, 149.02868096256574, 159.18593281709042, 168.84366672851263, 179.33395331048447]
cp_exp = [0.9245808357735688, 0.7297486303540259, 0.3275138038679759, -0.16270942100085994, -0.5900837898819916, -0.9734637447552554, -1.1682962723995107, -1.092877215581317, -0.9546088731424699, -0.8917600171547664, -0.8603355891609146, -0.8666200451267354, -0.8791898163242766, -0.8666200451267354, -0.8729049307255057, -0.8854749167395213, -0.8854749167395213, -0.8791898163242766]

# Initialize plot
plt.scatter(deg_exp, cp_exp, label='Experimental data', color='black')

# Loop through all simulations
for folder, label in simulations.items():
    cp_mean, theta_array, p_inf = compute_cp_average(
        folder=folder,
        theta_array=theta_array,
        z_array=z_array,
        cylinder_center=center,
        rho=rho,
        U_inf=U_inf,
    )

    # Align from stagnation point
    index_pt_stag = np.argmax(cp_mean)
    aligned_theta = np.rad2deg(theta_array[index_pt_stag]) - np.rad2deg(theta_array[:index_pt_stag])
    aligned_cp = cp_mean[:index_pt_stag]
    plt.plot(aligned_theta, aligned_cp, label=f"Simulation: {label}", linestyle='--')

# Finalize plot
plt.xlabel("Angle θ (degrees, aligned from stagnation point)")
plt.ylabel("Mean pressure coefficient ⟨Cp⟩")
plt.title("Pressure coefficient comparison across simulations")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
