from lethe_pyvista_tools import *
import matplotlib.pyplot as plt
from pyvista_utilities import *
from log_utilities import *

#############################
# Path and filenames
#############################
names = ["base", "asc", "lb", "asc-lb"]
dict_names = {"base": "Baseline", "asc": "Adaptive Sparse Contacts", "lb": "Load Balancing", "asc-lb": "ASC-LB"}
prm_filename = [f"discharge-plate_{name}.prm" for name in names]
pvd_name = "out"
log_filename = [f"log2_{name}.out" for name in names]

#############################
# Read data and process
#############################
# In this post-processing, we need the ID of the particles, the position, and the velocity, the other data is not
# necessary and takes a lot of memory
ignore_data = ['omega', 'type', 'diameter', 'mass', 'fem_force', 'fem_torque']

# Data we want to process
process_data = False
process_log = True

# All plots in one figure
fig_mass = plt.figure(figsize=(6, 5))
fig_angle = plt.figure(figsize=(10, 5))
fig_performance = plt.figure(figsize=(6, 5))

# Average fluid velocity and beta force
mass_plot = fig_mass.add_subplot(111)
angle_top_plot = fig_angle.add_subplot(121)
angle_bottom_plot = fig_angle.add_subplot(122)
performance_plot = fig_performance.add_subplot(111)

color = 0

for name in names:
# Create the particle object
    if process_data:
        particles = lethe_pyvista_tools(f"data", f"plate-discharge_{name}.prm", f"{pvd_name}.pvd",
                                        ignore_data=ignore_data, step=5)
        prm_dict = particles.prm_dict
        dataframe = []

        # Some plate informations
        plate_norm = 'y'
        plate_direction = 'x'

        for i in range(len(particles.list_vtu)):
            df = particles.get_df(i)
            data = pd.DataFrame()
            # Get the ID of the particles
            data['ID'] = df['ID'].astype(int)

            # Get the position of the particles
            data['x'] = df.points[:, 0]
            data['y'] = df.points[:, 1]
            data['z'] = df.points[:, 2]
            data.set_index('ID', inplace=True)

            dataframe.append(data)

        # Create the post-processing and log object
        postprocessing = PlateDischargeUtilities(dataframe, prm_dict, particles.time_list,
                                                 plate_norm, plate_direction)

        # Solid mass flow rate
        postprocessing.calculate_solid_mass_flow_rate()

        # Process the top left and right angle of repose
        postprocessing.calculate_angle_of_repose(top=True, start_time=5, x_min=-0.35, x_max=-0.15)
        postprocessing.calculate_angle_of_repose(top=True, start_time=5, x_min=0.15, x_max=0.35)
        postprocessing.calculate_angle_from_symmetry(top=True, start=5)

        # Process the bottom left and right angle of repose
        postprocessing.calculate_angle_of_repose(top=False, start_time=5, x_min=-0.35, x_max=-0.15)
        postprocessing.calculate_angle_of_repose(top=False, start_time=5, x_min=0.15, x_max=0.35)
        postprocessing.calculate_angle_from_symmetry(top=False, start=5)


        # Get the dataframes
        postprocessing_dataframe = postprocessing.get_data()


        mass_plot.plot(postprocessing_dataframe['time'], postprocessing_dataframe['mass'], color='C' + str(color),
                       label=f'{dict_names[name]}', linewidth=.75)

        angle_top_plot.plot(postprocessing_dataframe['time'], abs(postprocessing_dataframe['angle_top']),
                            color='C' + str(color), label=f'{dict_names[name]}', linewidth=0.75)
        angle_top_plot.fill_between(postprocessing_dataframe['time'].values,
                                    abs(pd.to_numeric(postprocessing_dataframe['left_angle_top'].values, errors='coerce')),
                                    abs(pd.to_numeric(postprocessing_dataframe['right_angle_top'].values, errors='coerce')),
                                    color='C' + str(color), alpha=0.5)

        angle_bottom_plot.plot(postprocessing_dataframe['time'], abs(postprocessing_dataframe['angle_bottom']),
                               color='C' + str(color), label=f'{dict_names[name]}', linewidth=0.75)
        angle_bottom_plot.fill_between(postprocessing_dataframe['time'],
                                       abs(pd.to_numeric(postprocessing_dataframe['left_angle_bottom'].values,
                                                         errors='coerce')),
                                       abs(pd.to_numeric(postprocessing_dataframe['right_angle_bottom'].values,
                                                         errors='coerce')), color='C' + str(color), alpha=0.5)


    if process_log:
        log = LogUtilities(f"performance/log_{name}.out")
        log.dem_performance_monitoring()
        log_dataframe = log.get_log_data()


        performance_plot.plot(log_dataframe['time'], log_dataframe['dem_walltime'], color='C' + str(color), label=f'{dict_names[name]}', linewidth=0.75)

    color += 1

if process_data:
    start = 10
    end = 15
    mass_plot.set(ylim=[0, 45], xlim=[0, end])
    mass_plot.set_xlabel('Time (s)')
    mass_plot.set_ylabel('Solid mass (kg)')
    mass_plot.grid(color='k', linestyle='-', linewidth=0.4)
    #mass_plot.xaxis.set_major_locator(plt.MaxNLocator(11))
    #mass_plot.axvspan(0, start, color='gray', alpha=0.15, lw=0)
    mass_plot.legend(edgecolor="k", fancybox=0,
                        facecolor="white", framealpha=1).get_frame().set_linewidth(0.75)

    #angle_top_plot.set_title('Top left angle of repose')
    angle_top_plot.set(ylim=[16, 26], xlim=[start , end])
    angle_top_plot.set_xlabel('Time (s)')
    angle_top_plot.set_ylabel('Top angle of repose (°)')
    angle_top_plot.grid(color='k', linestyle='-', linewidth=0.4)
    angle_top_plot.xaxis.set_major_locator(plt.MaxNLocator(6))
    angle_top_plot.axvspan(0, start, color='gray', alpha=0.15, lw=0)
    angle_top_plot.legend(edgecolor="k", fancybox=0,
                        facecolor="white", framealpha=1).get_frame().set_linewidth(0.75)

    #angle_bottom_plot.set_title('Bottom left angle of repose')
    angle_bottom_plot.set(ylim=[16, 26], xlim=[start, end])
    angle_bottom_plot.set_xlabel('Time (s)')
    angle_bottom_plot.set_ylabel('Bottom angle of repose (°)')
    angle_bottom_plot.grid(color='k', linestyle='-', linewidth=0.4)
    angle_bottom_plot.xaxis.set_major_locator(plt.MaxNLocator(6))
    angle_bottom_plot.axvspan(0, start, color='gray', alpha=0.15, lw=0)
    angle_bottom_plot.legend(edgecolor="k", fancybox=0,
                        facecolor="white", framealpha=1, loc='lower right').get_frame().set_linewidth(0.75)

    fig_mass.savefig(f'./plate-discharge-data.png', bbox_inches='tight', dpi=300)
    fig_angle.savefig(f'./plate-discharge-angle-data.png', bbox_inches='tight', dpi=300)

if process_log:
    performance_plot.set(ylim=[0, 3000], xlim=[0, 15])
    performance_plot.set_xlabel('Time (s)')
    performance_plot.set_ylabel('Walltime (s)')
    performance_plot.grid(color='k', linestyle='-', linewidth=0.4)
    #performance_plot.xaxis.set_major_locator(plt.MaxNLocator(11))
    performance_plot.legend(edgecolor="k", fancybox=0,
                        facecolor="white", framealpha=1).get_frame().set_linewidth(0.75)

    fig_performance.savefig(f'./plate-discharge-performance-data.png', bbox_inches='tight', dpi=300)

plt.tight_layout()
plt.show()