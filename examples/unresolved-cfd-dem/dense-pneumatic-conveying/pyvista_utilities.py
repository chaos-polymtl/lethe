import numpy as np
import pandas as pd

def position_int_to_string(coord_int):
    int_to_string = {0: 'x', 1: 'y', 2: 'z'}
    return int_to_string[coord_int]
def velocity_int_to_string(coord_int):
    int_to_string = {0: 'u', 1: 'v', 2: 'w'}
    return int_to_string[coord_int]

class PneumaticConveyingUtilities:
    def __init__(self, dataframe, prm_dict, time_list, L, vol_triangulation):
        self.dataframe = dataframe
        self.prm_dict = prm_dict
        self.time_list = time_list
        self.L = L
        self.vol_triangulation = vol_triangulation

        # Storing post-processing data
        postprocessing_data_type = ['time', 'front_position', 'rear_position', 'distance', 'number_of_particles',
                                    'mass', 'mass_flow_rate', 'average_velocity', 'slug_velocity']
        self.postprocessing_dataframe = pd.DataFrame(columns=postprocessing_data_type)
        self.postprocessing_dataframe['time'] = time_list

        # Useful variables for parameter file
        self.compute_info()

    def compute_info(self):
        self.dp = self.prm_dict['diameter']
        self.rhop = self.prm_dict['density particles']
        self.vp = 4. / 3. * np.pi * (self.dp / 2.) ** 3
        self.direction = position_int_to_string(self.prm_dict['periodic direction'])


    def get_data(self):
        return self.postprocessing_dataframe

    def calculate_slug_length(self, void_fraction_threshold=0.55):
        # Area of the cross-section of the pipe
        A = self.vol_triangulation/self.L

        # Get the number of slices and their length
        n_L = int(self.L / self.dp)
        dL = self.L / n_L

        for i, df in enumerate(self.dataframe):
            # Create a list to store the number of particles and the position limits of slices
            number_of_particles = np.zeros(n_L)
            position_limits = np.linspace(-self.L/2.0, self.L/2.0, n_L+1)

            # Get the number of particles in each slice of dL
            x_positions = df[self.direction]
            for j in range(0, n_L):
                # Get the number of particles per slice
                number_of_particles[j] = np.array([(x_positions >= position_limits[j]) &
                                                        (x_positions < position_limits[j+1])], dtype=bool).sum()

            # Get the void fraction of the slice
            void_fractions = 1 - number_of_particles * self.vp / (dL * A)

            # Smoothing the void fraction with a convolution operator (moving average)
            period = 5
            void_fractions_smoothed = np.convolve(np.ones(period)/period, void_fractions, mode='same')

            # Since the convolution operator do not take the periodic boundary into account, we add the contribution of
            # the slide at boundaries
            void_fractions_smoothed[0] += 1/period * void_fractions[-2:].sum()
            void_fractions_smoothed[1] += 1/period * void_fractions[-1]
            void_fractions_smoothed[-2] += 1/period * void_fractions[0]
            void_fractions_smoothed[-1] += 1/period * void_fractions[0:2].sum()

            # Find the state of the slide with the void fraction threshold (no slug = 0, slug = 1)
            slug_state = np.where(void_fractions_smoothed <= void_fraction_threshold, 1, 0)

            correction_term = 0.0

            # If the first and last slug state are 1, it means slug in the periodic boundaries
            if slug_state[0] == 1 and slug_state[-1] == 1:
                # Slug in PBC:     ‾‾‾‾\_________/‾‾‾‾
                # Discretized:     ‾|‾|\|_|_|_|_|/|‾|‾
                #                  1 1 0 0 0 0 0 0 1 1

                # Find the indices where the slug is not and get the front and rear indices
                no_slug_index = np.where(slug_state == 0)[0]
                front_index = no_slug_index[0] - 1
                rear_index =  no_slug_index[-1] + 1

                # Since the slug is in the periodic boundaries, the distance is corrected by the length of the pipe
                correction_term = self.L
            else:
                # Slug in middle:  ____/‾‾‾‾‾‾‾‾‾\____
                # Discretized:     _|_|/|‾|‾|‾|‾|\|_|_
                #                  0 0 0 1 1 1 1 0 0 0

                # Find the indices where the slug is and get the front and rear indices
                slug_index = np.where(slug_state == 1)[0]
                front_index = slug_index[-1]
                rear_index = slug_index[0]

            # Get the front and rear position of the slug
            # Note: the length tends to be overestimated, so we take the position limits towards the center
            front_position = position_limits[front_index]
            rear_position = position_limits[rear_index+1]

            # Compute the distance between the front and rear of the slug at middle of slice
            distance = front_position - rear_position + correction_term

            # Store information in the post-processing dataframe
            self.postprocessing_dataframe.loc[i, 'front_position'] = front_position
            self.postprocessing_dataframe.loc[i, 'rear_position'] = rear_position
            self.postprocessing_dataframe.loc[i, 'distance'] = distance

        print("Slug length computing done.")


    def calculate_solid_mass_flow_rate(self):
        # Virtual wall to limit the domain in the pipe for evaluation
        dL = 25 * self.dp
        wall = self.L / 2. - dL

        # Create a list to store the number of particles
        last_list = np.array([])

        for i, df in enumerate(self.dataframe):
            # Get the IDs and the positions of the particles
            ids = df['ID'].values
            positions = df[self.direction].values

            # Find the particles that crossed the wall
            particle_list = ids[positions > wall]

            # Remove the particles that already crossed the wall at last time step from the list
            gone_particles = np.setdiff1d(particle_list, last_list)
            last_list = particle_list

            # Number of particles below and mass discharge
            number_of_particles = len(gone_particles)
            mass = number_of_particles * self.vp * self.rhop

            self.postprocessing_dataframe.loc[i, 'number_of_particles'] = number_of_particles
            self.postprocessing_dataframe.loc[i, 'mass'] = mass

            if i == 0:
                self.postprocessing_dataframe.loc[i, 'mass_flow_rate'] = 0
            else:
                mass_list = self.postprocessing_dataframe['mass'].values
                self.postprocessing_dataframe.loc[i, 'mass_flow_rate'] = (mass_list[0:i].sum() /
                                                                          (self.time_list[i] - self.time_list[0]))

        print("Solid mass flow rate computing done.")

    def calculate_average_particle_velocity(self):
        for i, df in enumerate(self.dataframe):
            # ids = current_df['ID'].astype(int)
            # particle_id = np.where(ids == int(id))[0][0]
            velocities = df[['u', 'v', 'w']].values
            positions = df[self.direction]
            front_position = self.postprocessing_dataframe['front_position']
            rear_position = self.postprocessing_dataframe['rear_position']

            avg_velocities = np.zeros(3)

            # Get the particles index in the slug
            if front_position[i] < rear_position[i]:
                in_slug_index = np.where((positions >= rear_position[i]) | (positions < front_position[i]))[0]
            else:
                in_slug_index = np.where((positions >= rear_position[i]) & (positions < front_position[i]))[0]

            # Compute the average velocity of the particles in the slug
            for j in in_slug_index:
                avg_velocities += velocities[j, :]

            average_velocity = avg_velocities / len(in_slug_index)
            average_velocity = np.sqrt(average_velocity[0]**2 + average_velocity[1]**2 + average_velocity[2]**2)
            #print("Slug's particle velocity computing done.")

            self.postprocessing_dataframe.loc[i, 'average_velocity'] = average_velocity

        print("Particle velocity in slug computing done.")

    def calculate_slug_velocity(self, sample_step=1):
        # Periodic position handling
        front_positions = self.postprocessing_dataframe['front_position'].values
        continuous_front_positions = front_positions.copy()
        dt = self.time_list[1] - self.time_list[0]

        # Transform the slug front position as if the pipe was continuous, it helps for the numerical differentiation
        i_pass = 0
        for i in range(1, len(front_positions)):
            if front_positions[i] < front_positions[i-1]:
                i_pass += 1

            continuous_front_positions[i] += i_pass * self.L

        # First value is computed with forward difference
        velocity = (continuous_front_positions[1] - continuous_front_positions[0]) / dt
        self.postprocessing_dataframe.loc[0, 'slug_velocity'] = velocity

        # Central difference for the rest of the values
        for i in list(range(1, sample_step)) + list(range(len(continuous_front_positions)-sample_step, len(continuous_front_positions) - 1)):
            velocity = (continuous_front_positions[i+1] - continuous_front_positions[i-1]) / (2 * dt)
            self.postprocessing_dataframe.loc[i, 'slug_velocity'] = velocity

        # Central difference with skipping data in the sample step because there's a lot of fluctuation
        for i in range(sample_step, len(continuous_front_positions) - sample_step):
            velocity = (continuous_front_positions[i+sample_step] - continuous_front_positions[i-sample_step]) / (2 * sample_step * dt)
            self.postprocessing_dataframe.loc[i, 'slug_velocity'] = velocity

        # Backward difference for the last values
        velocity = (continuous_front_positions[-1] - continuous_front_positions[-2]) / dt
        self.postprocessing_dataframe.loc[-1, 'slug_velocity'] = velocity

        print("Slug velocity computing done.")