import numpy as np
import pandas as pd
import scipy as sc

def position_string_to_int(coord_str):
    string_to_int = {'x': 0, 'y': 1, 'z': 2}
    return string_to_int[coord_str]

class PlateDischargeUtilities:
    def __init__(self, dataframe, prm_dict, time_list,  plate_norm, plate_direction):
        self.dataframe = dataframe
        self.prm_dict = prm_dict
        self.time_list = time_list
        self.direction = plate_direction#position_string_to_int(plate_direction)
        self.direction_str = plate_direction
        self.norm = plate_norm
        self.norm_str = position_string_to_int(plate_norm)
        p0 = prm_dict['initial translation'].split(', ')
        self.H = float(p0[self.norm_str])
        self.plate_extremities = prm_dict[plate_direction]
        self.L = self.plate_extremities[1] - self.plate_extremities[0]

        self.particles_x_top = [pd.Series()] * len(time_list)
        self.particles_y_top = [pd.Series()] * len(time_list)
        self.particles_x_bottom = [pd.Series()] * len(time_list)
        self.particles_y_bottom = [pd.Series()] * len(time_list)

        # Storing post-processing data
        postprocessing_data_type = ['time',
                                    'number_of_particles', 'mass', 'mass_flow_rate',
                                    'left_angle_top', '_left_angle_top_rsqrt', 'left_angle_bottom', 'left_angle_bottom_rsqrt',
                                    'right_angle_top', '_right_angle_top_rsqrt', 'right_angle_bottom', 'right_angle_bottom_rsqrt',
                                    'angle_top', 'angle_bottom', 'angle_top_rsqrt', 'angle_bottom_rsqrt']
        self.postprocessing_dataframe = pd.DataFrame(columns=postprocessing_data_type)
        self.postprocessing_dataframe['time'] = time_list

        # Useful variables for parameter file
        self.compute_particle_info()

    def compute_particle_info(self):
        self.dp = self.prm_dict['diameter']
        self.rhop = self.prm_dict['density particles']
        self.vp = 4. / 3. * np.pi * (self.dp / 2.) ** 3

    def get_data(self):
        return self.postprocessing_dataframe

        print("Solid mass flow rate computing done.")
    def calculate_solid_mass_flow_rate(self):
        # Virtual wall to limit the domain in the pipe for evaluation
        wall = self.H

        for i, df in enumerate(self.dataframe):
            # Get the IDs and the positions of the particles
            positions = df[self.norm].values

            # Number of particles below and mass discharge
            number_of_particles = (positions < wall).sum()
            mass = number_of_particles * self.vp * self.rhop

            self.postprocessing_dataframe.loc[i, 'number_of_particles'] = number_of_particles
            self.postprocessing_dataframe.loc[i, 'mass'] = mass

            self.postprocessing_dataframe.loc[i, 'mass_flow_rate'] = (mass / self.time_list[i])

        print("Solid mass flow rate computing done.")

    def calculate_angle_of_repose(self, x_min, x_max, f=1,
                      ignore_particles=False,
                      top=True, start_time=0.0):

        def get_bands(x_min, x_max, size):
            nx = round((x_max - x_min) / size)
            bounds = np.linspace(x_min, x_max, nx + 1)
            bands = [(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]
            return bands

        for i, df in enumerate(self.dataframe):
            t = self.time_list[i]
            if t < start_time:
                if top:
                    if x_min > 0.0:
                        self.postprocessing_dataframe.loc[i, 'right_angle_top'] = 0.0
                        self.postprocessing_dataframe.loc[i, 'right_angle_top_rsqrt'] = 0.0
                    else:
                        self.postprocessing_dataframe.loc[i, 'left_angle_top'] = 0.0
                        self.postprocessing_dataframe.loc[i, 'left_angle_top_rsqrt'] = 0.0
                else:
                    if x_min > 0.0:
                        self.postprocessing_dataframe.loc[i, 'right_angle_bottom'] = 0.0
                        self.postprocessing_dataframe.loc[i, 'right_angle_bottom_rsqrt'] = 0.0
                    else:
                        self.postprocessing_dataframe.loc[i, 'left_angle_bottom'] = 0.0
                        self.postprocessing_dataframe.loc[i, 'left_angle_bottom_rsqrt'] = 0.0
                continue


            bands = get_bands(x_min, x_max, f * self.dp)

            # Get the particles in each band
            particles_in_band = []
            highest_particles = pd.DataFrame(columns=['x', 'y', 'z'])

            # Get the IDs and the positions of the particles
            ids = df.index
            positions_x = df[self.direction].values
            positions_y = df[self.norm].values
            if top:
                particle_list = df[(positions_x > x_min) &
                                    (positions_x < x_max) &
                                    (positions_y > self.H)]
            else:
                particle_list = df[(positions_x > x_min) &
                                    (positions_x < x_max) &
                                   (positions_y < self.H)]

            for band in bands:
                # Get the bool vector of particles in the band and extract them
                bool_vect = (band[0] <= particle_list[self.direction_str]) & (particle_list[self.direction_str] < band[1])
                if any(bool_vect):
                    particles_in_band.append(particle_list[bool_vect])

            # Get the highest particle in axis for each band
            if ignore_particles is False:
                max_particles = [particles[particles[self.norm] == particles[self.norm].max()] for particles in particles_in_band]
                if len(max_particles) > 0:
                    highest_particles = pd.concat(max_particles, axis=0, ignore_index=True)
                else:
                    continue
            else:
                max_particles = []
                for particles in particles_in_band:
                    particles = particles.sort_values(by=[f'{self.norm}'], ascending=False)

                    if len(particles) > 1:
                        if (particles.iloc[0][self.norm_str] - particles.iloc[1][self.norm_str]) < 1.01 * self.dp:
                            max_particles.append(particles.iloc[[0]])
                        else:
                            max_particles.append(particles.iloc[[1]])
                    else:
                        if len(particles) == 0:
                            continue
                        else:
                            max_particles.append(particles.iloc[[0]])

                    highest_particles = pd.concat(max_particles, axis=0, ignore_index=True)

            # Get the coordinates of the highest particles
            particle_x = highest_particles[self.direction_str]
            particle_y = highest_particles[self.norm]

            self.particle_x = particle_x
            self.particle_y = particle_y

            a, b, r_value = self.linear_regression(particle_x, particle_y)
            angle = np.arctan(a) * 180 / np.pi
            y_model = a * particle_x + b
            self.model = y_model
            R_sq = r_value ** 2

            if top:
                self.particles_y_top[i] = pd.concat([self.particles_y_top[i], particle_y])
                self.particles_x_top[i] = pd.concat([self.particles_x_top[i], particle_x])
                if x_min > 0.0:
                    self.postprocessing_dataframe.loc[i, 'right_angle_top'] = angle
                    self.postprocessing_dataframe.loc[i, 'right_angle_top_rsqrt'] = R_sq
                else:
                    self.postprocessing_dataframe.loc[i, 'left_angle_top'] = angle
                    self.postprocessing_dataframe.loc[i, 'left_angle_top_rsqrt'] = R_sq
            else:
                self.particles_y_bottom[i] = pd.concat([self.particles_y_bottom[i], particle_y])
                self.particles_x_bottom[i] = pd.concat([self.particles_x_bottom[i], particle_x])
                if x_min > 0.0:
                    self.postprocessing_dataframe.loc[i, 'right_angle_bottom'] = angle
                    self.postprocessing_dataframe.loc[i, 'right_angle_bottom_rsqrt'] = R_sq
                else:
                    self.postprocessing_dataframe.loc[i, 'left_angle_bottom'] = angle
                    self.postprocessing_dataframe.loc[i, 'left_angle_bottom_rsqrt'] = R_sq

        print("Angle of repose computing done.")


    def calculate_angle_from_symmetry(self, top, start):
        if top:
            for i in range(len(self.dataframe)):
                t = self.time_list[i]
                if t < start:
                    self.postprocessing_dataframe.loc[i, 'angle_top'] = 0.0
                    self.postprocessing_dataframe.loc[i, 'angle_top_rsqrt'] = 0.0
                    continue

                particle_x = abs(self.particles_x_top[i])
                particle_y = self.particles_y_top[i]

                a, b, r_value = self.linear_regression(particle_x, particle_y)
                angle = np.arctan(a) * 180 / np.pi
                R_sq = r_value ** 2

                self.postprocessing_dataframe.loc[i, 'angle_top'] = angle
                self.postprocessing_dataframe.loc[i, 'angle_top_rsqrt'] = R_sq
        else:
            for i in range(len(self.dataframe)):
                t = self.time_list[i]
                if t < start:
                    self.postprocessing_dataframe.loc[i, 'angle_bottom'] = 0.0
                    self.postprocessing_dataframe.loc[i, 'angle_bottom_rsqrt'] = 0.0
                    continue

                particle_x = abs(self.particles_x_bottom[i])
                particle_y = self.particles_y_bottom[i]

                a, b, r_value = self.linear_regression(particle_x, particle_y)
                angle = np.arctan(a) * 180 / np.pi
                R_sq = r_value ** 2

                self.postprocessing_dataframe.loc[i, 'angle_bottom'] = angle
                self.postprocessing_dataframe.loc[i, 'angle_bottom_rsqrt'] = R_sq


    def linear_regression(self, x, y):
        slope, intercept, r_value, p_value, std_err = sc.stats.linregress(x, y)
        return slope, intercept, r_value


    def theoretical_angle_of_repose(self):
        dp = self.dp * 1000 # mm
        mu_f_pp = self.prm_dict['friction coefficient particles']
        mu_f_pw = self.prm_dict['friction coefficient wall']
        mu_r_pp = self.prm_dict['rolling friction particles'] * dp
        mu_r_pw = self.prm_dict['rolling friction wall'] * dp

        self.theoritical_angle = (68.61 * mu_f_pp ** 0.27 * mu_f_pw ** 0.22 *
                                  mu_r_pp ** 0.06 * mu_r_pw ** 0.12 * dp ** -0.2)
        print(self.theoritical_angle)

        print(f"Theoretical angle of repose computing done.")
