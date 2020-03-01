#include <dem/nonuniform_insertion.h>

using namespace DEM;

template <int dim>
NonUniformInsertion<dim>::NonUniformInsertion(
    const DEMSolverParameters<dim> &dem_parameters) {
  // Getting properties as local parameters
  const auto physical_properties = dem_parameters.physicalProperties;
  const auto insertion_information = dem_parameters.insertionInfo;

  // This variable is used for calculation of the maximum number of particles
  // that can fit in the chosen insertion box
  int maximum_particle_number;

  // distance_threshold shows the ratio of the distance between the centers of
  // two adjacent particles to the diameter of particles
  maximum_particle_number =
      int((insertion_information.x_max - insertion_information.x_min) /
          (insertion_information.distance_threshold *
           physical_properties.diameter)) *
      int((insertion_information.y_max - insertion_information.y_min) /
          (insertion_information.distance_threshold *
           physical_properties.diameter)) *
      int((insertion_information.z_max - insertion_information.z_min) /
          (insertion_information.distance_threshold *
           physical_properties.diameter));

  // If the inserted number of particles at this step exceeds the maximum
  // number, a warning is printed
  if (insertion_information.inserted_this_step > maximum_particle_number) {
    std::cout << "The inserted number of particles ("
              << insertion_information.inserted_this_step
              << ") is higher than maximum expected number of particles ("
              << maximum_particle_number << ")" << std::endl;
    std::cout << "Inserting " << maximum_particle_number << " at this step"
              << std::endl;
  }
} // add error here

template <int dim>
void NonUniformInsertion<dim>::insert(
    Particles::ParticleHandler<dim> &particle_handler,
    const Triangulation<dim> &tr, Particles::PropertyPool &pool,
    const DEMSolverParameters<dim> &dem_parameters) {
  // Getting properties as local parameters
  const auto physical_properties = dem_parameters.physicalProperties;
  const auto insertion_information = dem_parameters.insertionInfo;

  // nx, ny and nz are the results of discretization of the insertion domain in
  // x, y and z directions
  int nx = int((insertion_information.x_max - insertion_information.x_min) /
               (insertion_information.distance_threshold *
                physical_properties.diameter));
  int ny = int((insertion_information.y_max - insertion_information.y_min) /
               (insertion_information.distance_threshold *
                physical_properties.diameter));
  int nz = int((insertion_information.z_max - insertion_information.z_min) /
               (insertion_information.distance_threshold *
                physical_properties.diameter));

  // inserted_sofar_step shows the number of inserted particles, so far, at
  // this step
  int inserted_sofar_step = 0;

  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j) {
      // Adapt the last index to the dimensionality of the problem
      int dim_nz = (dim == 3) ? nz : 1;
      for (int k = 0; k < dim_nz; ++k)

        // We need to check if the number of inserted particles so far at this
        // step reached the total desired number of inserted particles at this
        // step
        if (inserted_sofar_step < insertion_information.inserted_this_step) {
          Point<dim> position;
          Point<dim> reference_position;
          unsigned int id;

          // Obtaning position of the inserted particle
          // In non-uniform insertion, two random numbers are created and
          // added to the position of particles

          // This numbers 101 and 400 are hard-coded, I will fix them when
          // adding new paramters to the parameter handler file
          int randNum1 = rand() % 300;
          int randNum2 = rand() % 300;
          position[0] = insertion_information.x_min +
                        (physical_properties.diameter / 2) +
                        (k * insertion_information.distance_threshold *
                         physical_properties.diameter) +
                        randNum1 * (physical_properties.diameter / 400.0);
          position[1] = insertion_information.y_min +
                        (physical_properties.diameter / 2) +
                        (j * insertion_information.distance_threshold *
                         physical_properties.diameter) +
                        randNum2 * (physical_properties.diameter / 400.0);
          if (dim == 3)
            position[2] = insertion_information.z_min +
                          (physical_properties.diameter / 2) +
                          (i * insertion_information.distance_threshold *
                           physical_properties.diameter);

          // Since the id of each particle should be unique, we need to use
          // the total number of particles in the system to calculate the
          // ids of new particles
          if (dim == 3)
            id = i * ny * nz + j * nz + k +
                 particle_handler.n_global_particles();
          if (dim == 2)
            id = i * ny + j + particle_handler.n_global_particles();

          // Inserting the new particle using its location, id and
          // containing cell
          Particles::Particle<dim> particle(position, reference_position, id);
          typename Triangulation<dim>::active_cell_iterator cell =
              GridTools::find_active_cell_around_point(tr,
                                                       particle.get_location());
          Particles::ParticleIterator<dim> pit =
              particle_handler.insert_particle(particle, cell);

          // Setting property pool of inserted particle
          particle.set_property_pool(pool);

          // Initialization of the properties of the new particle
          pit->get_properties()[PropertiesIndex::id] = id;
          pit->get_properties()[PropertiesIndex::type] = 1;
          pit->get_properties()[PropertiesIndex::dp] =
              physical_properties.diameter;
          pit->get_properties()[PropertiesIndex::rho] =
              physical_properties.density;
          // Velocity
          pit->get_properties()[PropertiesIndex::v_x] = 0;
          pit->get_properties()[PropertiesIndex::v_y] = 0;
          if (dim == 3)
            pit->get_properties()[PropertiesIndex::v_z] = 0;
          // Acceleration
          pit->get_properties()[PropertiesIndex::acc_x] = 0;
          pit->get_properties()[PropertiesIndex::acc_y] = 0;
          if (dim == 3)
            pit->get_properties()[PropertiesIndex::acc_z] = 0;
          // Force
          pit->get_properties()[PropertiesIndex::force_x] = 0;
          pit->get_properties()[PropertiesIndex::force_y] = 0;
          if (dim == 3)
            pit->get_properties()[PropertiesIndex::force_z] = 0;
          // w
          pit->get_properties()[PropertiesIndex::omega_x] = 0;
          pit->get_properties()[PropertiesIndex::omega_y] = 0;
          if (dim == 3)
            pit->get_properties()[PropertiesIndex::omega_z] = 0;
          // mass and moi
          pit->get_properties()[PropertiesIndex::mass] =
              physical_properties.density *
              ((4.0 / 3.0) * 3.1415 *
               pow((pit->get_properties()[PropertiesIndex::dp] / 2.0), 3.0));

          pit->get_properties()[PropertiesIndex::mom_inertia] =
              (2.0 / 5.0) * (pit->get_properties()[PropertiesIndex::mass]) *
              pow((pit->get_properties()[PropertiesIndex::dp] / 2.0), 2.0);
          // Torque
          pit->get_properties()[PropertiesIndex::M_x] = 0;
          pit->get_properties()[PropertiesIndex::M_y] = 0;
          if (dim == 3)
            pit->get_properties()[PropertiesIndex::M_z] = 0;

          ++inserted_sofar_step;
        }
    }
}

template class NonUniformInsertion<2>;
template class NonUniformInsertion<3>;
