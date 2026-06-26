// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_insertion_volume_h
#define lethe_insertion_volume_h

#include <dem/insertion.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
class InsertionVolume : public Insertion<dim, PropertiesIndex>
{
public:
  /**
   * @brief The constructor investigates if the insertion box is large enough
   * to handle the desired number of particles with the specified insertion
   * parameters.
   * If the insertion box is not adequately large, the number of
   * inserted particles at each insertion step is updated. It also finds the
   * insertion points in each direction (number_of_particles_x_direction,
   * number_of_particles_y_direction and number_of_particles_z_direction).
   *
   * @param size_distribution_object_container Contains all distribution for each
   * particle type
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   * @param maximum_particle_diameter Maximum particle diameter based on values
   * defined in the parameter handler
   */
  InsertionVolume(
    const std::vector<std::shared_ptr<Distribution>>
      &size_distribution_object_container,
    const parallel::distributed::Triangulation<dim> &triangulation,
    const DEMSolverParameters<dim>                  &dem_parameters,
    const double                                     maximum_particle_diameter);

  /**
   * @brief Carries out the volume insertion of particles.
   *
   * @param particle_handler The particle handler of particles which are being
   * inserted
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   *
   */
  virtual void
  insert(Particles::ParticleHandler<dim>                 &particle_handler,
         const parallel::distributed::Triangulation<dim> &triangulation,
         const DEMSolverParameters<dim> &dem_parameters) override;



  /**
   * @brief Serialize the volume insertion object to an output archive.
   *
   * @param ar Output archive where the attributes are stored.
   */
  virtual void
  serialize(boost::archive::text_oarchive &ar) override
  {
    ar &particles_of_each_type_remaining &current_inserting_particle_type;
  }

  /**
   * @brief Deserialize an input archive to the plane insertion object.
   *
   * @param ar Input archive where the attributes are stored.
   *
   */
  virtual void
  deserialize(boost::archive::text_iarchive &ar) override
  {
    ar &particles_of_each_type_remaining &current_inserting_particle_type;
  }

private:
  /**
   * @brief Carries out finding the maximum number of inserted particles based on the
   * insertion box size. If the requested number of particles for insertion in
   * each insertion step is larger than this maximum, it is limited to this
   * value and a warning is printed.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   * @param pcout Printing in parallel
   * @param inserted_this_step Number of particle to insert at the current
   * insertion step. This value is updated if the requested number of particles
   * for insertion in the current  insertion step is larger than the maximum n
   * umber of inserted particles based on the insertion box size.
   */
  void
  calculate_insertion_domain_maximum_particle_number(
    const DEMSolverParameters<dim> &dem_parameters,
    const ConditionalOStream       &pcout,
    unsigned int                   &inserted_this_step);

  /**
   * @brief Converts id of particles to volume insertion location
   *
   * @param insertion_location Insertion location of the particle
   * @param id Particle_id
   * @param random_number1 A random number to create randomness in volume insertion
   * @param random_number2 A random number to create randomness in volume insertion
   * @param insertion_information DEM insertion parameters declared in the .prm
   * file
   */
  void
  find_insertion_location(
    Point<dim>                                       &insertion_location,
    const unsigned int                                id,
    const double                                      random_number1,
    const double                                      random_number2,
    const Parameters::Lagrangian::InsertionInfo<dim> &insertion_information);

  /**
   * @brief Fills the filted_id_to_box_id map.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
s   */
  void
  set_filtered_id_map(const DEMSolverParameters<dim> &dem_parameters);

  unsigned int current_inserting_particle_type;

  // Number of particles of each type that remain to be inserted in the
  // upcoming insertion steps
  unsigned int particles_of_each_type_remaining;

  // Number of insertion points in the x, y and z directions
  std::vector<int> number_of_particles_directions;

  // Minimum and maximum boundaries of the insertion box in the direction order
  // It means that axis 0 is not necessarily x, since it depends on the order
  // of the insertion direction.
  std::vector<double> axis_min, axis_max;

  // Function used to insert particle inside our outside an arbitrary shape.
  std::shared_ptr<Function<dim>> acceptance_fct;

  // Map used to transfer from the real particle ID to the id of the particle
  // if the sign_function was accepting every insertion point.
  std::map<unsigned int, unsigned int> filted_id_to_box_id;
};
#endif
