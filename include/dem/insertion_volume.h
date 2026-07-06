// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
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
   * @param[in] size_distribution_object_container Contains all distribution for
   * each particle type
   * @param[in] triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param[in] dem_parameters DEM parameters declared in the .prm file
   * @param[in] maximum_particle_diameter Maximum particle diameter based on
   * values defined in the parameter handler
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
   * @param[in] particle_handler The particle handler of particles which are
   * being inserted
   * @param[in] triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param[in] dem_parameters DEM parameters declared in the .prm file
   *
   */
  virtual void
  insert(Particles::ParticleHandler<dim>                 &particle_handler,
         const parallel::distributed::Triangulation<dim> &triangulation,
         const DEMSolverParameters<dim> &dem_parameters) override;



  /**
   * @brief Serialize the volume insertion object to an output archive.
   *
   * @param[out] ar Output archive where the attributes are stored.
   */
  virtual void
  serialize(boost::archive::text_oarchive &ar) override
  {
    ar &particles_of_each_type_remaining &current_inserting_particle_type;
  }

  /**
   * @brief Deserialize an input archive to the plane insertion object.
   *
   * @param[in] ar Input archive where the attributes are stored.
   *
   */
  virtual void
  deserialize(boost::archive::text_iarchive &ar) override
  {
    ar &particles_of_each_type_remaining &current_inserting_particle_type;
  }

private:
  /**
   * @brief Adjust the number of inserted particles for this insertion time step
   * if the number required by the used is larger than what is allowed by the
   * insertion box and the acceptance function.
   *
   * @param[in] insertion_information DEM insertion parameters declared in the
   * .prm
   * @param[in] pcout Printing in parallel
   * @param[in] inserted_this_step Number of particle to insert at the current
   * insertion step. This value is updated if the requested number of particles
   * for insertion in the current  insertion step is larger than the maximum
   * number of inserted particles based on the insertion box size.
   */
  void
  adjust_insertion_count_by_insertion_box_capacity(
    const InsertionInfo<dim> &insertion_information,
    const ConditionalOStream &pcout,
    unsigned int             &inserted_this_step);

  /**
   * @brief Converts id of particles to volume insertion location
   *
   * @param[out] insertion_location Insertion location of the particle
   * @param[in] id Particle_id
   * @param[in] random_number1 A random number to create randomness in volume
   * insertion
   * @param[in] random_number2 A random number to create randomness in volume
   * insertion
   * @param[in] insertion_information DEM insertion parameters declared in the
   * .prm file
   */
  void
  find_insertion_location(Point<dim>               &insertion_location,
                          const unsigned int        id,
                          const double              random_number1,
                          const double              random_number2,
                          const InsertionInfo<dim> &insertion_information);

  /**
   * @brief Fills the filted_box_index vector. Also finds how many insertion
   * points are valid inside the insertion box after applying the acceptance
   * function.
   *
   * @param[in] insertion_information DEM insertion parameters declared in the
   * .prm
   */
  void
  set_filtered_index(const InsertionInfo<dim> &insertion_information);

  unsigned int current_inserting_particle_type;

  /// Number of particles of each type that remain to be inserted in the
  /// upcoming insertion steps
  unsigned int particles_of_each_type_remaining;

  // Number of insertion points in the x, y and z directions
  std::vector<int> number_of_particles_directions;

  /// Minimum and maximum boundaries of the insertion box in the direction order
  /// It means that axis 0 is not necessarily x, since it depends on the order
  /// of the insertion direction.
  std::vector<double> axis_min, axis_max;

  /// Function used to insert particle inside our outside an arbitrary shape.
  std::shared_ptr<Function<dim>> acceptance_fct;

  /// Map used to transfer from the filtered to the pre-filted numbering of the
  /// insertion point inside the insertion box.
  std::vector<unsigned int> filted_box_index;

  /// The first and last index that will get inserted by this process following
  /// the filtered numeration.
  unsigned int first_index_this_proc, last_index_this_proc;

  /// The total number of valid points inside the box on all processes after
  /// applying the acceptance function.
  unsigned int number_of_valid_insertion_point_global;

  /// Maximum particle diameter based on the values defined in the parameter
  const double maximum_diameter;
};
#endif
