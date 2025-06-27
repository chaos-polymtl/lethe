// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_insertion_file_h
#define lethe_insertion_file_h

#include <dem/insertion.h>

template <int dim, typename PropertiesIndex>
class InsertionFile : public Insertion<dim, PropertiesIndex>
{
public:
  /**
   * @brief The InsertionFile class inserts particles using data stored in a file.
   * This allows the insertion of any number of particles at a well-controlled
   * location with any diameter value, translation, angular velocity and
   * temperature.
   *
   * @param size_distribution_object_container Contains all distribution for each
   * particle type
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  InsertionFile(const std::vector<std::shared_ptr<Distribution>>
                  &size_distribution_object_container,
                const parallel::distributed::Triangulation<dim> &triangulation,
                const DEMSolverParameters<dim> &dem_parameters);

  /**
   * @brief The InsertionFile class inserts particles using data stored in a file.
   * This allows the insertion of any number of particles at a well-controlled
   * location with any diameter value, translation, angular velocity and
   * temperature.
   *
   * @param particle_handler The particle handler of particles which are being
   * inserted
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void
  insert(Particles::ParticleHandler<dim>                 &particle_handler,
         const parallel::distributed::Triangulation<dim> &triangulation,
         const DEMSolverParameters<dim> &dem_parameters) override;


  /**
   * @brief Carries out assigning the properties of inserted particles specifically
   * for the file insertion method. In this method, the initial translation
   * and angular velocities, the temperature and the diameter of each particles
   * is set.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   *
   * @param inserted_this_step_this_proc Number of particles that are inserted
   * at each insertion step on each processor. This value can change in the last
   * insertion step to reach the desired number of particles
   *
   * @param particles_data Contains the particles
   *
   * @param particle_properties Properties of all inserted particles at each insertion step
   */
  void
  assign_particle_properties_for_file_insertion(
    const DEMSolverParameters<dim>             &dem_parameters,
    const unsigned int                         &inserted_this_step_this_proc,
    std::map<std::string, std::vector<double>> &particles_data,
    std::vector<std::vector<double>>           &particle_properties);



  /**
   * @brief Serialize the list insertion object to an output archive.
   *
   * @param ar Output archive where the attributes are stored.
   */
  virtual void
  serialize(boost::archive::text_oarchive &ar) override
  {
    ar &remaining_particles_of_each_type &current_inserting_particle_type
      &current_file_id;
  }

  /**
   * @brief Deserialize an input archive to the list insertion object.
   *
   * @param ar Input archive where the attributes are stored.
   *
   */
  virtual void
  deserialize(boost::archive::text_iarchive &ar) override
  {
    ar &remaining_particles_of_each_type &current_inserting_particle_type
      &current_file_id;
  }

  // Number of remaining particles of each type that should be inserted in the
  // upcoming insertion steps
  unsigned int remaining_particles_of_each_type;
  unsigned int current_inserting_particle_type;

  // File id iterator
  unsigned int       current_file_id;
  const unsigned int number_of_files;

  // Files where the particles properties are stored.
  const std::vector<std::string> insertion_files;
};

#endif /* lethe_insertion_file_h */
