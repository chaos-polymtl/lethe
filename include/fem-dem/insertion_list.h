// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_insertion_list_h
#define lethe_insertion_list_h

#include <dem/insertion.h>

template <int dim, typename PropertiesIndex>
class InsertionList : public Insertion<dim, PropertiesIndex>
{
public:
  /**
   * @brief Insert particles using a list specific positions. This allows the
   * insertion of any number of particles at a well-controlled location which is
   * especially useful from a testing perspective. The code ensures that the
   * number of positions provided in the x,y (and possibly z) direction is
   * coherent. If more particles than the number of positions in the list is
   * requested, the class will continue inserting particles at the insertion
   * frequency using the list of positions. There is no mechanism in place that
   * prevents the overlap of these new particles with previous ones.
   *
   * @param size_distribution_object_container Contains all distribution for each
   * particle type
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  InsertionList(const std::vector<std::shared_ptr<Distribution>>
                  &size_distribution_object_container,
                const parallel::distributed::Triangulation<dim> &triangulation,
                const DEMSolverParameters<dim> &dem_parameters);

  /**
   * @brief The InsertionList class inserts particles using a list specific position.
   * This allows the insertion of any number of particles at a well-controled
   * location which is especially useful from a testing perspective. The code
   * ensures that the number of positions provided in the x,y (and possibly z)
   * direction is coherent. If more particles than the number of position in the
   * list are requested, the class will continue inserting particles at the
   * insertion frequency using the list of position. There is no mechanism in
   * place that prevents the overlap of these new particles with previous ones.
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
   * @brief Carries out assigning the properties of inserted particles specificly
   * for the list insertion method. In this method, the initial translationnal
   * and angular velocities, the temperature and the diameter of each particles
   * is set.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   * @param inserted_this_step_this_proc Number of particles that are inserted
   * at each insertion step on each processor. This value can change in the last
   * insertion step to reach the desired number of particles
   * @param current_inserting_particle_type Type of inserting particles
   * @param particle_properties Properties of all inserted particles at each insertion step
   */
  void
  assign_particle_properties_for_list_insertion(
    const DEMSolverParameters<dim>   &dem_parameters,
    const unsigned int               &inserted_this_step_this_proc,
    const unsigned int               &current_inserting_particle_type,
    std::vector<std::vector<double>> &particle_properties);



  /**
   * @brief Serialize the list insertion object to an output archive.
   *
   * @param ar Output archive where the attributes are stored.
   */
  virtual void
  serialize(boost::archive::text_oarchive &ar) override
  {
    ar &remaining_particles_of_each_type &current_inserting_particle_type;
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
    ar &remaining_particles_of_each_type &current_inserting_particle_type;
  }

  // Number of remaining particles of each type that should be inserted in the
  // upcoming insertion steps
  unsigned int remaining_particles_of_each_type;
  unsigned int current_inserting_particle_type;

  // Vector of the location, velocity and angular velocity, diameter and
  // temperature used to insert the particles
  std::vector<Point<dim>>   insertion_points;
  std::vector<Tensor<1, 3>> velocities;
  std::vector<Tensor<1, 3>> angular_velocities;
  std::vector<double>       diameters;
  std::vector<double>       temperatures;
};

#endif /* lethe_insertion_list_h */
