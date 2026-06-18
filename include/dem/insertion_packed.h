// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_insertion_packed_h
#define lethe_insertion_packed_h

#include <dem/insertion.h>

#include <random>

using namespace dealii;

template <int dim, typename PropertiesIndex>
class InsertionPacked : public Insertion<dim, PropertiesIndex>
{
public:
  /**
   * @brief Insert particles using the packed insertion method.
   *
   * Particles are inserted by randomly generating points inside the insertion
   * region and assigning them particle properties (such as their diameter).
   * Because the generated configuration may contain particle-particle and
   * particle-wall overlaps, particles are iteratively displaced along the
   * corresponding contact normal directions until all overlaps are resolved.
   *
   * @param size_distribution_object_container Contains all distribution for each
   * particle type
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  InsertionPacked(
    const std::vector<std::shared_ptr<Distribution>>
      &size_distribution_object_container,
    const parallel::distributed::Triangulation<dim> &triangulation,
    const DEMSolverParameters<dim>                  &dem_parameters);

  /**
   * @brief Carries out the packed insertion of particles.
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
   * @brief Serialize the packed insertion object to an output archive.
   *
   * @param ar Output archive where the attributes are stored.
   */
  virtual void
  serialize(boost::archive::text_oarchive &ar) override
  {
    ar &particles_of_each_type_remaining &current_inserting_particle_type;
  }

  /**
   * @brief Deserialize an input archive to the packed insertion object.
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
   * @brief Generate an insertion point inside an insertion box.
   *
   * @param insertion_location Generated insertion location.
   * @param rng Random number generator.
   * @param insertion_information Insertion parameters read from the .prm file.
   */
  void
  generate_insertion_location(

    Point<dim>               &insertion_location,
    std::mt19937             &rng,
    const InsertionInfo<dim> &insertion_information);

  unsigned int current_inserting_particle_type = 0;

  // Number of particles of each type that remain to be inserted in the
  // upcoming insertion steps
  unsigned int particles_of_each_type_remaining;

  // Minimum and maximum boundaries of the insertion box in the direction order
  // It means that axis 0 is not necessarily x, since it depends on the order
  // of the insertion direction.
  std::vector<double> axis_min, axis_max;
};
#endif
