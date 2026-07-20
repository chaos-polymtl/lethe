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
   * For the packing insertion, nothing needs to be stored.
   *
   * @param ar Output archive where the attributes are stored.
   */
  virtual void
  serialize(boost::archive::text_oarchive &ar) override
  {
    (void)ar;
  }

  /**
   * @brief Deserialize an input archive to the packed insertion object.
   * For the packing insertion, nothing needs to be stored.
   *
   * @param ar Input archive where the attributes are stored.
   *
   */
  virtual void
  deserialize(boost::archive::text_iarchive &ar) override
  {
    (void)ar;
  }

  /**
   * @brief Updates the previous positions of particles in all locally owned cells.
   *
   * This method iterates through all active cells in the triangulation,
   * identifying locally owned cells. For each particle in these cells, it
   * retrieves the current particle position and stores it into the particle's
   * properties. The position is stored within certain indices that are
   * conventionally designated as velocity components (v_x, v_y, and v_z for 3D
   * configurations). This allows the tracking of previous particle positions
   * without additional data structures, which would drastically increase the
   * compilation time. This tracking is required for the packed insertion method.
   *
   * @param[out] particle_handler The particle handler of particles which are
   * being inserted.
   *
   */
  static void
  update_previous_position(Particles::ParticleHandler<dim> &particle_handler);

  /**
   * @brief Restricts the displacement of particles to a maximum allowable
   * value, preventing excessive movement during the simulation.
   *
   * This method iterates over all particles and checks their displacement
   * norm relative to their previous positions. If the displacement exceeds
   * the maximum particle diameter, their positions are clamped such that
   * the displacement does not exceed this threshold. The displacement
   * values are updated accordingly for each particle.
   *
   * @param[in,out] particle_handler The particle handler of particles which are
   * being inserted.
   * @param[in] max_disp Maximum displacement that a particle is able to move
   * during a single pseudo time-step.
   * @param[in] displacements Vector that stored the displacement of each
   * particle since the last contact detection time step.
   */
  static void
  clamp_displacement(Particles::ParticleHandler<dim> &particle_handler,
                     const double                     max_disp,
                     std::vector<double>             &displacements);

private:
  /**
   * @brief Generate an insertion point inside an insertion box.
   *
   * @param insertion_location Generated insertion location.
   * @param rng Random number generator.
   * @param insertion_information Insertion parameters read from the .prm file.
   */
  void
  generate_insertion_location(Point<dim>               &insertion_location,
                              std::mt19937             &rng,
                              const InsertionInfo<dim> &insertion_information);

  /// Current particle type being inserted
  unsigned int current_inserting_particle_type = 0;

  // Number of particles of each type that remain to be inserted in the
  // upcoming insertion steps
  unsigned int particles_of_each_type_remaining;

  // Minimum and maximum boundaries of the insertion box in the direction order
  // It means that axis 0 is not necessarily x, since it depends on the order
  // of the insertion direction.
  std::vector<double> axis_min, axis_max;

  /// Function used to insert particles inside or outside an arbitrary shape.
  std::shared_ptr<Function<dim>> acceptance_fct;
};
#endif
