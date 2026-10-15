// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_insertion_packed_h
#define lethe_insertion_packed_h

#include <dem/insertion.h>
#include <dem/particle_interaction_outcomes.h>

#include <random>
#include <vector>

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
   * @brief Displace the locally owned particles by the displacements
   * accumulated by the shift contact models.
   *
   * Each contact of a particle adds a displacement to the entry of that
   * particle in the force container of the contact outcomes. Each particle is
   * moved by the sum of the displacements of its contacts, limited to a
   * maximum norm, and the norm of the applied displacement is added to the
   * displacement used by the dynamic contact detection. The force container is
   * then reset to zero for the next pseudo-time step.
   *
   * @param[in,out] particle_handler The particle handler of particles which are
   * being inserted.
   * @param[in] max_displacement Maximum displacement of a particle during a
   * single pseudo-time step.
   * @param[in,out] contact_outcome Contact outcomes holding the accumulated
   * displacements of the particles.
   * @param[in,out] displacements Displacement of each particle since the last
   * contact detection step. The norm of the applied displacement is added.
   */
  static void
  apply_shift_displacement(
    Particles::ParticleHandler<dim>              &particle_handler,
    const double                                  max_displacement,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome,
    std::vector<double>                          &displacements);

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


  // Minimum and maximum boundaries of the insertion box in the direction order
  // It means that axis 0 is not necessarily x, since it depends on the order
  // of the insertion direction.
  std::vector<double> axis_min, axis_max;

  /// Function used to insert particles inside or outside an arbitrary shape.
  std::shared_ptr<Function<dim>> acceptance_fct;
};
#endif
