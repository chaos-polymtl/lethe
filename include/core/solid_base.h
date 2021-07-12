/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Carole-Anne Daunais, Valérie Bibeau, Polytechnique Montreal, 2020
 * Jeanne Joachim, Polytechnique Montreal, 2020-
 */

#ifndef lethe_solid_base_h
#define lethe_solid_base_h


// Lethe Includes
#include <core/nitsche.h>
#include <core/parameters.h>

// Dealii Includes

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/particle_handler.h>



using namespace dealii;

/**
 * A base class that generates a particle handler for the solid
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 * @tparam spacedim An integer that denotes the dimension of the space occupied
 * by the embedded solid
 *
 * @author Carole-Anne Daunais, Valerie Bibeau, 2020
 */

template <int dim, int spacedim = dim>
class SolidBase
{
public:
  // Member functions
  SolidBase(std::shared_ptr<Parameters::NitscheSolid<spacedim>> &param,
            std::shared_ptr<parallel::DistributedTriangulationBase<spacedim>>
                                               fluid_tria,
            std::shared_ptr<Mapping<spacedim>> fluid_mapping,
            const unsigned int                 degree_velocity);

  /**
   * @brief Manages solid triangulation and particles setup
   */
  void
  initial_setup();

  /**
   * @brief Generates a solid triangulation from a dealii or gmsh mesh
   */
  void
  setup_triangulation(const bool restart);

  /**
   * @brief Creates a particle handler in the fluid triangulation domain
   */
  void
  setup_particles_handler();

  /**
   * @brief Sets-up particles with particle handler, in the fluid triangulation domain that holds the particles of the solid
   * according to a specific quadrature
   */
  void
  setup_particles();

  /**
   * @brief Loads a solid triangulation from a restart file
   */
  void
  load_triangulation(const std::string filename_tria);

  /**
   * @brief Loads a particle handler in the fluid triangulation domain that holds the particles of the solid
   * according to a specific quadrature, and sets up dofs
   */
  void
  load_particles(const std::string filename_part);

  /**
   * @return the reference to the std::shared_ptr of a Particles::ParticleHandler<spacedim> that contains the solid particle handler
   */
  std::shared_ptr<Particles::ParticleHandler<spacedim>> &
  get_solid_particle_handler();

  /**
   * @return shared_ptr of the solid triangulation
   */
  std::shared_ptr<parallel::DistributedTriangulationBase<dim, spacedim>>
  get_solid_triangulation();

  /**
   * @return the reference to the solid dof handler
   */
  DoFHandler<dim, spacedim> &
  get_solid_dof_handler();

  /**
   * @return Function<spacedim> of the solid velocity
   */
  Function<spacedim> *
  get_solid_velocity();

  /**
   * @brief Updates particle positions in solid_particle_handler by integrating velocity using a Runge-Kutta method
   */
  void
  integrate_velocity(double time_step);

  /**
   * @brief Moves the dofs of the solid_triangulation
   */
  void
  move_solid_triangulation(double time_step);

  /**
   * @brief prints the positions of the particles
   */
  void
  print_particle_positions();


private:
  // Member variables
  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  std::shared_ptr<parallel::DistributedTriangulationBase<dim, spacedim>>
                                                                    solid_tria;
  std::shared_ptr<parallel::DistributedTriangulationBase<spacedim>> fluid_tria;
  DoFHandler<dim, spacedim>                                         solid_dh;
  std::shared_ptr<Particles::ParticleHandler<spacedim>> solid_particle_handler;

  std::shared_ptr<FiniteElement<dim, spacedim>> fe;

  std::shared_ptr<Mapping<dim, spacedim>> solid_mapping;
  std::shared_ptr<Mapping<spacedim>>      fluid_mapping;
  std::shared_ptr<Quadrature<dim>>        quadrature;

  std::shared_ptr<Parameters::NitscheSolid<spacedim>> &param;

  Function<spacedim> *velocity;

  const unsigned int degree_velocity;
  unsigned int       initial_number_of_particles;
};

#endif
