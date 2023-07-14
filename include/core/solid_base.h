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
#include <core/parameters.h>
#include <core/solid_objects_parameters.h>

// Dealii Includes

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/lac/trilinos_vector.h>

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
  SolidBase(std::shared_ptr<Parameters::NitscheObject<spacedim>> &param,
            std::shared_ptr<parallel::DistributedTriangulationBase<spacedim>>
                                               fluid_tria,
            std::shared_ptr<Mapping<spacedim>> fluid_mapping);

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
   * @brief Sets-up particles with particle handler, in the fluid triangulation domain that holds the particles of the solid
   * according to a specific quadrature
   */
  void
  setup_displacement();

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
   * @return the reference to the displacement dof handler
   */
  DoFHandler<dim, spacedim> &
  get_displacement_dof_handler();

  /**
   * @return the reference to the displacement vector
   */
  TrilinosWrappers::MPI::Vector &
  get_displacement_vector();

  /**
   * @return Function<spacedim> of the solid velocity
   */
  Function<spacedim> *
  get_solid_velocity();

  Function<spacedim> *
  get_solid_temperature();

  /**
   * @brief Updates particle positions in solid_particle_handler by integrating velocity using an explicit Runge-Kutta 4 method.
   *
   * @param time_step The time_step value for this iteration
   *
   * @param initial_time The initial time (time t) of the timestep. This is used to
   * set the time of the velocity function.
   */
  void
  integrate_velocity(double time_step, double initial_time);

  /**
   * @brief Moves the vertices of the solid triangulation. This function
   * uses an Runge-Kutta 4 explicit time integrator to displace the vertices
   * of the solid triangulation and stores the displacement in an array
   * in order to allow correct checkpointing of the triangulation.
   *               See
   https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
               for more details
   *
   * @param time_step The time_step value for this iteration
   *
   * @param initial_time The initial time (time t) of the timestep. This is used to
   * set the time of the velocity function.
   */
  void
  move_solid_triangulation(double time_step, double initial_time);

  /**
   * @brief Moves the dofs of the solid_triangulation by using the displacement vector.
   * This is only use to move the solid triangulation at the correct location
   * when the simulation is restarted.
   */
  void
  move_solid_triangulation_with_displacement();

  /**
   * @brief prints the positions of the particles
   */
  void
  print_particle_positions();

  /**
   * @brief Rotates the grid by a given angle around an axis
   *
   * @param angle The angle (in rad) with which the grid is rotated
   *
   * @param axis The axis over which the grid is rotated.
   *
   */
  void
  rotate_grid(const double angle, const Tensor<1, 3> axis);

  /**
   * @brief Translate the grid. In spacedim=2, the third component is ignore
   *
   * @param translate The vector with which the solid is translated.
   */
  void
  translate_grid(const Tensor<1, 3> translate);

  /**
   * @brief Updates the time in the function used to describe the solid temperature
   */
  void
  update_temperature_time(double time);

  /**
   * @brief read solid base triangulation checkpoint
   */
  void
  read_checkpoint(std::string prefix_name);

  /**
   * @brief write solid base triangulation checkpoint
   */
  void
  write_checkpoint(std::string prefix_name);


private:
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
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


  DoFHandler<dim, spacedim>                displacement_dh;
  std::shared_ptr<FESystem<dim, spacedim>> displacement_fe;
  TrilinosWrappers::MPI::Vector            displacement;
  TrilinosWrappers::MPI::Vector            displacement_relevant;


  std::shared_ptr<Parameters::NitscheObject<spacedim>> &param;

  Function<spacedim> *velocity;
  Function<spacedim> *temperature;

  unsigned int initial_number_of_particles;
};

#endif
