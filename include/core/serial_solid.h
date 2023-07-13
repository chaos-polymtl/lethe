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
 */

#ifndef lethe_serial_solid_h
#define lethe_serial_solid_h


// Lethe Includes
#include <core/parameters.h>
#include <core/pvd_handler.h>
#include <core/simulation_control.h>
#include <core/solid_objects_parameters.h>
#include <core/tensors_and_points_dimension_manipulation.h>

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
 * A class that generates an object to manage a serial triangulation that
 * represents a solid object. This class can represent solid objects for which
 * the dimension do not match the spatial dimension. For example, a 2D surface
 * defined in 3D or a line defined in 2D.
 *
 * @tparam dim An integer that denotes the dimensionality of the geometry
 * @tparam spacedim An integer that denotes the dimension of the space occupied
 * by the embedded solid
 *
 */

template <int dim, int spacedim = dim>
class SerialSolid
{
public:
  SerialSolid(std::shared_ptr<Parameters::RigidSolidObject<spacedim>> &param,
              unsigned int                                             id);

  /**
   * @brief Maps the solid object in the background triangulation
   */
  std::vector<
    std::pair<typename Triangulation<spacedim>::active_cell_iterator,
              typename Triangulation<dim, spacedim>::active_cell_iterator>>
  map_solid_in_background_triangulation(
    const parallel::TriangulationBase<spacedim> &background_tr);

  /**
   * @brief Manages solid triangulation setup
   */
  void
  initial_setup();

  /**
   * @brief Generates a solid triangulation from a dealii or gmsh mesh
   *
   * @param restart Boolean variable that indicates if the simulation is being restarted
   */
  void
  setup_triangulation(const bool restart);

  /**
   * @brief Set-ups the displacement vector to store the motion of the triangulation
   */
  void
  setup_dof_handler();

  /**
   * @brief Returns the shared pointer to the solid object triangulation
   *
   * @return shared_ptr of the solid triangulation
   */
  std::shared_ptr<Triangulation<dim, spacedim>>
  get_triangulation();

  /**
   * @brief Returns the dof_handler associated with the displacement of the solid object
   *
   * @return the reference to the displacement dof handler
   */
  DoFHandler<dim, spacedim> &
  get_dof_handler();

  /**
   * @brief Returns the displacement vector of all the vertices of the solid object
   *
   * @return the reference to the displacement vector
   */
  Vector<double> &
  get_displacement_vector();

  /**
   * @brief Calculates the max displacement of the solid object. The max displacement is defined as
   *  the L^\infty norm of the displacement components. This is a measure of the
   * maximal displacement in any direction.
   *
   * @return The reference to the vector of the displacement since the mapping
   */
  double
  get_max_displacement_since_mapped();

  /**
   * @brief Returns translational velocity of the center of mass of the solid object
   *
   * @return The translational velocity of the solid object
   */
  inline Tensor<1, 3>
  get_translational_velocity() const
  {
    if constexpr (spacedim == 3)
      return this->current_translational_velocity;

    if constexpr (spacedim == 2)
      return tensor_nd_to_3d(this->current_translational_velocity);
  }

  /**
   * @brief Returns the angular velocity of the center of mass of the solid object
   *
   * @return The angular velocity of the solid object
   */
  Tensor<1, 3>
  get_angular_velocity() const
  {
    return this->current_angular_velocity;
  }

  /**
   * @brief Returns the center of rotation of the center of mass of the solid object
   *
   * @return The center of rotation of the solid object
   */
  Point<3>
  get_center_of_rotation() const
  {
    if constexpr (spacedim == 3)
      return this->center_of_rotation;

    if constexpr (spacedim == 2)
      return point_nd_to_3d(this->center_of_rotation);
  }

  /**
   * @brief Returns the solid id of the solid object.
   * This ID is an unsigned integer which is given to the solid object at
   * compile time and is mostly used when writing files.
   *
   * @return id of the solid
   */
  unsigned int
  get_solid_id() const
  {
    return id;
  }

  /**
   * @brief Moves the vertices of the solid triangulation. This function
   * uses explicitEuler scheme time integrator to displace the vertices
   * of the solid triangulation and stores the displacement in an array
   * in order to allow correct checkpointing of the triangulation.
   *
   * @param time_step The time_step value for this iteration
   *
   * @param initial_time The initial time (time t) of the timestep. This is used to
   * set the time of the velocity function.
   */
  void
  move_solid_triangulation(double time_step, double initial_time);

  /**
   * @brief Write the output file of the solid object triangulation
   *
   * @param simulation_control The simulation control object
   */
  void
  write_output_results(std::shared_ptr<SimulationControl> simulation_control);

  /**
   * @brief read solid base triangulation checkpoint
   *
   * @param prefix_name The prefix of the checkpoint of the simulation
   *
   * @warning This is currently not implemented
   */
  void
  read_checkpoint(std::string prefix_name);

  /**
   * @brief write solid base triangulation checkpoint
   *
   * @param prefix_name The prefix of the checkpoint of the simulation
   *
   * @warning This is currently not implemented
   */
  void
  write_checkpoint(std::string prefix_name);


private:
  /**
   * @brief Moves the dofs of the solid_triangulation by using the displacement vector.
   * This is only used to move the solid triangulation at the correct location
   * when the simulation is restarted.
   */
  void
  displace_solid_triangulation();

  /**
   * @brief Reset displacements since intersection for contact detection
   */
  void
  reset_displacement_monitoring()
  {
    displacement_since_mapped = 0;
  }

  /**
   * @brief Rotates the grid by a given angle around an axis
   *
   * @param angle The angle (in rad) with which the solid is rotated
   *
   * @param axis The axis over which the solid is rotated.
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

  // Member variables

  // MPI
  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  // Parameters
  std::shared_ptr<Parameters::RigidSolidObject<spacedim>> &param;

  // Identifier of the solid
  unsigned int id;

  // Triangulation and DOFs management
  std::shared_ptr<Triangulation<dim, spacedim>> solid_tria;
  std::shared_ptr<FiniteElement<dim, spacedim>> fe;
  std::shared_ptr<Mapping<dim, spacedim>>       solid_mapping;
  DoFHandler<dim, spacedim>                     displacement_dh;
  std::shared_ptr<FESystem<dim, spacedim>>      displacement_fe;
  Vector<double>                                displacement;
  Vector<double>                                displacement_since_mapped;

  // Output management
  PVDHandler pvdhandler;

  // Elements used to store the velocity of the solid object
  Function<spacedim> *translational_velocity;
  Function<spacedim> *angular_velocity;
  Point<spacedim>     center_of_rotation;
  Tensor<1, spacedim> current_translational_velocity;
  Tensor<1, 3>        current_angular_velocity;
};

#endif
