// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_serial_solid_h
#define lethe_serial_solid_h

#include <core/pvd_handler.h>
#include <core/simulation_control.h>
#include <core/solid_objects_parameters.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_fe.h>

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
   * @brief Stores the neighboring cells information for each cell of the solid object.
   * These containers store for each cell which neighboring cells are sharing a
   * vertex or an edge (two vertex). Also stores which neighboring cells are
   * coplanar or non-coplanar.
   */
  void
  setup_containers();

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
   * @brief Calculates the max displacement of the solid object. The max
   * displacement is defined as the L^\infty norm of the displacement
   * components. This is a measure of the maximal displacement in any direction.
   *
   * @return The reference to the vector of the displacement since the mapping
   */
  double
  get_max_displacement_since_mapped() const;

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

    else if constexpr (spacedim == 2) // spacedim == 2
      return tensor_nd_to_3d(this->current_translational_velocity);
  }

  /**
   * @brief Returns the angular velocity of the center of mass of the solid object
   *
   * @return The angular velocity of the solid object
   */
  inline Tensor<1, 3>
  get_angular_velocity() const
  {
    return this->current_angular_velocity;
  }

  /**
   * @brief Returns the center of rotation of the center of mass of the solid object
   *
   * @return The center of rotation of the solid object
   */
  inline Point<3>
  get_center_of_rotation() const
  {
    if constexpr (spacedim == 3)
      return this->center_of_rotation;

    else
      return point_nd_to_3d(this->center_of_rotation);
  }

  /**
   * @brief Returns the temperature of the solid object
   *
   * @return The temperature of the solid object
   */
  inline double
  get_temperature() const
  {
    return this->current_solid_temperature;
  }

  /**
   * @brief Returns the thermal boundary type of the solid object
   *
   * @return The thermal boundary type of the solid object
   */
  inline Parameters::ThermalBoundaryType
  get_thermal_boundary_type() const
  {
    return this->thermal_boundary_type;
  }

  /**
   * @brief Returns the solid id of the solid object.
   * This ID is an unsigned integer which is given to the solid object at
   * compile time and is mostly used when writing files.
   *
   * @return id of the solid
   */
  inline unsigned int
  get_solid_id() const
  {
    return id;
  }

  struct cell_comparison
  {
    inline bool
    operator()(
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell_1,
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell_2)
      const
    {
      return cell_1->global_active_cell_index() <
             cell_2->global_active_cell_index();
    }
  };

  typedef std::map<
    typename Triangulation<dim, spacedim>::active_cell_iterator,
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
    cell_comparison>
    triangulation_cell_neighbors_map;

  inline std::tuple<triangulation_cell_neighbors_map,
                    triangulation_cell_neighbors_map,
                    triangulation_cell_neighbors_map,
                    triangulation_cell_neighbors_map>
  get_neighbors_maps() const
  {
    return std::make_tuple(cp_es_neighbors,
                           cp_vs_neighbors,
                           np_es_neighbors,
                           np_vs_neighbors);
  }

  /**
   * @brief Update the temperature of the solid object, which is a function of time.
   *
   * @param initial_time The initial time (time t) of the timestep. This is used to
   * set the time of the temperature function.
   */
  void
  update_solid_temperature(const double initial_time);

  /**
   * @brief Moves the vertices of the solid triangulation. This function
   * uses explicitEuler scheme time integrator to displace the vertices
   * of the solid triangulation and stores the displacement in an array
   * in order to allow correct checkpointing of the triangulation.
   *
   * @param time_step The time_step value for this iteration
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
  write_output_results(
    const std::shared_ptr<SimulationControl> &simulation_control);

  /**
   * @brief read solid base triangulation checkpoint and replaces the
   * triangulation at the location where it was when it was checkpointed
   *
   * @param prefix_name The prefix of the checkpoint of the simulation
   *
   */
  void
  read_checkpoint(const std::string &prefix_name);

  /**
   * @brief write solid base triangulation checkpoint
   *
   * @param prefix_name The prefix of the checkpoint of the simulation
   */
  void
  write_checkpoint(const std::string &prefix_name);


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
   * @brief Moves the vertices of the solid triangulation using the displacement
   * vector. This is mostly used when restarting simulations which contain
   * a serial solid in order to replace it at the position it was when the
   * simulation finished.
   */
  void
  move_solid_triangulation_with_displacement();

  /**
   * @brief Rotates the grid by a given angle around an axis
   *
   * @param angle The angle (in rad) with which the solid is rotated
   * @param axis The axis over which the solid is rotated.
   *
   */
  void
  rotate_grid(const double angle, [[maybe_unused]] const Tensor<1, 3> &axis);

  /**
   * @brief Translate the grid. In spacedim=2, the third component is ignore
   *
   * @param translate The vector with which the solid is translated.
   */
  void
  translate_grid(const Tensor<1, 3> &translate);

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
  const bool output_bool;

  // Elements used to store the velocity of the solid object
  std::shared_ptr<Function<spacedim>> translational_velocity;
  std::shared_ptr<Function<spacedim>> angular_velocity;
  Point<spacedim>                     center_of_rotation;
  Tensor<1, spacedim>                 current_translational_velocity;
  Tensor<1, 3>                        current_angular_velocity;
  std::shared_ptr<Function<spacedim>> solid_temperature;
  Parameters::ThermalBoundaryType     thermal_boundary_type;
  double                              current_solid_temperature;

  // CP : coplanar
  // NP : non-coplanar
  // ES : edge-sharing
  // VS : vertex-sharing
  triangulation_cell_neighbors_map cp_es_neighbors;
  triangulation_cell_neighbors_map np_es_neighbors;
  triangulation_cell_neighbors_map cp_vs_neighbors;
  triangulation_cell_neighbors_map np_vs_neighbors;
};

#endif
