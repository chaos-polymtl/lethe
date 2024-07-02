/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#ifndef lethe_navier_stokes_base_h
#define lethe_navier_stokes_base_h


// Lethe Includes
#include <core/mesh_controller.h>
#include <core/parameters.h>
#include <core/physics_solver.h>
#include <core/pvd_handler.h>
#include <core/simulation_control.h>

#include <solvers/flow_control.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/postprocessing_velocities.h>
#include <solvers/postprocessors.h>
#include <solvers/postprocessors_smoothing.h>
#include <solvers/simulation_parameters.h>

// Dealii Includes
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/data_out.h>

using namespace dealii;

/**
 * @brief Struct containing fluid id, temperature and phase fraction range
 * information, and flag containers for DOFs used in temperature-dependent
 * stasis constraints.
 *
 * @param[in] fluid_id Identifier of the fluid that is constrained.
 *
 * @param[in] min_solid_temperature Lower threshold value of the constraining
 * field (temperature).
 *
 * @param[in] max_solid_temperature Upper threshold values of the constraining
 * field (temperature).
 *
 * @param[in] filtered_phase_fraction_tolerance Tolerance applied on filtered
 * phase fraction.
 */
struct StasisConstraintWithTemperature
{
  StasisConstraintWithTemperature(
    const unsigned int fluid_id,
    const double       min_solid_temperature,
    const double       max_solid_temperature,
    const double       filtered_phase_fraction_tolerance)
    : fluid_id(fluid_id)
    , min_solid_temperature(min_solid_temperature)
    , max_solid_temperature(max_solid_temperature)
    , filtered_phase_fraction_tolerance(filtered_phase_fraction_tolerance)
  {}
  /// Identifier of the fluid that is constrained.
  const unsigned int fluid_id;
  /// Lower threshold values of the constraining field (temperature)
  const double min_solid_temperature;
  /// Upper threshold values of the constraining field (temperature)
  const double max_solid_temperature;
  /// Tolerance applied on filtered phase fraction
  const double filtered_phase_fraction_tolerance;
  /// Container of global DOF indices located in solid cells
  std::unordered_set<types::global_dof_index> dofs_are_in_solid;
  /// Container of global DOF indices connected to at least one fluid cell
  std::unordered_set<types::global_dof_index> dofs_are_connected_to_fluid;
};

/**
 * A base class for all the Navier-Stokes equation
 * This class regroups common facilities that are shared by all
 * the Navier-Stokes implementations to reduce code multiplicity
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @tparam VectorType  The Vector type used for the solvers
 *
 * @tparam DofsType the type of dof storage indices
 *
 * @ingroup solvers
 * @author Bruno Blais, 2019
 */

template <int dim, typename VectorType, typename DofsType>
class NavierStokesBase : public PhysicsSolver<VectorType>
{
protected:
  NavierStokesBase(SimulationParameters<dim> &nsparam);

  virtual ~NavierStokesBase()
  {}

  /**
   * @brief Getter methods to get the private attributes for the physic currently solved
   *
   * @param current_physics_id Indicates the number associated with the physic currently solved.
   * If the solver is only solving for a fluid dynamics problem, then the value
   * will always be PhysicsID::fluid_dynamics.
   *
   */
  virtual VectorType &
  get_evaluation_point() override
  {
    return evaluation_point;
  };
  virtual VectorType &
  get_local_evaluation_point() override
  {
    return local_evaluation_point;
  };
  virtual VectorType &
  get_newton_update() override
  {
    return newton_update;
  };
  virtual VectorType &
  get_present_solution() override
  {
    return present_solution;
  };
  virtual VectorType &
  get_system_rhs() override
  {
    return system_rhs;
  };
  virtual AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return nonzero_constraints;
  };

  /**
   * @brief Output the L2 and Linfty norms of the correction vector.
   *
   * @param[in] display_precision Number of outputted digits.
   */
  virtual void
  output_newton_update_norms(const unsigned int display_precision) override;

  /**
   *  Generic interface routine to allow the CFD solver
   *  to cooperate with the multiphysics modules
   **/

  /**
   * @brief finish_simulation
   * Finish the simulation by calling all
   * the post-processing elements that are required
   */
  void
  finish_simulation()
  {
    finish_simulation_fd();
    multiphysics->finish_simulation();
  }

  /**
   * @brief postprocess
   * Post-process simulation after an iteration
   *
   * @param first_iteration Indicator if the simulation is at it's first simulation or not.
   */
  virtual void
  postprocess(bool first_iteration)
  {
    postprocess_fd(first_iteration);
    multiphysics->postprocess(first_iteration);
  };

  /**
   * @brief setup_dofs
   *
   * Initialize the degree of freedom and the memory
   * associated with them for fluid dynamics and enabled auxiliary physics.
   */
  virtual void
  setup_dofs()
  {
    setup_dofs_fd();
    multiphysics->setup_dofs();
  };

  /**
   * @brief set_initial_conditions
   *
   * @param initial_condition_type Type of method  use to impose initial condition.
   *
   * @param restart Indicator if the simulation is being restarted or not.
   *
   **/

  virtual void
  set_initial_condition(Parameters::InitialConditionType initial_condition_type,
                        bool                             restart = false)
  {
    unsigned int ref_iter = 0;
    do
      {
        if (ref_iter > 0)
          this->refine_mesh();

        set_initial_condition_fd(initial_condition_type, restart);
        if (!restart)
          {
            multiphysics->set_initial_conditions();
            this->postprocess_fd(true);
            multiphysics->postprocess(true);
          }
        ref_iter++;
      }
    while (
      ref_iter <
        (this->simulation_parameters.mesh_adaptation.initial_refinement + 1) &&
      restart == false);
  }

  /**
   * Key physics component for fluid dynamics
   **/

  /**
   * @brief finish_time_step
   * Finishes the time step of the fluid dynamics
   * Post-processing and time stepping
   */
  virtual void
  finish_time_step();

  /**
   * @brief finish_time_step
   * Finishes the time step of the fluid dynamics
   * Post-processing and time stepping
   */
  virtual void
  percolate_time_vectors_fd();

  /**
   * @brief finish_simulation
   * Finishes the simulation for fluid dynamics by calling
   * the post-processing elements that are required
   */
  void
  finish_simulation_fd();

  /**
   * @brief postprocess
   * Post-process fluid dynamics after an iteration
   */
  virtual void
  postprocess_fd(bool first_iteration);

  /**
   * @brief setup_dofs
   *
   * Initialize the dofs for fluid dynamics
   */
  virtual void
  setup_dofs_fd() = 0;

  /**
   * @brief  Update the time average velocity field solution
   */
  virtual void
  update_multiphysics_time_average_solution() = 0;

  virtual void
  set_initial_condition_fd(
    Parameters::InitialConditionType initial_condition_type,
    bool                             restart = false) = 0;

  /**
   * End of key physics components for fluid dynamics
   **/

  /**
   * @brief postprocessing_forces
   * Post-processing function
   * Outputs the forces acting on each boundary condition
   */
  void
  postprocessing_forces(const VectorType &evaluation_point);

  /**
   * @brief postprocessing_torques
   * Post-processing function
   * Outputs the torque acting on each boundary condition
   */
  void
  postprocessing_torques(const VectorType &evaluation_point);

  /**
   * @brief dynamic_flow_control
   * If set to enable, dynamic_flow_control allows to control the flow by
   * executing space-average velocity and beta coefficient force calculation at
   * each time step.
   */
  virtual void
  dynamic_flow_control();

  /**
   * @brief iterate
   * Do a regular CFD iteration
   */
  virtual void
  iterate();

  /**
   * @brief Enable the use of dynamic zero constraints by initializing required
   * FEValues objects.
   *
   * @note At the moment, the only solution-dependent dynamic constraint depends
   * on the temperature field obtained from the Heat Transfer (HT) auxiliary
   * physic.
   */
  virtual void
  enable_dynamic_zero_constraints_fd();

  /**
   * @brief Allow the initial refinement of all cells of the principal mesh that are partially
   * contained in one of the cells of the box refinement mesh given in the
   * parameter.
   *
   * @param restart Indicator if the simulation is being restarted or not.
   */
  void
  box_refine_mesh(const bool restart);

  /**
   * @brief Allow the refinement of the mesh according to one of the 2 methods proposed
   */
  virtual void
  refine_mesh();

  /**
   * @brief Allow the refinement of the mesh based on the Kelly error estimator.
   * See :
   * https://www.dealii.org/current/doxygen/deal.II/classKellyErrorEstimator.html
   * for more information on the Kelly error estimator.
   */
  void
  refine_mesh_kelly();

  /**
   * @brief Allow the uniform refinement of all the mesh.
   */
  void
  refine_mesh_uniform();

  /**
   * @brief read_checkpoint
   */
  virtual void
  read_checkpoint();

  /**
   * @brief set_nodal_values
   * Set the nodal values of velocity and pressure
   */
  void
  set_nodal_values();

  /**
   * @brief Check if a specifique boundary condition exist
   * @param bc, the boundary type that we want to check if it exists
   */
  bool
  check_existance_of_bc(BoundaryConditions::BoundaryType bc);

  /**
   * @brief Turn regions of the mesh where the @p material_id>0 into a solid
   * block by injecting velocity and pressure DOFs into the zero constraints.
   *
   * It is achieved by imposing \f$\mathbf{u}=0\f$ within the cells which have a
   * @p material_id>0. In addition, solid cells which are not connected to the
   * fluid by any means also get a pressure Dirichlet boundary condition which
   * fixes the pressure to 0. It ensures that the linear system is well-posed.
   * Right now, this routine only supports the usage of 1 solid domain, but
   * eventually it could be extended to more than one. By default, the fluid
   * domain is assumed to have a @p material_id=0 and the rest of the domains
   * have a @p material_id>0.
   *
   * @param[in] non_zero_constraints If this parameter is true, it indicates
   * that non-zero constraints are being constrained for the solid domain. If
   * this is set to false, homogeneous constraints are constrained in the solid
   * domain.
   */
  void
  establish_solid_domain(const bool non_zero_constraints);

  /**
   * @brief Get cell's local temperature values at quadrature points.
   *
   * @param[in] cell Pointer to an active cell of the fluid dynamics DoFHandler.
   *
   * @param[in] dof_handler_ht DoFHandler of the Heat Transfer (HT) auxiliary
   * physic.
   *
   * @param[in] temperature_solution Temperature solution vector from the HT
   * auxiliary physic.
   *
   * @param[out] local_temperature_values Cell's local temperature values at
   * quadrature points.
   */
  inline void
  get_cell_temperature_values(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const DoFHandler<dim>                                *dof_handler_ht,
    const TrilinosWrappers::MPI::Vector                  &temperature_solution,
    std::vector<double> &local_temperature_values)
  {
    const typename DoFHandler<dim>::active_cell_iterator temperature_cell(
      &(*(this->triangulation)), cell->level(), cell->index(), dof_handler_ht);

    this->fe_values_temperature->reinit(temperature_cell);
    this->fe_values_temperature->get_function_values(temperature_solution,
                                                     local_temperature_values);
  }

  /**
   * @brief Get cell's local filtered phase fraction values at quadrature
   * points.
   *
   * @param[in] cell Pointer to an active cell of the fluid dynamics DoFHandler.
   *
   * @param[in] dof_handler_vof DoFHandler of the Volume of Fluid (VOF)
   * auxiliary physic.
   *
   * @param[in] filtered_phase_fraction_solution Filtered phase fraction
   * solution vector from the VOF auxiliary physic.
   *
   * @param[out] local_filtered_phase_fraction_values Cell's local filtered
   * phase fraction values at quadrature points.
   */
  inline void
  get_cell_filtered_phase_fraction_values(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const DoFHandler<dim>                                *dof_handler_vof,
    const TrilinosWrappers::MPI::Vector &filtered_phase_fraction_solution,
    std::vector<double>                 &local_filtered_phase_fraction_values)
  {
    const typename DoFHandler<dim>::active_cell_iterator vof_cell(
      &(*(this->triangulation)), cell->level(), cell->index(), dof_handler_vof);

    this->fe_values_vof->reinit(vof_cell);
    this->fe_values_vof->get_function_values(
      filtered_phase_fraction_solution, local_filtered_phase_fraction_values);
  }

  /**
   * @brief Flag cell DOFs (either "solid" or "connected to fluid") and
   * constrain velocity DOFs with homogeneous constraints if cell is considered
   * solid.
   *
   * @param[in] local_dof_indices Vector of a cell's local DOF indices.
   *
   * @param[in] local_temperature_values Cell's local temperature values at
   * quadrature points.
   *
   * @param[in,out] stasis_constraint_struct Struct containing flagged DOF
   * containers, temperature range information and fluid id.
   */
  void
  add_flags_and_constrain_velocity(
    const std::vector<types::global_dof_index> &local_dof_indices,
    const std::vector<double>                  &local_temperature_values,
    StasisConstraintWithTemperature            &stasis_constraint_struct);

  /**
   * @brief Check if solid cells are connected to fluid ones and constrain null
   * pressure DOFs if they are not.
   *
   * The check is done by looking if the global DOF indices are located in the
   * flag containers (@p dofs_are_in_solid and @p dofs_are_connected_to_fluid)
   * of @p stasis_constraint_struct.
   *
   * @param[in] stasis_constraint_struct Struct containing flagged DOF
   * containers, temperature range information and fluid id.
   *
   * @param[in,out] local_dof_indices Vector for storing a cell's local DOF
   * indices.
   */
  void
  check_and_constrain_pressure(
    const StasisConstraintWithTemperature &stasis_constraint_struct,
    std::vector<types::global_dof_index>  &local_dof_indices);

  /**
   * @brief Constrain velocity DOFs of a solid cell.
   *
   * @param[in] non_zero_constraints If this parameter is true, it indicates
   * that non-zero constraints are applied in the solid domain. If
   * this is set to false, homogeneous constraints are applied in the solid
   * domain.
   *
   * @param[in] local_dof_indices Vector of a cell's local DOF indices.
   *
   * @param[out] zero_constraints Homogeneous constraints holding object.
   */
  inline void
  constrain_solid_cell_velocity_dofs(
    const bool                                 &non_zero_constraints,
    const std::vector<types::global_dof_index> &local_dof_indices,
    AffineConstraints<double>                  &zero_constraints)
  {
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        const unsigned int component =
          this->fe->system_to_component_index(i).first;
        if (component < dim) // velocity DOFs
          {
            // We apply a constraint to all DOFs in the solid region, whether
            // they are locally owned or not.
            if (non_zero_constraints)
              {
                this->nonzero_constraints.add_line(local_dof_indices[i]);
                this->nonzero_constraints.set_inhomogeneity(
                  local_dof_indices[i], 0);
              }
            else
              zero_constraints.add_line(local_dof_indices[i]);
          }
      }
  }

  /**
   * @brief Flag DOFs in solid cells.
   *
   * @param[in] local_dof_indices Vector of a cell's local DOF indices.
   *
   * @param[out] dofs_are_in_solid Container of global DOF indices located in
   * solid cells.
   */
  inline void
  flag_dofs_in_solid(
    const std::vector<types::global_dof_index>  &local_dof_indices,
    std::unordered_set<types::global_dof_index> &dofs_are_in_solid)
  {
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        const unsigned int component =
          this->fe->system_to_component_index(i).first;
        if (component == dim)
          {
            dofs_are_in_solid.insert(local_dof_indices[i]);
          }
      }
  }

  /**
   * @brief Check if the cell is a solid.
   *
   * @param[in] dofs_are_in_solid Container of global DOF indices located in
   * solid cells.
   *
   * @param[in] local_dof_indices Vector of a cell's local DOF indices.
   *
   * @return Boolean indicating if the cell is in a solid (true) or not (false).
   */
  inline bool
  check_cell_is_in_solid(
    const std::unordered_set<types::global_dof_index> &dofs_are_in_solid,
    const std::vector<types::global_dof_index>        &local_dof_indices)
  {
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        auto search = dofs_are_in_solid.find(local_dof_indices[i]);
        if (search != dofs_are_in_solid.end())
          return true;
      }
    return false;
  }

  /**
   * @brief Flag DOFs connected to fluid cells.
   *
   * @param[in] local_dof_indices Vector of a cell's local DOF indices.
   *
   * @param[out] dofs_are_connected_to_fluid Container of global DOF indices
   * connected to at least one fluid cell.
   */
  inline void
  flag_dofs_connected_to_fluid(
    const std::vector<types::global_dof_index>  &local_dof_indices,
    std::unordered_set<types::global_dof_index> &dofs_are_connected_to_fluid)
  {
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        const unsigned int component =
          this->fe->system_to_component_index(i).first;
        if (component == dim)
          {
            dofs_are_connected_to_fluid.insert(local_dof_indices[i]);
          }
      }
  }

  /**
   * @brief Check if the cell is connected to a fluid cell.
   *
   * @param[in] dofs_are_connected_to_fluid Container of global DOF indices
   * connected to at least one fluid cell.
   *
   * @param[in] local_dof_indices Vector of a cell's local DOF indices.
   *
   * @return Boolean indicating if the cell is connected to a fluid (true) or
   * not (false).
   */
  inline bool
  check_cell_is_connected_to_fluid(
    const std::unordered_set<types::global_dof_index>
                                               &dofs_are_connected_to_fluid,
    const std::vector<types::global_dof_index> &local_dof_indices)
  {
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        auto search = dofs_are_connected_to_fluid.find(local_dof_indices[i]);
        if (search != dofs_are_connected_to_fluid.end())
          return true;
      }
    return false;
  }

  /**
   * @brief Constrain pressure DOFs if cells are not connected to fluid and DOFs
   * are locally owned.
   *
   * @param[in] non_zero_constraints If this parameter is true, it indicates
   * that non-zero constraints are applied in the solid domain. If
   * this is set to false, homogeneous constraints are applied in the solid
   * domain.
   *
   * @param[in] local_dof_indices Vector of a cell's local DOF indices.
   *
   * @param[out] zero_constraints Homogeneous constraints holding object.
   */
  inline void
  constrain_pressure(
    const bool                                 &non_zero_constraints,
    const std::vector<types::global_dof_index> &local_dof_indices,
    AffineConstraints<double>                  &zero_constraints)
  {
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        const unsigned int component =
          this->fe->system_to_component_index(i).first;

        // Only pressure DOFs have an additional Dirichlet condition
        if (component == dim) // pressure DOFs
          {
            // We only apply the constraint on the locally owned pressure DOFs
            // since we have no way of verifying if the locally relevant DOFs
            // are connected to a fluid cell.
            bool dof_is_locally_owned = false;

            // For the GLS-family of solvers, we only have a single index set
            if constexpr (std::is_same_v<DofsType, IndexSet>)
              {
                dof_is_locally_owned =
                  this->locally_owned_dofs.is_element(local_dof_indices[i]);
              }

            // For the GD-family of solvers, we have two index sets. One for
            // velocities and one for pressure.
            if constexpr (std::is_same_v<DofsType, std::vector<IndexSet>>)
              {
                dof_is_locally_owned =
                  this->locally_owned_dofs[1].is_element(local_dof_indices[i]);
              }

            if (dof_is_locally_owned)
              {
                if (non_zero_constraints)
                  {
                    this->nonzero_constraints.add_line(local_dof_indices[i]);
                    this->nonzero_constraints.set_inhomogeneity(
                      local_dof_indices[i], 0);
                  }
                else
                  {
                    zero_constraints.add_line(local_dof_indices[i]);
                  }
              }
          }
      }
  }

  /**
   * @brief Constrain a fluid's subdomains according to the temperature field to null
   * velocity and pressure fields to model solid subdomains.
   *
   * @note Its equivalent for Volume of Fluid (VOF) simulations is
   * NavierStokesBase<dim, VectorType, DofsType>::constrain_solid_domain_vof.
   *
   * @param[in] dof_handler_ht DoFHandler of the Heat Transfer (HT) auxiliary
   * physic.
   */
  void
  constrain_stasis_with_temperature(const DoFHandler<dim> *dof_handler_ht);

  /**
   * @brief Constrain fluids' subdomains according to the temperature field to
   * null velocity and pressure fields to model solid subdomains in Volume of
   * Fluid (VOF) simulations.
   *
   * @param[in] dof_handler_vof DoFHandler of the VOF auxiliary physic.
   *
   * @param[in] dof_handler_ht DoFHandler of the Heat Transfer (HT) auxiliary
   * physic.
   */
  void
  constrain_stasis_with_temperature_vof(const DoFHandler<dim> *dof_handler_vof,
                                        const DoFHandler<dim> *dof_handler_ht);

  /**
   * @brief write_checkpoint
   */
  virtual void
  write_checkpoint();

  /**
   * @brief write_output_results
   * Post-processing as parallel VTU files
   */
  void
  write_output_results(const VectorType &solution);

  /**
   * @brief output_field_hook
   * This function is to be redefined in specialized classes to adapt the output
   * to each solver.
   */
  virtual void
  output_field_hook(DataOut<dim> &);

  /**
   * @brief write_output_forces
   * Writes the forces per boundary condition to a text file output
   */
  void
  write_output_forces();

  /**
   * @brief write_output_torques
   * Writes the torques per boundary condition to a text file output
   */
  void
  write_output_torques();

  /**
   * @brief rescale_pressure_dofs_in_newton_update
   * This function is used to rescale pressure DOFs in the newton correction
   */
  void
  rescale_pressure_dofs_in_newton_update();

  /**
   * @brief init_temporary_vector
   * This function initializes correctly a temporary vector depending on the
   * vector type
   */
  inline VectorType
  init_temporary_vector();

  // Member variables
protected:
  DofsType locally_owned_dofs;
  DofsType locally_relevant_dofs;

  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  DoFHandler<dim>                                              dof_handler;
  std::shared_ptr<FESystem<dim>>                               fe;

  TimerOutput computing_timer;

  SimulationParameters<dim> simulation_parameters;
  PVDHandler                pvdhandler;

  // Functions used for source term and error analysis
  Function<dim>                 *exact_solution;
  std::shared_ptr<Function<dim>> forcing_function;

  // Dynamic flow control
  FlowControl<dim> flow_control;

  // Constraints for Dirichlet boundary conditions
  AffineConstraints<double> zero_constraints;
  AffineConstraints<double> nonzero_constraints;

  // Present solution and non-linear solution components
  VectorType evaluation_point;
  VectorType local_evaluation_point;
  VectorType newton_update;
  VectorType present_solution;
  VectorType system_rhs;

  // Previous solutions vectors
  std::vector<VectorType> previous_solutions;

  // Finite element order used
  const unsigned int velocity_fem_degree;
  const unsigned int pressure_fem_degree;
  unsigned int       number_quadrature_points;

  // Mappings and Quadratures
  std::shared_ptr<Mapping<dim>>        mapping;
  std::shared_ptr<Quadrature<dim>>     cell_quadrature;
  std::shared_ptr<Quadrature<dim - 1>> face_quadrature;

  // Assemblers for the matrix and rhs
  std::vector<std::shared_ptr<NavierStokesAssemblerBase<dim>>> assemblers;

  // Multiphysics interface
  std::shared_ptr<MultiphysicsInterface<dim>> multiphysics;

  // Simulation control for time stepping and I/Os
  std::shared_ptr<SimulationControl> simulation_control;

  // Post-processing variables
  TableHandler enstrophy_table;
  TableHandler pressure_power_table;
  TableHandler viscous_dissipation_table;
  TableHandler kinetic_energy_table;
  TableHandler apparent_viscosity_table;
  TableHandler pressure_drop_table;
  TableHandler flow_rate_table;
  std::shared_ptr<AverageVelocities<dim, VectorType, DofsType>>
             average_velocities;
  VectorType average_solution;

  // Refinement control
  MeshController mesh_controller;

  // Convergence Analysis
  ConvergenceTable error_table;

  // Force analysis
  std::vector<std::vector<Tensor<1, dim>>> forces_on_boundaries;
  std::vector<Tensor<1, 3>>                torques_on_boundaries;
  std::vector<TableHandler>                forces_tables;
  std::vector<TableHandler>                torques_tables;

  /// FEValues object used for temperature-dependent solid domain constraints
  std::shared_ptr<FEValues<dim>> fe_values_temperature;
  /// FEValues object used for temperature-dependent solid domain constraints in
  /// VOF simulations
  std::shared_ptr<FEValues<dim>> fe_values_vof;
  /// Vector containing solid domain constraint structs for
  /// temperature-dependent solid domain constraints in VOF simulations
  std::vector<StasisConstraintWithTemperature> stasis_constraint_structs;
  /// Dynamic homogeneous constraints used for temperature-dependent solid
  /// domain constraints
  AffineConstraints<double> dynamic_zero_constraints;
};

#endif
