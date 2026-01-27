// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_navier_stokes_base_h
#define lethe_navier_stokes_base_h

#include <core/mesh_controller.h>
#include <core/mortar_coupling_manager.h>
#include <core/output_struct.h>
#include <core/parameters.h>
#include <core/physics_solver.h>
#include <core/pvd_handler.h>
#include <core/sdirk_stage_data.h>
#include <core/simulation_control.h>
#include <core/solutions_output.h>

#include <solvers/flow_control.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/postprocessing_scalar.h>
#include <solvers/postprocessing_velocities.h>
#include <solvers/postprocessors.h>
#include <solvers/simulation_parameters.h>

// Dealii Includes
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/data_out.h>

using namespace dealii;



DeclException1(
  FluidDynamicsBoundaryConditionMissing,
  types::boundary_id,
  << "The boundary id: " << arg1
  << " is defined in the triangulation, but not as a boundary condition for the fluid dynamics physics. Lethe does not assign a    default boundary condition to boundary ids. Every boundary id defined within the triangulation must have a corresponding boundary condition defined in the input file.");

/**
 * @brief Struct containing fluid id, temperature and phase fraction range
 * information, and flag containers for DOFs used in temperature-dependent
 * stasis constraints.
 */
struct StasisConstraintWithTemperature
{
  /**
   * @brief Default constructor of the struct.
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
    return *present_solution;
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
   * @brief Return the metric for residual rescaling. By default, should return 1.
   * If the rescale_residual_by_volume is set to true, the method
   * returns the global volume of the triangulation.
   *
   * @return Rescale metric.
   */
  double
  get_residual_rescale_metric() const override
  {
    return simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
               .rescale_residual_by_volume ?
             std::sqrt(
               GridTools::volume(*this->triangulation, *this->mapping)) :
             1.;
  }

  /**
   *  Generic interface routine to allow the CFD solver
   *  to cooperate with the multiphysics modules
   **/

  /**
   * @brief Finish the simulation by calling all the post-processing elements that are required
   */
  void
  finish_simulation()
  {
    finish_simulation_fd();
    multiphysics->finish_simulation();
  }

  /**
   * @brief Post-process simulation after an iteration
   *
   * @param first_iteration Indicator if the simulation is at its first simulation or not.
   */
  virtual void
  postprocess(bool first_iteration)
  {
    postprocess_fd(first_iteration);
    multiphysics->postprocess(first_iteration);

    if (this->simulation_control->is_output_iteration())
      this->write_output_results(*this->present_solution);
  };

  /**
   * @brief Initialize the degree of freedom and the memory associated with them for fluid dynamics and enabled auxiliary physics.
   */
  virtual void
  setup_dofs()
  {
    verify_consistency_of_boundary_conditions();
    setup_dofs_fd();
    multiphysics->setup_dofs();
  };

  /**
   * @brief Set the initial condition
   *
   * @param initial_condition_type Type of method  use to impose initial condition.
   *
   * @param restart Indicator if the simulation is being restarted or not.
   *
   **/

  virtual void
  set_initial_condition(
    const Parameters::FluidDynamicsInitialConditionType initial_condition_type,
    const bool                                          restart = false)
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
          }
        ref_iter++;
      }
    while (
      ref_iter <
        (this->simulation_parameters.mesh_adaptation.initial_refinement + 1) &&
      restart == false);

    if (!restart)
      {
        this->postprocess_fd(true);
        multiphysics->postprocess(true);
        if (this->simulation_control->is_output_iteration())
          this->write_output_results(*this->present_solution);
      }
  }

  /**
   * Key physics component for fluid dynamics
   **/

  /**
   * @brief Finishes the time step of the fluid dynamics. Post-processing and time stepping
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
   * @brief Finishes the simulation for fluid dynamics by calling the post-processing elements that are required
   */
  void
  finish_simulation_fd();

  /**
   * @brief Post-process fluid dynamics after an iteration
   */
  virtual void
  postprocess_fd(bool first_iteration);

  /**
   * @brief Initialize the dofs for fluid dynamics
   */
  virtual void
  setup_dofs_fd() = 0;

  /**
   * @brief Update mortar configuration.
   *
   * When the rotor domain is rotated, the mortar cells need to be reinitialized
   * according to the new rotor-stator interface configuration.
   */
  virtual void
  update_mortar_configuration()
  {}

  /**
   * @brief  Update the time average velocity field solution
   */
  virtual void
  update_multiphysics_time_average_solution() = 0;

  /**
   * @brief  Required only for the matrix-free solver to update solutions within iterate function.
   */
  virtual void
  update_solutions_for_multiphysics()
  {}

  /**
   * @brief  Required only for the matrix-free solver to update multiphysics solutions within iterate function.
   */
  virtual void
  update_solutions_for_fluid_dynamics()
  {}

  virtual void
  set_initial_condition_fd(
    Parameters::FluidDynamicsInitialConditionType initial_condition_type,
    bool                                          restart = false) = 0;

  /**
   * End of key physics components for fluid dynamics
   **/

  /**
   * @brief Post-processing function
   * Outputs the forces acting on each boundary condition
   */
  void
  postprocessing_forces(const VectorType &evaluation_point);

  /**
   * @brief Post-processing function
   * Outputs the torque acting on each boundary condition
   */
  void
  postprocessing_torques(const VectorType &evaluation_point);

  /**
   * @brief If set to enable, dynamic_flow_control allows to control the flow by executing space-average velocity and beta coefficient force calculation at each time step.
   */
  virtual void
  dynamic_flow_control();

  /**
   * @brief Update the sum_over_previous_stages variable to use it in the solver.
   *
   * @param stage An unsigned integer which gives the index of the current stage
   * @param method SDIRK method for now (BDF methods are single stage)
   */
  virtual void
  multi_stage_preresolution(
    unsigned int                                      stage,
    Parameters::SimulationControl::TimeSteppingMethod method);

  /**
   * @brief Calculate \f$ k_i \f$ from \f$ u*_i : k_i = \frac{u*_i - u_n}{time_step*a_ii} - \sum{j=1}^{i-1} a_ij * k_j \f$
   * Update \f$ \sum_{i=1}^{N_{stages}} b_i * k_i \f$
   *
   * @param stage An unsigned integer which gives the index of the current stage
   * @param method SDIRK method for now (BDF methods are single stage)
   * @param time_step The current time step to build the solution at the next time step
   */
  virtual void
  multi_stage_postresolution(
    unsigned int                                      stage,
    Parameters::SimulationControl::TimeSteppingMethod method,
    double                                            time_step);

  /**
   * @brief \f$ u_{n+1} = u_n + time_step * \sum_{i=1}^{N_{stages}} b_i * k_i \f$
   *
   * @param time_step The current time step to build the solution at the next time step
   */
  virtual void
  update_multi_stage_solution(double time_step);

  /**
   * @brief Do a regular CFD iteration
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
   * @brief Transfer solution after mesh refinement
   */
  void
  transfer_solution(SolutionTransfer<dim, VectorType> &solution_transfer,
                    std::vector<SolutionTransfer<dim, VectorType>>
                      &previous_solutions_transfer);

  /**
   * @brief Restart a previous simulation from a checkpoint file.
   */
  virtual void
  read_checkpoint();

  /**
   * @brief Read the solution from a checkpoint file and set it as the current solution.
   */
  virtual void
  set_solution_from_checkpoint(std::string checkpoint_file_prefix);

  /**
   * @brief Set the nodal values of velocity and pressure
   */
  void
  set_nodal_values();

  /**
   * @brief Define the non-zero constraints used to solve the problem.
   */
  void
  define_non_zero_constraints();

  /**
   * @brief Define the zero constraints used to solve the problem.
   */
  void
  define_zero_constraints();

  /**
   * @brief Initialize mortar coupling manager, operator, and evaluator
   */
  void
  reinit_mortar_operators();

  /**
   * @brief Returns the mapping shared pointer. A MappingQCache is
   * necessary for prescribed rotation in rotor-stator configurations
   */
  inline std::shared_ptr<Mapping<dim>>
  get_mapping()
  {
    if (!this->simulation_parameters.mortar_parameters.enable)
      return this->mapping;
    else
      return this->mapping_cache;
  }

  /**
   * @brief Rotate rotor mapping in mortar method
   *
   * @param[in] is_first Whether this is the first mapping rotation. This
   * boolean is used only for printing output purposes; since the mapping needs
   * to be rotated before setting up dofs, but also after the constraints are
   * first defined, this function is called twice before the iterate() loop. The
   * parameter is_first just prevents the mortar verbosity from being printed
   * twice
   */
  void
  rotate_rotor_mapping(const bool is_first);

  /**
   * @brief Update non-zero constraints if the boundary is time dependent.
   * Note: not implemented for the block fluid dynamics application.
   */
  void
  update_boundary_conditions();

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
   * @brief Checks and identifies if the cell is located in the constraining
   * domain defined by a plane with its outward-pointing normal vector. The
   * check is done through a scalar product between a vector formed by a vertex
   * of the cell to @plane_point and @plane_normal_vector. If one of the
   * vertices of the cell results in a negative scalar product result, the cell
   * is excluded from the constrained domain.
   *
   * @param[in] cell Pointer to an active cell of the fluid dynamics DoFHandler.
   *
   * @param[in] plane_point Coordinates of a point on the restriction plane for
   * the stasis constraint application domain.
   *
   * @param[in] plane_normal_vector Outward-pointing normal vector to define the
   * restriction plane for the stasis constraint application domain.
   *
   * @return Boolean indicating if the cell is in the domain of interest
   * (@p true) or not (@p false).
   */
  inline bool
  cell_in_constraining_domain(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const Point<dim>                                     &plane_point,
    const Tensor<1, dim>                                 &plane_normal_vector)
  {
    for (unsigned int v = 0; v < cell->n_vertices(); ++v)
      {
        Point<dim>     cell_vertex = cell->vertex(v);
        Tensor<1, dim> cell_vertex_to_plane_point_vector =
          plane_point - cell_vertex;
        double scalar_product_result =
          scalar_product(cell_vertex_to_plane_point_vector,
                         plane_normal_vector);
        if (scalar_product_result <= 0)
          return false;
      }
    return true;
  }

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
    const GlobalVectorType                               &temperature_solution,
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
    const GlobalVectorType &filtered_phase_fraction_solution,
    std::vector<double>    &local_filtered_phase_fraction_values)
  {
    const typename DoFHandler<dim>::active_cell_iterator vof_cell(
      &(*(this->triangulation)), cell->level(), cell->index(), dof_handler_vof);

    this->fe_values_vof->reinit(vof_cell);
    this->fe_values_vof->get_function_values(
      filtered_phase_fraction_solution, local_filtered_phase_fraction_values);
  }

  /**
   * @brief Identify if temperature DOFs of the cell are within the
   * constraining range. If they are, velocity DOFs of the cell are constrained.
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
  identify_cell_and_constrain_velocity(
    const std::vector<types::global_dof_index> &local_dof_indices,
    const std::vector<double>                  &local_temperature_values,
    StasisConstraintWithTemperature            &stasis_constraint_struct);

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
   * @brief Returns a vector of references to TableHandler objects that needs to be serialized/
   * deserialized for a given solver.
   *
   * @return Structure containing a vector of references to TableHandler objects that needs to be
   * serialized/deserialized for a given solver, and their corresponding file
   * names.
   */
  virtual std::vector<OutputStructTableHandler>
  gather_tables();

  /**
   * @brief Write the checkpoint
   */
  virtual void
  write_checkpoint();

  /**
   * @brief Gather and return vector of output structs that are particular to some applications.
   *
   * @return Vector of OutputStructs that will be used to write the output results as VTU files.
   */
  virtual std::vector<OutputStruct<dim, VectorType>>
  gather_output_hook();

  /**
   * @brief Gather solution information to generate output results
   *
   * @param[in] solution Vector of the present solution
   * @param[in,out] solution_output_structs Vector of OutputStructs that will be
   * used to write the output results as VTU files
   */
  void
  gather_output_results(
    const VectorType                           &solution,
    std::vector<OutputStruct<dim, VectorType>> &solution_output_structs);

  /**
   * @brief Generate post-processing parallel VTU files from vector of
   * solution_output_struct
   *
   * @param[in] solution Vector of present solution
   */
  void
  write_output_results(const VectorType &solution);

  /**
   * @brief Writes the forces per boundary condition to a text file output
   */
  void
  write_output_forces();

  /**
   * @brief Writes the torques per boundary condition to a text file output
   */
  void
  write_output_torques();

  /**
   * @brief This function is used to rescale pressure DOFs in the newton correction
   */
  void
  rescale_pressure_dofs_in_newton_update();

  /**
   * @brief This function initializes correctly a temporary vector depending on the
   * vector type
   */
  inline VectorType
  init_temporary_vector();

  // Member variables
protected:
  /**
   * @brief Verify consistency of the input parameters for boundary
   * conditions to ensure that for every boundary condition within the
   * triangulation, a boundary condition has been specified in the input file.
   */
  void
  verify_consistency_of_boundary_conditions()
  {
    // Sanity check all the boundary conditions of the triangulation to ensure
    // that they have a type.
    std::vector<types::boundary_id> boundary_ids_in_triangulation =
      this->triangulation->get_boundary_ids();
    for (auto const &boundary_id_in_tria : boundary_ids_in_triangulation)
      {
        AssertThrow(simulation_parameters.boundary_conditions.type.find(
                      boundary_id_in_tria) !=
                      simulation_parameters.boundary_conditions.type.end(),
                    FluidDynamicsBoundaryConditionMissing(boundary_id_in_tria));
      }
  }


  DofsType locally_owned_dofs;
  DofsType locally_relevant_dofs;

  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  std::shared_ptr<DoFHandler<dim>>                             dof_handler;
  std::shared_ptr<FESystem<dim>>                               fe;

  TimerOutput computing_timer;

  SimulationParameters<dim> simulation_parameters;
  PVDHandler                pvdhandler;
  PVDHandler                pvdhandler_boundary;

  // Functions used for source term and error analysis
  Function<dim>                 *exact_solution;
  std::shared_ptr<Function<dim>> forcing_function;

  // Dynamic flow control
  FlowControl<dim> flow_control;

  // Constraints for Dirichlet boundary conditions
  AffineConstraints<double> zero_constraints;
  AffineConstraints<double> nonzero_constraints;

  // Define whether manifolds are used for the normal vectors computation in
  // slip BCs
  bool use_manifold_for_normal;

  // Present solution and non-linear solution components
  VectorType                  evaluation_point;
  VectorType                  local_evaluation_point;
  VectorType                  newton_update;
  std::shared_ptr<VectorType> present_solution;
  VectorType                  system_rhs;

  // Previous solutions vectors
  std::shared_ptr<std::vector<VectorType>> previous_solutions;

  /**
   * @brief Structure that stores all SDIRK-related vectors used during the time integration process.
   */
  struct SDIRKVectors
  {
    /// Vector to hold the locally non-relevant part of the solution (for
    /// calculation purposes)
    VectorType locally_owned_for_calculation;

    /// Stores the previous k_j stage vectors (one per stage j < i)
    std::vector<VectorType> previous_k_j_solutions;

    /// Stores the sum of b_i * k_i across all stages
    VectorType sum_bi_ki;
    VectorType local_sum_bi_ki;

    /// Stores the sum of a_ij * k_j for j < i
    VectorType sum_over_previous_stages;
    VectorType local_sum_over_previous_stages;
  };

  /**
   * @brief Instance of SDIRK-related vectors.
   */
  SDIRKVectors sdirk_vectors;

  // Finite element order used
  const unsigned int velocity_fem_degree;
  const unsigned int pressure_fem_degree;
  unsigned int       number_quadrature_points;

  // Mappings and Quadratures
  std::shared_ptr<Mapping<dim>>        mapping;
  std::shared_ptr<Quadrature<dim>>     cell_quadrature;
  std::shared_ptr<Quadrature<dim - 1>> face_quadrature;

  // Mortar coupling manager and operator
  double                                         mortar_interface_radius;
  std::shared_ptr<MortarManagerCircle<dim>>      mortar_manager;
  std::shared_ptr<CouplingOperator<dim, double>> mortar_coupling_operator;
  std::shared_ptr<NavierStokesCouplingEvaluation<dim, double>>
    mortar_coupling_evaluator;

  // Mapping cache used in rotor mesh rotation in mortar method
  std::shared_ptr<MappingQCache<dim>> mapping_cache;

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
  std::vector<std::map<types::boundary_id, Tensor<1, dim>>>
                                             forces_on_boundaries;
  std::map<types::boundary_id, TableHandler> forces_tables;
  std::map<types::boundary_id, TableHandler> torques_tables;

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
