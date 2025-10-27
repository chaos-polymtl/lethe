// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_heat_transfer_h
#define lethe_heat_transfer_h

#include <core/bdf.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/heat_transfer_assemblers.h>
#include <solvers/heat_transfer_scratch_data.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/postprocessing_scalar.h>
#include <solvers/postprocessors.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>


DeclException1(
  HeatTransferBoundaryConditionMissing,
  types::boundary_id,
  << "The boundary id: " << arg1
  << " is defined in the triangulation, but not as a boundary condition for the heat transfer physics. Lethe does not assign a default boundary condition to boundary ids. Every boundary id defined within the triangulation must have a corresponding boundary condition defined in the input file.");


/**
 * @brief Implementation of heat transfer as an auxiliary physics. The heat
 * equation is weakly coupled to the velocity field. Equation solved:
 * \f$ \rho C_p \cdot \left( \frac{dT}{dt} + u \cdot \nabla T \right) =
 * k \cdot \nabla^2 T + \tau : \nabla u \f$
 *
 */

template <int dim>
class HeatTransfer : public AuxiliaryPhysics<dim, GlobalVectorType>
{
public:
  /**
   * @brief Constructor of the HeatTransfer object.
   *
   * @param multiphysics_interface Map of the auxiliary physics that will be
   * solved on top of a computational fluid dynamic simulation.
   *
   * @param p_simulation_parameters Contain the simulation parameter file
   * information.
   *
   * @param p_triangulation Contain the mesh information. In a
   * parallel::DistributedTriangulationBase<dim> not every detail may be known
   * on each processor. The mesh is distributed between the processors.
   *
   * @param p_simulation_control Object responsible for the control of
   * steady-state and transient simulations. Contains all the information
   * related to time stepping and the stopping criteria.
   *
   */
  HeatTransfer(MultiphysicsInterface<dim>      *multiphysics_interface,
               const SimulationParameters<dim> &p_simulation_parameters,
               std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                                  p_triangulation,
               std::shared_ptr<SimulationControl> p_simulation_control)
    : AuxiliaryPhysics<dim, GlobalVectorType>(
        p_simulation_parameters.non_linear_solver.at(PhysicsID::heat_transfer))
    , multiphysics(multiphysics_interface)
    , computing_timer(p_triangulation->get_mpi_communicator(),
                      this->pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , simulation_parameters(p_simulation_parameters)
    , triangulation(p_triangulation)
    , simulation_control(p_simulation_control)
    , dof_handler(*triangulation)
    , thermal_conductivity_models(
        p_simulation_parameters.physical_properties_manager
          .get_thermal_conductivity_vector())

  {
    if (simulation_parameters.mesh.simplex)
      {
        // for simplex meshes
        fe = std::make_shared<FE_SimplexP<dim>>(
          simulation_parameters.fem_parameters.temperature_order);
        temperature_mapping = std::make_shared<MappingFE<dim>>(*fe);
        cell_quadrature = std::make_shared<QGaussSimplex<dim>>(fe->degree + 1);
        face_quadrature =
          std::make_shared<QGaussSimplex<dim - 1>>(fe->degree + 1);
      }
    else
      {
        // Usual case, for quad/hex meshes
        fe = std::make_shared<FE_Q<dim>>(
          simulation_parameters.fem_parameters.temperature_order);
        temperature_mapping = std::make_shared<MappingQ<dim>>(fe->degree);
        cell_quadrature     = std::make_shared<QGauss<dim>>(fe->degree + 1);
        face_quadrature     = std::make_shared<QGauss<dim - 1>>(fe->degree + 1);
      }

    // Allocate solution transfer
    solution_transfer =
      std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(dof_handler);

    // Set size of previous solutions using BDF schemes information
    previous_solutions.resize(maximum_number_of_previous_solutions());

    // Prepare previous solutions transfer
    previous_solutions_transfer.reserve(previous_solutions.size());
    for (unsigned int i = 0; i < previous_solutions.size(); ++i)
      {
        previous_solutions_transfer.emplace_back(
          SolutionTransfer<dim, GlobalVectorType>(this->dof_handler));
      }

    // Change the behavior of the timer for situations when you don't want
    // outputs
    if (simulation_parameters.timer.type == Parameters::Timer::Type::none)
      this->computing_timer.disable_output();

    if (simulation_parameters.post_processing.calculate_average_temp_and_hf)
      {
        average_temperature =
          std::make_shared<AverageScalar<dim>>(this->dof_handler);
      }
  }

  /**
   * @brief Call for the assembly of the matrix and the right-hand side.
   *
   * @deprecated This function is to be deprecated when the new assembly mechanism
   * is integrated to this solver.
   */
  void
  assemble_matrix_and_rhs();

  /**
   * @brief Call for the assembly of the right-hand side.
   *
   * @deprecated This function is to be deprecated when the new assembly mechanism
   * is integrated to this solver.
   */
  void
  assemble_rhs();

  /**
   * @brief Call for the assembly of the matrix and the right-hand side of the Nitsche restriction for the heat transfert equation.
   *
   * @param assemble_matrix Boolean that is true for matrix assembly, and false for rhs assembly.
   */
  void
  assemble_nitsche_heat_restriction(const bool assemble_matrix);

  /**
   * @brief Gather and return vector of output structs that are particular to some applications.
   *
   * @return Vector of OutputStructs that will be used to write the output results as VTU files.
   */
  std::vector<OutputStruct<dim, GlobalVectorType>>
  gather_output_hook() override;

  /**
   * @brief Calculate delta_T_ref for the DCDD shock capture mechanism. delta_T_ref = T_max - T_min.
   *
   * @param minimum_delta_T_ref Minimum temperature value acceptable as reference to calculate the DCDD shock capture stabilization term.
   * DCDD shock capture elements are divided by delta_T_ref. Limiting
   * delta_T_ref prevents overestimation of the virtual diffusivity at
   * simulations with low differences between minimum and maximum temperatures.
   *
   * @return The difference between the maximal temperature and the minimal temperature in the whole domain.
   */
  double
  calculate_delta_T_ref(double minimum_delta_T_ref = 1.);

  /**
   * @brief Calculate the L2 error of the solution.
   */
  double
  calculate_L2_error();


  /**
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  void
  finish_simulation() override;

  /**
   * @brief Rearrange vector solution correctly for transient simulations.
   */
  void
  percolate_time_vectors() override;

  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  void
  postprocess(bool first_iteration) override;


  /**
   * @brief Prepare the auxiliary physics variables for a
   * mesh refinement/coarsening.
   */
  void
  pre_mesh_adaptation() override;

  /**
   * @brief Interpolate the auxiliary physics variables to the new mesh.
   */
  void
  post_mesh_adaptation() override;

  /**
   * @brief Compute the Kelly error estimator for mesh refinement.
   *
   * @param ivar The current element of the map simulation_parameters.mesh_adaptation.variables.
   *
   * @param estimated_error_per_cell The deal.II vector of estimated_error_per_cell.
   */
  void
  compute_kelly(const std::pair<const Variable,
                                Parameters::MultipleAdaptationParameters> &ivar,
                dealii::Vector<float> &estimated_error_per_cell) override;

  /**
   * @brief Prepare Heat Transfer to write checkpoint.
   */
  void
  write_checkpoint() override;

  /**
   * @brief Allow Heat Transfer to set up solution vector from checkpoint file.
   */
  void
  read_checkpoint() override;

  /**
   * @brief Returns a vector of references to TableHandler objects that needs to
   * be serialized/deserialized for the Heat Transfer solver.
   *
   * @return Structure containing a vector of references to TableHandler objects
   * that needs to be serialized/deserialized for the Heat Transfer solver, and
   * their corresponding file names.
   */
  std::vector<OutputStructTableHandler>
  gather_tables() override;

  /**
   * @brief Set up the DofHandler and the degrees of freedom associated with the physics.
   */
  void
  setup_dofs() override;

  /**
   * @brief Set up the initial conditions associated with the physics.
   * heat_transfer allows a temperature function initial condition over the
   * the domain.
   */
  void
  set_initial_conditions() override;

  /**
   * @brief Update non zero constraints if the boundary is time-dependent.
   */
  void
  update_boundary_conditions() override;

  /**
   * @brief Call for the solution of the linear system of equation using ILU
   * or GMRES.
   *
   * @param initial_step Provide the linear solver with indication if this
   * solution is the first one for the system of equation or not.
   *
   * @param renewed_matrix Indicate to the linear solve if the system matrix
   * has been recalculated or not.
   */
  void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) override;

  /**
   * @brief Getter method to access the private attribute dof_handler for the
   * physic currently solved. NB : dof_handler is now passed to the
   * multiphysics interface at the end of the setup_dofs method.
   *
   * @return A list of the indices of the degrees of freedom.
   */
  const DoFHandler<dim> &
  get_dof_handler() override
  {
    return dof_handler;
  }

  /**
   * @brief Getter method to access the private attribute evaluation_point for
   * the physic currently solved.
   *
   * @return The vector at which the evaluation is performed.
   */
  GlobalVectorType &
  get_evaluation_point() override
  {
    return evaluation_point;
  }

  /**
   * @brief Getter method to access the private attribute
   * local_evaluation_point for the physic currently solved.
   *
   * @return The local evaluation point. Ghosts cells are not considered in
   * this evaluation.
   */
  GlobalVectorType &
  get_local_evaluation_point() override
  {
    return local_evaluation_point;
  }

  /**
   * @brief Getter method to access the private attribute
   * newton_update for the physic currently solved.
   *
   * @return The direction used to perform the newton iteration.
   */
  GlobalVectorType &
  get_newton_update() override
  {
    return newton_update;
  }

  /**
   * @brief Getter method to access the private attribute
   * present_solution for the physic currently solved. NB : present_solution is
   * now passed to the multiphysics interface at the end of the setup_dofs
   * method.
   *
   * @return A vector containing all the values of the solution.
   */
  GlobalVectorType &
  get_present_solution() override
  {
    return present_solution;
  }

  /**
   * @brief Getter method to access the private attribute
   * system_rhs for the physic currently solved.
   *
   * @return Right hand side vector.
   */
  GlobalVectorType &
  get_system_rhs() override
  {
    return system_rhs;
  }

  /**
   * @brief Getter method to access the private attribute
   * nonzero_constraints for the physic currently solved.
   *
   * @return Store the nonzero constraints that arise from several sources such
   * as boundary conditions and hanging nodes in the mesh. See the deal.II
   * documentation on constraints on degrees of freedom for more information.
   */
  AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return nonzero_constraints;
  }

  /**
   * @brief Output the L2 and Linfty norms of the correction vector.
   *
   * @param[in] display_precision Number of outputted digits.
   */
  void
  output_newton_update_norms(const unsigned int display_precision) override
  {
    this->pcout << std::setprecision(display_precision)
                << "\t||dT||_L2 = " << std::setw(6) << newton_update.l2_norm()
                << std::setw(6)
                << "\t||dT||_Linfty = " << std::setprecision(display_precision)
                << newton_update.linfty_norm() << std::endl;
  }

  /**
   * @brief Get volume for residual normalization. By default, should return 1. In solvers, if normalize by volume is activated, the overriden method should return the global volume of the triangulation.
   *
   * @return Normalization volume.
   */
  virtual double
  get_residual_normalize_volume() const override
  {
    return simulation_parameters.non_linear_solver.at(PhysicsID::fluid_dynamics)
               .normalize_residual_by_volume ?
             GridTools::volume(*this->triangulation,
                               *this->temperature_mapping) :
             1.;
  }

private:
  /**
   * @brief Verify consistency of the input parameters for boundary
   * conditions to ensure that for every boundary condition within the
   * triangulation, a boundary condition has been specified in the input file.
   */
  void
  verify_consistency_of_boundary_conditions()
  {
    // Sanity check all of the boundary conditions of the triangulation to
    // ensure that they have a type.
    std::vector<types::boundary_id> boundary_ids_in_triangulation =
      this->triangulation->get_boundary_ids();
    for (auto const &boundary_id_in_tria : boundary_ids_in_triangulation)
      {
        AssertThrow(simulation_parameters.boundary_conditions_ht.type.find(
                      boundary_id_in_tria) !=
                      simulation_parameters.boundary_conditions_ht.type.end(),
                    HeatTransferBoundaryConditionMissing(boundary_id_in_tria));
      }
  }

  /**
   *  @brief Assemble the matrix associated with the solver
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the rhs associated with the solver
   */
  void
  assemble_system_rhs() override;


  /**
   * @brief Assemble the local matrix for a given cell.
   *
   * Use the WorkStream class to assemble the system matrix. It is
   * a thread safe function.
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for HeatTransferScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  virtual void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    HeatTransferScratchData<dim>                         &scratch_data,
    StabilizedMethodsCopyData                            &copy_data);

  /**
   * @brief Assemble the local rhs for a given cell
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for HeatTransferScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  virtual void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    HeatTransferScratchData<dim>                         &scratch_data,
    StabilizedMethodsCopyData                            &copy_data);

  /**
   * @brief Set up the vector of assembler functions
   */
  virtual void
  setup_assemblers();


  /**
   * @brief Copy local cell information to global matrix
   */

  virtual void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsCopyData &copy_data);

  /**
   * @brief Copy local cell rhs information to global rhs
   */

  virtual void
  copy_local_rhs_to_global_rhs(const StabilizedMethodsCopyData &copy_data);

  /**
   * @brief Post-processing.
   * Calculate temperature statistics on the domain : Max, min, average and
   * standard-deviation.
   *
   * @param gather_vof Boolean true when VOF=true (multiphase flow), used to gather
   * VOF information.
   *
   * @param monitored_fluid Fluid indicator (fluid0 or fluid1 or both) corresponding
   * to the phase of interest.
   *
   * @param domain_name String indicating the monitored_fluid in the output filename.
   *
   * @param time_average Boolean true when calculating the spacial average of the time-averaged temperature solution.
   */

  void
  postprocess_temperature_statistics(
    const bool                       gather_vof,
    const Parameters::FluidIndicator monitored_fluid,
    const std::string                domain_name,
    const bool                       time_average);

  /**
   * @brief Post-processing. Write the temperature statistics to an output file.
   *
   * @param domain_name string indicating the postprocessed_fluid in the
   * console output, table and filename.
   */

  void
  write_temperature_statistics(const std::string domain_name);

  /**
   * @brief Post-processing.
   * Calculate liquid fraction on the domain.
   *
   * @param gather_vof boolean true when VOF=true (multiphase flow), used to gather
   * VOF information
   */

  void
  postprocess_liquid_fraction(const bool gather_vof);

  /**
   * @brief Post-processing. Write the liquid fraction to an output file.
   */

  void
  write_liquid_fraction();

  /**
   * Post-processing. Calculate the heat flux at heat transfer boundary
   * conditions.
   *
   * @param gather_vof boolean true when VOF=true (multiphase flow), used to gather
   * VOF information
   *
   * @param current_solution_fd current solution for the fluid dynamics, parsed
   * by postprocess
   */

  template <typename VectorType>
  void
  postprocess_heat_flux_on_bc(const bool        gather_vof,
                              const VectorType &current_solution_fd);

  /**
   * Post-processing. Calculate the heat flux in the Nitsche immersed boundary.
   *
   */
  void
  postprocess_heat_flux_on_nitsche_ib();

  /**
   * Post-processing. Calculate the thermal energy \f$ (\rho \cdot C_p \cdot T)
   * \f$ in a fluid domain.
   *
   * @param gather_vof Boolean true when VOF=true (multiphase flow), used to gather
   * VOF information.
   *
   * @param monitored_fluid Fluid indicator (fluid0 or fluid1 or both) corresponding
   * to the phase of interest.
   *
   * @param domain_name String indicating the monitored_fluid in the console output,
   * table and filename.
   *
   * @param current_solution_fd Current solution for the fluid dynamics, parsed
   * by postprocess.
   */

  template <typename VectorType>
  void
  postprocess_thermal_energy_in_fluid(
    const bool                       gather_vof,
    const Parameters::FluidIndicator monitored_fluid,
    const std::string                domain_name,
    const VectorType                &current_solution_fd);

  /**
   * @brief Post-processing. Write the heat transfer values to an output file.
   *
   * @param domain_name string indicating the postprocessed_fluid in the
   * console output, table and filename.
   */

  void
  write_heat_flux(const std::string domain_name);

  /**
   * @brief Set a phase_coefficient that is used in the postprocessing, to
   * account for multiphase flow.
   *
   * Returns
   * - a double: phase_coefficient from 0 to 1, see the documentation on VOF for
   * more information,
   * - a boolean: point_is_in_postprocessed_fluid, true if the quadrature point
   * is in the postprocessed_fluid (used for min_max_temperature calculation)
   *
   * @param gather_vof boolean true when VOF=true (multiphase flow), used to gather
   * VOF information
   *
   * @param monitored_fluid Fluid indicator (fluid0 or fluid1 or both) corresponding
   * to the phase of interest.
   *
   * @param phase_value_q double corresponding to the phase value at this quadrature point
   */

  std::pair<double, bool>
  set_phase_coefficient(const bool                       gather_vof,
                        const Parameters::FluidIndicator monitored_fluid,
                        const double                     phase_value_q);


  MultiphysicsInterface<dim> *multiphysics;

  /**
   * @brief Store information related to the computing time such as CPU times or
   * wall time.
   */
  TimerOutput computing_timer;

  /**
   * @brief Contain the simulation parameter file information.
   */
  const SimulationParameters<dim> &simulation_parameters;


  // Core elements for the heat transfer simulation

  /**
   * @brief Collection of cells that cover the domain on which one wants to
   * solve a partial differential equation.
   */
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  /**
   * @brief Responsible for the control of steady-state and transient
   * simulations. Contains all the information related to time stepping and the
   * stopping criteria. See simulation_control abstract class for more
   * information.
   */
  std::shared_ptr<SimulationControl> simulation_control;
  /**
   * @brief Given a triangulation and a description of a finite element, this
   * class enumerates degrees of freedom on all vertices, edges, faces, and
   * cells of the triangulation.
   */
  DoFHandler<dim> dof_handler;
  /**
   * @brief The base class for finite element.
   */
  std::shared_ptr<FiniteElement<dim>> fe;

  /**
   * @brief Store some convergence data, such as residuals of the cg-method,
   * or some evaluated <i>L<sup>2</sup></i>-errors of discrete solutions.
   * Evaluate convergence rates or orders.
   */
  ConvergenceTable error_table;

  // Mapping and Quadrature

  /**
   * @brief Transformation which maps point in the reference cell to
   * points in the actual grid cell.
   */
  std::shared_ptr<Mapping<dim>> temperature_mapping;
  /**
   * @brief Approximate an integral by evaluating the integrand at specific
   * points and summing the point values with specific weights.
   */
  std::shared_ptr<Quadrature<dim>> cell_quadrature;
  /**
   * @brief Approximate an integral by evaluating the integrand at specific
   * points and summing the point values with specific weights.
   */
  std::shared_ptr<Quadrature<dim - 1>> face_quadrature;


  // Solution storage:

  /**
   * @brief Store a subset of the DoF indices that are stored on a
   * particular processor on a distributed architecture.
   */
  IndexSet locally_owned_dofs;
  /**
   * @brief Store a subset of the DoF indices that contains both the
   * locally_owned_dofs and the DoF indices on ghost cells.
   */
  IndexSet locally_relevant_dofs;

  /**
   * @brief The point at which the evaluation is performed.
   */
  GlobalVectorType evaluation_point;
  /**
   * @brief The local evaluation point. Ghosts cells are not considered in
   * this evaluation.
   */
  GlobalVectorType local_evaluation_point;
  /**
   * @brief The direction \f$ \delta u^n \f$ used to perform the newton iteration
   * in \f$ u^{n+1} = u^n + \alpha ^n \delta u ^n \f$. Where \f$ \alpha ^n \f$
   * is the size of the step.
   */
  GlobalVectorType newton_update;
  /**
   * @brief A vector containing all the values of the solution.
   */
  GlobalVectorType present_solution;
  /**
   * @brief The right hand side vector.
   */
  GlobalVectorType system_rhs;
  /**
   * @brief Store the nonzero constraints that arise from several sources such
   * as boundary conditions and hanging nodes in the mesh. See the deal.II
   * documentation on constraints on degrees of freedom for more information.
   */
  AffineConstraints<double> nonzero_constraints;
  AffineConstraints<double> zero_constraints;
  /**
   * @brief The system matrix.
   */
  TrilinosWrappers::SparseMatrix system_matrix;


  /**
   * @brief Previous solution vector.
   */
  std::vector<GlobalVectorType> previous_solutions;

  /**
   * @brief SolutionTransfer<dim, GlobalVectorType>> is
   * used to implement the transfer of a discrete FE function
   * (e.g. a solution vector) from one mesh to another. This Deal.ii class is
   * used for mesh_refinement and simulation restarts.
   */
  std::shared_ptr<SolutionTransfer<dim, GlobalVectorType>> solution_transfer;

  /**
   * @brief previous_solutions_transfer occupies the same role as
   * solution_transfer but for solutions at previous time steps.
   */
  std::vector<SolutionTransfer<dim, GlobalVectorType>>
    previous_solutions_transfer;

  /**
   * @brief Reference for GGLS https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.2324
   * Warning, this GGLS implementation is valid only for Linear elements
   * Quad elements will be lacking the third derivative of the diffusion
   * operator Whether this affects or not the final result is unclear to me at
   * the moment. Additionnaly, this formulation does not use the gradient of the
   * source term. The same applies, I have no clue if this is detrimental or not
   * to the solution since anyway the GGLS term scales as h^(order+1).
   */
  const bool GGLS = true;

  /**
   * @brief Assemblers for the matrix and the right hand side (RHS).
   */
  std::vector<std::shared_ptr<HeatTransferAssemblerBase<dim>>> assemblers;

  /**
   * @brief Temperature statistics table used for post-processing. It contains
   * the minimum , maximum, average and standard deviation of the temperature.
   */
  TableHandler statistics_table;

  /**
   * @brief Contain heat flux information used for post-processing:
   * - the total fluxes on a boundary: \f$(-k \nabla T + u \rho C_p) \cdot n \f$
   * - the convective heat flux on a boundary: \f$ h(T-T_{inf}) \f$
   * - the total fluxes on the nitsche immersed boundaries (if active)
   */
  TableHandler heat_flux_table;

  // Heat flux postprocessing

  /**
   * @brief Abstract class that allows to calculate the thermal conductivity on
   * each quadrature point using the temperature of the fluid.
   */
  std::vector<std::shared_ptr<ThermalConductivityModel>>
    thermal_conductivity_models;

  /**
   * @brief Compute the average temperature in time.
   */
  std::shared_ptr<AverageScalar<dim>> average_temperature;

  /**
   * @brief Locally owned average temperature calculated using the AverageScalar object.
   */
  GlobalVectorType average_temperature_to_output;

  /**
   * @brief Compute local heat flux quantities from temperature field and
   * material properties.
   */
  std::vector<HeatFluxPostprocessor<dim>> heat_flux_postprocessors;

  /**
   * @brief Compute local average heat flux quantities from temperature field and
   * material properties.
   */
  std::vector<HeatFluxPostprocessor<dim>> average_heat_flux_postprocessors;

  /*
   * Phase change post-processing. These parameters track the presence of a
   * phase change physical property and the associated post-processing
   * information
   */

  /**
   * @brief Liquid fraction in the domain.
   */
  TableHandler liquid_fraction_table;
};


#endif
