// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_algebraic_interface_reinitialization_h
#define lethe_vof_algebraic_interface_reinitialization_h

#include <core/pvd_handler.h> // for debugging

#include <solvers/physics_subequations_solver.h>
#include <solvers/vof_assemblers.h>
#include <solvers/vof_filter.h> // for debugging
#include <solvers/vof_scratch_data.h>
#include <solvers/vof_subequations_interface.h>

#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

/**
 * @brief VOF algebraic reinitialization solver.
 *
 * @tparam dim Number of dimensions of the problem.
 */
template <int dim>
class VOFAlgebraicInterfaceReinitialization
  : public PhysicsNonlinearSubequationsSolver<dim, GlobalVectorType>
{
public:
  /**
   * TODO AA
   * @brief Constructor of the VOF algebraic interface reinitialization.
   *
   * @param[in] p_simulation_parameters Simulation parameters.
   *
   * @param[in] p_pcout Parallel cout used to print the information.
   *
   * @param[in] p_simulation_control Object containing simulation time-stepping
   * information.
   *
   * @param[in] p_triangulation Distributed mesh information.
   *
   * @param[in] p_multiphysics_interface Multiphysics interface object used to
   * get information from physics.
   *
   * @param[in, out] p_subequations_interface Subequations interface object used
   * to get information from other subequations and store information from the
   * current one.
   */
  VOFAlgebraicInterfaceReinitialization(
    const SimulationParameters<dim> &p_simulation_parameters,
    const ConditionalOStream        &p_pcout,
    const SimulationControl         &p_simulation_control,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                  &p_triangulation,
    MultiphysicsInterface<dim>    *p_multiphysics_interface,
    VOFSubequationsInterface<dim> *p_subequations_interface)
    : PhysicsNonlinearSubequationsSolver<dim, GlobalVectorType>(
        p_simulation_parameters.non_linear_solver.at(PhysicsID::VOF),
        p_pcout)
    , subequation_id(VOFSubequationsID::algebraic_interface_reinitialization)
    , subequations_interface(p_subequations_interface)
    , multiphysics_interface(p_multiphysics_interface)
    , simulation_parameters(p_simulation_parameters)
    , simulation_control(p_simulation_control)
    , triangulation(p_triangulation)
    , dof_handler(*this->triangulation)
    , linear_solver_verbosity(
        p_simulation_parameters.linear_solver.at(PhysicsID::VOF).verbosity)
    , subequation_verbosity(p_simulation_parameters.multiphysics.vof_parameters
                              .algebraic_interface_reinitialization.verbosity)
  {
    if (this->simulation_parameters.mesh.simplex)
      {
        // For simplex meshes
        this->fe = std::make_shared<FE_SimplexP<dim>>(
          this->simulation_parameters.fem_parameters.VOF_order);
        this->mapping = std::make_shared<MappingFE<dim>>(*this->fe);
        this->cell_quadrature =
          std::make_shared<QGaussSimplex<dim>>(this->fe->degree + 1);
      }
    else
      {
        // Usual case, for quad/hex meshes
        this->fe = std::make_shared<FE_Q<dim>>(
          this->simulation_parameters.fem_parameters.VOF_order);
        this->mapping = std::make_shared<MappingQ<dim>>(this->fe->degree);
        this->cell_quadrature =
          std::make_shared<QGauss<dim>>(this->fe->degree + 1);
      }
  }

  /**
   * @brief Default destructor.
   */
  ~VOFAlgebraicInterfaceReinitialization() = default;


  /**
   * @brief Set up the DofHandler and the degree of freedom associated with
   * the subequation.
   */
  void
  setup_dofs() override;


  /**
   * TODO AA make something with this
   * @brief Assemble and solve linear system when the equation to solve is
   * linear without using the non-linear solver interface.
   *
   * @param[in] is_post_mesh_adaptation Indicates if the equation is being
   * solved during post_mesh_adaptation(), for verbosity.
   */
  void
  solve(const bool & /*is_post_mesh_adaptation = false*/) override;

  /**
   * @brief Getter methods to get the private attributes for the physic currently solved
   * NB : dof_handler and present_solution are passed to the multiphysics
   * interface at the end of the setup_dofs method
   */

  /**
   * TODO AA for all get methods
   * @return
   */
  GlobalVectorType &
  get_evaluation_point() override
  {
    return this->evaluation_point;
  }

  GlobalVectorType &
  get_local_evaluation_point() override
  {
    return this->local_evaluation_point;
  }

  GlobalVectorType &
  get_newton_update() override
  {
    return this->newton_update;
  }

  GlobalVectorType &
  get_present_solution() override
  {
    return this->present_solution;
  }

  GlobalVectorType &
  get_system_rhs() override
  {
    return this->system_rhs;
  }

  AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return this->nonzero_constraints;
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
                << "\t||dphi_reinit||_L2 = " << std::setw(6)
                << newton_update.l2_norm() << std::setw(6)
                << "\t||dphi_reinit||_Linfty = "
                << std::setprecision(display_precision)
                << newton_update.linfty_norm() << std::endl;
  }

private:
  // functions for condition number calculation


  /**
   * TODO AA
   * @brief write_output_results
   * Post-processing as parallel VTU files
   */
  void
  write_output_results(const unsigned int it);

  /**
   * TODO AA
   * @return
   */
  inline double
  identify_minimum_cell_size() const
  {
    // Initialize FEValues for interface algebraic reinitialization
    FEValues<dim> fe_values(*this->mapping,
                            this->dof_handler.get_fe(),
                            *this->cell_quadrature,
                            update_JxW_values);

    // Initialize cell diameter
    double h = DBL_MAX;

    // Element degree
    double degree = double(fe_values.get_fe().degree);

    for (const auto &cell : this->dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            // Compute cell diameter
            double cell_measure =
              compute_cell_measure_with_JxW(fe_values.get_JxW_values());
            double h_local = compute_cell_diameter<dim>(cell_measure, degree);

            // Update cell diameter to minimum value
            h = std::min(h, h_local);
          }
      }
    return h;
  }

  /**
   * TODO AA do model casting?
   * @param[in] cell_size
   * @return
   */
  inline double
  compute_diffusivity(const double cell_size) const
  {
    if (this->simulation_parameters.multiphysics.vof_parameters
          .algebraic_interface_reinitialization.diffusivity_type ==
        Parameters::ReinitializationDiffusivityType::constant)
      {
        return this->simulation_parameters.multiphysics.vof_parameters
          .algebraic_interface_reinitialization.diffusivity_value;
      }
    else // Mesh-dependant diffusivity model
      {
        const double multiplier =
          this->simulation_parameters.multiphysics.vof_parameters
            .algebraic_interface_reinitialization.diffusivity_multiplier;
        const double power =
          this->simulation_parameters.multiphysics.vof_parameters
            .algebraic_interface_reinitialization.diffusivity_power;
        return multiplier * std::pow(cell_size, power);
      }
  }

  /**
   * TODO AA
   * @return
   */
  inline double
  compute_time_step()
  {
    // Get CFL value
    const double cfl =
      this->simulation_parameters.multiphysics.vof_parameters
        .algebraic_interface_reinitialization.reinitialization_cfl;

    const double h_min = identify_minimum_cell_size();
    //    h = std::min(h_min, vof_time_step);
    return h_min * cfl;
  }

  /**
   * @brief Define the zero constraints used to solve the problem.
   */
  void
  define_zero_constraints();

  /**
   * @brief  Define the non zero constraints used to solve the problem.
   */
  void
  define_non_zero_constraints();

  /**
   * @brief Sets-up the initial conditions associated with the subequation.
   */
  void
  set_initial_conditions();

  /**
   * @brief Assemble the matrix.
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the right-hand side (rhs).
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief Solve the linear system.
   *
   * @param initial_step Provides the linear solver with indication if this
   * solution is the first one for the system of equation or not.
   *
   * @param renewed_matrix Indicates to the linear solve if the system matrix
   * has been recalculated or not.
   */
  void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) override;

  /**
   * TODO AA
   * @param time_step_inv
   * @param tolerance
   * @return
   */
  inline bool
  continue_iterating(const double time_step_inv, const double tolerance)
  {
    // TODO AA cleanup
    auto solution_diff = local_evaluation_point;
    solution_diff -= previous_local_evaluation_point;

    double stop_criterion = time_step_inv * solution_diff.l2_norm();

    this->pcout << "TIME-STEP INV = " << time_step_inv << std::endl;
    this->pcout << "SOLUTION DIFF NORM = " << solution_diff.l2_norm()
                << std::endl;
    this->pcout << "STOP CRITERION = " << stop_criterion << std::endl;
    this->pcout << "TOLERANCE = " << tolerance << std::endl;

    return (stop_criterion > tolerance);
  }

  const VOFSubequationsID        subequation_id;
  VOFSubequationsInterface<dim> *subequations_interface;
  MultiphysicsInterface<dim>
    *multiphysics_interface; // to get VOF DoFHandler and solution

  // Parameters
  const SimulationParameters<dim> &simulation_parameters;
  PVDHandler                       pvdhandler;

  // VOF Simulation control
  const SimulationControl
        &simulation_control; // TODO AA Only needed for CFL value... work on it
  double current_time_step;
  std::vector<double> time_step_vector = {
    0}; // TODO AA resee how it is initialized

  // Core elements
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  DoFHandler<dim>                                              dof_handler;
  std::shared_ptr<FiniteElement<dim>>                          fe;

  // Mapping and Quadrature
  std::shared_ptr<Mapping<dim>>    mapping;
  std::shared_ptr<Quadrature<dim>> cell_quadrature;

  // Solution storage
  IndexSet                       locally_owned_dofs;
  IndexSet                       locally_relevant_dofs;
  GlobalVectorType               evaluation_point;
  GlobalVectorType               local_evaluation_point;
  GlobalVectorType               previous_local_evaluation_point;
  GlobalVectorType               newton_update;
  GlobalVectorType               present_solution;
  GlobalVectorType               previous_solution; // Only used with bdf1
  GlobalVectorType               system_rhs;
  AffineConstraints<double>      zero_constraints;
  AffineConstraints<double>      nonzero_constraints;
  TrilinosWrappers::SparseMatrix system_matrix;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  // Diffusivity storage TODO AA for non constant diffusivity
  //  GlobalVectorType diffusivity_values;

  // Verbosity
  const Parameters::Verbosity linear_solver_verbosity;
  const Parameters::Verbosity subequation_verbosity;
};

#endif
