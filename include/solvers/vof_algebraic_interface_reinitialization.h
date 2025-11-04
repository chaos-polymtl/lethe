// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_algebraic_interface_reinitialization_h
#define lethe_vof_algebraic_interface_reinitialization_h

#include <core/pvd_handler.h>

#include <solvers/physics_subequations_solver.h>
#include <solvers/vof_scratch_data.h>
#include <solvers/vof_subequations_interface.h>

#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>

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
   * @brief Constructor of the VOF algebraic interface reinitialization.
   *
   * @param[in] p_simulation_parameters Simulation parameters.
   *
   * @param[in] p_pcout Parallel cout used to print the information.
   *
   * @param[in] p_triangulation Distributed mesh information.
   *
   * @param[in] p_simulation_control SimulationControl object.
   *
   * @param[in, out] p_subequations_interface Subequations interface object used
   * to get information from other subequations and store information from the
   * current one.
   */
  VOFAlgebraicInterfaceReinitialization(
    const SimulationParameters<dim> &p_simulation_parameters,
    const ConditionalOStream        &p_pcout,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                              p_triangulation,
    const std::shared_ptr<SimulationControl> &p_simulation_control,
    VOFSubequationsInterface<dim>            &p_subequations_interface)
    : PhysicsNonlinearSubequationsSolver<dim, GlobalVectorType>(
        p_simulation_parameters.vof_subequations_non_linear_solvers.at(
          VOFSubequationsID::algebraic_interface_reinitialization),
        p_pcout)
    , subequation_id(VOFSubequationsID::algebraic_interface_reinitialization)
    , subequations_interface(p_subequations_interface)
    , simulation_parameters(p_simulation_parameters)
    , simulation_control(p_simulation_control)
    , triangulation(p_triangulation)
    , dof_handler(std::make_shared<DoFHandler<dim>>(*this->triangulation))
    , subequation_verbosity(p_simulation_parameters.multiphysics.vof_parameters
                              .regularization_method.verbosity)
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

    // Ensure that the shared pointer is properly allocated
    this->present_solution = std::make_shared<GlobalVectorType>();
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
   * @brief Solve interface algebraic reinitialization process until one of the
   * stop criteria is met.
   */
  void
  solve() override;

  /**
   * @brief Getter method to access the private attribute evaluation_point for
   * the subequation currently solved.
   *
   * @return The vector in which the evaluation is performed.
   */
  GlobalVectorType &
  get_evaluation_point() override
  {
    return this->evaluation_point;
  }

  /**
   * @brief Getter method to access the private attribute
   * local_evaluation_point for the subequation currently solved.
   *
   * @return The local evaluation point. Ghosts cells are not considered in
   * this evaluation.
   */
  GlobalVectorType &
  get_local_evaluation_point() override
  {
    return this->local_evaluation_point;
  }

  /**
   * @brief Getter method to access the private attribute
   * newton_update for the subequation currently solved.
   *
   * @return The direction used to perform the newton iteration.
   */
  GlobalVectorType &
  get_newton_update() override
  {
    return this->newton_update;
  }

  /**
   * @brief Getter method to access the private attribute
   * present_solution for the subequation currently solved.
   *
   * @return Vector containing all the values of the solution.
   */
  GlobalVectorType &
  get_present_solution() override
  {
    return *this->present_solution;
  }

  /**
   * @brief Getter method to access the private attribute
   * system_rhs for the subequation currently solved.
   *
   * @return Right-hand side vector.
   */
  GlobalVectorType &
  get_system_rhs() override
  {
    return this->system_rhs;
  }

  /**
   * @brief Getter method to access the private attribute
   * nonzero_constraints for the subequation currently solved.
   *
   * @return Container of nonzero constraints that arise from several sources such
   * as boundary conditions and hanging nodes in the mesh. See the deal.II
   * documentation on constraints on degrees of freedom for more information.
   */
  AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return this->nonzero_constraints;
  }


  /**
   * @brief Output the \f$L_2\f$ and \f$L_\infty\f$ norms of the correction vector.
   *
   * @param[in] display_precision Number of outputted digits.
   */
  void
  output_newton_update_norms(const unsigned int display_precision) override
  {
    this->pcout << std::setprecision(display_precision)
                << "\t||dphi_reinit||_L2 = " << std::setw(6)
                << this->newton_update.l2_norm() << std::setw(6)
                << "\t||dphi_reinit||_Linfty = "
                << std::setprecision(display_precision)
                << this->newton_update.linfty_norm() << std::endl;
  }

  /**
   * @brief Return the metric for residual normalization. By default, should return 1.
   * If the normalize_residual_by_volume is set to true, the method
   * returns the global volume of the triangulation.
   *
   * @return Normalization metric.
   */
  double
  get_residual_normalization_metric() const override
  {
    return simulation_parameters.non_linear_solver.at(PhysicsID::VOF)
               .normalize_residual_by_volume ?
             GridTools::volume(*this->triangulation, *this->mapping) :
             1.;
  }

private:
  /**
   * @brief Write parallel VTU files of quantities of interest for
   * the algebraic interface reinitialization process.
   *
   * @param[in] step Integer indicating the reinitialization step number.
   *
   * @note Only the reinitialization steps of the last global time-step
   * (simulation time-step) ran are outputted.
   */
  void
  write_output_results(const unsigned int step);

  /**
   * @brief Computes mesh-dependant diffusivity coefficient value.
   *
   * @param[in] min_cell_size Smallest cell's measure.
   *
   * @return Constant diffusivity coefficient value for a given mesh.
   */
  inline double
  compute_diffusivity(const double min_cell_size) const
  {
    const double multiplier =
      this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.algebraic_interface_reinitialization
        .diffusivity_multiplier;
    const double power =
      this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.algebraic_interface_reinitialization
        .diffusivity_power;
    return multiplier * std::pow(min_cell_size, power);
  }

  /**
   * @brief Computes algebraic reinitialization time-step value with
   * user-imposed CFL and minimum cell size.
   *
   * @return Algebraic reinitialization time-step.
   */
  inline double
  compute_time_step()
  {
    // Get CFL value
    const double cfl =
      this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.algebraic_interface_reinitialization
        .reinitialization_cfl;

    // Get the minimum cell size
    const double h_min =
      identify_minimum_cell_size(*this->mapping,
                                 *this->dof_handler,
                                 *this->cell_quadrature,
                                 this->triangulation->get_mpi_communicator());

    return h_min * cfl;
  }

  /**
   * @brief Define the zero constraints used to solve the problem.
   */
  void
  define_zero_constraints();

  /**
   * @brief Define the non zero constraints used to solve the problem.
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
   * @param[in] initial_step Provides the linear solver with indication if this
   * solution is the first one for the system of equation or not.
   *
   * @param[in] renewed_matrix Indicates to the linear solve if the system
   * matrix has been recalculated or not.
   */
  void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) override;

  /**
   * @brief Indicate if the algebraic reinitialization should continue by
   * checking if at least one of the two stop criteria (steady-state criterion
   * or maximum number of reinitialization steps) is met.
   *
   * @param[in] time_step_inv Inverse of the current reinitialization time-step.
   *
   * @param[in] step_number Algebraic interface reinitialization step number.
   *
   * @return Boolean indicating if the algebraic reinitialization should
   * continue
   */
  inline bool
  continue_iterating(const double time_step_inv, const unsigned int step_number)
  {
    if (step_number == 1)
      {
        return (step_number <
                this->simulation_parameters.multiphysics.vof_parameters
                    .regularization_method.algebraic_interface_reinitialization
                    .max_steps_number +
                  1);
      }
    else
      {
        // Get the stop criterion of the pseudo-time-stepping scheme
        double steady_state_criterion =
          this->simulation_parameters.multiphysics.vof_parameters
            .regularization_method.algebraic_interface_reinitialization
            .steady_state_criterion;

        // Evaluate the solution difference between the 2 last solutions
        auto solution_diff = local_evaluation_point;
        solution_diff -= previous_local_evaluation_point;

        // Evaluate the current steady-state criterion value
        double stop_criterion = time_step_inv * solution_diff.l2_norm();

        if (this->subequation_verbosity == Parameters::Verbosity::extra_verbose)
          {
            this->pcout
              << "Algebraic reinitialization solution norm difference = "
              << solution_diff.l2_norm() << std::endl;
            this->pcout
              << "Algebraic reinitialization steady-state criterion value = "
              << stop_criterion << std::endl;
            this->pcout
              << "Algebraic reinitialization fixed steady-state criterion = "
              << steady_state_criterion << std::endl;
          }

        return ((stop_criterion > steady_state_criterion) &&
                (step_number <
                 this->simulation_parameters.multiphysics.vof_parameters
                     .regularization_method.algebraic_interface_reinitialization
                     .max_steps_number +
                   1));
      }
  }

  /**
   * @brief Check if the phase gradient and curvature L2 projections have been
   * solved.
   */
  void
  check_dependencies_validity();

  const VOFSubequationsID        subequation_id;
  VOFSubequationsInterface<dim> &subequations_interface;

  // Parameters
  const SimulationParameters<dim> &simulation_parameters;

  // Simulation control object for simulation iteration number
  std::shared_ptr<SimulationControl> simulation_control;

  // Handler to output algebraic reinitialization steps
  PVDHandler pvdhandler;

  // Time-stepping with BDF1
  double              current_time_step;
  std::vector<double> time_step_vector = {0.0};

  // Core elements
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  std::shared_ptr<DoFHandler<dim>>                             dof_handler;
  std::shared_ptr<FiniteElement<dim>>                          fe;

  // Mapping and Quadrature
  std::shared_ptr<Mapping<dim>>    mapping;
  std::shared_ptr<Quadrature<dim>> cell_quadrature;

  // Solution storage
  IndexSet                          locally_owned_dofs;
  IndexSet                          locally_relevant_dofs;
  GlobalVectorType                  evaluation_point;
  GlobalVectorType                  local_evaluation_point;
  GlobalVectorType                  previous_local_evaluation_point;
  GlobalVectorType                  newton_update;
  std::shared_ptr<GlobalVectorType> present_solution;
  GlobalVectorType                  previous_solution; // Only used with bdf1
  GlobalVectorType                  system_rhs;
  AffineConstraints<double>         zero_constraints;
  AffineConstraints<double>         nonzero_constraints;
  TrilinosWrappers::SparseMatrix    system_matrix;

  // Verbosity
  const Parameters::Verbosity subequation_verbosity;
};

#endif
