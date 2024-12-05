// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_algebraic_interface_reinitialization_h
#define lethe_vof_algebraic_interface_reinitialization_h

#include <solvers/physics_subequations_solver.h>
#include <solvers/vof_assemblers.h>
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
   * TODO AMISHGA
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
    , dof_handler(*triangulation)
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
   * TODO AMISHGA make something with this
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
   * TODO AMISHGA for all get methods
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
   * TODO AMISHGA
   * @param time_step_inv
   * @param tolerance
   * @return
   */
  inline bool
  continue_iterating(const double time_step_inv, const double tolerance)
  {
    auto solution_diff = present_solution;
    solution_diff -= previous_solution;

    double stop_criterion = time_step_inv * solution_diff.l2_norm();

    return (stop_criterion < tolerance);
  }



  const VOFSubequationsID        subequation_id;
  VOFSubequationsInterface<dim> *subequations_interface;
  MultiphysicsInterface<dim>
    *multiphysics_interface; // to get VOF DoFHandler and solution

  // Parameters
  const SimulationParameters<dim> &simulation_parameters;

  // VOF Simulation control
  const SimulationControl &simulation_control;

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
  GlobalVectorType               newton_update;
  GlobalVectorType               present_solution;
  GlobalVectorType               previous_solution; // Only used with bdf1
  GlobalVectorType               system_rhs;
  AffineConstraints<double>      zero_constraints;
  AffineConstraints<double>      nonzero_constraints;
  TrilinosWrappers::SparseMatrix system_matrix;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  // Verbosity
  const Parameters::Verbosity linear_solver_verbosity;
  const Parameters::Verbosity subequation_verbosity;
};

#endif
