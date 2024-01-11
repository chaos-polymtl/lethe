
#ifndef lethe_post_processors_smoothing_h
#define lethe_post_processors_smoothing_h

#include <core/vector.h>

#include <solvers/simulation_parameters.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_postprocessor.h>

#include <vector>

using namespace dealii;


/**
 * A base class that holds the dof_handler, fe_values and the simulation output
 * to calculate post-processed parameters
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the variable being smoothed is solved
 * @param triangulation Flow triangulation
 * @param simulation_parameters The simulation parameters
 * @param number_quadrature_points The number of quadrature points
 * @param mpi_communicator Allowing communication between cores
 */

template <int dim, typename VectorType>
class PostProcessorSmoothing
{
public:
  // Member functions
  PostProcessorSmoothing(
    const parallel::DistributedTriangulationBase<dim> &triangulation,
    const SimulationParameters<dim>                   &simulation_parameters,
    const unsigned int &number_quadrature_points);

  /**
   * @brief Generates the mass matrix, that is independent of the physics.
   */
  void
  generate_mass_matrix();

  /**
   * @brief Generates the right hand side based on the fluid's solution.
   */
  virtual void
  generate_rhs(const VectorType &,
               const DoFHandler<dim> &,
               std::shared_ptr<Mapping<dim>>)
  {}

  /**
   * @brief Solves the matrix system and outputs the smoothed field solution
   */
  const GlobalVectorType &
  solve_L2_projection();

  /**
   * @brief Returns the smoothed field solution.
   */
  const GlobalVectorType &
  calculate_smoothed_field(const VectorType             &solution,
                           const DoFHandler<dim>        &dof_handler_velocity,
                           std::shared_ptr<Mapping<dim>> mapping_velocity);

  /**
   * @brief Returns a constant reference to the dof_handler
   */
  const DoFHandler<dim> &
  get_dof_handler() const;

protected:
  FE_Q<dim>                                       fe_q;
  DoFHandler<dim>                                 dof_handler;
  SimulationParameters<dim>                       simulation_parameters;
  unsigned int                                    number_quadrature_points;
  std::shared_ptr<Mapping<dim>>                   mapping;
  std::shared_ptr<TrilinosWrappers::SparseMatrix> system_matrix;
  GlobalVectorType                                system_rhs;
  MPI_Comm                                        mpi_communicator;
  AffineConstraints<double>                       constraints;
  IndexSet                                        locally_owned_dofs;
  IndexSet                                        locally_active_dofs;
  IndexSet                                        locally_relevant_dofs;
  GlobalVectorType completely_distributed_solution;

private:
};

/**
 * A class that assembles the rhs and solves the L2projection of the Qcriterion
 * on nodes
 */
template <int dim, typename VectorType>
class QcriterionPostProcessorSmoothing
  : public PostProcessorSmoothing<dim, VectorType>
{
public:
  // Member functions
  QcriterionPostProcessorSmoothing(
    const parallel::DistributedTriangulationBase<dim> &triangulation,
    const SimulationParameters<dim>                   &simulation_parameters,
    const unsigned int &number_quadrature_points);

  /**
   * @brief Generates the right hand side based on the fluid's solution.
   *
   * @param solution The fluid's solution
   * @param dof_handler_fluid The dof_handler of the fluid solution
   * @param mapping_fluid The mapping of the fluid
   */
  void
  generate_rhs(const VectorType             &solution,
               const DoFHandler<dim>        &dof_handler_fluid,
               std::shared_ptr<Mapping<dim>> mapping_fluid) override;

private:
};


/**
 * A class that assembles the rhs and solves the L2projection of the continuity
 * equation (div u) at the nodes
 */
template <int dim, typename VectorType>
class ContinuityPostProcessorSmoothing
  : public PostProcessorSmoothing<dim, VectorType>
{
public:
  // Member functions
  ContinuityPostProcessorSmoothing(
    const parallel::DistributedTriangulationBase<dim> &triangulation,
    const SimulationParameters<dim>                   &simulation_parameters,
    const unsigned int &number_quadrature_points);

  /**
   * @brief Generates the right hand side based on the fluid dynamics solution.
   *
   * @param solution The fluid's solution
   * @param dof_handler_fluid The dof_handler of the fluid dynamics solution
   * @param mapping_fluid The mapping of the fluid
   */
  void
  generate_rhs(const VectorType             &solution,
               const DoFHandler<dim>        &dof_handler_fluid,
               std::shared_ptr<Mapping<dim>> mapping_fluid) override;

private:
};

#endif
