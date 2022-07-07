
#ifndef lethe_post_processors_smoothing_h
#define lethe_post_processors_smoothing_h

// DEALII INCLUDES
// Base
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>

// Numerics
#include <deal.II/numerics/data_postprocessor.h>

// Rheological models
#include <core/rheological_model.h>

#include <solvers/simulation_parameters.h>


// standard library includes includes
#include <vector>

using namespace dealii;


/**
 * A base class that holds the dof_handler, fe_values and the simulation output
 * to calculate post-processed parameters
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @tparam triangulation Flow triangulation
 * @tparam simulation_parameters The simulation parameters
 * @tparam number_quadrature_points The number of quadrature points
 */

template <int dim>
class PostProcessorSmoothing
{
public:
  // Member functions
  PostProcessorSmoothing(
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
    SimulationParameters<dim> simulation_parameters,
    unsigned int              number_quadrature_points);

  void
  generate_mass_matrix();
  /**
   * @brief Outputs a solution for the field on the nodes.
   */
  virtual void
  evaluate_smoothed_field()
  {}

private:
  FE_Q<dim>                                       fe_q;
  DoFHandler<dim>                                 dof_handler;
  SimulationParameters<dim>                       simulation_parameters;
  unsigned int                                    number_quadrature_points;
  std::shared_ptr<TrilinosWrappers::SparseMatrix> system_matrix;
  TrilinosWrappers::MPI::Vector                   system_rhs;
};

/**
 * A class that assembles the rhs and solves the L2projection of the Vorticity
 * on nodes
 */
template <int dim>
class VorticitySmoothing : public PostProcessorSmoothing<dim>
{
public:
  // Member functions
  VorticitySmoothing(
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
    SimulationParameters<dim> simulation_parameters,
    unsigned int              number_quadrature_points);

  /**
   * @brief Outputs a solution for the field on the nodes.
   */
  void
  evaluate_smoothed_field();
};

/**
 * A class that assembles the rhs and solves the L2projection of the Qcriterion
 * on nodes
 */
template <int dim>
class QcriterionSmoothing : public PostProcessorSmoothing<dim>
{
public:
  // Member functions
  QcriterionSmoothing(
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
    SimulationParameters<dim> simulation_parameters,
    unsigned int              number_quadrature_points);

  /**
   * @brief Outputs a solution for the field on the nodes.
   */
  void
  evaluate_smoothed_field();
};

#endif
