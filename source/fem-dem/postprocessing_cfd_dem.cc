// Base
#include <deal.II/base/quadrature_lib.h>

// Lac
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/vector.h>

// grid
#include <deal.II/grid/grid_tools.h>

// Dofs
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

// Lac - Trilinos includes
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

// Fe
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

// Lethe includes
#include <core/boundary_conditions.h>
#include <core/parameters.h>
#include <core/rheological_model.h>
#include <core/vector.h>

#include <fem-dem/postprocessing_cfd_dem.h>


using namespace dealii;


template <int dim, typename VectorType>
std::pair<double, double>
calculate_total_volume(const DoFHandler<dim> &void_fraction_dof_handler,
                           const VectorType   &present_void_fraction_solution,
                           const Quadrature<dim> &quadrature_formula,
                           const Mapping<dim>    &mapping)
{
  // Set up for void fraction fe values
  const unsigned int               n_q_points = quadrature_formula.size();
  const FESystem<dim, dim> fe_void_fraction =
    void_fraction_dof_handler.get_fe();
  const FEValuesExtractors::Scalar void_fraction(0);
  std::vector<double>              void_fraction_values(n_q_points);

  FEValues<dim> fe_vf_values(mapping,
                             fe_void_fraction,
                             quadrature_formula,
                             update_values | update_quadrature_points);

  // Initialize variables for summation
  double total_volume_fluid = 0;
  double total_volume_solid = 0;

  for (const auto &cell : void_fraction_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_vf_values.reinit(cell);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              // Get the void fraction at the quadrature point
              fe_vf_values[void_fraction].get_function_values(
                present_void_fraction_solution, void_fraction_values);

              // Calculate total volume
              total_volume_fluid += void_fraction_values[q] * fe_vf_values.JxW(q);
              total_volume_solid += (1 - void_fraction_values[q]) * fe_vf_values.JxW(q); 
            }
        }
    }

  const MPI_Comm mpi_communicator = void_fraction_dof_handler.get_communicator();
  total_volume_fluid = Utilities::MPI::sum(total_volume_fluid, mpi_communicator);
  total_volume_solid = Utilities::MPI::sum(total_volume_solid, mpi_communicator);


  return {total_volume_fluid, total_volume_solid};
}

template std::pair<double, double>
calculate_total_volume<2, GlobalVectorType>(
  const DoFHandler<2>    &void_fraction_dof_handler,
  const GlobalVectorType &present_void_fraction_solution,
  const Quadrature<2>    &quadrature_formula,
  const Mapping<2>       &mapping);

template std::pair<double, double>
calculate_total_volume<3, GlobalVectorType>(
  const DoFHandler<3>    &void_fraction_dof_handler,
  const GlobalVectorType &present_void_fraction_solution,
  const Quadrature<3>    &quadrature_formula,
  const Mapping<3>       &mapping);

template std::pair<double, double>
calculate_total_volume<2, GlobalBlockVectorType>(
  const DoFHandler<2>    &void_fraction_dof_handler,
  const GlobalBlockVectorType &present_void_fraction_solution,
  const Quadrature<2>    &quadrature_formula,
  const Mapping<2>       &mapping);

template std::pair<double, double>
calculate_total_volume<3, GlobalBlockVectorType>(
  const DoFHandler<3>    &void_fraction_dof_handler,
  const GlobalBlockVectorType &present_void_fraction_solution,
  const Quadrature<3>    &quadrature_formula,
  const Mapping<3>       &mapping);

#ifndef LETHE_USE_LDV
template std::pair<double, double>
calculate_total_volume<2, LinearAlgebra::distributed::Vector<double>>(
  const DoFHandler<2>    &void_fraction_dof_handler,
  const LinearAlgebra::distributed::Vector<double> &present_void_fraction_solution,
  const Quadrature<2>    &quadrature_formula,
  const Mapping<2>       &mapping);

template std::pair<double, double>
calculate_total_volume<3, LinearAlgebra::distributed::Vector<double>>(
  const DoFHandler<3>    &void_fraction_dof_handler,
  const LinearAlgebra::distributed::Vector<double> &present_void_fraction_solution,
  const Quadrature<3>    &quadrature_formula,
  const Mapping<3>       &mapping);
#endif
