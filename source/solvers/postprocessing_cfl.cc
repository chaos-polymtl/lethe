// Base
#include <deal.II/base/quadrature_lib.h>

// Lac
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/vector.h>

// Dofs
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

// Grid
#include <deal.II/grid/grid_tools.h>

// Lac - Trilinos includes
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

// Fe
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

// Lethe includes
#include <core/parameters.h>
#include <solvers/postprocessing_cfl.h>


using namespace dealii;

// This is a primitive first implementation that could be greatly improved by
// doing a single pass instead of N boundary passes
template <int dim, typename VectorType>
double
calculate_CFL(const DoFHandler<dim> &dof_handler,
              const VectorType &     evaluation_point,
              const Parameters::FEM &fem_parameters,
              const double time_step,
              const MPI_Comm &       mpi_communicator)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  QGauss<dim>         quadrature_formula(1);
  const MappingQ<dim> mapping(fe.degree, fem_parameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                                fe,
                                quadrature_formula,
                                update_values | update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const unsigned int                   n_q_points = quadrature_formula.size();


  std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);

  // Element size
  double h;

  // Element degree
  double degree = double(fe.degree);

  // CFL
  double CFL = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          if (dim == 2)
            h = std::sqrt(4. * cell->measure() / M_PI) / degree;
          else if (dim == 3)
            h =
              pow(6 * cell->measure() / M_PI, 1. / 3.) / degree;
          fe_values.reinit(cell);
          fe_values[velocities].get_function_values(evaluation_point,
                                                    present_velocity_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double localCFL =
                present_velocity_values[q].norm() / h * time_step;
              CFL = std::max(CFL, localCFL);
            }
        }
    }
  CFL = Utilities::MPI::max(CFL, mpi_communicator);
  return (CFL);
}

template double
calculate_CFL<2, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<2> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Parameters::FEM &              fem_parameters,
  const double time_step,
  const MPI_Comm &                     mpi_communicator);

template double
calculate_CFL<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Parameters::FEM &              fem_parameters,
  const double time_step,
  const MPI_Comm &                     mpi_communicator);

template double
calculate_CFL<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Parameters::FEM &                   fem_parameters,
  const double time_step,
  const MPI_Comm &                          mpi_communicator);

template double
calculate_CFL<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Parameters::FEM &                   fem_parameters,
  const double time_step,
  const MPI_Comm &                          mpi_communicator);
