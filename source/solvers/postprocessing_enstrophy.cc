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
#include <solvers/postprocessing_enstrophy.h>


using namespace dealii;

// This is a primitive first implementation that could be greatly improved by
// doing a single pass instead of N boundary passes
template <int dim, typename VectorType>
double
calculate_enstrophy(const FESystem<dim> &  fe,
                    const DoFHandler<dim> &dof_handler,
                    const VectorType &     evaluation_point,
                    const Parameters::FEM &fem_parameters,
                    const MPI_Comm &       mpi_communicator)
{
  QGauss<dim>         quadrature_formula(fe.degree + 1);
  const MappingQ<dim> mapping(fe.degree, fem_parameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);

  const unsigned int n_q_points = quadrature_formula.size();

  std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);
  double                      en = 0.0;
  double domain_volume = GridTools::volume(dof_handler.get_triangulation());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          fe_values[velocities].get_function_gradients(
            evaluation_point, present_velocity_gradients);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              // Find the values of gradient of ux and uy (the finite element
              // solution) at the quadrature points
              double ux_y = present_velocity_gradients[q][0][1];
              double uy_x = present_velocity_gradients[q][1][0];

              if (dim == 2)
                {
                  en += 0.5 * (uy_x - ux_y) * (uy_x - ux_y) * fe_values.JxW(q) /
                        domain_volume;
                }
              else
                {
                  double uz_y = present_velocity_gradients[q][2][1];
                  double uy_z = present_velocity_gradients[q][1][2];
                  double ux_z = present_velocity_gradients[q][0][2];
                  double uz_x = present_velocity_gradients[q][2][0];
                  en += 0.5 * (uz_y - uy_z) * (uz_y - uy_z) * fe_values.JxW(q) /
                        domain_volume;
                  en += 0.5 * (ux_z - uz_x) * (ux_z - uz_x) * fe_values.JxW(q) /
                        domain_volume;
                  en += 0.5 * (uy_x - ux_y) * (uy_x - ux_y) * fe_values.JxW(q) /
                        domain_volume;
                }
            }
        }
    }
  en = Utilities::MPI::sum(en, mpi_communicator);
  return (en);
}

template double
calculate_enstrophy<2, TrilinosWrappers::MPI::Vector>(
  const FESystem<2> &                  fe,
  const DoFHandler<2> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Parameters::FEM &              fem_parameters,
  const MPI_Comm &                     mpi_communicator);

template double
calculate_enstrophy<3, TrilinosWrappers::MPI::Vector>(
  const FESystem<3> &                  fe,
  const DoFHandler<3> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Parameters::FEM &              fem_parameters,
  const MPI_Comm &                     mpi_communicator);

template double
calculate_enstrophy<2, TrilinosWrappers::MPI::BlockVector>(
  const FESystem<2> &                       fe,
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Parameters::FEM &                   fem_parameters,
  const MPI_Comm &                          mpi_communicator);

template double
calculate_enstrophy<3, TrilinosWrappers::MPI::BlockVector>(
  const FESystem<3> &                       fe,
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Parameters::FEM &                   fem_parameters,
  const MPI_Comm &                          mpi_communicator);
