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
#include <solvers/postprocessing_kinetic_energy.h>


using namespace dealii;

// This is a primitive first implementation that could be greatly improved by
// doing a single pass instead of N boundary passes
template <int dim, typename VectorType>
double
calculate_kinetic_energy(const DoFHandler<dim> &dof_handler,
                         const VectorType &     evaluation_point,
                         const Parameters::FEM &fem_parameters,
                         const MPI_Comm &       mpi_communicator)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();

  QGauss<dim>         quadrature_formula(fe.degree + 1);
  const MappingQ<dim> mapping(fe.degree, fem_parameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const unsigned int               n_q_points = quadrature_formula.size();

  std::vector<Tensor<1, dim>> local_velocity_values(n_q_points);
  double domain_volume = GridTools::volume(dof_handler.get_triangulation());

  double KEU = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[velocities].get_function_values(evaluation_point,
                                                    local_velocity_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              double ux_sim = local_velocity_values[q][0];
              double uy_sim = local_velocity_values[q][1];

              KEU += 0.5 * ((ux_sim) * (ux_sim)*fe_values.JxW(q));
              KEU += 0.5 * ((uy_sim) * (uy_sim)*fe_values.JxW(q));
              if (dim == 3)
                {
                  double uz_sim = local_velocity_values[q][2];
                  KEU += 0.5 * ((uz_sim) * (uz_sim)*fe_values.JxW(q));
                }
            }
        }
    }
  KEU = Utilities::MPI::sum(KEU / domain_volume, mpi_communicator);
  return (KEU);
}

template double
calculate_kinetic_energy<2, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<2> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Parameters::FEM &              fem_parameters,
  const MPI_Comm &                     mpi_communicator);

template double
calculate_kinetic_energy<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Parameters::FEM &              fem_parameters,
  const MPI_Comm &                     mpi_communicator);

template double
calculate_kinetic_energy<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Parameters::FEM &                   fem_parameters,
  const MPI_Comm &                          mpi_communicator);

template double
calculate_kinetic_energy<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Parameters::FEM &                   fem_parameters,
  const MPI_Comm &                          mpi_communicator);
