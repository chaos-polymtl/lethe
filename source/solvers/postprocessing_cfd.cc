// Lac - Trilinos includes
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <solvers/postprocessing_cfd.h>



template <int dim, typename VectorType>
std::pair<double, double>
calculate_flow_rate(const DoFHandler<dim> &dof_handler,
                    const VectorType &     present_solution,
                    const unsigned int &   boundary_id,
                    const Parameters::FEM &fem_parameters,
                    const MPI_Comm &       mpi_communicator)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  const MappingQ<dim>       mapping(fe.degree, fem_parameters.qmapping_all);
  QGauss<dim - 1>           face_quadrature_formula(fe.degree + 1);
  const unsigned int        n_q_points = face_quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      velocity_values(n_q_points);
  Tensor<1, dim>                   normal_vector;

  FEFaceValues<dim> fe_face_values(mapping,
                                   fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_JxW_values | update_normal_vectors);

  double flow_rate = 0;
  double area      = 0;

  // Calculating area and volumetric flow rate at the inlet flow
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() && cell->at_boundary())
        {
          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               face++)
            {
              if (cell->face(face)->at_boundary())
                {
                  fe_face_values.reinit(cell, face);
                  if (cell->face(face)->boundary_id() == boundary_id)
                    {
                      for (unsigned int q = 0; q < n_q_points; q++)
                        {
                          area += fe_face_values.JxW(q);
                          normal_vector = fe_face_values.normal_vector(q);
                          fe_face_values[velocities].get_function_values(
                            present_solution, velocity_values);
                          flow_rate += velocity_values[q] * normal_vector *
                                       fe_face_values.JxW(q);
                        }
                    }
                }
            }
        }
    }

  area      = Utilities::MPI::sum(area, mpi_communicator);
  flow_rate = Utilities::MPI::sum(flow_rate, mpi_communicator);

  return std::make_pair(flow_rate, area);
}

template std::pair<double, double>
calculate_flow_rate(const DoFHandler<2> &                dof_handler,
                    const TrilinosWrappers::MPI::Vector &present_solution,
                    const unsigned int &                 boundary_id,
                    const Parameters::FEM &              fem_parameters,
                    const MPI_Comm &                     mpi_communicator);

template std::pair<double, double>
calculate_flow_rate(const DoFHandler<3> &                dof_handler,
                    const TrilinosWrappers::MPI::Vector &present_solution,
                    const unsigned int &                 boundary_id,
                    const Parameters::FEM &              fem_parameters,
                    const MPI_Comm &                     mpi_communicator);

template std::pair<double, double>
calculate_flow_rate(const DoFHandler<2> &                     dof_handler,
                    const TrilinosWrappers::MPI::BlockVector &present_solution,
                    const unsigned int &                      boundary_id,
                    const Parameters::FEM &                   fem_parameters,
                    const MPI_Comm &                          mpi_communicator);

template std::pair<double, double>
calculate_flow_rate(const DoFHandler<3> &                     dof_handler,
                    const TrilinosWrappers::MPI::BlockVector &present_solution,
                    const unsigned int &                      boundary_id,
                    const Parameters::FEM &                   fem_parameters,
                    const MPI_Comm &                          mpi_communicator);
