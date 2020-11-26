// Base
#include <deal.II/base/quadrature_lib.h>

// Lac
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/vector.h>

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
#include <solvers/postprocessing_cfd.h>


using namespace dealii;

// This is a primitive first implementation that could be greatly improved by
// doing a single pass instead of N boundary passes
template <int dim, typename VectorType>
double
calculate_CFL(const DoFHandler<dim> &dof_handler,
              const VectorType &     evaluation_point,
              const double           time_step,
              const MPI_Comm &       mpi_communicator)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  QGauss<dim>               quadrature_formula(1);
  const MappingQ<dim>       mapping(fe.degree, false);
  FEValues<dim>             fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const unsigned int               n_q_points = quadrature_formula.size();


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
            h = pow(6 * cell->measure() / M_PI, 1. / 3.) / degree;
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
  const double                         time_step,
  const MPI_Comm &                     mpi_communicator);

template double
calculate_CFL<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const double                         time_step,
  const MPI_Comm &                     mpi_communicator);

template double
calculate_CFL<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const double                              time_step,
  const MPI_Comm &                          mpi_communicator);

template double
calculate_CFL<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const double                              time_step,
  const MPI_Comm &                          mpi_communicator);


// This is a primitive first implementation that could be greatly improved by
// doing a single pass instead of N boundary passes
template <int dim, typename VectorType>
std::vector<Tensor<1, dim>>
calculate_forces(
  const DoFHandler<dim> &                              dof_handler,
  const VectorType &                                   evaluation_point,
  const Parameters::PhysicalProperties &               physical_properties,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const MPI_Comm &                                     mpi_communicator)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();

  double viscosity = physical_properties.viscosity;

  QGauss<dim - 1>                  face_quadrature_formula(fe.degree + 1);
  const MappingQ<dim>              mapping(fe.degree, true);
  const unsigned int               n_q_points = face_quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  std::vector<double>              pressure_values(n_q_points);
  std::vector<Tensor<2, dim>>      velocity_gradients(n_q_points);
  Tensor<1, dim>                   normal_vector;
  Tensor<2, dim>                   fluid_stress;
  Tensor<2, dim>                   fluid_pressure;
  Tensor<1, dim>                   force;

  std::vector<Tensor<1, dim>> force_vector(boundary_conditions.size);

  FEFaceValues<dim> fe_face_values(mapping,
                                   fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_gradients | update_JxW_values |
                                     update_normal_vectors);

  for (unsigned int i_bc = 0; i_bc < boundary_conditions.size; ++i_bc)
    {
      unsigned int boundary_id = boundary_conditions.id[i_bc];
      force                    = 0;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              if (cell->at_boundary())
                {
                  for (unsigned int face = 0;
                       face < GeometryInfo<dim>::faces_per_cell;
                       face++)
                    {
                      if (cell->face(face)->at_boundary())
                        {
                          fe_face_values.reinit(cell, face);
                          if (cell->face(face)->boundary_id() == boundary_id)
                            {
                              std::vector<Point<dim>> q_points =
                                fe_face_values.get_quadrature_points();
                              fe_face_values[velocities].get_function_gradients(
                                evaluation_point, velocity_gradients);
                              fe_face_values[pressure].get_function_values(
                                evaluation_point, pressure_values);
                              for (unsigned int q = 0; q < n_q_points; q++)
                                {
                                  normal_vector =
                                    -fe_face_values.normal_vector(q);
                                  for (int d = 0; d < dim; ++d)
                                    {
                                      fluid_pressure[d][d] = pressure_values[q];
                                    }
                                  fluid_stress =
                                    viscosity *
                                      (velocity_gradients[q] +
                                       transpose(velocity_gradients[q])) -
                                    fluid_pressure;
                                  force += fluid_stress * normal_vector *
                                           fe_face_values.JxW(q);
                                }
                            }
                        }
                    }
                }
            }
        }
      force_vector[i_bc] = Utilities::MPI::sum(force, mpi_communicator);
    }
  return force_vector;
}

template std::vector<Tensor<1, 2>>
calculate_forces<2, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<2> &                              dof_handler,
  const TrilinosWrappers::MPI::Vector &              evaluation_point,
  const Parameters::PhysicalProperties &             physical_properties,
  const BoundaryConditions::NSBoundaryConditions<2> &boundary_conditions,
  const MPI_Comm &                                   mpi_communicator);
template std::vector<Tensor<1, 3>>
calculate_forces<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                              dof_handler,
  const TrilinosWrappers::MPI::Vector &              evaluation_point,
  const Parameters::PhysicalProperties &             physical_properties,
  const BoundaryConditions::NSBoundaryConditions<3> &boundary_conditions,
  const MPI_Comm &                                   mpi_communicator);

template std::vector<Tensor<1, 2>>
calculate_forces<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                              dof_handler,
  const TrilinosWrappers::MPI::BlockVector &         evaluation_point,
  const Parameters::PhysicalProperties &             physical_properties,
  const BoundaryConditions::NSBoundaryConditions<2> &boundary_conditions,
  const MPI_Comm &                                   mpi_communicator);

template std::vector<Tensor<1, 3>>
calculate_forces<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                              dof_handler,
  const TrilinosWrappers::MPI::BlockVector &         evaluation_point,
  const Parameters::PhysicalProperties &             physical_properties,
  const BoundaryConditions::NSBoundaryConditions<3> &boundary_conditions,
  const MPI_Comm &                                   mpi_communicator);



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
