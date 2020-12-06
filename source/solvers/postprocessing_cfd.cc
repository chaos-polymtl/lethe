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
#include <solvers/postprocessing_cfd.h>


using namespace dealii;

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



template <int dim, typename VectorType>
double
calculate_enstrophy(const DoFHandler<dim> &dof_handler,
                    const VectorType &     evaluation_point,
                    const Parameters::FEM &fem_parameters,
                    const MPI_Comm &       mpi_communicator)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  QGauss<dim>               quadrature_formula(fe.degree + 1);
  const MappingQ<dim>       mapping(fe.degree, fem_parameters.qmapping_all);
  FEValues<dim>             fe_values(mapping,
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
  const DoFHandler<2> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Parameters::FEM &              fem_parameters,
  const MPI_Comm &                     mpi_communicator);

template double
calculate_enstrophy<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Parameters::FEM &              fem_parameters,
  const MPI_Comm &                     mpi_communicator);

template double
calculate_enstrophy<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Parameters::FEM &                   fem_parameters,
  const MPI_Comm &                          mpi_communicator);

template double
calculate_enstrophy<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Parameters::FEM &                   fem_parameters,
  const MPI_Comm &                          mpi_communicator);


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
std::vector<Tensor<1, 3>>
calculate_torques(
  const DoFHandler<dim> &                              dof_handler,
  const VectorType &                                   evaluation_point,
  const Parameters::PhysicalProperties &               physical_properties,
  const Parameters::FEM &                              fem_parameters,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const MPI_Comm &                                     mpi_communicator)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();

  double viscosity = physical_properties.viscosity;

  QGauss<dim - 1>     face_quadrature_formula(fe.degree + 1);
  const MappingQ<dim> mapping(fe.degree, fem_parameters.qmapping_all);
  const unsigned int  n_q_points = face_quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  std::vector<double>              pressure_values(n_q_points);
  std::vector<Tensor<2, dim>>      velocity_gradients(n_q_points);
  Tensor<1, dim>                   normal_vector;
  Tensor<2, dim>                   fluid_stress;
  Tensor<2, dim>                   fluid_pressure;
  // torque tensor had to be considered in 3D at all time...
  Tensor<1, 3> torque;

  std::vector<Tensor<1, 3>> torque_vector(boundary_conditions.size);

  FEFaceValues<dim> fe_face_values(mapping,
                                   fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_gradients | update_JxW_values |
                                     update_normal_vectors);

  for (unsigned int i_bc = 0; i_bc < boundary_conditions.size; ++i_bc)
    {
      unsigned int boundary_id = boundary_conditions.id[i_bc];
      torque                   = 0;
      Point<dim> center_of_rotation =
        boundary_conditions.bcFunctions[boundary_id].cor;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
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
                              normal_vector = -fe_face_values.normal_vector(q);
                              for (int d = 0; d < dim; ++d)
                                {
                                  fluid_pressure[d][d] = pressure_values[q];
                                }
                              fluid_stress =
                                viscosity * (velocity_gradients[q] +
                                             transpose(velocity_gradients[q])) -
                                fluid_pressure;
                              auto force = fluid_stress * normal_vector *
                                           fe_face_values.JxW(q);

                              auto distance = q_points[q] - center_of_rotation;
                              if (dim == 2)
                                {
                                  torque[0] = 0.;
                                  torque[1] = 0.;
                                  torque[2] += distance[0] * force[1] -
                                               distance[1] * force[0];
                                }
                              else if (dim == 3)
                                {
                                  torque[0] += distance[1] * force[2] -
                                               distance[2] * force[1];
                                  torque[1] += distance[2] * force[0] -
                                               distance[0] * force[2];
                                  torque[2] += distance[0] * force[1] -
                                               distance[1] * force[0];
                                }
                            }
                        }
                    }
                }
            }
        }
      torque_vector[i_bc] = Utilities::MPI::sum(torque, mpi_communicator);
    }
  return torque_vector;
}

template std::vector<Tensor<1, 3>>
calculate_torques<2, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<2> &                              dof_handler,
  const TrilinosWrappers::MPI::Vector &              evaluation_point,
  const Parameters::PhysicalProperties &             physical_properties,
  const Parameters::FEM &                            fem_parameters,
  const BoundaryConditions::NSBoundaryConditions<2> &boundary_conditions,
  const MPI_Comm &                                   mpi_communicator);
template std::vector<Tensor<1, 3>>
calculate_torques<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                              dof_handler,
  const TrilinosWrappers::MPI::Vector &              evaluation_point,
  const Parameters::PhysicalProperties &             physical_properties,
  const Parameters::FEM &                            fem_parameters,
  const BoundaryConditions::NSBoundaryConditions<3> &boundary_conditions,
  const MPI_Comm &                                   mpi_communicator);

template std::vector<Tensor<1, 3>>
calculate_torques<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                              dof_handler,
  const TrilinosWrappers::MPI::BlockVector &         evaluation_point,
  const Parameters::PhysicalProperties &             physical_properties,
  const Parameters::FEM &                            fem_parameters,
  const BoundaryConditions::NSBoundaryConditions<2> &boundary_conditions,
  const MPI_Comm &                                   mpi_communicator);

template std::vector<Tensor<1, 3>>
calculate_torques<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                              dof_handler,
  const TrilinosWrappers::MPI::BlockVector &         evaluation_point,
  const Parameters::PhysicalProperties &             physical_properties,
  const Parameters::FEM &                            fem_parameters,
  const BoundaryConditions::NSBoundaryConditions<3> &boundary_conditions,
  const MPI_Comm &                                   mpi_communicator);


// Find the l2 norm of the error between the finite element sol'n and the exact
// sol'n for both the velocity and the pressure
// Mean pressure is removed from both the analytical and the simulation solution
template <int dim, typename VectorType>
std::pair<double, double>
calculate_L2_error(const DoFHandler<dim> &dof_handler,
                   const VectorType &     evaluation_point,
                   const Function<dim> *  exact_solution,
                   const Parameters::FEM &fem_parameters,
                   const MPI_Comm &       mpi_communicator)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();


  QGauss<dim>         quadrature_formula(fe.degree + 2);
  const MappingQ<dim> mapping(fe.degree, fem_parameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  const unsigned int dofs_per_cell =
    fe.dofs_per_cell; // This gives you dofs per cell
  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points = quadrature_formula.size();

  std::vector<Vector<double>> q_exactSol(n_q_points, Vector<double>(dim + 1));

  std::vector<Tensor<1, dim>> local_velocity_values(n_q_points);
  std::vector<double>         local_pressure_values(n_q_points);

  double pressure_integral       = 0;
  double exact_pressure_integral = 0;

  // loop over elements to calculate average pressure
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          fe_values[pressure].get_function_values(evaluation_point,
                                                  local_pressure_values);
          // Get the exact solution at all gauss points
          exact_solution->vector_value_list(fe_values.get_quadrature_points(),
                                            q_exactSol);


          // Retrieve the effective "connectivity matrix" for this element
          cell->get_dof_indices(local_dof_indices);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              pressure_integral += local_pressure_values[q] * fe_values.JxW(q);
              exact_pressure_integral += q_exactSol[q][dim] * fe_values.JxW(q);
            }
        }
    }

  pressure_integral = Utilities::MPI::sum(pressure_integral, mpi_communicator);
  exact_pressure_integral =
    Utilities::MPI::sum(exact_pressure_integral, mpi_communicator);

  double global_volume    = GridTools::volume(dof_handler.get_triangulation());
  double average_pressure = pressure_integral / global_volume;
  double average_exact_pressure = exact_pressure_integral / global_volume;


  double l2errorU = 0.;
  double l2errorP = 0.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[velocities].get_function_values(evaluation_point,
                                                    local_velocity_values);
          fe_values[pressure].get_function_values(evaluation_point,
                                                  local_pressure_values);

          // Retrieve the effective "connectivity matrix" for this element
          cell->get_dof_indices(local_dof_indices);

          // Get the exact solution at all gauss points
          exact_solution->vector_value_list(fe_values.get_quadrature_points(),
                                            q_exactSol);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              // Find the values of x and u_h (the finite element solution) at
              // the quadrature points
              double ux_sim   = local_velocity_values[q][0];
              double ux_exact = q_exactSol[q][0];

              double uy_sim   = local_velocity_values[q][1];
              double uy_exact = q_exactSol[q][1];

              l2errorU +=
                (ux_sim - ux_exact) * (ux_sim - ux_exact) * fe_values.JxW(q);
              l2errorU +=
                (uy_sim - uy_exact) * (uy_sim - uy_exact) * fe_values.JxW(q);

              if (dim == 3)
                {
                  double uz_sim   = local_velocity_values[q][2];
                  double uz_exact = q_exactSol[q][2];
                  l2errorU += (uz_sim - uz_exact) * (uz_sim - uz_exact) *
                              fe_values.JxW(q);
                }

              double p_sim   = local_pressure_values[q] - average_pressure;
              double p_exact = q_exactSol[q][dim] - average_exact_pressure;
              l2errorP +=
                (p_sim - p_exact) * (p_sim - p_exact) * fe_values.JxW(q);
            }
        }
    }
  l2errorU = Utilities::MPI::sum(l2errorU, mpi_communicator);
  l2errorP = Utilities::MPI::sum(l2errorP, mpi_communicator);

  return std::make_pair(std::sqrt(l2errorU), std::sqrt(l2errorP));
}

template std::pair<double, double>
calculate_L2_error(const DoFHandler<2> &                dof_handler,
                   const TrilinosWrappers::MPI::Vector &present_solution,
                   const Function<2> *                  l_exact_solution,
                   const Parameters::FEM &              fem_parameters,
                   const MPI_Comm &                     mpi_communicator);

template std::pair<double, double>
calculate_L2_error(const DoFHandler<3> &                dof_handler,
                   const TrilinosWrappers::MPI::Vector &present_solution,
                   const Function<3> *                  l_exact_solution,
                   const Parameters::FEM &              fem_parameters,
                   const MPI_Comm &                     mpi_communicator);

template std::pair<double, double>
calculate_L2_error(const DoFHandler<2> &                     dof_handler,
                   const TrilinosWrappers::MPI::BlockVector &present_solution,
                   const Function<2> *                       l_exact_solution,
                   const Parameters::FEM &                   fem_parameters,
                   const MPI_Comm &                          mpi_communicator);

template std::pair<double, double>
calculate_L2_error(const DoFHandler<3> &                     dof_handler,
                   const TrilinosWrappers::MPI::BlockVector &present_solution,
                   const Function<3> *                       l_exact_solution,
                   const Parameters::FEM &                   fem_parameters,
                   const MPI_Comm &                          mpi_communicator);

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
