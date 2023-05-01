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

#include <solvers/postprocessing_cfd.h>


using namespace dealii;

template <int dim, typename VectorType>
std::pair<double, double>
calculate_pressure_drop(const DoFHandler<dim> &       dof_handler,
                        std::shared_ptr<Mapping<dim>> mapping,
                        const VectorType &            evaluation_point,
                        const Quadrature<dim> &       cell_quadrature_formula,
                        const Quadrature<dim - 1> &   face_quadrature_formula,
                        const unsigned int            inlet_boundary_id,
                        const unsigned int            outlet_boundary_id)
{
  FESystem<dim, dim> fe = dof_handler.get_fe();
  FEValues<dim>      fe_values(*mapping,
                          fe,
                          cell_quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients |
                            update_hessians);
  FEFaceValues<dim>  fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  const unsigned int face_n_q_points = face_quadrature_formula.size();

  std::vector<Tensor<1, dim>> velocity_values(face_n_q_points);
  std::vector<double>         pressure_values(face_n_q_points);

  double static_pressure_outlet_boundary = 0.;
  double static_pressure_inlet_boundary  = 0.;
  double static_pressure_drop            = 0.;

  double dynamic_pressure_outlet_boundary = 0.;
  double dynamic_pressure_inlet_boundary  = 0.;
  double dynamic_pressure_drop            = 0.;

  double outlet_surface = 0.;
  double inlet_surface  = 0.;

  double temp_surface                      = 0.;
  double temp_static_pressure_at_boundary  = 0.;
  double temp_dynamic_pressure_at_boundary = 0.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int face_id = 0;
               face_id < GeometryInfo<dim>::faces_per_cell;
               face_id++)
            {
              if (cell->face(face_id)->at_boundary())
                {
                  unsigned int boundary_id = cell->face(face_id)->boundary_id();
                  if (boundary_id != inlet_boundary_id &&
                      boundary_id != outlet_boundary_id)
                    continue;
                  fe_face_values.reinit(cell, face_id);
                  fe_face_values[velocities].get_function_values(
                    evaluation_point, velocity_values);
                  fe_face_values[pressure].get_function_values(evaluation_point,
                                                               pressure_values);
                  temp_surface                      = 0.;
                  temp_static_pressure_at_boundary  = 0.;
                  temp_dynamic_pressure_at_boundary = 0.;
                  for (unsigned int q = 0; q < face_n_q_points; q++)
                    {
                      temp_surface += fe_face_values.JxW(q);
                      // Integration of velocity^2 and
                      // pressure, since the pressure in Lethe has units of
                      // Length^2/Time^2
                      temp_static_pressure_at_boundary +=
                        fe_face_values.JxW(q) * pressure_values[q];
                      temp_dynamic_pressure_at_boundary +=
                        0.5 * fe_face_values.JxW(q) * velocity_values[q] *
                        velocity_values[q];
                    }
                  if (boundary_id == inlet_boundary_id)
                    {
                      inlet_surface += temp_surface;
                      dynamic_pressure_inlet_boundary +=
                        temp_dynamic_pressure_at_boundary;
                      static_pressure_inlet_boundary +=
                        temp_static_pressure_at_boundary;
                    }
                  if (boundary_id == outlet_boundary_id)
                    {
                      outlet_surface += temp_surface;
                      dynamic_pressure_outlet_boundary +=
                        temp_dynamic_pressure_at_boundary;
                      static_pressure_outlet_boundary +=
                        temp_static_pressure_at_boundary;
                    }
                }
            }
        }
    }

  const MPI_Comm mpi_communicator = dof_handler.get_communicator();
  static_pressure_inlet_boundary =
    Utilities::MPI::sum(static_pressure_inlet_boundary, mpi_communicator);
  static_pressure_outlet_boundary =
    Utilities::MPI::sum(static_pressure_outlet_boundary, mpi_communicator);

  dynamic_pressure_inlet_boundary =
    Utilities::MPI::sum(dynamic_pressure_inlet_boundary, mpi_communicator);
  dynamic_pressure_outlet_boundary =
    Utilities::MPI::sum(dynamic_pressure_outlet_boundary, mpi_communicator);

  inlet_surface  = Utilities::MPI::sum(inlet_surface, mpi_communicator);
  outlet_surface = Utilities::MPI::sum(outlet_surface, mpi_communicator);

  static_pressure_outlet_boundary =
    static_pressure_outlet_boundary / outlet_surface;
  static_pressure_inlet_boundary =
    static_pressure_inlet_boundary / inlet_surface;
  static_pressure_drop =
    static_pressure_inlet_boundary - static_pressure_outlet_boundary;

  dynamic_pressure_outlet_boundary =
    dynamic_pressure_outlet_boundary / outlet_surface;
  dynamic_pressure_inlet_boundary =
    dynamic_pressure_inlet_boundary / inlet_surface;
  dynamic_pressure_drop =
    dynamic_pressure_inlet_boundary - dynamic_pressure_outlet_boundary;

  return {static_pressure_drop, static_pressure_drop + dynamic_pressure_drop};
}

template std::pair<double, double>
calculate_pressure_drop<2, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<2> &                dof_handler,
  std::shared_ptr<Mapping<2>>          mapping,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Quadrature<2> &                cell_quadrature_formula,
  const Quadrature<1> &                face_quadrature_formula,
  const unsigned int                   inlet_boundary_id,
  const unsigned int                   outlet_boundary_id);

template std::pair<double, double>
calculate_pressure_drop<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                dof_handler,
  std::shared_ptr<Mapping<3>>          mapping,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Quadrature<3> &                cell_quadrature_formula,
  const Quadrature<2> &                face_quadrature_formula,
  const unsigned int                   inlet_boundary_id,
  const unsigned int                   outlet_boundary_id);

template std::pair<double, double>
calculate_pressure_drop<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                     dof_handler,
  std::shared_ptr<Mapping<2>>               mapping,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Quadrature<2> &                     cell_quadrature_formula,
  const Quadrature<1> &                     face_quadrature_formula,
  const unsigned int                        inlet_boundary_id,
  const unsigned int                        outlet_boundary_id);

template std::pair<double, double>
calculate_pressure_drop<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                     dof_handler,
  std::shared_ptr<Mapping<3>>               mapping,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Quadrature<3> &                     cell_quadrature_formula,
  const Quadrature<2> &                     face_quadrature_formula,
  const unsigned int                        inlet_boundary_id,
  const unsigned int                        outlet_boundary_id);

template <int dim, typename VectorType>
double
calculate_CFL(const DoFHandler<dim> &dof_handler,
              const VectorType &     evaluation_point,
              const double           time_step,
              const Quadrature<dim> &quadrature_formula,
              const Mapping<dim> &   mapping)
{
  FESystem<dim, dim> fe = dof_handler.get_fe();
  FEValues<dim>      fe_values(mapping,
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
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();
  CFL                             = Utilities::MPI::max(CFL, mpi_communicator);
  return (CFL);
}

template double
calculate_CFL<2, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<2> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const double                         time_step,
  const Quadrature<2> &                quadrature_formula,
  const Mapping<2> &                   mapping);

template double
calculate_CFL<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const double                         time_step,
  const Quadrature<3> &                quadrature_formula,
  const Mapping<3> &                   mapping);

template double
calculate_CFL<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const double                              time_step,
  const Quadrature<2> &                     quadrature_formula,
  const Mapping<2> &                        mapping);

template double
calculate_CFL<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const double                              time_step,
  const Quadrature<3> &                     quadrature_formula,
  const Mapping<3> &                        mapping);



template <int dim, typename VectorType>
double
calculate_enstrophy(const DoFHandler<dim> &dof_handler,
                    const VectorType &     evaluation_point,
                    const Quadrature<dim> &quadrature_formula,
                    const Mapping<dim> &   mapping)
{
  const FESystem<dim, dim> fe = dof_handler.get_fe();
  FEValues<dim>            fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);

  const unsigned int n_q_points = quadrature_formula.size();

  std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);
  double                      en = 0.0;
  double                      domain_volume =
    GridTools::volume(dof_handler.get_triangulation(), mapping);

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
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();
  en                              = Utilities::MPI::sum(en, mpi_communicator);
  return (en);
}

template double
calculate_enstrophy<2, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<2> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Quadrature<2> &                quadrature_formula,
  const Mapping<2> &                   mapping);

template double
calculate_enstrophy<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Quadrature<3> &                quadrature_formula,
  const Mapping<3> &                   mapping);

template double
calculate_enstrophy<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Quadrature<2> &                     quadrature_formula,
  const Mapping<2> &                        mapping);

template double
calculate_enstrophy<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Quadrature<3> &                     quadrature_formula,
  const Mapping<3> &                        mapping);


template <int dim, typename VectorType>
double
calculate_kinetic_energy(const DoFHandler<dim> &dof_handler,
                         const VectorType &     evaluation_point,
                         const Quadrature<dim> &quadrature_formula,
                         const Mapping<dim> &   mapping)
{
  const FESystem<dim, dim> fe = dof_handler.get_fe();
  FEValues<dim>            fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const unsigned int               n_q_points = quadrature_formula.size();

  std::vector<Tensor<1, dim>> local_velocity_values(n_q_points);
  double                      domain_volume =
    GridTools::volume(dof_handler.get_triangulation(), mapping);
  // double domain_volume =
  // GridTools::volume(dof_handler.get_triangulation(),*mapping);

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
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();
  KEU = Utilities::MPI::sum(KEU / domain_volume, mpi_communicator);
  return (KEU);
}

template double
calculate_kinetic_energy<2, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<2> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Quadrature<2> &                quadrature_formula,
  const Mapping<2> &                   mapping);

template double
calculate_kinetic_energy<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Quadrature<3> &                quadrature_formula,
  const Mapping<3> &                   mapping);

template double
calculate_kinetic_energy<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Quadrature<2> &                     quadrature_formula,
  const Mapping<2> &                        mapping);

template double
calculate_kinetic_energy<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Quadrature<3> &                     quadrature_formula,
  const Mapping<3> &                        mapping);


template <int dim, typename VectorType>
double
calculate_apparent_viscosity(const DoFHandler<dim> &    dof_handler,
                             const VectorType &         evaluation_point,
                             const Quadrature<dim> &    quadrature_formula,
                             const Mapping<dim> &       mapping,
                             PhysicalPropertiesManager &properties_manager)
{
  double         integral_viscosity_x_shear_rate = 0;
  double         integral_shear_rate             = 0;
  double         shear_rate_magnitude;
  Tensor<2, dim> shear_rate;
  double         viscosity;
  const auto     rheological_model = properties_manager.get_rheology();

  const FESystem<dim, dim> fe = dof_handler.get_fe();
  FEValues<dim>            fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);

  const unsigned int n_q_points = quadrature_formula.size();

  std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          fe_values[velocities].get_function_gradients(
            evaluation_point, present_velocity_gradients);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              shear_rate = present_velocity_gradients[q] +
                           transpose(present_velocity_gradients[q]);

              double shear_rate_x_velocity_gradient = 0;
              for (int i = 0; i < dim; ++i)
                {
                  for (int j = 0; j < dim; ++j)
                    {
                      shear_rate_x_velocity_gradient +=
                        shear_rate[i][j] * present_velocity_gradients[q][j][i];
                    }
                }

              shear_rate_magnitude = calculate_shear_rate_magnitude(shear_rate);

              std::map<field, double> field_values;
              field_values[field::shear_rate] = shear_rate_magnitude;

              viscosity = rheological_model->value(field_values);

              integral_viscosity_x_shear_rate +=
                viscosity * shear_rate_x_velocity_gradient * fe_values.JxW(q);
              integral_shear_rate +=
                shear_rate_x_velocity_gradient * fe_values.JxW(q);
            }
        }
    }
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();
  integral_viscosity_x_shear_rate =
    Utilities::MPI::sum(integral_viscosity_x_shear_rate, mpi_communicator);
  integral_shear_rate =
    Utilities::MPI::sum(integral_shear_rate, mpi_communicator);
  return integral_viscosity_x_shear_rate / integral_shear_rate;
}

template double
calculate_apparent_viscosity<2, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<2> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Quadrature<2> &                quadrature_formula,
  const Mapping<2> &                   mapping,
  PhysicalPropertiesManager &          properties_manager);

template double
calculate_apparent_viscosity<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Quadrature<3> &                quadrature_formula,
  const Mapping<3> &                   mapping,
  PhysicalPropertiesManager &          properties_manager);

template double
calculate_apparent_viscosity<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Quadrature<2> &                     quadrature_formula,
  const Mapping<2> &                        mapping,
  PhysicalPropertiesManager &               properties_manager);

template double
calculate_apparent_viscosity<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &evaluation_point,
  const Quadrature<3> &                     quadrature_formula,
  const Mapping<3> &                        mapping,
  PhysicalPropertiesManager &               properties_manager);

template <int dim, typename VectorType>
std::vector<std::vector<Tensor<1, dim>>>
calculate_forces(
  const DoFHandler<dim> &                              dof_handler,
  const VectorType &                                   evaluation_point,
  PhysicalPropertiesManager &                          properties_manager,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const Quadrature<dim - 1> &                          face_quadrature_formula,
  const Mapping<dim> &                                 mapping)
{
  const FESystem<dim, dim> fe = dof_handler.get_fe();

  // Rheological model for viscosity properties
  double     viscosity;
  const auto rheological_model = properties_manager.get_rheology();


  const unsigned int               n_q_points = face_quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  std::vector<double>              pressure_values(n_q_points);
  std::vector<Tensor<2, dim>>      velocity_gradients(n_q_points);
  Tensor<1, dim>                   normal_vector;
  Tensor<2, dim>                   shear_rate;
  Tensor<2, dim>                   fluid_stress;
  Tensor<2, dim>                   fluid_viscous_stress;
  Tensor<2, dim>                   fluid_pressure;

  std::vector<Tensor<1, dim>> viscous_force_vector(boundary_conditions.size);
  std::vector<Tensor<1, dim>> pressure_force_vector(boundary_conditions.size);
  std::vector<Tensor<1, dim>> force_vector(boundary_conditions.size);


  FEFaceValues<dim> fe_face_values(mapping,
                                   fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_gradients | update_JxW_values |
                                     update_normal_vectors);
  const MPI_Comm    mpi_communicator = dof_handler.get_communicator();

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          if (cell->at_boundary())
            {
              for (const auto face : cell->face_indices())
                {
                  if (cell->face(face)->at_boundary())
                    {
                      const auto boundary_id =
                        std::find(begin(boundary_conditions.id),
                                  end(boundary_conditions.id),
                                  cell->face(face)->boundary_id());


                      if (boundary_id != end(boundary_conditions.id))
                        {
                          unsigned int vector_index =
                            boundary_id - begin(boundary_conditions.id);
                          fe_face_values.reinit(cell, face);

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
                              shear_rate = velocity_gradients[q] +
                                           transpose(velocity_gradients[q]);

                              const double shear_rate_magnitude =
                                calculate_shear_rate_magnitude(shear_rate);

                              std::map<field, double> field_values;
                              field_values[field::shear_rate] =
                                shear_rate_magnitude;

                              viscosity =
                                rheological_model->value(field_values);
                              fluid_viscous_stress = -viscosity * shear_rate;
                              fluid_stress =
                                -fluid_viscous_stress - fluid_pressure;

                              viscous_force_vector[vector_index] -=
                                fluid_viscous_stress * normal_vector *
                                fe_face_values.JxW(q);
                              pressure_force_vector[vector_index] -=
                                fluid_pressure * normal_vector *
                                fe_face_values.JxW(q);
                              force_vector[vector_index] +=
                                fluid_stress * normal_vector *
                                fe_face_values.JxW(q);
                            }
                        }
                    }
                }
            }
        }
    }

  for (unsigned int i_bc = 0; i_bc < boundary_conditions.size; ++i_bc)
    {
      viscous_force_vector[i_bc] =
        Utilities::MPI::sum(viscous_force_vector[i_bc], mpi_communicator);
      pressure_force_vector[i_bc] =
        Utilities::MPI::sum(pressure_force_vector[i_bc], mpi_communicator);
      force_vector[i_bc] =
        Utilities::MPI::sum(force_vector[i_bc], mpi_communicator);
    }
  std::vector<std::vector<Tensor<1, dim>>> forces{force_vector,
                                                  viscous_force_vector,
                                                  pressure_force_vector};
  return forces;
}

template std::vector<std::vector<Tensor<1, 2>>>
calculate_forces<2, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<2> &                              dof_handler,
  const TrilinosWrappers::MPI::Vector &              evaluation_point,
  PhysicalPropertiesManager &                        properties_manager,
  const BoundaryConditions::NSBoundaryConditions<2> &boundary_conditions,
  const Quadrature<1> &                              face_quadrature_formula,
  const Mapping<2> &                                 mapping);
template std::vector<std::vector<Tensor<1, 3>>>
calculate_forces<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                              dof_handler,
  const TrilinosWrappers::MPI::Vector &              evaluation_point,
  PhysicalPropertiesManager &                        properties_manager,
  const BoundaryConditions::NSBoundaryConditions<3> &boundary_conditions,
  const Quadrature<2> &                              face_quadrature_formula,
  const Mapping<3> &                                 mapping);

template std::vector<std::vector<Tensor<1, 2>>>
calculate_forces<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                              dof_handler,
  const TrilinosWrappers::MPI::BlockVector &         evaluation_point,
  PhysicalPropertiesManager &                        properties_manager,
  const BoundaryConditions::NSBoundaryConditions<2> &boundary_conditions,
  const Quadrature<1> &                              face_quadrature_formula,
  const Mapping<2> &                                 mapping);

template std::vector<std::vector<Tensor<1, 3>>>
calculate_forces<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                              dof_handler,
  const TrilinosWrappers::MPI::BlockVector &         evaluation_point,
  PhysicalPropertiesManager &                        properties_manager,
  const BoundaryConditions::NSBoundaryConditions<3> &boundary_conditions,
  const Quadrature<2> &                              face_quadrature_formula,
  const Mapping<3> &                                 mapping);


template <int dim, typename VectorType>
std::vector<Tensor<1, 3>>
calculate_torques(
  const DoFHandler<dim> &                              dof_handler,
  const VectorType &                                   evaluation_point,
  PhysicalPropertiesManager &                          properties_manager,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const Quadrature<dim - 1> &                          face_quadrature_formula,
  const Mapping<dim> &                                 mapping)
{
  const FESystem<dim, dim> fe = dof_handler.get_fe();

  // Rheological model for viscosity properties
  double     viscosity;
  const auto rheological_model = properties_manager.get_rheology();


  const unsigned int               n_q_points = face_quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  std::vector<double>              pressure_values(n_q_points);
  std::vector<Tensor<2, dim>>      velocity_gradients(n_q_points);
  Tensor<1, dim>                   normal_vector;
  Tensor<2, dim>                   shear_rate;
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
  const MPI_Comm    mpi_communicator = dof_handler.get_communicator();
  for (unsigned int i_bc = 0; i_bc < boundary_conditions.size; ++i_bc)
    {
      unsigned int boundary_id = boundary_conditions.id[i_bc];
      torque                   = 0;
      Point<dim> center_of_rotation =
        boundary_conditions.bcFunctions[boundary_id].center_of_rotation;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              for (const auto face : cell->face_indices())
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
                              shear_rate = velocity_gradients[q] +
                                           transpose(velocity_gradients[q]);
                              const double shear_rate_magnitude =
                                calculate_shear_rate_magnitude(shear_rate);

                              std::map<field, double> field_values;
                              field_values[field::shear_rate] =
                                shear_rate_magnitude;

                              viscosity =
                                rheological_model->value(field_values);

                              fluid_stress =
                                viscosity * shear_rate - fluid_pressure;
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
  PhysicalPropertiesManager &                        properties_manager,
  const BoundaryConditions::NSBoundaryConditions<2> &boundary_conditions,
  const Quadrature<1> &                              face_quadrature_formula,
  const Mapping<2> &                                 mapping);
template std::vector<Tensor<1, 3>>
calculate_torques<3, TrilinosWrappers::MPI::Vector>(
  const DoFHandler<3> &                              dof_handler,
  const TrilinosWrappers::MPI::Vector &              evaluation_point,
  PhysicalPropertiesManager &                        properties_manager,
  const BoundaryConditions::NSBoundaryConditions<3> &boundary_conditions,
  const Quadrature<2> &                              face_quadrature_formula,
  const Mapping<3> &                                 mapping);

template std::vector<Tensor<1, 3>>
calculate_torques<2, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<2> &                              dof_handler,
  const TrilinosWrappers::MPI::BlockVector &         evaluation_point,
  PhysicalPropertiesManager &                        properties_manager,
  const BoundaryConditions::NSBoundaryConditions<2> &boundary_conditions,
  const Quadrature<1> &                              face_quadrature_formula,
  const Mapping<2> &                                 mapping);

template std::vector<Tensor<1, 3>>
calculate_torques<3, TrilinosWrappers::MPI::BlockVector>(
  const DoFHandler<3> &                              dof_handler,
  const TrilinosWrappers::MPI::BlockVector &         evaluation_point,
  PhysicalPropertiesManager &                        properties_manager,
  const BoundaryConditions::NSBoundaryConditions<3> &boundary_conditions,
  const Quadrature<2> &                              face_quadrature_formula,
  const Mapping<3> &                                 mapping);


// Find the l2 norm of the error between the finite element sol'n and the exact
// sol'n for both the velocity and the pressure
// Mean pressure is removed from both the analytical and the simulation solution
template <int dim, typename VectorType>
std::pair<double, double>
calculate_L2_error(const DoFHandler<dim> &dof_handler,
                   const VectorType &     evaluation_point,
                   const Function<dim> *  exact_solution,
                   const Quadrature<dim> &quadrature_formula,
                   const Mapping<dim> &   mapping)
{
  const FESystem<dim, dim> fe = dof_handler.get_fe();
  FEValues<dim>            fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  const unsigned int               n_q_points = quadrature_formula.size();

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

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              pressure_integral += local_pressure_values[q] * fe_values.JxW(q);
              exact_pressure_integral += q_exactSol[q][dim] * fe_values.JxW(q);
            }
        }
    }
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();
  pressure_integral = Utilities::MPI::sum(pressure_integral, mpi_communicator);
  exact_pressure_integral =
    Utilities::MPI::sum(exact_pressure_integral, mpi_communicator);

  // double global_volume    =
  // GridTools::volume(dof_handler.get_triangulation(),*mapping);
  double global_volume =
    GridTools::volume(dof_handler.get_triangulation(), mapping);
  double average_pressure       = pressure_integral / global_volume;
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
                   const Quadrature<2> &                quadrature_formula,
                   const Mapping<2> &                   mapping);

template std::pair<double, double>
calculate_L2_error(const DoFHandler<3> &                dof_handler,
                   const TrilinosWrappers::MPI::Vector &present_solution,
                   const Function<3> *                  l_exact_solution,
                   const Quadrature<3> &                quadrature_formula,
                   const Mapping<3> &                   mapping);

template std::pair<double, double>
calculate_L2_error(const DoFHandler<2> &                     dof_handler,
                   const TrilinosWrappers::MPI::BlockVector &present_solution,
                   const Function<2> *                       l_exact_solution,
                   const Quadrature<2> &                     quadrature_formula,
                   const Mapping<2> &                        mapping);

template std::pair<double, double>
calculate_L2_error(const DoFHandler<3> &                     dof_handler,
                   const TrilinosWrappers::MPI::BlockVector &present_solution,
                   const Function<3> *                       l_exact_solution,
                   const Quadrature<3> &                     quadrature_formula,
                   const Mapping<3> &                        mapping);

template <int dim, typename VectorType>
std::pair<double, double>
calculate_flow_rate(const DoFHandler<dim> &    dof_handler,
                    const VectorType &         present_solution,
                    const unsigned int &       boundary_id,
                    const Quadrature<dim - 1> &face_quadrature_formula,
                    const Mapping<dim> &       mapping)
{
  const FESystem<dim, dim> fe = dof_handler.get_fe();

  const unsigned int               n_q_points = face_quadrature_formula.size();
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

  // Calculating area and volumetric flow rate at the boundary
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

  const MPI_Comm mpi_communicator = dof_handler.get_communicator();
  area                            = Utilities::MPI::sum(area, mpi_communicator);
  flow_rate = Utilities::MPI::sum(flow_rate, mpi_communicator);

  return std::make_pair(flow_rate, area);
}

template std::pair<double, double>
calculate_flow_rate(const DoFHandler<2> &                dof_handler,
                    const TrilinosWrappers::MPI::Vector &present_solution,
                    const unsigned int &                 boundary_id,
                    const Quadrature<1> &face_quadrature_formula,
                    const Mapping<2> &   mapping);

template std::pair<double, double>
calculate_flow_rate(const DoFHandler<3> &                dof_handler,
                    const TrilinosWrappers::MPI::Vector &present_solution,
                    const unsigned int &                 boundary_id,
                    const Quadrature<2> &face_quadrature_formula,
                    const Mapping<3> &   mapping);

template std::pair<double, double>
calculate_flow_rate(const DoFHandler<2> &                     dof_handler,
                    const TrilinosWrappers::MPI::BlockVector &present_solution,
                    const unsigned int &                      boundary_id,
                    const Quadrature<1> &face_quadrature_formula,
                    const Mapping<2> &   mapping);

template std::pair<double, double>
calculate_flow_rate(const DoFHandler<3> &                     dof_handler,
                    const TrilinosWrappers::MPI::BlockVector &present_solution,
                    const unsigned int &                      boundary_id,
                    const Quadrature<2> &face_quadrature_formula,
                    const Mapping<3> &   mapping);

template <int dim, typename VectorType>
double
calculate_average_velocity(const DoFHandler<dim> &    dof_handler,
                           const VectorType &         present_solution,
                           const unsigned int &       boundary_id,
                           const Quadrature<dim - 1> &face_quadrature_formula,
                           const Mapping<dim> &       mapping)
{
  std::pair<double, double> flow_rate_and_area =
    calculate_flow_rate(dof_handler,
                        present_solution,
                        boundary_id,
                        face_quadrature_formula,
                        mapping);

  // In dynamic flow control, the flow rate is negative if the velocity is
  // positive when calculated at inlet boundary because of the normal vector
  // direction. Hence, we need to invert the sign to get the average velocity.
  return -flow_rate_and_area.first / flow_rate_and_area.second;
}

template double
calculate_average_velocity(
  const DoFHandler<2> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &present_solution,
  const unsigned int &                 boundary_id,
  const Quadrature<1> &                face_quadrature_formula,
  const Mapping<2> &                   mapping);

template double
calculate_average_velocity(
  const DoFHandler<3> &                dof_handler,
  const TrilinosWrappers::MPI::Vector &present_solution,
  const unsigned int &                 boundary_id,
  const Quadrature<2> &                face_quadrature_formula,
  const Mapping<3> &                   mapping);

template double
calculate_average_velocity(
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &present_solution,
  const unsigned int &                      boundary_id,
  const Quadrature<1> &                     face_quadrature_formula,
  const Mapping<2> &                        mapping);

template double
calculate_average_velocity(
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &present_solution,
  const unsigned int &                      boundary_id,
  const Quadrature<2> &                     face_quadrature_formula,
  const Mapping<3> &                        mapping);

template <int dim, typename VectorType>
double
calculate_average_velocity(const DoFHandler<dim> &dof_handler,
                           const DoFHandler<dim> &void_fraction_dof_handler,
                           const VectorType &     present_solution,
                           const VectorType &  present_void_fraction_solution,
                           const unsigned int &flow_direction,
                           const Quadrature<dim> &quadrature_formula,
                           const Mapping<dim> &   mapping)
{
  // Set up for velocity fe values
  const auto &                     tria       = dof_handler.get_triangulation();
  const FESystem<dim, dim>         fe         = dof_handler.get_fe();
  const unsigned int               n_q_points = quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      velocity_values(n_q_points);
  Tensor<1, dim>                   normal_vector;

  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  // Set up for void fraction fe values
  const FESystem<dim, dim> fe_void_fraction =
    void_fraction_dof_handler.get_fe();
  const FEValuesExtractors::Scalar void_fraction(0);
  std::vector<double>              void_fraction_values(n_q_points);

  FEValues<dim> fe_vf_values(mapping,
                             fe_void_fraction,
                             quadrature_formula,
                             update_values | update_quadrature_points);

  // Initialize variables for summation
  double average_velocity = 0;
  double volume           = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          typename DoFHandler<dim>::active_cell_iterator vf_cell(
            &tria, cell->level(), cell->index(), &void_fraction_dof_handler);
          fe_vf_values.reinit(vf_cell);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              // Get the void fraction at the quadrature point
              fe_vf_values[void_fraction].get_function_values(
                present_void_fraction_solution, void_fraction_values);

              // Add the area of the face on boundary
              // weighted by void fraction
              volume += void_fraction_values[q] * fe_values.JxW(q);

              // Add the flow rate at the face on boundary weighted
              // by void fraction
              fe_values[velocities].get_function_values(present_solution,
                                                        velocity_values);
              average_velocity += void_fraction_values[q] *
                                  velocity_values[q][flow_direction] *
                                  fe_values.JxW(q);
            }
        }
    }

  const MPI_Comm mpi_communicator = dof_handler.get_communicator();
  average_velocity = Utilities::MPI::sum(average_velocity, mpi_communicator);
  volume           = Utilities::MPI::sum(volume, mpi_communicator);

  // Calculate the average velocity
  average_velocity /= volume;

  return average_velocity;
}

template double
calculate_average_velocity(
  const DoFHandler<2> &                dof_handler,
  const DoFHandler<2> &                void_fraction_dof_handler,
  const TrilinosWrappers::MPI::Vector &present_solution,
  const TrilinosWrappers::MPI::Vector &present_void_fraction_solution,
  const unsigned int &                 flow_direction,
  const Quadrature<2> &                quadrature_formula,
  const Mapping<2> &                   mapping);

template double
calculate_average_velocity(
  const DoFHandler<3> &                dof_handler,
  const DoFHandler<3> &                void_fraction_dof_handler,
  const TrilinosWrappers::MPI::Vector &present_solution,
  const TrilinosWrappers::MPI::Vector &present_void_fraction_solution,
  const unsigned int &                 flow_direction,
  const Quadrature<3> &                quadrature_formula,
  const Mapping<3> &                   mapping);
