/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#include <core/bdf.h>
#include <core/grids.h>
#include <core/lethegridtools.h>
#include <core/sdirk.h>
#include <core/solutions_output.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/navier_stokes_vof_assemblers.h>
#include <solvers/postprocessing_cfd.h>

#include <fem-dem/gls_sharp_navier_stokes.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/full_matrix.h>

// Constructor for class GLSNavierStokesSolver
template <int dim>
GLSSharpNavierStokesSolver<dim>::GLSSharpNavierStokesSolver(
  SimulationParameters<dim> &p_nsparam)
  : GLSNavierStokesSolver<dim>(p_nsparam)
{}

template <int dim>
GLSSharpNavierStokesSolver<dim>::~GLSSharpNavierStokesSolver()
{}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::vertices_cell_mapping()
{
  // Find all the cells around each vertices
  TimerOutput::Scope t(this->computing_timer, "vertices_to_cell_map");

  LetheGridTools::vertices_cell_mapping(this->dof_handler, vertices_to_cell);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::generate_cut_cells_map()
{
  // check all the cells if they are cut or not. Put the information in a map
  // with the key being the cell.
  TimerOutput::Scope t(this->computing_timer, "cut_cells_mapping");
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(*this->mapping,
                                       this->dof_handler,
                                       support_points);
  cut_cells_map.clear();
  const auto        &cell_iterator = this->dof_handler.active_cell_iterators();
  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // // Loop on all the cells and check if they are cut.
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          bool cell_is_cut;
          // is the particle index that cut the cell if it's cut.  If the cell
          // is not cut the default value is stored (0). If the cell is not cut
          // this value will never be used.
          unsigned int p;

          std::tie(cell_is_cut, p, local_dof_indices) =
            cell_cut(cell, local_dof_indices, support_points);

          // Add information about if the cell is cut "cell_is_cut" and the
          // particle id that cuts it "p" in the map.
          cut_cells_map[cell] = {cell_is_cut, p};
          std::tie(cell_is_cut, p, local_dof_indices) =
            cell_inside(cell, local_dof_indices, support_points);
          cells_inside_map[cell] = {cell_is_cut, p};
        }
    }
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::define_particles()
{
  // initialized the particles

  particles.resize(this->simulation_parameters.particlesParameters->nb);
  for (unsigned int i = 0;
       i < this->simulation_parameters.particlesParameters->nb;
       ++i)
    {
      particles[i] =
        this->simulation_parameters.particlesParameters->particles[i];
    }

  table_p.resize(particles.size());
  ib_dem.initialize(this->simulation_parameters,this->mpi_communicator);
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::refine_ib()
{
  Point<dim>                                    center_immersed;
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(*this->mapping,
                                       this->dof_handler,
                                       support_points);

  const unsigned int                   dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  const auto &cell_iterator = this->dof_handler.active_cell_iterators();
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int p = 0; p < particles.size(); ++p)
            {
              unsigned int   count_small     = 0;
              Point<dim>     center_immersed = particles[p].position;
              Tensor<1, dim> r;
              r[0] = particles[p].radius;

              bool cell_as_ib_inside =
                cell->point_inside(particles[p].position + r);
              for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
                {
                  // Count the number of dofs that are smaller or larger than
                  // the radius of the particles if all the dof are on one side
                  // the cell is not cut by the boundary meaning we don’t have
                  // to do anything

                  if (((support_points[local_dof_indices[j]] - center_immersed)
                          .norm() <= particles[p].radius *
                                       this->simulation_parameters
                                         .particlesParameters->outside_radius &&
                      (support_points[local_dof_indices[j]] - center_immersed)
                          .norm() >= particles[p].radius *
                                       this->simulation_parameters
                                         .particlesParameters->inside_radius)
                    {
                      ++count_small;
                    }
                }
              if (count_small > 0 || cell_as_ib_inside)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }
    }
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::force_on_ib()
{
  // This function defines the force and torque applied on an Immersed Boundary
  // based on the sharp edge method on a hyper_sphere of dim=2 or dim=3

  TimerOutput::Scope t(this->computing_timer, "new force_eval");

  const FEValuesExtractors::Scalar pressure(dim);
  const FEValuesExtractors::Vector velocities(0);

  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  // Define a map to all dofs and their support points
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(*this->mapping,
                                       this->dof_handler,
                                       support_points);

  // Initalize fe value objects in order to do calculation with it later
  QGauss<dim>            q_formula(this->number_quadrature_points);
  MappingQ<dim - 1, dim> local_face_map(
    this->velocity_fem_degree,
    this->simulation_parameters.fem_parameters.qmapping_all);

  FESystem<dim - 1, dim> local_face_fe(
    FE_Q<dim - 1, dim>(
      this->simulation_parameters.fem_parameters.velocity_order),
    dim,
    FE_Q<dim - 1, dim>(
      this->simulation_parameters.fem_parameters.pressure_order),
    1);
  FEValues<dim - 1, dim> fe_face_projection_values(local_face_map,
                                                   local_face_fe,
                                                   *this->face_quadrature,
                                                   update_values |
                                                     update_quadrature_points |
                                                     update_gradients |
                                                     update_JxW_values |
                                                     update_normal_vectors);



  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int dofs_per_face = this->fe->dofs_per_face;

  int    order   = this->simulation_parameters.particlesParameters->order;
  double density = this->simulation_parameters.particlesParameters->density;
  double length_ratio =
    this->simulation_parameters.particlesParameters->length_ratio;
  IBStencil<dim>      stencil;
  std::vector<double> ib_coef = stencil.coefficients(order, length_ratio);

  // Rheological model for viscosity properties
  double viscosity;
  // Cast rheological model to either a Newtonian model or one of the
  // non Newtonian models according to the physical properties
  std::shared_ptr<RheologicalModel> rheological_model =
    RheologicalModel::model_cast(
      this->simulation_parameters.physical_properties);

  const unsigned int vertices_per_face = GeometryInfo<dim>::vertices_per_face;
  const unsigned int n_q_points_face   = this->face_quadrature->size();

  std::vector<double>                      pressure_values(ib_coef.size());
  Tensor<1, dim>                           normal_vector;
  std::vector<Tensor<2, dim>>              velocity_gradients(ib_coef.size());
  std::vector<std::vector<Tensor<1, dim>>> velocity_gradients_component(dim +
                                                                        1);
  for (unsigned int i = 0; i < dim + 1; ++i)
    velocity_gradients_component[i].resize(ib_coef.size());
  Tensor<2, dim>              fluid_stress;
  Tensor<2, dim>              fluid_pressure;
  Tensor<2, dim>              fluid_stress_at_ib;
  Tensor<2, dim>              shear_rate;
  DoFHandler<dim - 1, dim>    local_face_dof_handler;
  Triangulation<dim - 1, dim> local_face_projection_triangulation;

  std::vector<Point<dim>>        vertices_of_face_projection(vertices_per_face);
  std::vector<CellData<dim - 1>> local_face_cell_data(1);

  typename Mapping<dim>::InternalDataBase mapping_data;
  // Define multiple local_dof_indices one for the cell iterator one for the
  // cell with the second point for the sharp edge stencil and one for
  // manipulation on the neighbour’s cell.

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_2(dofs_per_cell);

  std::vector<types::global_dof_index> local_face_dof_indices(dofs_per_face);

  std::vector<Point<dim>> unite_cell_interpolation_points(ib_coef.size());
  std::vector<Point<dim>> cell_interpolation_points(ib_coef.size());
  std::vector<double>     local_interp_sol(ib_coef.size());

  std::map<unsigned int, std::pair<bool, Tensor<2, dim>>> force_eval_done;

  // Define cell iterator
  const auto &cell_iterator = this->dof_handler.active_cell_iterators();

  // Clear particle force and torque
  for (unsigned int i = 0; i < particles.size(); ++i)
    {
      particles[i].forces  = 0;
      particles[i].torques = 0;
    }

  double       total_area    = 0;
  unsigned int nb_evaluation = 0;
  // Loop over all the cell
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_locally_owned())
        {
          // Particle id that cut the cell.
          unsigned int p;
          bool         cell_is_cut;
          std::tie(cell_is_cut, p) = cut_cells_map[cell];
          // If the cell is cut
          if (cell_is_cut)
            {
              // Loop over all the face of the cell that is cut.
              for (const auto face : cell->face_indices())
                {
                  auto local_face = cell->face(face);
                  cell->face(face)->get_dof_indices(local_face_dof_indices);

                  // Check if the face is cut
                  unsigned int nb_dof_inside = 0;
                  for (unsigned int j = 0; j < local_face_dof_indices.size();
                       ++j)
                    {
                      // Count the number of DOFs that are inside
                      // of the particles. If all the DOfs are on one side
                      // the cell is not cut by the boundary.
                      if ((support_points[local_face_dof_indices[j]] -
                           particles[p].position)
                            .norm() <= particles[p].radius)
                        ++nb_dof_inside;
                    }

                  // If the face is not cut and the face is outside of the IB,
                  // the face of the cell is on the boundary of the
                  // computational domain.
                  if (nb_dof_inside == 0)
                    {
                      // Projects the face on the surface of the IB. This
                      // creates a surface cell where we can evaluate the
                      // solution. Define the triangulation of the surface cell.
                      for (unsigned int i = 0; i < vertices_per_face; ++i)
                        {
                          // Define the vertices of the surface cell.
                          Tensor<1, dim, double> vertex_projection =
                            particles[p].position +
                            particles[p].radius *
                              (local_face->vertex(i) - particles[p].position) /
                              (local_face->vertex(i) - particles[p].position)
                                .norm();
                          // Create the list of vertices
                          for (unsigned int j = 0; j < dim; ++j)
                            {
                              vertices_of_face_projection[i][j] =
                                vertex_projection[j];
                            }
                          // Create the connectivity of the vertices of the cell
                          local_face_cell_data[0].vertices[i] = i;
                        }

                      local_face_cell_data[0].material_id = 0;

                      SphericalManifold<dim - 1, dim> sphere_manifold(
                        particles[p].position);
                      local_face_projection_triangulation =
                        Triangulation<dim - 1, dim>();
                      // Create a dof handler that contains the triangulation of
                      // the projection of the face on the IB. This create the
                      // surface cell on the IB
                      local_face_projection_triangulation.create_triangulation(
                        vertices_of_face_projection,
                        local_face_cell_data,
                        SubCellData());
                      local_face_projection_triangulation.set_all_manifold_ids(
                        0);
                      local_face_projection_triangulation.set_manifold(
                        0, sphere_manifold);

                      local_face_dof_handler.reinit(
                        local_face_projection_triangulation);
                      local_face_dof_handler.distribute_dofs(local_face_fe);

                      // Defined the solution on  IB surface cell by using the
                      // IB stencil to extrapolate the fluid stress tensor.

                      std::vector<Tensor<2, dim>> local_face_tensor(
                        dofs_per_face);
                      for (unsigned int i = 0;
                           i < local_face_dof_indices.size();
                           ++i)
                        {
                          const unsigned int component_i =
                            this->fe->system_to_component_index(i).first;
                          // Check if that dof already have been used to
                          // extrapolate the fluid stress tensor on the IB
                          // surface.
                          if (force_eval_done[local_face_dof_indices[i]]
                                .first == false)
                            {
                              // Only need one extrapolation by dof location;
                              if (component_i == 0)
                                {
                                  // Count the number of evaluation

                                  auto [point, interpolation_points] =
                                    stencil.points(
                                      order,
                                      length_ratio,
                                      particles[p],
                                      support_points
                                        [local_face_dof_indices[i]]);

                                  auto cell_2 =
                                    ib_done[local_face_dof_indices[i]].second;
                                  // Check if we already have the cell used to
                                  // defined the IB constraint of that dof. We
                                  // always have that information except if the
                                  // dof is not owned.
                                  if (ib_done[local_face_dof_indices[i]]
                                        .first == false)
                                    {
                                      // Get the cell use for the extrapolation
                                      auto point_to_find_cell =
                                        stencil.point_for_cell_detection(
                                          particles[p],
                                          support_points
                                            [local_face_dof_indices[i]]);

                                      try
                                        {
                                          cell_2 = LetheGridTools::
                                            find_cell_around_point_with_neighbors<
                                              dim>(this->dof_handler,
                                                   vertices_to_cell,
                                                   cell,
                                                   point_to_find_cell);
                                        }
                                      catch (...)
                                        {
                                          cell_2 = cell;
                                        }
                                    }

                                  cell_2->get_dof_indices(local_dof_indices_2);

                                  unite_cell_interpolation_points[0] =
                                    this->mapping->transform_real_to_unit_cell(
                                      cell_2, point);
                                  for (unsigned int j = 1; j < ib_coef.size();
                                       ++j)
                                    {
                                      unite_cell_interpolation_points[j] =
                                        this->mapping
                                          ->transform_real_to_unit_cell(
                                            cell_2,
                                            interpolation_points[j - 1]);
                                    }

                                  fluid_stress_at_ib = 0;

                                  // Create a quadrature that is based on the IB
                                  // stencil
                                  Quadrature<dim> q_local(
                                    unite_cell_interpolation_points, ib_coef);
                                  FEValues<dim> fe_values_cell2(
                                    *this->fe,
                                    q_local,
                                    update_quadrature_points |
                                      update_gradients | update_values);

                                  // Evaluate the relevant information at the
                                  // quadrature points to do the extrapolation.
                                  fe_values_cell2.reinit(cell_2);
                                  fe_values_cell2[velocities]
                                    .get_function_gradients(
                                      this->evaluation_point,
                                      velocity_gradients);
                                  fe_values_cell2[pressure].get_function_values(
                                    this->evaluation_point, pressure_values);

                                  // Extrapolate the fluid stress tensor on the
                                  // surface of the IB.
                                  for (unsigned int k = 0; k < ib_coef.size();
                                       ++k)
                                    {
                                      fluid_pressure = 0;
                                      for (int d = 0; d < dim; ++d)
                                        {
                                          fluid_pressure[d][d] =
                                            pressure_values[k];
                                        }
                                      shear_rate =
                                        velocity_gradients[k] +
                                        transpose(velocity_gradients[k]);


                                      const double shear_rate_magnitude =
                                        calculate_shear_rate_magnitude(
                                          shear_rate);

                                      std::map<field, double> field_values;
                                      field_values[field::shear_rate] =
                                        shear_rate_magnitude;

                                      viscosity =
                                        rheological_model->value(field_values);

                                      fluid_stress =
                                        viscosity * shear_rate - fluid_pressure;

                                      fluid_stress_at_ib +=
                                        fluid_stress * ib_coef[k];
                                    }
                                  // Store the stress tensor that results from
                                  // the extrapolation in the local evaluation
                                  // vector of the IB surface cell and in a map
                                  // that is used if the same extrapolation is
                                  // needed in another face.
                                  local_face_tensor[i] = fluid_stress_at_ib;
                                  force_eval_done[local_face_dof_indices[i]] =
                                    std::make_pair(true, fluid_stress_at_ib);
                                }
                            }
                          else
                            {
                              // Use the results from a previously evaluated
                              // extrapolation. This step comes with an error du
                              // to the curvature of the surface in Q2 and
                              // higher order elements.
                              local_face_tensor[i] =
                                force_eval_done[local_face_dof_indices[i]]
                                  .second;
                            }
                        }
                      // Use the extrapolation of fluid stress tensor at the
                      // dof location of the IB surface cell to integrate the
                      // stress tensor on the surface of the IB
                      auto local_face_tensor_old = local_face_tensor;
                      local_face_tensor.clear();
                      local_face_tensor.resize(local_face_dof_indices.size());
                      for (const auto &projection_cell_face :
                           local_face_dof_handler.active_cell_iterators())
                        {
                          fe_face_projection_values.reinit(
                            projection_cell_face);
                          std::vector<Point<dim>> q_points =
                            fe_face_projection_values.get_quadrature_points();
                          if (this->simulation_parameters.fem_parameters
                                .velocity_order > 1)
                            {
                              FullMatrix<double> interpolation_matrix(
                                local_face_dof_indices.size(),
                                local_face_dof_indices.size());
                              FullMatrix<double> inv_interpolation_matrix(
                                local_face_dof_indices.size(),
                                local_face_dof_indices.size());

                              // Define the interpolation matrix of the surface
                              // cell
                              for (unsigned int i = 0;
                                   i < local_face_dof_indices.size();
                                   ++i)
                                {
                                  Point<dim> point_projection =
                                    particles[p].position +
                                    particles[p].radius *
                                      (support_points
                                         [local_face_dof_indices[i]] -
                                       particles[p].position) /
                                      (support_points
                                         [local_face_dof_indices[i]] -
                                       particles[p].position)
                                        .norm();

                                  auto projected_point_unit =
                                    local_face_map.transform_real_to_unit_cell(
                                      projection_cell_face, point_projection);
                                  for (unsigned int j = 0;
                                       j < local_face_dof_indices.size();
                                       ++j)
                                    {
                                      interpolation_matrix[i][j] = 0;
                                      if (this->fe
                                            ->face_system_to_component_index(j)
                                            .first ==
                                          this->fe
                                            ->face_system_to_component_index(i)
                                            .first)
                                        interpolation_matrix[i][j] +=
                                          local_face_fe.shape_value(
                                            j, projected_point_unit);
                                    }
                                }
                              inv_interpolation_matrix.invert(
                                interpolation_matrix);
                              // Define the value of the fluid stress tensor on
                              // the surface cell at the DOF support points
                              // location.
                              for (unsigned int i = 0;
                                   i < local_face_dof_indices.size();
                                   ++i)
                                {
                                  for (unsigned int j = 0;
                                       j < local_face_dof_indices.size();
                                       ++j)
                                    {
                                      local_face_tensor[i] +=
                                        inv_interpolation_matrix[i][j] *
                                        local_face_tensor_old[j];
                                    }
                                }
                            }
                          else
                            {
                              local_face_tensor = local_face_tensor_old;
                            }
                          for (unsigned int q = 0; q < n_q_points_face; q++)
                            {
                              // Evaluate the total surface
                              // Redefined the normal at the quadrature point
                              // since we dont control the orientation of the
                              // cell.
                              normal_vector =
                                (q_points[q] - particles[p].position) /
                                (q_points[q] - particles[p].position).norm();
                              fluid_stress        = 0;
                              double local_weight = 0;
                              // Integrate
                              for (unsigned int i = 0;
                                   i < local_face_dof_indices.size();
                                   ++i)
                                {
                                  const unsigned int component_i =
                                    local_face_fe.system_to_component_index(i)
                                      .first;
                                  if (component_i == 0)
                                    {
                                      fluid_stress +=
                                        fe_face_projection_values.shape_value(
                                          i, q) *
                                        local_face_tensor[i];
                                      total_area +=
                                        fe_face_projection_values.JxW(q) *
                                        fe_face_projection_values.shape_value(
                                          i, q);
                                      local_weight +=
                                        fe_face_projection_values.shape_value(
                                          i, q);
                                    }
                                }

                              auto force = fluid_stress * normal_vector *
                                           fe_face_projection_values.JxW(q);
                              if (force.norm() > 0)
                                {
                                  nb_evaluation += local_weight;
                                }



                              // Add the local contribution of this surface
                              // cell.
                              particles[p].forces += force;

                              auto distance =
                                q_points[q] - particles[p].position;
                              if (dim == 2)
                                {
                                  particles[p].torques[0] += 0.;
                                  particles[p].torques[1] += 0.;
                                  particles[p].torques[2] +=
                                    distance[0] * force[1] -
                                    distance[1] * force[0];
                                }
                              else if (dim == 3)
                                {
                                  particles[p].torques[0] +=
                                    distance[1] * force[2] -
                                    distance[2] * force[1];
                                  particles[p].torques[1] +=
                                    distance[2] * force[0] -
                                    distance[0] * force[2];
                                  particles[p].torques[2] +=
                                    distance[0] * force[1] -
                                    distance[1] * force[0];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  // Sums the force evaluation on each of the processor.
  for (unsigned int i = 0; i < particles.size(); ++i)
    {
      particles[i].forces =
        Utilities::MPI::sum(particles[i].forces, this->mpi_communicator) *
        density;
      particles[i].torques =
        Utilities::MPI::sum(particles[i].torques, this->mpi_communicator) *
        density;
    }
  // total_area = Utilities::MPI::sum(total_area, this->mpi_communicator);
  // std::cout << "total area " << total_area << std::endl;

}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::write_force_ib()
{
  TimerOutput::Scope t(this->computing_timer, "output_forces_ib");
  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      {
        if (this->this_mpi_process == 0)
          {
            std::string filename =
              this->simulation_parameters.simulation_control.output_folder +
              this->simulation_parameters.particlesParameters
                ->ib_force_output_file +
              "." + Utilities::int_to_string(p, 2) + ".dat";
            std::ofstream output(filename.c_str());

            table_p[p].write_text(output);
            std::string filename_residual =
              this->simulation_parameters.simulation_control.output_folder +
              "residual" + ".dat";
            std::ofstream output_residual(filename_residual.c_str());
            table_residual.write_text(output_residual);
          }
        MPI_Barrier(this->mpi_communicator);
      }
    }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::postprocess_fd(bool firstIter)
{
  if (this->simulation_control->is_output_iteration())
    {
      this->write_output_results(this->present_solution);
    }

  bool enable =
    this->simulation_parameters.analytical_solution->calculate_error();
  this->simulation_parameters.analytical_solution->set_enable(false);
  NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
    postprocess_fd(firstIter);
  this->simulation_parameters.analytical_solution->set_enable(enable);
  // Calculate the error with respect to the analytical solution
  if (!firstIter &&
      this->simulation_parameters.analytical_solution->calculate_error())
    {
      // Update the time of the exact solution to the actual time
      this->exact_solution->set_time(
        this->simulation_control->get_current_time());
      const std::pair<double, double> error =
        this->calculate_L2_error_particles();

      if (this->simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        {
          this->error_table.add_value(
            "cells", this->triangulation->n_global_active_cells());
          this->error_table.add_value("error_velocity", error.first);
          this->error_table.add_value("error_pressure", error.second);

          auto summary = this->computing_timer.get_summary_data(
            this->computing_timer.total_wall_time);
          double total_time = 0;
          for (auto it = summary.begin(); it != summary.end(); ++it)
            {
              total_time += summary[it->first];
            }
          this->error_table.add_value("total_time", total_time);
        }
      else
        {
          this->error_table.add_value(
            "time", this->simulation_control->get_current_time());
          this->error_table.add_value("error_velocity", error.first);
          this->error_table.add_value("error_pressure", error.second);

          if (this->simulation_parameters.timer.write_time_in_error_table)
            {
              auto summary = this->computing_timer.get_summary_data(
                this->computing_timer.total_wall_time);
              double total_time = 0;
              for (auto it = summary.begin(); it != summary.end(); ++it)
                {
                  total_time += summary[it->first];
                }
              this->error_table.add_value("total_time", total_time);
            }
        }
      if (this->simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error velocity : " << error.first
                      << " L2 error pressure: " << error.second << std::endl;
        }
    }
}


template <int dim>
std::pair<double, double>
GLSSharpNavierStokesSolver<dim>::calculate_L2_error_particles()
{
  TimerOutput::Scope t(this->computing_timer, "error");
  QGauss<dim>        quadrature_formula(this->number_quadrature_points + 1);
  FEValues<dim>      fe_values(*this->mapping,
                          *this->fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
  FEFaceValues<dim>  fe_face_values(*this->mapping,
                                   *this->fe,
                                   *this->face_quadrature,
                                   update_values | update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  const unsigned int dofs_per_cell =
    this->fe->dofs_per_cell; // This gives you dofs per cell
  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points      = quadrature_formula.size();
  const unsigned int n_q_points_face = this->face_quadrature->size();

  std::vector<Vector<double>> q_exactSol(n_q_points, Vector<double>(dim + 1));

  std::vector<Tensor<1, dim>> local_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> local_face_velocity_values(n_q_points_face);
  std::vector<double>         local_pressure_values(n_q_points);
  std::vector<double>         div_phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);

  Function<dim> *l_exact_solution = this->exact_solution;

  Point<dim>                                    center_immersed;
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(*this->mapping,
                                       this->dof_handler,
                                       support_points);

  double l2errorU                  = 0.;
  double l2errorU_boundary         = 0.;
  double l2errorP                  = 0.;
  double total_velocity_divergence = 0.;
  double pressure_integral         = 0;
  double exact_pressure_integral   = 0;
  double volume                    = 0;

  // loop over elements
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();

  // loop over elements to calculate average pressure
  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          bool cell_is_cut;
          // std::ignore is used because we don't care about what particle cut
          // the cell.
          std::tie(cell_is_cut, std::ignore) = cut_cells_map[cell];
          bool cell_is_inside;
          std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];

          if (cell_is_cut == false && cell_is_inside == false)
            {
              auto &evaluation_point = this->evaluation_point;
              fe_values.reinit(cell);

              fe_values[pressure].get_function_values(evaluation_point,
                                                      local_pressure_values);
              // Get the exact solution at all gauss points
              l_exact_solution->vector_value_list(
                fe_values.get_quadrature_points(), q_exactSol);


              // Retrieve the effective "connectivity matrix" for this element
              cell->get_dof_indices(local_dof_indices);

              for (unsigned int q = 0; q < n_q_points; q++)
                {
                  pressure_integral +=
                    local_pressure_values[q] * fe_values.JxW(q);
                  exact_pressure_integral +=
                    q_exactSol[q][dim] * fe_values.JxW(q);
                  volume += fe_values.JxW(q);
                }
            }
        }
    }

  pressure_integral =
    Utilities::MPI::sum(pressure_integral, this->mpi_communicator);
  exact_pressure_integral =
    Utilities::MPI::sum(exact_pressure_integral, this->mpi_communicator);
  volume = Utilities::MPI::sum(volume, this->mpi_communicator);


  double average_pressure       = pressure_integral / volume;
  double average_exact_pressure = exact_pressure_integral / volume;
  cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();

  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);

          bool cell_is_cut;
          // std::ignore is used because we don't care about what particle cut
          // the cell.
          std::tie(cell_is_cut, std::ignore) = cut_cells_map[cell];

          bool cell_is_inside;
          std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];
          if (cell->at_boundary() &&
              this->check_existance_of_bc(
                BoundaryConditions::BoundaryType::function_weak))
            {
              for (unsigned int i_bc = 0;
                   i_bc < this->simulation_parameters.boundary_conditions.size;
                   ++i_bc)
                {
                  if (this->simulation_parameters.boundary_conditions
                        .type[i_bc] ==
                      BoundaryConditions::BoundaryType::function_weak)
                    {
                      for (const auto face : cell->face_indices())
                        {
                          if (cell->face(face)->at_boundary())
                            {
                              unsigned int boundary_id =
                                cell->face(face)->boundary_id();
                              if (boundary_id ==
                                  this->simulation_parameters
                                    .boundary_conditions.id[i_bc])
                                {
                                  NavierStokesFunctionDefined<dim> function_v(
                                    &this->simulation_parameters
                                       .boundary_conditions
                                       .bcFunctions[boundary_id]
                                       .u,
                                    &this->simulation_parameters
                                       .boundary_conditions
                                       .bcFunctions[boundary_id]
                                       .v,
                                    &this->simulation_parameters
                                       .boundary_conditions
                                       .bcFunctions[boundary_id]
                                       .w);

                                  fe_face_values.reinit(cell, face);
                                  fe_face_values[velocities]
                                    .get_function_values(
                                      this->present_solution,
                                      local_face_velocity_values);
                                  for (unsigned int q = 0; q < n_q_points_face;
                                       q++)
                                    {
                                      double u_x =
                                        local_face_velocity_values[q][0];
                                      double u_y =
                                        local_face_velocity_values[q][1];

                                      double u_x_a = function_v.value(
                                        fe_face_values.quadrature_point(q), 0);
                                      double u_y_a = function_v.value(
                                        fe_face_values.quadrature_point(q), 1);

                                      l2errorU_boundary +=
                                        ((u_x - u_x_a) * (u_x - u_x_a) +
                                         (u_y - u_y_a) * (u_y - u_y_a)) *
                                        fe_face_values.JxW(q);
                                      if (dim == 3)
                                        {
                                          double u_z =
                                            local_face_velocity_values[q][2];
                                          double u_z_a = function_v.value(
                                            fe_face_values.quadrature_point(q),
                                            2);
                                          l2errorU_boundary +=
                                            (u_z - u_z_a) * (u_z - u_z_a) *
                                            fe_face_values.JxW(q);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

          if (cell_is_cut == false)
            {
              auto &evaluation_point = this->evaluation_point;
              auto &present_solution = this->present_solution;
              fe_values.reinit(cell);
              fe_values[velocities].get_function_values(present_solution,
                                                        local_velocity_values);
              fe_values[pressure].get_function_values(present_solution,
                                                      local_pressure_values);
              fe_values[velocities].get_function_gradients(
                evaluation_point, present_velocity_gradients);


              // Retrieve the effective "connectivity matrix" for this element
              cell->get_dof_indices(local_dof_indices);

              // Get the exact solution at all gauss points
              l_exact_solution->vector_value_list(
                fe_values.get_quadrature_points(), q_exactSol);
              for (unsigned int q = 0; q < n_q_points; q++)
                {
                  double present_velocity_divergence =
                    trace(present_velocity_gradients[q]);
                  double mass_source =
                    this->simulation_parameters.source_term
                      ->navier_stokes_source.value(
                        fe_values.get_quadrature_points()[q], dim);


                  // Find the values of x and u_h (the finite element solution)
                  // at the quadrature points
                  double ux_sim   = local_velocity_values[q][0];
                  double ux_exact = q_exactSol[q][0];

                  double uy_sim   = local_velocity_values[q][1];
                  double uy_exact = q_exactSol[q][1];
                  l2errorU += (ux_sim - ux_exact) * (ux_sim - ux_exact) *
                              fe_values.JxW(q);
                  l2errorU += (uy_sim - uy_exact) * (uy_sim - uy_exact) *
                              fe_values.JxW(q);

                  if (dim == 3)
                    {
                      double uz_sim   = local_velocity_values[q][2];
                      double uz_exact = q_exactSol[q][2];
                      l2errorU += (uz_sim - uz_exact) * (uz_sim - uz_exact) *
                                  fe_values.JxW(q);
                    }
                  if (cell_is_inside == false)
                    {
                      total_velocity_divergence +=
                        (present_velocity_divergence - mass_source) *
                        fe_values.JxW(q);
                      double p_sim =
                        local_pressure_values[q] - average_pressure;
                      double p_exact =
                        q_exactSol[q][dim] - average_exact_pressure;
                      l2errorP += (p_sim - p_exact) * (p_sim - p_exact) *
                                  fe_values.JxW(q);
                    }
                }
            }
        }
    }
  l2errorU = Utilities::MPI::sum(l2errorU, this->mpi_communicator);
  l2errorU_boundary =
    Utilities::MPI::sum(l2errorU_boundary, this->mpi_communicator);
  l2errorP = Utilities::MPI::sum(l2errorP, this->mpi_communicator);
  total_velocity_divergence =
    Utilities::MPI::sum(total_velocity_divergence, this->mpi_communicator);

  this->pcout << "div u : " << total_velocity_divergence << std::endl;
  if (this->check_existance_of_bc(
        BoundaryConditions::BoundaryType::function_weak))
    {
      this->pcout << "error at the weak BC : " << std::sqrt(l2errorU_boundary)
                  << std::endl;
    }

  return std::make_pair(std::sqrt(l2errorU), std::sqrt(l2errorP));
}




template <int dim>
void
GLSSharpNavierStokesSolver<dim>::particles_dem()
{ // add refilling containers
  using numbers::PI;
  Tensor<1, dim> g   = this->simulation_parameters.particlesParameters.gravity;
  double         rho = this->simulation_parameters.particlesParameters.density;
  double         dt  = this->simulation_control->get_time_steps_vector()[0];
  double         dt_dem             = dt / 10000;
  // local time for the dem step
  std::vector<Tensor<1, dim>> contact_force(particles.size());
  std::vector<Tensor<1, 3>>   contact_torque(particles.size());
  std::vector<Tensor<1, 3>>   contact_wall_torque(particles.size());
  std::vector<Tensor<1, dim>> contact_wall_force(particles.size());

  std::vector<Tensor<1, dim>> current_fluid_force(particles.size());
  std::vector<Tensor<1, 3>>   current_fluid_torque(particles.size());

  std::vector<Tensor<1, dim>> velocity(particles.size());
  std::vector<Point<dim>>     position(particles.size());

  double         t = 0;
  Tensor<1, dim> gravity;

  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      dem_particles[p_i].position    = particles[p_i].last_position;
      dem_particles[p_i].velocity    = particles[p_i].last_velocity;
      dem_particles[p_i].omega       = particles[p_i].last_omega;
      particles[p_i].impulsion       = 0;
      particles[p_i].omega_impulsion = 0;
    }

  while (t < dt)
    {
      current_fluid_force.clear();
      current_fluid_force.resize(particles.size());
      current_fluid_torque.clear();
      current_fluid_torque.resize(particles.size());



      contact_torque.clear();
      contact_torque.resize(particles.size());
      contact_force.clear();
      contact_force.resize(particles.size());
      contact_wall_force.clear();
      contact_wall_force.resize(particles.size());
      contact_wall_torque.clear();
      contact_wall_torque.resize(particles.size());



      // Calculate particle-particle and particle-wall contact force
      calculate_pp_contact_force(dt_dem, contact_force, contact_torque);
      calculate_pw_contact_force(dt_dem,
                                 contact_wall_force,
                                 contact_wall_torque);


      for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
        {
          auto inv_inertia = invert(particles[p_i].inertia);
          if (dim == 2)
            gravity =
              g * (particles[p_i].mass -
                   particles[p_i].radius * particles[p_i].radius * PI * rho);
          if (dim == 3)
            {
              gravity =
                g * (particles[p_i].mass - 4.0 / 3.0 * particles[p_i].radius *
                                             particles[p_i].radius *
                                             particles[p_i].radius * PI * rho);
            }

          current_fluid_force[p_i] =
            particles[p_i].last_forces +
            (particles[p_i].forces - particles[p_i].last_forces) * t / dt;
          current_fluid_torque[p_i] =
            particles[p_i].last_torques +
            (particles[p_i].torques - particles[p_i].last_torques) * t / dt;



          dem_particles[p_i].velocity =
            dem_particles[p_i].velocity +
            (current_fluid_force[p_i] + contact_force[p_i] +
             contact_wall_force[p_i] + gravity) *
              dt_dem / particles[p_i].mass;
          dem_particles[p_i].position =
            dem_particles[p_i].position + dem_particles[p_i].velocity * dt_dem;

          dem_particles[p_i].omega =
            dem_particles[p_i].omega +
            inv_inertia *
              (current_fluid_torque[p_i] + contact_torque[p_i] +
               contact_wall_torque[p_i]) *
              dt_dem;

          particles[p_i].impulsion +=
            (current_fluid_force[p_i] + contact_force[p_i] +
             contact_wall_force[p_i] + gravity) *
            dt_dem;
          particles[p_i].omega_impulsion +=
            (current_fluid_torque[p_i] + contact_torque[p_i] +
             contact_wall_torque[p_i]) *
            dt_dem;
        }
      t += dt_dem;
    }
  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      particles[p_i].position = dem_particles[p_i].position;
    }
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::calculate_pp_contact_force(
  const double                &dt_dem,
  std::vector<Tensor<1, dim>> &contact_force,
  std::vector<Tensor<1, 3>>   &contact_torque)
{
  for (auto &particle_one : dem_particles)
    {
      for (auto &particle_two : dem_particles)
        {
          if (particle_one.particle_id != particle_two.particle_id and
              particle_one.particle_id < particle_two.particle_id)
            {
              const Point<dim> particle_one_location = particle_one.position;
              const Point<dim> particle_two_location = particle_two.position;
              ContactTangentialHistory contact_history;
              ContactTangentialHistory contact_info;
              try
                {
                  contact_info = pp_contact_map[particle_one.particle_id]
                                               [particle_two.particle_id];
                }
              catch (...)
                {
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_history.tangential_overlap[d]           = 0;
                      contact_history.tangential_relative_velocity[d] = 0;
                    }
                  pp_contact_map[particle_one.particle_id]
                                [particle_two.particle_id] = contact_history;
                }

              double effective_mass = (particle_one.mass * particle_two.mass) /
                                      (particle_one.mass + particle_two.mass);
              double effective_radius =
                (particle_one.radius * particle_two.radius) /
                (particle_one.radius + particle_two.radius);


              // Calculation of normal overlap
              double normal_overlap =
                (particle_one.radius + particle_two.radius) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0)
                // This means that the adjacent particles are in contact
                {
                  // Calculation of the contact vector (vector from particle one
                  // to particle two
                  auto contact_vector =
                    particle_two_location - particle_one_location;

                  // Using contact_vector, the contact normal vector is obtained
                  auto normal_unit_vector =
                    contact_vector / contact_vector.norm();
                  Tensor<1, 3> normal;
                  if (dim == 2)
                    {
                      normal[0] = normal_unit_vector[0];
                      normal[1] = normal_unit_vector[1];
                    }
                  if (dim == 3)
                    {
                      normal[0] = normal_unit_vector[0];
                      normal[1] = normal_unit_vector[1];
                      normal[2] = normal_unit_vector[2];
                    }

                  // Defining relative contact velocity
                  Tensor<1, dim> contact_relative_velocity;

                  if (dim == 3)
                    {
                      // Calculation of contact relative velocity
                      Tensor<1, 3> rotational_velocity = cross_product_3d(
                        (particle_one.radius * particle_one.omega +
                         particle_two.radius * particle_two.omega),
                        normal);

                      contact_relative_velocity[0] =
                        (particle_one.velocity - particle_two.velocity)[0] +
                        rotational_velocity[0];
                      contact_relative_velocity[1] =
                        (particle_one.velocity - particle_two.velocity)[1] +
                        rotational_velocity[1];
                      contact_relative_velocity[2] =
                        (particle_one.velocity - particle_two.velocity)[2] +
                        rotational_velocity[2];
                    }
                  else
                    {
                      // TODO: Correct this after correcting dem_2d
                      contact_relative_velocity =
                        particle_one.velocity - particle_two.velocity;
                    }

                  auto normal_relative_velocity_value =
                    contact_relative_velocity * normal_unit_vector;
                  Tensor<1, dim> normal_relative_velocity =
                    normal_relative_velocity_value * normal_unit_vector;

                  // Calculation of tangential relative velocity
                  Tensor<1, dim> tangential_relative_velocity =
                    contact_relative_velocity - normal_relative_velocity;

                  Tensor<1, dim> modified_tangential_overlap =
                    contact_info.tangential_overlap +
                    contact_info.tangential_relative_velocity * dt_dem;

                  // Updating the contact_info container based on the new
                  // calculated values
                  contact_history.tangential_overlap =
                    modified_tangential_overlap;
                  contact_history.tangential_relative_velocity =
                    tangential_relative_velocity;
                  pp_contact_map[particle_one.particle_id]
                                [particle_two.particle_id] = contact_history;


                  const double effective_mass =
                    (particle_one.mass * particle_two.mass) /
                    (particle_one.mass + particle_two.mass);
                  const double effective_radius =
                    (particle_one.radius * particle_two.radius) /
                    (particle_one.radius + particle_two.radius);
                  const double effective_youngs_modulus =
                    (particle_one.youngs_modulus *
                     particle_two.youngs_modulus) /
                    ((particle_two.youngs_modulus *
                      (1 - particle_one.poisson_ratio *
                             particle_one.poisson_ratio)) +
                     (particle_one.youngs_modulus *
                      (1 - particle_two.poisson_ratio *
                             particle_two.poisson_ratio)) +
                     DBL_MIN);

                  const double effective_shear_modulus =
                    (particle_one.youngs_modulus *
                     particle_two.youngs_modulus) /
                    (2 * ((particle_two.youngs_modulus *
                           (2 - particle_one.poisson_ratio) *
                           (1 + particle_one.poisson_ratio)) +
                          (particle_one.youngs_modulus *
                           (2 - particle_two.poisson_ratio) *
                           (1 + particle_two.poisson_ratio))) +
                     DBL_MIN);


                  const double effective_coefficient_of_restitution =
                    2 * particle_one.restitution_coefficient *
                    particle_two.restitution_coefficient /
                    (particle_one.restitution_coefficient +
                     particle_two.restitution_coefficient + DBL_MIN);

                  const double effective_coefficient_of_friction =
                    2 * particle_one.friction_coefficient *
                    particle_two.friction_coefficient /
                    (particle_one.friction_coefficient +
                     particle_two.friction_coefficient + DBL_MIN);


                  const double effective_coefficient_of_rolling_friction =
                    2 * particle_one.rolling_friction_coefficient *
                    particle_two.rolling_friction_coefficient /
                    (particle_one.rolling_friction_coefficient +
                     particle_two.rolling_friction_coefficient + DBL_MIN);

                  const double restitution_coefficient_particle_log =
                    std::log(effective_coefficient_of_restitution);

                  const double model_parameter_beta =
                    restitution_coefficient_particle_log /
                    sqrt(restitution_coefficient_particle_log *
                           restitution_coefficient_particle_log +
                         9.8696);

                  const double radius_times_overlap_sqrt =
                    sqrt(effective_radius * normal_overlap);

                  const double model_parameter_sn =
                    2 * effective_youngs_modulus * radius_times_overlap_sqrt;

                  const double model_parameter_st =
                    8 * effective_youngs_modulus * radius_times_overlap_sqrt;

                  // Calculation of normal and tangential spring and dashpot
                  // constants using particle properties
                  const double normal_spring_constant =
                    0.66665 * model_parameter_sn;
                  const double normal_damping_constant =
                    -1.8257 * model_parameter_beta *
                    sqrt(model_parameter_sn * effective_mass);
                  const double tangential_spring_constant =
                    8 * effective_shear_modulus * radius_times_overlap_sqrt +
                    DBL_MIN;
                  const double tangential_damping_constant =
                    normal_damping_constant *
                    sqrt(model_parameter_st / model_parameter_sn);

                  // Calculation of normal force using spring and dashpot normal
                  // forces
                  Tensor<1, dim> normal_force =
                    ((normal_spring_constant * normal_overlap) *
                     normal_unit_vector) +
                    ((normal_damping_constant *
                      normal_relative_velocity_value) *
                     normal_unit_vector);


                  // Calculation of tangential force using spring and dashpot
                  // tangential forces. Since we need dashpot tangential force
                  // in the gross sliding again, we define it as a separate
                  // variable
                  Tensor<1, dim> damping_tangential_force =
                    tangential_damping_constant *
                    contact_info.tangential_relative_velocity;
                  Tensor<1, dim> tangential_force =
                    (tangential_spring_constant *
                     contact_info.tangential_overlap) +
                    damping_tangential_force;

                  const double coulomb_threshold =
                    effective_coefficient_of_friction * normal_force.norm();

                  // Check for gross sliding
                  if (tangential_force.norm() > coulomb_threshold)
                    {
                      // Gross sliding occurs and the tangential overlap and
                      // tangnetial force are limited to Coulumb's criterion
                      contact_info.tangential_overlap =
                        (coulomb_threshold *
                           (tangential_force /
                            (tangential_force.norm() + DBL_MIN)) -
                         damping_tangential_force) /
                        (tangential_spring_constant + DBL_MIN);

                      tangential_force = (tangential_spring_constant *
                                          contact_info.tangential_overlap) +
                                         damping_tangential_force;
                    }

                  // Calculation of torque
                  // Torque caused by tangential force (tangential_torque)
                  Tensor<1, 3> tangential_torque;
                  if (dim == 2)
                    {
                      Tensor<1, 3> tangential_force_3d;
                      tangential_force_3d[0] = tangential_force[0];
                      tangential_force_3d[1] = tangential_force[1];
                      tangential_torque =
                        cross_product_3d((effective_radius * normal),
                                         tangential_force_3d);
                    }
                  if (dim == 3)
                    {
                      Tensor<1, 3> tangential_force_3d;
                      tangential_force_3d[0] = tangential_force[0];
                      tangential_force_3d[1] = tangential_force[1];
                      tangential_force_3d[2] = tangential_force[2];
                      tangential_torque =
                        cross_product_3d((effective_radius * normal),
                                         tangential_force_3d);
                    }
                  // Rolling resistance torque using viscous rolling resistance
                  // model
                  auto omega_ij = particle_one.omega - particle_two.omega;
                  auto omega_ij_direction =
                    omega_ij / (omega_ij.norm() + DBL_MIN);

                  Tensor<1, 3> v_omega =
                    cross_product_3d(particle_one.omega,
                                     particle_one.radius * normal) -
                    cross_product_3d(particle_two.omega,
                                     particle_two.radius * -normal);

                  // Calculation of rolling resistance torque
                  Tensor<1, 3> rolling_resistance_torque =
                    -effective_coefficient_of_rolling_friction *
                    effective_radius * normal_force.norm() * v_omega.norm() *
                    omega_ij_direction;



                  contact_force[particle_one.particle_id] -=
                    (normal_force + tangential_force);
                  contact_force[particle_two.particle_id] +=
                    (normal_force + tangential_force);

                  contact_torque[particle_one.particle_id] =
                    contact_torque[particle_one.particle_id] -
                    tangential_torque + rolling_resistance_torque;
                  contact_torque[particle_two.particle_id] =
                    contact_torque[particle_two.particle_id] -
                    tangential_torque - rolling_resistance_torque;
                }

              else
                {
                  // if the adjacent pair is not in contact anymore
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_history.tangential_overlap[d]           = 0;
                      contact_history.tangential_relative_velocity[d] = 0;
                    }
                  pp_contact_map[particle_one.particle_id].erase(
                    particle_two.particle_id);
                }
            }
        }
    }
}
template <int dim>
void
GLSSharpNavierStokesSolver<dim>::update_particles_boundary_contact()
{
  TimerOutput::Scope t(this->computing_timer,
                       "update_particles_boundary_contact");
  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      boundary_cells[p_i].clear();
      auto cells_at_boundary =
        LetheGridTools::find_boundary_cell_in_sphere(this->dof_handler,
                                                     particles[p_i].position,
                                                     particles[p_i].radius *
                                                       1.5);
      // Initialize a simple quadrature for on the system. This will be used to
      // obtain a single sample point on the boundary faces
      const FE_Q<dim> fe(1);
      for (unsigned int i = 0; i < cells_at_boundary.size(); ++i)
        {
          QGauss<dim - 1>   face_quadrature_formula(1);
          unsigned int      n_face_q_points = face_quadrature_formula.size();
          FEFaceValues<dim> fe_face_values(fe,
                                           face_quadrature_formula,
                                           update_values |
                                             update_quadrature_points |
                                             update_normal_vectors);
          for (int face_id = 0;
               face_id < int(GeometryInfo<dim>::faces_per_cell);
               ++face_id)
            {
              if (cells_at_boundary[i]->face(face_id)->at_boundary())
                {
                  fe_face_values.reinit(cells_at_boundary[i], face_id);

                  for (unsigned int f_q_point = 0; f_q_point < n_face_q_points;
                       ++f_q_point)
                    {
                      Tensor<1, dim> normal_vector =
                        -fe_face_values.normal_vector(f_q_point);
                      BoundaryCellsInfo boundary_information;
                      boundary_information.normal_vector = normal_vector;
                      boundary_information.point_on_boundary =
                        fe_face_values.quadrature_point(f_q_point);
                      boundary_cells[p_i].push_back(boundary_information);
                    }
                }
            }
        }


            auto
         global_boundary_cell=Utilities::MPI::all_gather(this->mpi_communicator,boundary_cells[p_i]);
            boundary_cells[p_i].clear();
            for(unsigned int i=0;i<global_boundary_cell.size();++i){
                boundary_cells[p_i].insert(boundary_cells[p_i].end(),
         global_boundary_cell[i].begin(), global_boundary_cell[i].end());
              }

    }
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::calculate_pw_contact_force(
  const double                &dt_dem,
  std::vector<Tensor<1, dim>> &contact_force,
  std::vector<Tensor<1, 3>>   &contact_torque)
{
  double wall_youngs_modulus =
    this->simulation_parameters.particlesParameters.wall_youngs_modulus;
  double wall_poisson_ratio =
    this->simulation_parameters.particlesParameters.wall_poisson_ratio;
  double wall_rolling_friction_coefficient =
    this->simulation_parameters.particlesParameters
      .wall_rolling_friction_coefficient;
  double wall_friction_coefficient =
    this->simulation_parameters.particlesParameters.wall_friction_coefficient;
  double wall_restitution_coefficient =
    this->simulation_parameters.particlesParameters
      .wall_restitution_coefficient;


  for (auto &particle : dem_particles)
    {
      unsigned int boundary_index = 0;
      double       best_dist      = DBL_MAX;
      unsigned int best_index;

      for (auto &boundary_cell_iter : boundary_cells[particle.particle_id])
        {
          double dist =
            (boundary_cell_iter.point_on_boundary - particle.position).norm();
          if (dist < best_dist)
            {
              best_dist  = dist;
              best_index = boundary_index;
            }
          boundary_index += 1;
        }
      if (boundary_cells[particle.particle_id].size() > 0)
        {
          auto &boundary_cell =
            boundary_cells[particle.particle_id][best_index];


          auto                     boundary_cell_information = boundary_cell;
          ContactTangentialHistory contact_history;
          ContactTangentialHistory contact_info;
          try
            {
              contact_info =
                pw_contact_map[particle.particle_id][boundary_index];
            }
          catch (...)
            {
              for (int d = 0; d < dim; ++d)
                {
                  contact_history.tangential_overlap[d]           = 0;
                  contact_history.tangential_relative_velocity[d] = 0;
                }
              pw_contact_map[particle.particle_id][boundary_index] =
                contact_history;
            }

          auto normal_vector     = boundary_cell_information.normal_vector;
          auto point_on_boundary = boundary_cell_information.point_on_boundary;

          Tensor<1, 3> normal;
          if (dim == 2)
            {
              normal[0] = normal_vector[0];
              normal[1] = normal_vector[1];
            }
          if (dim == 3)
            {
              normal[0] = normal_vector[0];
              normal[1] = normal_vector[1];
              normal[2] = normal_vector[2];
            }

          // A vector (point_to_particle_vector) is defined which connects the
          // center of particle to the point_on_boundary. This vector will then
          // be projected on the normal vector of the boundary to obtain the
          // particle-wall distance
          Tensor<1, dim> point_to_particle_vector =
            particle.position - point_on_boundary;

          // Finding the projected vector on the normal vector of the boundary.
          // Here we have used the private function find_projection. Using this
          // projected vector, the particle-wall distance is calculated
          Tensor<1, dim> projected_vector =
            find_projection(point_to_particle_vector, normal_vector);

          double normal_overlap = particle.radius - (projected_vector.norm());

          if (normal_overlap > 0)
            {
              // Defining relative contact velocity
              Tensor<1, dim> contact_relative_velocity;
              if (dim == 3)
                {
                  Tensor<1, 3> rotational_velocity =
                    cross_product_3d((particle.radius * particle.omega),
                                     normal);
                  contact_relative_velocity[0] =
                    particle.velocity[0] + rotational_velocity[0];
                  contact_relative_velocity[1] =
                    particle.velocity[1] + rotational_velocity[1];
                  contact_relative_velocity[2] =
                    particle.velocity[2] + rotational_velocity[2];
                }
              if (dim == 2)
                {
                  // TODO correct this after dem_2d correction
                  contact_relative_velocity = particle.velocity;
                }

              // Calculation of normal relative velocity
              double normal_relative_velocity_value =
                contact_relative_velocity * normal_vector;
              Tensor<1, dim> normal_relative_velocity =
                normal_relative_velocity_value * normal_vector;

              // Calculation of tangential relative velocity
              Tensor<1, dim> tangential_relative_velocity =
                contact_relative_velocity - normal_relative_velocity;

              Tensor<1, dim> modified_tangential_overlap =
                contact_info.tangential_overlap +
                tangential_relative_velocity * dt_dem;


              const double effective_youngs_modulus =
                (particle.youngs_modulus * wall_youngs_modulus) /
                (wall_youngs_modulus *
                   (1 - particle.poisson_ratio * particle.poisson_ratio) +
                 particle.youngs_modulus *
                   (1 - wall_poisson_ratio * wall_poisson_ratio) +
                 DBL_MIN);

              const double effective_shear_modulus =
                (particle.youngs_modulus * wall_youngs_modulus) /
                ((2 * wall_youngs_modulus * (2 - particle.poisson_ratio) *
                  (1 + particle.poisson_ratio)) +
                 (2 * particle.youngs_modulus * (2 - wall_poisson_ratio) *
                  (1 + wall_poisson_ratio)) +
                 DBL_MIN);

              const double effective_coefficient_of_restitution =
                2 * particle.restitution_coefficient *
                wall_restitution_coefficient /
                (particle.restitution_coefficient +
                 wall_restitution_coefficient + DBL_MIN);

              const double effective_coefficient_of_friction =
                2 * particle.friction_coefficient * wall_friction_coefficient /
                (particle.friction_coefficient + wall_friction_coefficient +
                 DBL_MIN);

              const double effective_coefficient_of_rolling_friction =
                2 * particle.rolling_friction_coefficient *
                wall_rolling_friction_coefficient /
                (particle.rolling_friction_coefficient +
                 wall_rolling_friction_coefficient + DBL_MIN);

              const double radius_times_overlap_sqrt =
                sqrt(particle.radius * normal_overlap);
              const double log_coeff_restitution =
                log(effective_coefficient_of_restitution);
              const double model_parameter_beta =
                log_coeff_restitution /
                sqrt((log_coeff_restitution * log_coeff_restitution) + 9.8696);
              const double model_parameter_sn =
                2 * effective_youngs_modulus * radius_times_overlap_sqrt;

              // Calculation of normal and tangential spring and dashpot
              // constants using particle and wall properties
              const double normal_spring_constant =
                1.3333 * effective_youngs_modulus * radius_times_overlap_sqrt;
              const double normal_damping_constant =
                1.8257 * model_parameter_beta *
                sqrt(model_parameter_sn * particle.mass);
              const double tangential_spring_constant =
                -8 * effective_shear_modulus * radius_times_overlap_sqrt +
                DBL_MIN;

              // Calculation of normal force using spring and dashpot normal
              // forces

              Tensor<1, dim> normal_force =
                (normal_spring_constant * normal_overlap +
                 normal_damping_constant * normal_relative_velocity_value) *
                normal_vector;

              // Calculation of tangential force
              Tensor<1, dim> tangential_force =
                tangential_spring_constant * contact_info.tangential_overlap;

              const double coulomb_threshold =
                effective_coefficient_of_friction * normal_force.norm();

              // Check for gross sliding
              if (tangential_force.norm() > coulomb_threshold)
                {
                  // Gross sliding occurs and the tangential overlap and
                  // tangnetial force are limited to Coulumb's criterion
                  tangential_force =
                    coulomb_threshold *
                    (tangential_force / tangential_force.norm());

                  contact_info.tangential_overlap =
                    tangential_force / (tangential_spring_constant + DBL_MIN);
                }

              // Calculation of torque
              // Torque caused by tangential force (tangential_torque)
              Tensor<1, 3> tangential_torque;
              if (dim == 2)
                {
                  Tensor<1, 3> tangential_force_3d;
                  tangential_force_3d[0] = tangential_force[0];
                  tangential_force_3d[1] = tangential_force[1];
                  tangential_torque =
                    cross_product_3d((particle.radius * normal),
                                     tangential_force_3d);
                }
              if (dim == 3)
                {
                  Tensor<1, 3> tangential_force_3d;
                  tangential_force_3d[0] = tangential_force[0];
                  tangential_force_3d[1] = tangential_force[1];
                  tangential_force_3d[2] = tangential_force[2];
                  tangential_torque =
                    cross_product_3d((particle.radius * normal),
                                     tangential_force_3d);
                }

              // Rolling resistance torque
              Tensor<1, 3> angular_velocity = particle.omega;
              Tensor<1, 3> pw_angular_velocity;

              double omega_value = angular_velocity.norm();
              if (omega_value != 0)
                {
                  pw_angular_velocity = angular_velocity / omega_value;
                }

              Tensor<1, 3> v_omega =
                cross_product_3d(angular_velocity, particle.radius * normal);

              // Calculation of rolling resistance torque
              Tensor<1, 3> rolling_resistance_torque =
                -effective_coefficient_of_rolling_friction * particle.radius *
                normal_force.norm() * v_omega.norm() * pw_angular_velocity;


              // Updating the force of particles in the particle handler
              contact_force[particle.particle_id] +=
                normal_force + tangential_force;

              // Updating the torque acting on particles

              contact_torque[particle.particle_id] +=
                tangential_torque + rolling_resistance_torque;
            }
          else
            {
              for (int d = 0; d < dim; ++d)
                {
                  contact_info.tangential_overlap[d] = 0;
                }
            }
        }
    }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::integrate_particles()
{
  dem_particles = particles;
  particle_residual=0;

  particle_residual=0;
  ib_dem.update_particles(particles);
  TimerOutput::Scope t(this->computing_timer, "integrate particles");
  // Integrate the velocity of the particle. If integrate motion is defined as
  // true in the parameter this function will also integrate the force to update
  // the velocity. Otherwise the velocity is kept constant

  // To integrate the forces and update the velocity, this function uses the
  // implicit Euler algorithm. To find the force at t+dt the function use the
  // fix point algorithm in parallel to the newton iteration used for the fluid
  // resolution.
  using numbers::PI;
  double dt    = this->simulation_control->get_time_steps_vector()[0];
  double time  = this->simulation_control->get_current_time();
  double alpha = this->simulation_parameters.particlesParameters->alpha;
  this->simulation_parameters.particlesParameters->f_gravity->set_time(time);

  double density    = this->simulation_parameters.particlesParameters->density;
  particle_residual = 0;
  if (this->simulation_parameters.particlesParameters->integrate_motion)
    {
      ib_dem.particles_dem( dt, this->simulation_control->is_at_start());
      for (unsigned int p = 0; p < particles.size(); ++p)
        {
          Tensor<1, dim> g;
          // Translation
          // Define the gravity force applied on the particle based on his masse
          // and the density of fluid applied on it.

          double volume=0;
          if (dim == 2)
            {
              g[0] = this->simulation_parameters.particlesParameters->f_gravity
                       ->value(particles[p].position, 0);
              g[1] = this->simulation_parameters.particlesParameters->f_gravity
                       ->value(particles[p].position, 1);
              gravity =
                g * (particles[p].mass -
                     particles[p].radius * particles[p].radius * PI * rho);
            volume=particles[p].radius * particles[p].radius * PI * rho;
            volume = particles[p].radius * particles[p].radius * PI * rho;
            }
          if (dim == 3)
            {
              volume = 4.0 / 3.0 * particles[p].radius * particles[p].radius *
                       particles[p].radius * PI;
              g[0] = this->simulation_parameters.particlesParameters->f_gravity
                       ->value(particles[p].position, 0);
              g[1] = this->simulation_parameters.particlesParameters->f_gravity
                       ->value(particles[p].position, 1);
              g[2] = this->simulation_parameters.particlesParameters->f_gravity
                       ->value(particles[p].position, 2);
              gravity =
                g * (particles[p].mass - 4.0 / 3.0 * particles[p].radius *
                                           particles[p].radius *
                                           particles[p].radius * PI * density);
              volume= 4.0 / 3.0 * particles[p].radius *
                       particles[p].radius *particles[p].radius * PI;
            }
          // Translation
          // and the density of fluide applied on it.
          // Define the gravity force applied on the particle based on his masse
          // Evaluate the velocity of the particle

          Tensor<1, dim> residual_velocity =
            particles[p].last_velocity +
            (ib_dem.dem_particles[p].impulsion+ib_dem.dem_particles[p].contact_impulsion )/ particles[p].mass - particles[p].velocity;

          Tensor<2, dim> jac_velocity;
          if (particles[p].impulsion_iter.norm() == 0)
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  jac_velocity[d][d] =
                    -1 - 0.25 * volume * rho / particles[p].mass;
                }
            }
          else
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  if ((particles[p].velocity[d] -
                       particles[p].velocity_iter[d]) != 0)
                    jac_velocity[d][d] =
                      -1 + (ib_dem.dem_particles[p].impulsion[d] -
                            ib_dem.dem_particles[p].impulsion_iter[d]) /
                             (ib_dem.dem_particles[p].velocity[d] -
                              ib_dem.dem_particles[p].velocity_iter[d]) /
                             particles[p].mass;
                  else
                    jac_velocity[d][d] =
                      -1 - 0.25 * volume * rho / particles[p].mass;
                }
            }
          particles[p].velocity_iter = particles[p].velocity;
          particles[p].velocity =
            particles[p].velocity_iter -
            residual_velocity * invert(jac_velocity) * alpha  ;


          particles[p].impulsion_iter = particles[p].impulsion;
          if(particles[p].contact_impulsion.norm()<1e-12)
            particles[p].position=particles[p].last_position+(particles[p].last_velocity+particles[p].velocity)*0.5*dt;
          else
            particles[p].position=ib_dem.dem_particles[p].position;


          if (this->simulation_parameters.non_linear_solver.verbosity !=
              Parameters::Verbosity::quiet)
            {
              this->pcout<< " contact_impulsion " <<particles[p].contact_impulsion<<std::endl;
              this->pcout << "particle " << p << " residual "
                          << residual_velocity.norm() << " particle velocity"
                          << particles[p].velocity
                          << " residual "
                          << (particles[p].position-ib_dem.dem_particles[p].position).norm()
                          << " particle position "
                          << particles[p].position << std::endl;
            }


          // For the rotation velocity : same logic as the velocity.
          auto         inv_inertia = invert(particles[p].inertia);
          Tensor<1, 3> residual_omega =
            particles[p].last_omega +
            inv_inertia * particles[p].omega_impulsion - particles[p].omega;

          Tensor<2, 3> jac_omega;
          if (particles[p].omega_impulsion.norm() == 0)
            {
              for (unsigned int d = 0; d < 3; ++d)
                {
                  jac_omega[d][d] =
                    -1 - 0.5 * 2. / 5 * volume * rho * particles[p].radius *
                           particles[p].radius * inv_inertia[d][d];
                }
            }
          else
            {
              for (unsigned int d = 0; d < 3; ++d)
                {
                  if ((particles[p].omega[d] - particles[p].omega_iter[d]) != 0)
                    jac_omega[d][d] =
                      -1 +
                      (particles[p].omega_impulsion[d] -
                       particles[p].omega_impulsion_iter[d]) /
                        (particles[p].omega[d] - particles[p].omega_iter[d]) *
                        inv_inertia[d][d];
                  else
                    jac_omega[d][d] =
                      -1 - 0.5 * 2. / 5 * volume * rho * particles[p].radius *
                             particles[p].radius * inv_inertia[d][d];
                }
            }
          particles[p].omega_iter = particles[p].omega;
          particles[p].omega      = particles[p].omega_iter -
                               residual_omega * invert(jac_omega) * alpha;


          particles[p].omega_impulsion_iter = particles[p].omega_impulsion;
          double this_particle_residual=sqrt(residual_velocity.norm()*residual_velocity.norm()+residual_omega.norm()*residual_omega.norm());
          particles[p].residual=this_particle_residual;
          if(this_particle_residual>particle_residual)
            particle_residual=this_particle_residual;
        }

    }
  else
    {
      // direct integration of the movement of the particle if it's velocity is
      // predefined.
      for (unsigned int p = 0; p < particles.size(); ++p)
        {
          particles[p].f_position->set_time(time);
          particles[p].f_velocity->set_time(time);
          particles[p].f_omega->set_time(time);
          particles[p].last_position = particles[p].position;
          particles[p].position[0] =
            particles[p].f_position->value(particles[p].position, 0);
          particles[p].position[1] =
            particles[p].f_position->value(particles[p].position, 1);
          particles[p].velocity[0] =
            particles[p].f_velocity->value(particles[p].position, 0);
          particles[p].velocity[1] =
            particles[p].f_velocity->value(particles[p].position, 1);
          particles[p].omega[0] =
            particles[p].f_omega->value(particles[p].position, 0);
          particles[p].omega[1] =
            particles[p].f_omega->value(particles[p].position, 1);
          particles[p].omega[2] =
            particles[p].f_omega->value(particles[p].position, 2);
          if (dim == 3)
            {
              particles[p].position[2] =
                particles[p].f_position->value(particles[p].position, 2);
              particles[p].velocity[2] =
                particles[p].f_velocity->value(particles[p].position, 2);
            }
        }
      particle_residual=0;
    }
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::Visualization_IB::build_patches(
  std::vector<IBParticle<dim>> particles)
{
  properties_to_write = particles[0].get_properties_name();
  /**
   * A list of field names for all data components stored in patches.
   */

  // Defining property field position
  int field_position = 0;
  // Iterating over properties
  for (auto properties_iterator = properties_to_write.begin();
       properties_iterator != properties_to_write.end();
       ++properties_iterator, ++field_position)
    {
      // Get the property field name
      const std::string field_name = properties_iterator->first;

      // Number of components of the corresponding property
      const unsigned components_number = properties_iterator->second;

      // Check to see if the property is a vector
      if (components_number == dim)
        {
          vector_datasets.push_back(std::make_tuple(
            field_position,
            field_position + components_number - 1,
            field_name,
            DataComponentInterpretation::component_is_part_of_vector));
        }
      dataset_names.push_back(field_name);
    }

  // Building the patch data
  patches.resize(particles.size());


  // Looping over particle to get the properties from the particle_handler
  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      // Particle location
      patches[p].vertices[0] = particles[p].position;
      patches[p].patch_index = p;
      patches[p].data.reinit(particles[p].get_number_properties(), 1);

      // ID and other properties

      // Calculating force for visualization
      auto particle_properties = particles[p].get_properties();

      for (unsigned int property_index = 0;
           property_index < particles[p].get_number_properties();
           ++property_index)
        patches[p].data(property_index, 0) =
          particle_properties[property_index];
    }
}

template <int dim>
const std::vector<DataOutBase::Patch<0, dim>> &
GLSSharpNavierStokesSolver<dim>::Visualization_IB::get_patches() const
{
  return patches;
}

template <int dim>
std::vector<std::string>
GLSSharpNavierStokesSolver<dim>::Visualization_IB::get_dataset_names() const
{
  return dataset_names;
}

template <int dim>
std::vector<
  std::tuple<unsigned int,
             unsigned int,
             std::string,
             DataComponentInterpretation::DataComponentInterpretation>>
GLSSharpNavierStokesSolver<dim>::Visualization_IB::get_nonscalar_data_ranges()
  const
{
  return vector_datasets;
}

template <int dim>
GLSSharpNavierStokesSolver<dim>::Visualization_IB::~Visualization_IB()
{}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::finish_time_step_particles()
{
  // Store information about the particle used for the integration and print the
  // results if requested.

  const std::string folder =
    this->simulation_parameters.simulation_control.output_folder;
  const std::string particles_solution_name =
    this->simulation_parameters.particlesParameters->ib_particles_pvd_file;
  const unsigned int iter = this->simulation_control->get_step_number();
  const double       time = this->simulation_control->get_current_time();
  const unsigned int group_files =
    this->simulation_parameters.simulation_control.group_files;

  Visualization_IB ib_particles_data;
  ib_particles_data.build_patches(particles);


  write_vtu_and_pvd<0, dim>(ib_particles_pvdhandler,
                            ib_particles_data,
                            folder,
                            particles_solution_name,
                            time,
                            iter,
                            group_files,
                            this->mpi_communicator);
  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      particles[p].last_position      = particles[p].position;
      particles[p].last_velocity      = particles[p].velocity;
      particles[p].last_forces        = particles[p].forces;
      particles[p].last_omega         = particles[p].omega;
      particles[p].last_torques       = particles[p].torques;
      particles[p].local_alpha_torque = 1;
      particles[p].local_alpha_force  = 1;

      particles[p].velocity_iter=particles[p].velocity;
      particles[p].impulsion_iter=particles[p].impulsion;
      particles[p].omega_iter=particles[p].omega;
      particles[p].omega_impulsion_iter=particles[p].omega_impulsion;



      if (this->simulation_parameters.particlesParameters.integrate_motion)
        {
          this->pcout << "particule " << p << " position "
                      << particles[p].position << std::endl;
          this->pcout << "particule " << p << " velocity "
                      << particles[p].velocity << std::endl;
        }
      table_p[p].add_value("particle_ID", p);
      if (this->simulation_parameters.simulation_control.method !=
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        table_p[p].add_value("time",
                             this->simulation_control->get_current_time());
      if (dim == 3)
        {
          table_p[p].add_value("T_x", particles[p].torques[0]);
          table_p[p].set_precision(
            "T_x",
            this->simulation_parameters.simulation_control.log_precision);
          if (this->simulation_parameters.particlesParameters->integrate_motion)
            {
              table_p[p].add_value("omega_x", particles[p].omega[0]);
              table_p[p].set_precision(
                "omega_x",
                this->simulation_parameters.simulation_control.log_precision);
            }

          table_p[p].add_value("T_y", particles[p].torques[1]);
          table_p[p].set_precision(
            "T_y",
            this->simulation_parameters.simulation_control.log_precision);
          if (this->simulation_parameters.particlesParameters->integrate_motion)
            {
              table_p[p].add_value("omega_y", particles[p].omega[1]);
              table_p[p].set_precision(
                "omega_y",
                this->simulation_parameters.simulation_control.log_precision);
            }
        }

      table_p[p].add_value("T_z", particles[p].torques[2]);
      table_p[p].set_precision(
        "T_z", this->simulation_parameters.simulation_control.log_precision);
      if (this->simulation_parameters.particlesParameters->integrate_motion)
        {
          table_p[p].add_value("omega_z", particles[p].omega[2]);
          table_p[p].set_precision(
            "omega_z",
            this->simulation_parameters.simulation_control.log_precision);
        }



      table_p[p].add_value("f_x", particles[p].forces[0]);
      if (this->simulation_parameters.particlesParameters->integrate_motion)
        {
          table_p[p].add_value("v_x", particles[p].velocity[0]);
          table_p[p].add_value("p_x", particles[p].position[0]);
        }
      table_p[p].add_value("f_y", particles[p].forces[1]);
      if (this->simulation_parameters.particlesParameters->integrate_motion)
        {
          table_p[p].add_value("v_y", particles[p].velocity[1]);
          table_p[p].add_value("p_y", particles[p].position[1]);
        }
      table_p[p].set_precision(
        "f_x", this->simulation_parameters.simulation_control.log_precision);
      table_p[p].set_precision(
        "f_y", this->simulation_parameters.simulation_control.log_precision);
      if (this->simulation_parameters.particlesParameters->integrate_motion)
        {
          table_p[p].set_precision(
            "v_x",
            this->simulation_parameters.simulation_control.log_precision);
          table_p[p].set_precision(
            "v_y",
            this->simulation_parameters.simulation_control.log_precision);
          table_p[p].set_precision(
            "p_x",
            this->simulation_parameters.simulation_control.log_precision);
          table_p[p].set_precision(
            "p_y",
            this->simulation_parameters.simulation_control.log_precision);
        }
      if (dim == 3)
        {
          table_p[p].add_value("f_z", particles[p].forces[2]);
          table_p[p].set_precision(
            "f_z",
            this->simulation_parameters.simulation_control.log_precision);
          if (this->simulation_parameters.particlesParameters->integrate_motion)
            {
              table_p[p].add_value("v_z", particles[p].velocity[2]);
              table_p[p].add_value("p_z", particles[p].position[2]);
            }
        }
    }
  if (this->this_mpi_process == 0)
    {
      if (this->simulation_parameters.forces_parameters.verbosity ==
          Parameters::Verbosity::verbose)
        {
          for (unsigned int p = 0; p < particles.size(); ++p)
            {
              std::cout << "+------------------------------------------+"
                        << std::endl;
              std::cout << "|  Force  summary particle " << p
                        << "               |" << std::endl;
              std::cout << "+------------------------------------------+"
                        << std::endl;
              table_p[p].write_text(std::cout);
            }
        }
    }
}

template <int dim>
bool
GLSSharpNavierStokesSolver<dim>::cell_cut_by_p(
  std::vector<types::global_dof_index>          &local_dof_indices,
  std::map<types::global_dof_index, Point<dim>> &support_points,
  unsigned int                                   p)
{
  // Check if a cell is cut and if it's rerun the particle by which it's cut and
  // the local DOFs index. The check is done by counting the number of DOFs that
  // is on either side of the boundary define by a particle.

  unsigned int nb_dof_inside = 0;
  for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
    {
      // Count the number of DOFs that are inside
      // of the particles. If all the DOfs are on one side
      // the cell is not cut by the boundary.
      if ((support_points[local_dof_indices[j]] - particles[p].position)
            .norm() <= particles[p].radius)
        ++nb_dof_inside;
    }
  if (nb_dof_inside != 0 && nb_dof_inside != local_dof_indices.size())
    {
      // Some of the DOFs are inside the boundary, some are outside.
      // This mean that the cell is cut so we return that information and
      // the index of the particle that cut the cell as well as the
      // container containing local DOF of the cell.
      return true;
    }

  return false;
}

template <int dim>
std::tuple<bool, unsigned int, std::vector<types::global_dof_index>>
GLSSharpNavierStokesSolver<dim>::cell_cut(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  std::vector<types::global_dof_index>                 &local_dof_indices,
  std::map<types::global_dof_index, Point<dim>>        &support_points)
{
  // Check if a cell is cut and if it's rerun the particle by which it's cut and
  // the local DOFs index. The check is done by counting the number of DOFs that
  // is on either side of the boundary define by a particle.

  cell->get_dof_indices(local_dof_indices);

  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      if (cell_cut_by_p(local_dof_indices, support_points, p))
        {
          return {true, p, local_dof_indices};
        }
    }
  return {false, 0, local_dof_indices};
}

template <int dim>
std::tuple<bool, unsigned int, std::vector<types::global_dof_index>>
GLSSharpNavierStokesSolver<dim>::cell_inside(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  std::vector<types::global_dof_index>                 &local_dof_indices,
  std::map<types::global_dof_index, Point<dim>>        &support_points)
{
  // Check if a cell is cut and if it's rerun the particle by which it's cut and
  // the local DOFs index. The check is done by counting the number of DOFs that
  // is on either side of the boundary define by a particle.

  cell->get_dof_indices(local_dof_indices);

  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      unsigned int nb_dof_inside = 0;
      for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
        {
          // Count the number of DOFs that are inside
          // of the particles. If all the DOfs are on one side
          // the cell is not cut by the boundary.
          if ((support_points[local_dof_indices[j]] - particles[p].position)
                .norm() <= particles[p].radius)
            ++nb_dof_inside;
        }
      if (nb_dof_inside == local_dof_indices.size())
        {
          // Some of the DOFs are inside the boundary, some are outside.
          // This mean that the cell is cut so we return that information and
          // the index of the particle that cut the cell as well as the
          // container containing local DOF of the cell.
          return {true, p, local_dof_indices};
        }
    }
  return {false, 0, local_dof_indices};
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::sharp_edge()
{
  // This function defines an Immersed Boundary based on the sharp edge method
  // on a hyper_sphere of dim=2 or dim=3

  TimerOutput::Scope t(this->computing_timer, "assemble_sharp");
  using numbers::PI;
  Point<dim>                                                  center_immersed;
  Point<dim>                                                  pressure_bridge;
  std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
    active_neighbors_set;
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
                                   active_neighbors_2;
  const FEValuesExtractors::Scalar pressure(dim);

  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  // Define a map to all dofs and their support points
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(*this->mapping,
                                       this->dof_handler,
                                       support_points);

  // Initalize fe value objects in order to do calculation with it later
  QGauss<dim>        q_formula(this->number_quadrature_points);
  FEValues<dim>      fe_values(*this->fe,
                               q_formula,
                               update_quadrature_points | update_JxW_values);
  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

  int    order = this->simulation_parameters.particlesParameters->order;
  double length_ratio =
    this->simulation_parameters.particlesParameters->length_ratio;



  IBStencil<dim>      stencil;
  std::vector<double> ib_coef = stencil.coefficients(order, length_ratio);

  unsigned int n_q_points = q_formula.size();

  // Define multiple local_dof_indices one for the cell iterator one for the
  // cell with the second point for the sharp edge stencil and one for
  // manipulation on the neighbour’s cell.

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_2(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_3(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_4(dofs_per_cell);
  std::set<unsigned int>               clear_line;

  // Define minimal cell length
  double dr = GridTools::minimal_cell_diameter(*this->triangulation) / sqrt(2);

  // Define cell iterator
  const auto &cell_iterator = this->dof_handler.active_cell_iterators();
  double      dt            = time_steps_vector[0];
  if (Parameters::SimulationControl::TimeSteppingMethod::steady ==
      this->simulation_parameters.simulation_control.method)
    dt = 1;

  // impose pressure reference in each of the particle
  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      Point<dim> pressure_reference_location =
        particles[p].pressure_location + particles[p].position;


      const auto &cell = LetheGridTools::find_cell_around_point_with_tree(
        this->dof_handler, pressure_reference_location);

      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);
          double sum_line = 0;
          fe_values.reinit(cell);
          std::vector<int> set_pressure_cell;
          set_pressure_cell.resize(particles.size());

          // Define the order of magnitude for the stencil.
          for (unsigned int qf = 0; qf < n_q_points; ++qf)
            sum_line += fe_values.JxW(qf);

          sum_line = 1;
          // Clear the line in the matrix
          unsigned int inside_index = local_dof_indices[dim];
          // Check on which DOF of the cell to impose the pressure. If the dof
          // is on a hanging node, it is already constrained and the pressure
          // cannot be imposed there. So we just go to the next pressure DOF of
          // the cell.

          for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
            {
              const unsigned int component_i =
                this->fe->system_to_component_index(i).first;
              if (this->zero_constraints.is_constrained(local_dof_indices[i]) ==
                    false &&
                  this->locally_owned_dofs.is_element(local_dof_indices[i]) &&
                  component_i == dim)
                {
                  inside_index = local_dof_indices[i];
                  break;
                }
            }

          this->system_matrix.clear_row(inside_index);
          // this->system_matrix.clear_row(inside_index)
          // is not reliable on edge case

          // Set the new equation for the first pressure dofs of the
          // cell. this is the new reference pressure inside a
          // particle

          this->system_matrix.set(inside_index, inside_index, sum_line);
          this->system_rhs(inside_index) =
            0 - this->evaluation_point(inside_index) * sum_line;
        }
    }

  ib_done.clear();
  // Loop on all the cell to define if the sharp edge cut them
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          double sum_line = 0;
          fe_values.reinit(cell);

          double volume = 0;
          // Define the order of magnitude for the stencil.
          for (unsigned int qf = 0; qf < n_q_points; ++qf)
            volume += fe_values.JxW(qf);

          sum_line = volume / dt;



          cell->get_dof_indices(local_dof_indices);

          // Check if the cell is cut or not by the IB and what the particle the
          // cut the cell. If the particle is cut
          bool cell_is_cut;
          // The id of the particle that cut the cell. Returns 0 if the cell is
          // not cut.
          unsigned int ib_particle_id;
          std::tie(cell_is_cut, ib_particle_id) = cut_cells_map[cell];

          if (cell_is_cut)
            {
              // If we are here, the cell is cut by the IB.
              // Loops on the dof that represents the velocity  component
              // and pressure separately
              for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
                {
                  const unsigned int component_i =
                    this->fe->system_to_component_index(i).first;
                  unsigned int global_index_overwrite = local_dof_indices[i];
                  bool         dof_is_inside =
                    (support_points[local_dof_indices[i]] -
                     particles[ib_particle_id].position)
                      .norm() < particles[ib_particle_id].radius;
                  bool use_ib_for_pressure =
                    (dof_is_inside) && (component_i == dim) &&
                    (this->simulation_parameters.particlesParameters
                       ->assemble_navier_stokes_inside == false);

                  // Check if the DOfs is owned and if it's not a hanging node.
                  if (((component_i < dim) || use_ib_for_pressure) &&
                      this->locally_owned_dofs.is_element(
                        global_index_overwrite) &&
                      ib_done[global_index_overwrite].first == false)
                    {
                      // We are working on the velocity of the cell cut
                      // loops on the dof that are for vx or vy separately
                      // loops on all the dof of the cell that represent
                      // a specific component



                      // Define which dof is going to be redefined

                      // Clear the current line of this dof
                      this->system_matrix.clear_row(global_index_overwrite);

                      // Define the points for the IB stencil, based on the
                      // order and the particle position as well as the DOF
                      // position. Depending on the order, the output variable
                      // "point" change definition. In the case of stencil
                      // orders 1 to 4 the variable point returns the position
                      // of the DOF directly. In the case of high order stencil,
                      // it returns the position of the point that is on the IB.
                      // The variable "interpolation points" return the points
                      // used to define the cell used for the stencil definition
                      // and the locations of the points use in the stencil
                      // calculation.

                      auto [point, interpolation_points] =
                        stencil.points(order,
                                       length_ratio,
                                       particles[ib_particle_id],
                                       support_points[local_dof_indices[i]]);

                      // Find the cell used for the stencil definition.
                      auto point_to_find_cell =
                        stencil.point_for_cell_detection(
                          particles[ib_particle_id],
                          support_points[local_dof_indices[i]]);
                      typename DoFHandler<dim>::active_cell_iterator cell_2;
                      bool particle_close_to_wall = false;
                      (void)particle_close_to_wall;
                      try
                        {
                          cell_2 = LetheGridTools::
                            find_cell_around_point_with_neighbors<dim>(
                              this->dof_handler,
                              vertices_to_cell,
                              cell,
                              point_to_find_cell);
                        }
                      catch (...)
                        {
                          particle_close_to_wall = true;
                          cell_2                 = cell;
                        }

                      cell_2->get_dof_indices(local_dof_indices_2);
                      ib_done[global_index_overwrite] =
                        std::make_pair(true, cell_2);

                      bool skip_stencil = false;

                      // Check if the DOF intersect the IB
                      bool dof_on_ib = false;

                      // Check if this dof is a dummy dof or directly on IB and
                      // Check if the point used to define the cell used for the
                      // definition of the stencil ("cell_2") is on a face
                      // between the cell that is cut ("cell") and the "cell_2".
                      bool point_in_cell = cell->point_inside(
                        interpolation_points[stencil.nb_points(order) - 1]);

                      if (cell_2 == cell || point_in_cell ||
                          use_ib_for_pressure || particle_close_to_wall)
                        {
                          // Give the DOF an approximated value. This help
                          // with pressure shock when the DOF passe from part of
                          // the boundary to the fluid.

                          this->system_matrix.set(global_index_overwrite,
                                                  global_index_overwrite,
                                                  sum_line);
                          skip_stencil = true;

                          // Tolerence to define a intersection of
                          // the DOF and IB
                          if (abs((support_points[local_dof_indices[i]] -
                                   particles[ib_particle_id].position)
                                    .norm() -
                                  particles[ib_particle_id].radius) <=
                              1e-12 * dr)
                            {
                              dof_on_ib = true;
                            }
                        }
                      // Define the variable used for the
                      // extrapolation of the actual solution at the
                      // boundaries in order to define the correction

                      // Define the unit cell points for the points
                      // used in the stencil.
                      std::vector<Point<dim>> unite_cell_interpolation_points(
                        ib_coef.size());
                      unite_cell_interpolation_points[0] =
                        this->mapping->transform_real_to_unit_cell(cell_2,
                                                                   point);
                      for (unsigned int j = 1; j < ib_coef.size(); ++j)
                        {
                          unite_cell_interpolation_points[j] =
                            this->mapping->transform_real_to_unit_cell(
                              cell_2, interpolation_points[j - 1]);
                        }

                      std::vector<double> local_interp_sol(ib_coef.size());

                      // Define the new matrix entry for this dof
                      if (skip_stencil == false)
                        {
                          for (unsigned int j = 0;
                               j < local_dof_indices_2.size();
                               ++j)
                            {
                              const unsigned int component_j =
                                this->fe->system_to_component_index(j).first;
                              if (component_j == component_i)
                                {
                                  //  Define the solution at each point used for
                                  //  the stencil and applied the stencil for
                                  //  the specific DOF. For stencils of order 4
                                  //  or higher, the stencil is defined through
                                  //  direct extrapolation of the cell. This can
                                  //  only be done when using a structured mesh
                                  //  as this required a mapping of a point
                                  //  outside of a cell.

                                  // Define the local matrix entries of this DOF
                                  // based on its contribution of each of the
                                  // points used in the stencil definition and
                                  // the coefficient associated with this point.
                                  // This loop defined the current solution at
                                  // the boundary using the same stencil. This
                                  // is needed to define the residual.
                                  double local_matrix_entry = 0;
                                  for (unsigned int k = 0; k < ib_coef.size();
                                       ++k)
                                    {
                                      local_matrix_entry +=
                                        this->fe->shape_value(
                                          j,
                                          unite_cell_interpolation_points[k]) *
                                        ib_coef[k];
                                      local_interp_sol[k] +=
                                        this->fe->shape_value(
                                          j,
                                          unite_cell_interpolation_points[k]) *
                                        this->evaluation_point(
                                          local_dof_indices_2[j]);
                                    }
                                  // update the matrix.
                                  this->system_matrix.set(
                                    global_index_overwrite,
                                    local_dof_indices_2[j],
                                    local_matrix_entry * sum_line);
                                }
                            }
                        }

                      // Define the RHS of the equation.

                      double rhs_add = 0;
                      // Different boundary conditions depending
                      // on the component index of the DOF and
                      // the dimension.
                      double v_ib = stencil.ib_velocity(
                        particles[ib_particle_id],
                        support_points[local_dof_indices[i]],
                        component_i);

                      for (unsigned int k = 0; k < ib_coef.size(); ++k)
                        {
                          rhs_add +=
                            -local_interp_sol[k] * ib_coef[k] * sum_line;
                        }
                      this->system_rhs(global_index_overwrite) =
                        v_ib * sum_line + rhs_add;

                      if (dof_on_ib)
                        // Dof is on the immersed boundary
                        this->system_rhs(global_index_overwrite) =
                          v_ib * sum_line -
                          this->evaluation_point(global_index_overwrite) *
                            sum_line;

                      if (skip_stencil && dof_on_ib == false)
                        // Impose the value of the dummy dof. This help
                        // with pressure variation when the IB is
                        // moving.
                        this->system_rhs(global_index_overwrite) =
                          sum_line * v_ib -
                          this->evaluation_point(global_index_overwrite) *
                            sum_line;
                    }

                  // If the DOFs is hanging put back the equations of the
                  // hanging nodes
                  if (this->zero_constraints.is_constrained(
                        local_dof_indices[i]) &&
                      this->locally_owned_dofs.is_element(
                        global_index_overwrite))
                    {
                      // Clear the line if there is something on it
                      this->system_matrix.clear_row(global_index_overwrite);
                      // Get the constraint equations
                      auto local_entries =
                        *this->zero_constraints.get_constraint_entries(
                          local_dof_indices[i]);

                      double interpolation = 0;
                      // Write the equation
                      for (unsigned int j = 0; j < local_entries.size(); ++j)
                        {
                          unsigned int col     = local_entries[j].first;
                          double       entries = local_entries[j].second;

                          interpolation +=
                            this->evaluation_point(col) * entries;
                          this->system_matrix.set(local_dof_indices[i],
                                                  col,
                                                  entries * sum_line);
                        }
                      this->system_matrix.set(local_dof_indices[i],
                                              local_dof_indices[i],
                                              sum_line);
                      // Write the RHS
                      this->system_rhs(local_dof_indices[i]) =
                        -this->evaluation_point(local_dof_indices[i]) *
                          sum_line +
                        interpolation * sum_line +
                        this->zero_constraints.get_inhomogeneity(
                          local_dof_indices[i]) *
                          sum_line;
                    }

                  if (component_i == dim && this->locally_owned_dofs.is_element(
                                              global_index_overwrite))
                    {
                      // Applied equation on dof that have no equation
                      // defined for them. those DOF become Dummy dof. This
                      // is usefull for high order cells or when a dof is
                      // only element of cells that are cut.
                      unsigned int global_index_overwrite =
                        local_dof_indices[i];
                      bool dummy_dof = true;

                      // To check if the pressure dof is a dummy. first check if
                      // the matrix entry is close to 0.
                      if (abs(this->system_matrix.el(global_index_overwrite,
                                                     global_index_overwrite)) <=
                          1e-16 * dr)
                        {
                          // If the matrix entry on the diagonal of this DOF is
                          // close to zero, check if all the cells close are
                          // cut. If it's the case, the DOF is a dummy DOF.
                          active_neighbors_set =
                            LetheGridTools::find_cells_around_cell<dim>(
                              vertices_to_cell, cell);
                          for (unsigned int m = 0;
                               m < active_neighbors_set.size();
                               m++)
                            {
                              const auto &cell_3 = active_neighbors_set[m];
                              cell_3->get_dof_indices(local_dof_indices_3);
                              for (unsigned int o = 0;
                                   o < local_dof_indices_3.size();
                                   ++o)
                                {
                                  if (global_index_overwrite ==
                                      local_dof_indices_3[o])
                                    {
                                      // cell_3 contain the same dof
                                      // check if this cell is cut if
                                      // it's not cut this dof must not
                                      // be overwritten
                                      bool cell_is_cut;
                                      std::tie(cell_is_cut, std::ignore) =
                                        cut_cells_map[cell_3];


                                      if (cell_is_cut == false)
                                        {
                                          dummy_dof = false;
                                          break;
                                        }
                                    }
                                }
                              if (dummy_dof == false)
                                break;
                            }

                          if (dummy_dof)
                            {
                              // The DOF is dummy
                              this->system_matrix.set(global_index_overwrite,
                                                      global_index_overwrite,
                                                      sum_line);
                              auto &system_rhs = this->system_rhs;
                              system_rhs(global_index_overwrite) = 0;
                            }
                        }
                    }
                }
            }
        }
    }

  this->system_rhs.compress(VectorOperation::insert);
  this->system_matrix.compress(VectorOperation::insert);

  table_residual.add_value("matrix_residual",  this->system_rhs.l2_norm());
  table_residual.set_precision("matrix_residual", 12);
  for( unsigned int p=0;p<particles.size();++p){
      table_residual.add_value("particles_residual"+ Utilities::int_to_string(p, 2),  particles[p].residual);
      table_residual.set_precision("particles_residual"+ Utilities::int_to_string(p, 2), 12);
    }
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::setup_assemblers()
{
  this->assemblers.clear();
  assemblers_inside_ib.clear();

  if (this->check_existance_of_bc(
        BoundaryConditions::BoundaryType::function_weak))
    {
      this->assemblers.push_back(
        std::make_shared<WeakDirichletBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.physical_properties,
          this->simulation_parameters.boundary_conditions));
    }
  if (this->check_existance_of_bc(BoundaryConditions::BoundaryType::pressure))
    {
      this->assemblers.push_back(
        std::make_shared<PressureBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.physical_properties,
          this->simulation_parameters.boundary_conditions));
    }
  if (this->simulation_parameters.multiphysics.VOF)
    {
      // Time-stepping schemes
      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesVOFAssemblerBDF<dim>>(
              this->simulation_control,
              this->simulation_parameters.physical_properties));
        }
      // Core assembler
      this->assemblers.push_back(
        std::make_shared<GLSNavierStokesVOFAssemblerCore<dim>>(
          this->simulation_control,
          this->simulation_parameters.physical_properties));
    }
  else
    {
      // Time-stepping schemes
      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerBDF<dim>>(
              this->simulation_control));
        }
      else if (is_sdirk(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerSDIRK<dim>>(
              this->simulation_control));
        }

      // Velocity sources term
      if (this->simulation_parameters.velocity_sources.type ==
          Parameters::VelocitySource::VelocitySourceType::srf)
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerSRF<dim>>(
              this->simulation_parameters.velocity_sources));
        }

      // Core assemblers
      if (this->simulation_parameters.physical_properties.non_newtonian_flow)
        {
          // Core assembler with Non newtonian viscosity
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerNonNewtonianCore<dim>>(
              this->simulation_control,
              this->simulation_parameters.physical_properties));
        }
      else
        {
          // Core assembler
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerCore<dim>>(
              this->simulation_control,
              this->simulation_parameters.physical_properties));
        }
    }

  assemblers_inside_ib.push_back(std::make_shared<LaplaceAssembly<dim>>(
    this->simulation_control, this->simulation_parameters.physical_properties));
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim>                         &scratch_data,
  StabilizedMethodsTensorCopyData<dim>                 &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();


  if (!cell->is_locally_owned())
    return;

  // Check if the cell is cut or not by the IB and what the particle the
  // cut the cell. If the particle is cut
  bool cell_is_cut = false;
  // The id of the particle that cut the cell. Returns 0 if the cell is
  // not cut.
  unsigned int ib_particle_id;
  std::tie(cell_is_cut, ib_particle_id) = cut_cells_map[cell];
  copy_data.cell_is_cut                 = cell_is_cut;

  if (cell_is_cut)
    return;
  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      this->solution_stages,
                      this->forcing_function,
                      this->beta);
  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_fs =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_fs);

      std::vector<TrilinosWrappers::MPI::Vector> previous_solutions;
      previous_solutions.push_back(
        *this->multiphysics->get_solution_m1(PhysicsID::VOF));

      scratch_data.reinit_VOF(phase_cell,
                              *this->multiphysics->get_solution(PhysicsID::VOF),
                              previous_solutions,
                              std::vector<TrilinosWrappers::MPI::Vector>());
    }

  copy_data.reset();

  // check if we assemble the NS eqaution inside the particle or the Laplacien
  // of the variables
  bool cell_is_inside;
  std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];
  if (cell_is_inside && this->simulation_parameters.particlesParameters
                            ->assemble_navier_stokes_inside == false)
    {
      for (auto &assembler : this->assemblers_inside_ib)
        {
          assembler->assemble_matrix(scratch_data, copy_data);
        }
    }
  else
    {
      for (auto &assembler : this->assemblers)
        {
          assembler->assemble_matrix(scratch_data, copy_data);
        }
    }


  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::copy_local_matrix_to_global_matrix(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local || copy_data.cell_is_cut)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_matrix,
                                              copy_data.local_dof_indices,
                                              this->system_matrix);
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim>                         &scratch_data,
  StabilizedMethodsTensorCopyData<dim>                 &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();

  if (!cell->is_locally_owned())
    return;

  // Check if the cell is cut or not by the IB and what the particle the
  // cut the cell. If the particle is cut
  bool cell_is_cut = false;
  // The id of the particle that cut the cell. Returns 0 if the cell is
  // not cut.
  unsigned int ib_particle_id;
  std::tie(cell_is_cut, ib_particle_id) = cut_cells_map[cell];
  copy_data.cell_is_cut                 = cell_is_cut;

  if (cell_is_cut)
    return;
  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      this->solution_stages,
                      this->forcing_function,
                      this->beta);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_fs =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_fs);

      std::vector<TrilinosWrappers::MPI::Vector> previous_solutions;
      previous_solutions.push_back(
        *this->multiphysics->get_solution_m1(PhysicsID::VOF));


      scratch_data.reinit_VOF(phase_cell,
                              *this->multiphysics->get_solution(PhysicsID::VOF),
                              previous_solutions,
                              std::vector<TrilinosWrappers::MPI::Vector>());
    }

  copy_data.reset();

  // check if we assemble the NS eqaution inside the particle or the Laplacien
  // of the variables
  bool cell_is_inside;
  std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];
  if (cell_is_inside && this->simulation_parameters.particlesParameters
                            ->assemble_navier_stokes_inside == false)
    {
      for (auto &assembler : this->assemblers_inside_ib)
        {
          assembler->assemble_rhs(scratch_data, copy_data);
        }
    }
  else
    {
      for (auto &assembler : this->assemblers)
        {
          assembler->assemble_rhs(scratch_data, copy_data);
        }
    }



  cell->get_dof_indices(copy_data.local_dof_indices);
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::copy_local_rhs_to_global_rhs(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local || copy_data.cell_is_cut)
    return;


  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_rhs,
                                              copy_data.local_dof_indices,
                                              this->system_rhs);
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::write_checkpoint()
{
  this->GLSNavierStokesSolver<dim>::write_checkpoint();

  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;

  TableHandler particles_information_table;
  std::string  filename =
    this->simulation_parameters.simulation_control.output_folder + prefix +
    ".ib_particles";
  std::ofstream output(filename.c_str());
  // Write a table with all the relevant properties of the particle in a table.
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      this->simulation_control->save(prefix);

      this->pvdhandler.save(prefix);
      for (unsigned int i_particle = 0; i_particle < particles.size();
           ++i_particle)
        {
          particles_information_table.add_value("ID", i_particle);
          particles_information_table.add_value(
            "p_x", particles[i_particle].position[0]);
          particles_information_table.set_precision("p_x", 12);
          particles_information_table.add_value(
            "p_y", particles[i_particle].position[1]);
          particles_information_table.set_precision("p_y", 12);
          if (dim == 3)
            {
              particles_information_table.add_value(
                "p_z", particles[i_particle].position[2]);
              particles_information_table.set_precision("p_z", 12);
            }

          particles_information_table.add_value(
            "v_x", particles[i_particle].velocity[0]);
          particles_information_table.set_precision("v_x", 12);
          particles_information_table.add_value(
            "v_y", particles[i_particle].velocity[1]);
          particles_information_table.set_precision("v_y", 12);

          if (dim == 3)
            {
              particles_information_table.add_value(
                "v_z", particles[i_particle].velocity[2]);
              particles_information_table.set_precision("v_z", 12);
            }

          particles_information_table.add_value(
            "f_x", particles[i_particle].forces[0]);
          particles_information_table.set_precision("f_x", 12);
          particles_information_table.add_value(
            "f_y", particles[i_particle].forces[1]);
          particles_information_table.set_precision("f_y", 12);

          if (dim == 3)
            {
              particles_information_table.add_value(
                "f_z", particles[i_particle].forces[2]);
              particles_information_table.set_precision("f_z", 12);
            }

          if (dim == 3)
            {
              particles_information_table.add_value(
                "omega_x", particles[i_particle].omega[0]);
              particles_information_table.set_precision("omega_x", 12);
              particles_information_table.add_value(
                "omega_y", particles[i_particle].omega[1]);
              particles_information_table.set_precision("omega_y", 12);
            }
          particles_information_table.add_value("omega_z",
                                                particles[i_particle].omega[2]);
          particles_information_table.set_precision("omega_z", 12);
          if (dim == 3)
            {
              particles_information_table.add_value(
                "T_x", particles[i_particle].torques[0]);
              particles_information_table.set_precision("T_x", 12);
              particles_information_table.add_value(
                "T_y", particles[i_particle].torques[1]);
              particles_information_table.set_precision("T_y", 12);
            }
          particles_information_table.add_value(
            "T_z", particles[i_particle].torques[2]);
          particles_information_table.set_precision("T_z", 12);
        }
      // Write the table in the checkpoint file.
      particles_information_table.write_text(output);
    }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::read_checkpoint()
{
  this->GLSNavierStokesSolver<dim>::read_checkpoint();

  TimerOutput::Scope t(this->computing_timer,
                       "Reset Sharp-Edge particle information");

  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;

  std::string filename =
    this->simulation_parameters.simulation_control.output_folder + prefix +
    ".ib_particles";

  // refill the table from checkpoint
  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      std::string filename_table =
        this->simulation_parameters.simulation_control.output_folder +
        this->simulation_parameters.particlesParameters->ib_force_output_file +
        "." + Utilities::int_to_string(p_i, 2) + ".dat";
      fill_table_from_file(table_p[p_i], filename_table);
    }

  // Read the data of each particle and put the relevant information in a
  // vector.
  std::map<std::string, std::vector<double>> restart_data;
  fill_vectors_from_file(restart_data, filename);


  // Implement the data  in the particles.
  if (dim == 2)
    {
      for (unsigned int p_i = 0; p_i < restart_data.size(); ++p_i)
        {
          particles[p_i].position[0] = restart_data["p_x"][p_i];
          particles[p_i].position[1] = restart_data["p_y"][p_i];
          particles[p_i].velocity[0] = restart_data["v_x"][p_i];
          particles[p_i].velocity[1] = restart_data["v_y"][p_i];
          particles[p_i].forces[0]   = restart_data["f_x"][p_i];
          particles[p_i].forces[1]   = restart_data["f_y"][p_i];
          particles[p_i].omega[2]    = restart_data["omega_z"][p_i];
          particles[p_i].torques[2]  = restart_data["T_z"][p_i];
        }
    }
  if (dim == 3)
    {
      for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
        {
          particles[p_i].position[0] = restart_data["p_x"][p_i];
          particles[p_i].position[1] = restart_data["p_y"][p_i];
          particles[p_i].position[2] = restart_data["p_z"][p_i];
          particles[p_i].velocity[0] = restart_data["v_x"][p_i];
          particles[p_i].velocity[1] = restart_data["v_y"][p_i];
          particles[p_i].velocity[2] = restart_data["v_z"][p_i];
          particles[p_i].forces[0]   = restart_data["f_x"][p_i];
          particles[p_i].forces[1]   = restart_data["f_y"][p_i];
          particles[p_i].forces[2]   = restart_data["f_z"][p_i];
          particles[p_i].omega[0]    = restart_data["omega_x"][p_i];
          particles[p_i].omega[1]    = restart_data["omega_y"][p_i];
          particles[p_i].omega[2]    = restart_data["omega_z"][p_i];
          particles[p_i].torques[0]  = restart_data["T_x"][p_i];
          particles[p_i].torques[1]  = restart_data["T_y"][p_i];
          particles[p_i].torques[2]  = restart_data["T_z"][p_i];

          particles[p_i].last_position      = particles[p_i].position;
          particles[p_i].last_velocity      = particles[p_i].velocity;
          particles[p_i].last_forces        = particles[p_i].forces;
          particles[p_i].last_omega         = particles[p_i].omega;
          particles[p_i].local_alpha_torque = 1;
          particles[p_i].local_alpha_force  = 1;
        }
    }
  // Finish the time step of the particle.
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::solve()
{
  read_mesh_and_manifolds(
    this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  define_particles();
  this->setup_dofs();
  this->box_refine_mesh();

  if (this->simulation_parameters.restart_parameters.restart == false)
    {
      // To change once refinement is split into two function
      double temp_refine =
        this->simulation_parameters.mesh_adaptation.refinement_fraction;
      double temp_coarse =
        this->simulation_parameters.mesh_adaptation.coarsening_fraction;
      this->simulation_parameters.mesh_adaptation.refinement_fraction = 0;
      this->simulation_parameters.mesh_adaptation.coarsening_fraction = 0;

      for (unsigned int i = 0;
           i <
           this->simulation_parameters.particlesParameters->initial_refinement;
           ++i)
        {
          refine_ib();
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
        }
      this->simulation_parameters.mesh_adaptation.refinement_fraction =
        temp_refine;
      this->simulation_parameters.mesh_adaptation.coarsening_fraction =
        temp_coarse;

      vertices_cell_mapping();
      generate_cut_cells_map();
    }
  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      this->forcing_function->set_time(
        this->simulation_control->get_current_time());
      if ((this->simulation_control->get_step_number() %
               this->simulation_parameters.mesh_adaptation.frequency !=
             0 ||
           this->simulation_parameters.mesh_adaptation.type ==
             Parameters::MeshAdaptation::Type::none ||
           this->simulation_control->is_at_start()) &&
          this->simulation_parameters.boundary_conditions.time_dependent)
        {
          this->update_boundary_conditions();
        }

      if (this->simulation_parameters.particlesParameters->integrate_motion ==
          false)
        integrate_particles();



      if (this->simulation_control->is_at_start())
        {
          vertices_cell_mapping();
          generate_cut_cells_map();
          ib_dem.update_particles_boundary_contact(this->particles,this->dof_handler);
          this->iterate();
        }
      else
        {
          ib_done.clear();
          refine_ib();
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          vertices_cell_mapping();
          generate_cut_cells_map();
          ib_dem.update_particles_boundary_contact(this->particles,this->dof_handler);
          // add initialization
          this->iterate();
        }

      this->postprocess_fd(false);



      if (this->simulation_parameters.particlesParameters->calculate_force_ib)
        force_on_ib();
      finish_time_step_particles();
      write_force_ib();
      this->finish_time_step();
    }

  if (this->simulation_parameters.particlesParameters->calculate_force_ib)
    this->finish_simulation();
}

// Pre-compile the 2D and 3D versopm solver to ensure that the library is
// valid before we actually compile the final solver
template class GLSSharpNavierStokesSolver<2>;
template class GLSSharpNavierStokesSolver<3>;
