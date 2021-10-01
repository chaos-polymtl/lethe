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
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/gls_sharp_navier_stokes.h>
#include <solvers/navier_stokes_vof_assemblers.h>

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
  const auto &       cell_iterator = this->dof_handler.active_cell_iterators();
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
  particles = this->simulation_parameters.particlesParameters.particles;
  table_f.resize(particles.size());
  table_t.resize(particles.size());
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
              unsigned int count_small     = 0;
              Point<dim>   center_immersed = particles[p].position;

              for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
                {
                  // Count the number of dofs that are smaller or larger than
                  // the radius of the particles if all the dof are on one side
                  // the cell is not cut by the boundary meaning we don’t have
                  // to do anything
                  if ((support_points[local_dof_indices[j]] - center_immersed)
                          .norm() <= particles[p].radius *
                                       this->simulation_parameters
                                         .particlesParameters.outside_radius &&
                      (support_points[local_dof_indices[j]] - center_immersed)
                          .norm() >= particles[p].radius *
                                       this->simulation_parameters
                                         .particlesParameters.inside_radius)
                    {
                      ++count_small;
                    }
                }
              if (count_small > 0)
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

  int    order = this->simulation_parameters.particlesParameters.order;
  double mu    = this->simulation_parameters.physical_properties.viscosity;
  double rho   = this->simulation_parameters.particlesParameters.density;

  IBStencil<dim>      stencil;
  std::vector<double> ib_coef = stencil.coefficients(order);

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
  Tensor<2, dim>              fluide_stress_at_ib;
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
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          // Particle id that cut the cell.
          unsigned int p;
          bool         cell_is_cut;
          std::tie(cell_is_cut, p) = cut_cells_map[cell];
          // If the cell is cute
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
                                  .first == false &&
                              this->locally_owned_dofs.is_element(
                                local_face_dof_indices[i]))
                            {
                              // Only need one extrapolation by dof location;
                              if (component_i == 0)
                                {
                                  // Count the number of evaluation

                                  auto [point, interpolation_points] =
                                    stencil.points(
                                      order,
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
                                      cell_2 = LetheGridTools::
                                        find_cell_around_point_with_neighbors(
                                          this->dof_handler,
                                          vertices_to_cell,
                                          cell,
                                          interpolation_points
                                            [stencil.nb_points(order) - 1]);
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

                                  fluide_stress_at_ib = 0;

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
                                      fluid_stress =
                                        mu *
                                          (velocity_gradients[k] +
                                           transpose(velocity_gradients[k])) -
                                        fluid_pressure;

                                      fluide_stress_at_ib +=
                                        fluid_stress * ib_coef[k];
                                    }
                                  // Store the stress tensor that results from
                                  // the extrapolation in the local evaluation
                                  // vector of the IB surface cell and in a map
                                  // that is used if the same extrapolation is
                                  // needed in another face.
                                  local_face_tensor[i] = fluide_stress_at_ib;
                                  force_eval_done[local_face_dof_indices[i]] =
                                    std::make_pair(true, fluide_stress_at_ib);
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
                      // Use the extrapolation of fluide stress tensor at the
                      // dof location of the IB surface cell to integrate the
                      // stress tensor on the surface of the IB
                      for (const auto &projection_cell_face :
                           local_face_dof_handler.active_cell_iterators())
                        {
                          fe_face_projection_values.reinit(
                            projection_cell_face);
                          std::vector<Point<dim>> q_points =
                            fe_face_projection_values.get_quadrature_points();
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
                                    this->fe->system_to_component_index(i)
                                      .first;
                                  if (component_i == 0 &&
                                      this->locally_owned_dofs.is_element(
                                        local_face_dof_indices[i]))
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
        Utilities::MPI::sum(particles[i].forces, this->mpi_communicator) * rho;
      particles[i].torques =
        Utilities::MPI::sum(particles[i].torques, this->mpi_communicator) * rho;
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
                .ib_force_output_file +
              "." + Utilities::int_to_string(p, 2) + ".dat";
            std::ofstream output(filename.c_str());

            table_f[p].write_text(output);
            table_t[p].write_text(output);
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

  // Calculate the error with respect to the analytical solution
  if (!firstIter &&
      this->simulation_parameters.analytical_solution->calculate_error())
    {
      // Update the time of the exact solution to the actual time
      this->exact_solution->set_time(
        this->simulation_control->get_current_time());
      const double error = this->calculate_L2_error_particles();

      if (this->simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        {
          this->error_table.add_value(
            "cells", this->triangulation->n_global_active_cells());
          this->error_table.add_value("error_velocity", error);
          this->error_table.add_value("error_pressure", 0);

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
          this->error_table.add_value("error_velocity", error);
        }
      if (this->simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error velocity : " << error << std::endl;
        }
    }
}

template <int dim>
double
GLSSharpNavierStokesSolver<dim>::calculate_L2_error_particles()
{
  TimerOutput::Scope t(this->computing_timer, "error");

  QGauss<dim>   quadrature_formula(this->number_quadrature_points + 1);
  FEValues<dim> fe_values(*this->mapping,
                          *this->fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  const unsigned int dofs_per_cell =
    this->fe->dofs_per_cell; // This gives you dofs per cell
  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points = quadrature_formula.size();

  std::vector<Vector<double>> q_exactSol(n_q_points, Vector<double>(dim + 1));

  std::vector<Tensor<1, dim>> local_velocity_values(n_q_points);
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
  double total_velocity_divergence = 0.;

  // loop over elements
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);

          bool cell_is_cut;
          // std::ignore is used because we don't care about what particle cut
          // the cell.
          std::tie(cell_is_cut, std::ignore) = cut_cells_map[cell];

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
                  total_velocity_divergence +=
                    present_velocity_divergence * fe_values.JxW(q);

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
                }
            }
        }
    }
  l2errorU = Utilities::MPI::sum(l2errorU, this->mpi_communicator);
  total_velocity_divergence =
    Utilities::MPI::sum(total_velocity_divergence, this->mpi_communicator);

  this->pcout << "div u : " << total_velocity_divergence << std::endl;

  return std::sqrt(l2errorU);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::integrate_particles()
{
  // Integrate the velocity of the particle. If integrate motion is defined as
  // true in the parameter this function will also integrate the force to update
  // the velocity. Otherwise the velocity is kept constant

  // To integrate the forces and update the velocity, this function uses the
  // implicit Euler algorithm. To find the force at t+dt the function use the
  // fix point algorithm in parallel to the newton iteration used for the fluid
  // resolution.
  using numbers::PI;
  double         dt    = this->simulation_control->get_time_steps_vector()[0];
  double         alpha = this->simulation_parameters.particlesParameters.alpha;
  Tensor<1, dim> g   = this->simulation_parameters.particlesParameters.gravity;
  double         rho = this->simulation_parameters.particlesParameters.density;
  particle_residual  = 0;
  if (this->simulation_parameters.particlesParameters.integrate_motion)
    {
      Tensor<1, dim> gravity;

      for (unsigned int p = 0; p < particles.size(); ++p)
        {
          // Translation
          // Define the gravity force applied on the particle based on his masse
          // and the density of fluide applied on it.

          if (dim == 2)
            gravity =
              g * (particles[p].mass -
                   particles[p].radius * particles[p].radius * PI * rho);
          if (dim == 3)
            {
              gravity =
                g * (particles[p].mass - 4.0 / 3.0 * particles[p].radius *
                                           particles[p].radius *
                                           particles[p].radius * PI * rho);
            }
          // Evaluate the velocity of the particle

          Tensor<1, dim> velocity_iter;
          velocity_iter =
            particles[p].last_velocity +
            (particles[p].forces + gravity) * dt / particles[p].mass;

          // This section is used to check if the fix point iteration is
          // diverging. If, between 2 iterations, the correction changes its
          // direction the relaxation parameter alpha is divided by 2. A change
          // of direction is defined as a negative cross product of the
          // correction vector and the last correction vector will the norm of
          // the new correction vector is larger than the last one.
          Tensor<1, dim> last_variation_v =
            particles[p].velocity_iter - particles[p].last_velocity;
          Tensor<1, dim> variation_v =
            velocity_iter - particles[p].last_velocity;
          double cross_product_v;
          if (dim == 2)
            cross_product_v = (last_variation_v[0] * variation_v[0] +
                               last_variation_v[1] * variation_v[1]) /
                              (last_variation_v.norm() * variation_v.norm());
          if (dim == 3)
            cross_product_v = (last_variation_v[0] * variation_v[0] +
                               last_variation_v[1] * variation_v[1] +
                               last_variation_v[2] * variation_v[2]) /
                              (last_variation_v.norm() * variation_v.norm());

          // Evaluate the velocity of the particle with the relaxation parameter
          if (last_variation_v.norm() < 1e-10)
            {
              particles[p].velocity =
                particles[p].velocity +
                alpha * (velocity_iter - particles[p].velocity);
              ;
            }
          else
            {
              if (variation_v.norm() * cross_product_v >
                  -last_variation_v.norm())
                {
                  particles[p].velocity =
                    particles[p].velocity +
                    alpha * particles[p].local_alpha_force *
                      (velocity_iter - particles[p].velocity);
                }
              else
                {
                  // If a potential divergence is observed the norm of the
                  // correction vector is adjusted to be half of the last
                  // correction vector norm and alpha are divided by 2.
                  particles[p].velocity =
                    particles[p].velocity + variation_v *
                                              last_variation_v.norm() /
                                              variation_v.norm() / 2;
                  particles[p].local_alpha_force =
                    particles[p].local_alpha_force / 2;
                }
            }
          particles[p].velocity_iter = particles[p].velocity;

          particles[p].position =
            particles[p].last_position +
            (particles[p].velocity * 0.5 + particles[p].last_velocity * 0.5) *
              dt;

          // For the rotation velocity : same logic as the velocity.
          if (dim == 2)
            {
              double i_inverse;
              i_inverse = 1.0 / particles[p].inertia[2][2];
              Tensor<1, 3> omega_iter;
              omega_iter = particles[p].last_omega +
                           (i_inverse * particles[p].torques) * dt;
              Tensor<1, 3> last_variation =
                particles[p].omega_iter - particles[p].last_omega;
              Tensor<1, 3> variation = omega_iter - particles[p].last_omega;

              double cross_product = (last_variation[0] * variation[0] +
                                      last_variation[1] * variation[1] +
                                      last_variation[2] * variation[2]) /
                                     (last_variation.norm() * variation.norm());

              if (last_variation.norm() < 1e-10)
                {
                  particles[p].omega =
                    particles[p].omega +
                    alpha * (omega_iter - particles[p].omega);
                  ;
                }
              else
                {
                  if (variation.norm() * cross_product > -last_variation.norm())
                    {
                      particles[p].omega = particles[p].omega +
                                           alpha *
                                             particles[p].local_alpha_torque *
                                             (omega_iter - particles[p].omega);
                    }
                  else
                    {
                      particles[p].omega =
                        particles[p].omega + variation * last_variation.norm() /
                                               variation.norm() / 2;
                      particles[p].local_alpha_torque =
                        particles[p].local_alpha_torque / 2;
                    }
                }


              particles[p].omega_iter    = particles[p].omega;
              particles[p].velocity_iter = particles[p].velocity;

              particles[p].angular_position +
                (particles[p].omega * 0.5 + particles[p].last_omega * 0.5) * dt;
            }


          if (dim == 3)
            {
              Tensor<2, 3, double> i_inverse;
              i_inverse = invert(particles[p].inertia);
              Tensor<1, 3, double> omega_iter;
              omega_iter[0] = particles[p].last_omega[0] +
                              (i_inverse[0][0] * particles[p].torques[0] +
                               i_inverse[0][1] * particles[p].torques[1] +
                               i_inverse[0][2] * particles[p].torques[2]) *
                                dt;
              omega_iter[1] = particles[p].last_omega[1] +
                              (i_inverse[1][0] * particles[p].torques[0] +
                               i_inverse[1][1] * particles[p].torques[1] +
                               i_inverse[1][2] * particles[p].torques[2]) *
                                dt;
              omega_iter[2] = particles[p].last_omega[2] +
                              (i_inverse[2][0] * particles[p].torques[0] +
                               i_inverse[2][1] * particles[p].torques[1] +
                               i_inverse[2][2] * particles[p].torques[2]) *
                                dt;
              Tensor<1, 3> last_variation =
                particles[p].omega_iter - particles[p].last_omega;
              Tensor<1, 3> variation = omega_iter - particles[p].last_omega;

              double cross_product = (last_variation[0] * variation[0] +
                                      last_variation[1] * variation[1] +
                                      last_variation[2] * variation[2]) /
                                     (last_variation.norm() * variation.norm());


              if (last_variation.norm() < 1e-10)
                {
                  particles[p].omega =
                    particles[p].omega +
                    alpha * (omega_iter - particles[p].omega);
                  ;
                }
              else
                {
                  if (variation.norm() * cross_product > -last_variation.norm())
                    {
                      particles[p].omega = particles[p].omega +
                                           alpha *
                                             particles[p].local_alpha_torque *
                                             (omega_iter - particles[p].omega);
                    }
                  else
                    {
                      particles[p].omega =
                        particles[p].omega + variation * last_variation.norm() /
                                               variation.norm() / 2;
                      particles[p].local_alpha_torque =
                        particles[p].local_alpha_torque / 2;
                    }
                }
              particles[p].omega_iter = particles[p].omega;

              particles[p].angular_position +
                (particles[p].omega * 0.5 + particles[p].last_omega * 0.5) * dt;
            }
        }
    }
  else
    {
      // direct integration of the movement of the particle if it's velocity is
      // predefined.
      for (unsigned int p = 0; p < particles.size(); ++p)
        {
          particles[p].last_position = particles[p].position;
          particles[p].position[0] =
            particles[p].position[0] + dt * particles[p].velocity[0];
          particles[p].position[1] =
            particles[p].position[1] + dt * particles[p].velocity[1];
          if (dim == 3)
            particles[p].position[2] =
              particles[p].position[2] + dt * particles[p].velocity[2];
        }
    }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::finish_time_step_particles()
{
  // Store information about the particle used for the integration and print the
  // results if requested.

  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      particles[p].last_position      = particles[p].position;
      particles[p].last_velocity      = particles[p].velocity;
      particles[p].last_forces        = particles[p].forces;
      particles[p].last_omega         = particles[p].omega;
      particles[p].local_alpha_torque = 1;
      particles[p].local_alpha_force  = 1;

      if (this->simulation_parameters.particlesParameters.integrate_motion)
        {
          this->pcout << "particule " << p << " position "
                      << particles[p].position << std::endl;
          this->pcout << "particule " << p << " velocity "
                      << particles[p].velocity << std::endl;
        }
      table_t[p].add_value("particle ID", p);
      if (this->simulation_parameters.simulation_control.method !=
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        table_t[p].add_value("time",
                             this->simulation_control->get_current_time());
      if (dim == 3)
        {
          table_t[p].add_value("T_x", particles[p].torques[0]);
          table_t[p].set_precision(
            "T_x",
            this->simulation_parameters.simulation_control.log_precision);
          if (this->simulation_parameters.particlesParameters.integrate_motion)
            {
              table_t[p].add_value("omega_x", particles[p].omega[0]);
              table_t[p].set_precision(
                "omega_x",
                this->simulation_parameters.simulation_control.log_precision);
            }

          table_t[p].add_value("T_y", particles[p].torques[1]);
          table_t[p].set_precision(
            "T_y",
            this->simulation_parameters.simulation_control.log_precision);
          if (this->simulation_parameters.particlesParameters.integrate_motion)
            {
              table_t[p].add_value("omega_y", particles[p].omega[1]);
              table_t[p].set_precision(
                "omega_y",
                this->simulation_parameters.simulation_control.log_precision);
            }
        }

      table_t[p].add_value("T_z", particles[p].torques[2]);
      table_t[p].set_precision(
        "T_z", this->simulation_parameters.simulation_control.log_precision);
      if (this->simulation_parameters.particlesParameters.integrate_motion)
        {
          table_t[p].add_value("omega_z", particles[p].omega[2]);
          table_t[p].set_precision(
            "omega_z",
            this->simulation_parameters.simulation_control.log_precision);
        }



      table_f[p].add_value("particle ID", p);
      if (this->simulation_parameters.simulation_control.method !=
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        table_f[p].add_value("time",
                             this->simulation_control->get_current_time());

      table_f[p].add_value("f_x", particles[p].forces[0]);
      if (this->simulation_parameters.particlesParameters.integrate_motion)
        {
          table_f[p].add_value("v_x", particles[p].velocity[0]);
          table_f[p].add_value("p_x", particles[p].position[0]);
        }
      table_f[p].add_value("f_y", particles[p].forces[1]);
      if (this->simulation_parameters.particlesParameters.integrate_motion)
        {
          table_f[p].add_value("v_y", particles[p].velocity[1]);
          table_f[p].add_value("p_y", particles[p].position[1]);
        }
      table_f[p].set_precision(
        "f_x", this->simulation_parameters.simulation_control.log_precision);
      table_f[p].set_precision(
        "f_y", this->simulation_parameters.simulation_control.log_precision);
      if (this->simulation_parameters.particlesParameters.integrate_motion)
        {
          table_f[p].set_precision(
            "v_x",
            this->simulation_parameters.simulation_control.log_precision);
          table_f[p].set_precision(
            "v_y",
            this->simulation_parameters.simulation_control.log_precision);
          table_f[p].set_precision(
            "p_x",
            this->simulation_parameters.simulation_control.log_precision);
          table_f[p].set_precision(
            "p_y",
            this->simulation_parameters.simulation_control.log_precision);
        }
      if (dim == 3)
        {
          table_f[p].add_value("f_z", particles[p].forces[2]);
          table_f[p].set_precision(
            "f_z",
            this->simulation_parameters.simulation_control.log_precision);
          if (this->simulation_parameters.particlesParameters.integrate_motion)
            {
              table_f[p].add_value("v_z", particles[p].velocity[2]);
              table_f[p].add_value("p_z", particles[p].position[2]);
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
              table_f[p].write_text(std::cout);
              table_t[p].write_text(std::cout);
            }
        }
    }
}

template <int dim>
bool
GLSSharpNavierStokesSolver<dim>::cell_cut_by_p(
  std::vector<types::global_dof_index> &         local_dof_indices,
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
  std::vector<types::global_dof_index> &                local_dof_indices,
  std::map<types::global_dof_index, Point<dim>> &       support_points)
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
  std::vector<types::global_dof_index> &                local_dof_indices,
  std::map<types::global_dof_index, Point<dim>> &       support_points)
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

  int order = this->simulation_parameters.particlesParameters.order;



  IBStencil<dim>      stencil;
  std::vector<double> ib_coef = stencil.coefficients(order);

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

          sum_line = sum_line / dt;
          // Clear the line in the matrix
          unsigned int inside_index = local_dof_indices[dim];
          // Check on which DOF of the cell to impose the pressure. If the dof
          // is on a hanging node, it is already constrained and
          // the pressure cannot be imposed there. So we just go to the next
          // pressure DOF of the cell.

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

          // Define the order of magnitude for the stencil.
          for (unsigned int qf = 0; qf < n_q_points; ++qf)
            sum_line += fe_values.JxW(qf);

          sum_line = sum_line / dt;


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


                  // Check if the DOfs is owned and if it's not a hanging node.
                  if (component_i < dim &&
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
                                       particles[ib_particle_id],
                                       support_points[local_dof_indices[i]]);

                      // Find the cell used for the stencil definition.
                      auto cell_2 =
                        LetheGridTools::find_cell_around_point_with_neighbors(
                          this->dof_handler,
                          vertices_to_cell,
                          cell,
                          interpolation_points[stencil.nb_points(order) - 1]);
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

                      if (cell_2 == cell || point_in_cell)
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

  this->system_matrix.compress(VectorOperation::insert);
  this->system_rhs.compress(VectorOperation::insert);
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::setup_assemblers()
{
  this->assemblers.clear();
  assemblers_inside_ib.clear();
  if (this->simulation_parameters.multiphysics.free_surface)
    {
      // Time-stepping schemes
      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesFreeSurfaceAssemblerBDF<dim>>(
              this->simulation_control,
              this->simulation_parameters.physical_properties));
        }
      // Core assembler
      this->assemblers.push_back(
        std::make_shared<GLSNavierStokesFreeSurfaceAssemblerCore<dim>>(
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


      // Core assembler
      this->assemblers.push_back(
        std::make_shared<GLSNavierStokesAssemblerCore<dim>>(
          this->simulation_control,
          this->simulation_parameters.physical_properties));
    }

  assemblers_inside_ib.push_back(std::make_shared<LaplaceAssembly<dim>>(
    this->simulation_control, this->simulation_parameters.physical_properties));
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim> &                        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &                copy_data)
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
  if (this->simulation_parameters.multiphysics.free_surface)
    {
      const DoFHandler<dim> *dof_handler_fs =
        this->multiphysics->get_dof_handler(PhysicsID::free_surface);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_fs);

      std::vector<TrilinosWrappers::MPI::Vector> previous_solutions;
      previous_solutions.push_back(
        *this->multiphysics->get_solution_m1(PhysicsID::free_surface));

      scratch_data.reinit_free_surface(
        phase_cell,
        *this->multiphysics->get_solution(PhysicsID::free_surface),
        previous_solutions,
        std::vector<TrilinosWrappers::MPI::Vector>());
    }

  copy_data.reset();

  // check if we assemble the NS eqaution inside the particle or the Laplacien
  // of the variables
  bool cell_is_inside;
  std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];
  if (cell_is_inside && this->simulation_parameters.particlesParameters
                            .assemble_navier_stokes_inside == false)
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
  NavierStokesScratchData<dim> &                        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &                copy_data)
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

  if (this->simulation_parameters.multiphysics.free_surface)
    {
      const DoFHandler<dim> *dof_handler_fs =
        this->multiphysics->get_dof_handler(PhysicsID::free_surface);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_fs);

      std::vector<TrilinosWrappers::MPI::Vector> previous_solutions;
      previous_solutions.push_back(
        *this->multiphysics->get_solution_m1(PhysicsID::free_surface));


      scratch_data.reinit_free_surface(
        phase_cell,
        *this->multiphysics->get_solution(PhysicsID::free_surface),
        previous_solutions,
        std::vector<TrilinosWrappers::MPI::Vector>());
    }

  copy_data.reset();

  // check if we assemble the NS eqaution inside the particle or the Laplacien
  // of the variables
  bool cell_is_inside;
  std::tie(cell_is_inside, std::ignore) = cells_inside_map[cell];
  if (cell_is_inside && this->simulation_parameters.particlesParameters
                            .assemble_navier_stokes_inside == false)
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

  // To change once refinement is split into two function
  double temp_refine =
    this->simulation_parameters.mesh_adaptation.refinement_fraction;
  double temp_coarse =
    this->simulation_parameters.mesh_adaptation.coarsening_fraction;
  this->simulation_parameters.mesh_adaptation.refinement_fraction = 0;
  this->simulation_parameters.mesh_adaptation.coarsening_fraction = 0;

  for (unsigned int i = 0;
       i < this->simulation_parameters.particlesParameters.initial_refinement;
       ++i)
    {
      refine_ib();
      NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
        refine_mesh();
    }
  this->simulation_parameters.mesh_adaptation.refinement_fraction = temp_refine;
  this->simulation_parameters.mesh_adaptation.coarsening_fraction = temp_coarse;

  vertices_cell_mapping();
  generate_cut_cells_map();

  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);

  while (this->simulation_control->integrate())
    {
      if (this->simulation_parameters.particlesParameters.integrate_motion ==
          false)
        integrate_particles();

      this->simulation_control->print_progression(this->pcout);
      if (this->simulation_control->is_at_start())
        {
          vertices_cell_mapping();
          generate_cut_cells_map();
          this->first_iteration();
        }
      else
        {
          refine_ib();
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          vertices_cell_mapping();
          generate_cut_cells_map();
          this->iterate();
        }

      this->postprocess_fd(false);

      this->finish_time_step();

      if (this->simulation_parameters.particlesParameters.calculate_force_ib)
        force_on_ib();
      finish_time_step_particles();
      write_force_ib();
    }

  if (this->simulation_parameters.particlesParameters.calculate_force_ib)


    this->finish_simulation();
}


// Pre-compile the 2D and 3D versopm solver to ensure that the library is
// valid before we actually compile the final solver
template class GLSSharpNavierStokesSolver<2>;
template class GLSSharpNavierStokesSolver<3>;
