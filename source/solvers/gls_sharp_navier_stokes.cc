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

#include "solvers/gls_sharp_navier_stokes.h"

#include "core/bdf.h"
#include "core/grids.h"
#include "core/sdirk.h"
#include "core/time_integration_utilities.h"
#include "core/utilities.h"

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
  vertices_to_cell.clear();
  vertices_to_cell.resize(this->dof_handler.n_dofs() / (dim + 1));
  const auto &cell_iterator = this->dof_handler.active_cell_iterators();

  // loop on all the cell and
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_locally_owned() | cell->is_ghost())
        {
          const unsigned int vertices_per_cell =
            GeometryInfo<dim>::vertices_per_cell;
          for (unsigned int i = 0; i < vertices_per_cell; i++)
            {
              // First obtain vertex index
              unsigned int v_index = cell->vertex_index(i);

              // Get the vector of active cell linked to that vertex
              std::vector<typename DoFHandler<dim>::active_cell_iterator>
                adjacent = vertices_to_cell[v_index];

              // Insert the cell found within that vector using a set
              // to prevent key duplication
              std::set<typename DoFHandler<dim>::active_cell_iterator>
                adjacent_2(adjacent.begin(), adjacent.end());
              adjacent_2.insert(cell);

              // Convert back the set to a vector and add it in the
              // vertices_to_cell;
              std::vector<typename DoFHandler<dim>::active_cell_iterator>
                adjacent_3(adjacent_2.begin(), adjacent_2.end());
              vertices_to_cell[v_index] = adjacent_3;
            }
        }
    }
}

// TO REFACTOR
template <int dim>
void
GLSSharpNavierStokesSolver<dim>::define_particles()
{
  particles = this->simulation_parameters.particlesParameters.particles;
  table_f.resize(particles.size());
  table_t.resize(particles.size());
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::refine_ib()
{
  Point<dim>                                    center_immersed;
  MappingQ1<dim>                                immersed_map;
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(immersed_map,
                                       this->dof_handler,
                                       support_points);

  const unsigned int                   dofs_per_cell = this->fe.dofs_per_cell;
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
                  // Count the number of dof that are smaller or larger then the
                  // radius of the particles if all the dof are on one side the
                  // cell is not cut by the boundary meaning we dont have to do
                  // anything
                  if ((support_points[local_dof_indices[j]] - center_immersed)
                          .norm() <=
                        particles[p].radius *
                          this->simulation_parameters.particlesParameters.outside_radius &&
                      (support_points[local_dof_indices[j]] - center_immersed)
                          .norm() >=
                        particles[p].radius *
                          this->simulation_parameters.particlesParameters.inside_radius)
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
  // Calculate the torque and force on a immersed boundary
  // The boundary is a circle in 2D or a sphere in 3D

  std::vector<typename DoFHandler<dim>::active_cell_iterator>
                                                              active_neighbors_set;
  std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;
  MappingQ1<dim>                                              map;
  std::pair<unsigned int, unsigned int>                       cell_vertex_map;
  const double min_cell_diameter =
    GridTools::minimal_cell_diameter(*this->triangulation);

  double dr = (min_cell_diameter) / std::sqrt(2);

  // Define stuff for later use
  using numbers::PI;
  Point<dim> center_immersed;

  if (dim == 2)
    {
      // Define general stuff useful for the evaluation of force with stencil
      QGauss<dim>   q_formula(this->fe.degree + 1);
      FEValues<dim> fe_values(this->fe, q_formula, update_quadrature_points);

      double mu = this->simulation_parameters.physical_properties.viscosity;

      MappingQ1<dim>                                immersed_map;
      std::map<types::global_dof_index, Point<dim>> support_points;
      DoFTools::map_dofs_to_support_points(immersed_map,
                                           this->dof_handler,
                                           support_points);


      std::vector<types::global_dof_index> local_dof_indices(
        this->fe.dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices_2(
        this->fe.dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices_3(
        this->fe.dofs_per_cell);
      const unsigned int nb_evaluation =
        this->simulation_parameters.particlesParameters.nb_force_eval;

      // Loop on all particles
      for (unsigned int p = 0; p < particles.size(); ++p)
        {
          // Define the center
          center_immersed = particles[p].position;

          // Initialise the output variable for this particle

          double t_torque = 0;
          // unsigned int nb_eval=0;
          double fx_v   = 0;
          double fy_v   = 0;
          double fx_p_2 = 0;
          double fy_p_2 = 0;


          // loop on all the evaluation point

          for (unsigned int i = 0; i < nb_evaluation; ++i)
            {
              // define the normal to the surface evaluated and the vector that
              // is along the surface.
              Tensor<1, dim, double> surf_normal;
              Tensor<1, dim, double> surf_vect;


              surf_normal[0] = dr * cos(i * 2 * PI / (nb_evaluation));
              surf_normal[1] = dr * sin(i * 2 * PI / (nb_evaluation));
              surf_vect[0]   = -surf_normal[1];
              surf_vect[1]   = surf_normal[0];
              double da      = 2 * PI * particles[p].radius / (nb_evaluation);
              // define the reference point for the surface evaluated.
              // the step ratio is the proportion constant used to multiplie the
              // size of the smallest cell (dr) the step done are then:
              // step_ratio * dr
              double           step_ratio = 0.5;
              const Point<dim> eval_point(
                particles[p].radius * cos(i * 2 * PI / (nb_evaluation)) +
                  center_immersed(0),
                particles[p].radius * sin(i * 2 * PI / (nb_evaluation)) +
                  center_immersed(1));

              // step in the normal direction of the surface until we find a
              // cell that is not cut by the immersed boundary of the particule
              // p.
              unsigned int     nb_step    = 0;
              bool             cell_found = false;
              const Point<dim> eval_point_2(
                eval_point[0] + surf_normal[0] * (nb_step + 1) * step_ratio,
                eval_point[1] + surf_normal[1] * (nb_step + 1) * step_ratio);

              // step in the normal direction to the surface until the point
              // used for the ib stencil is not in a cell that is cut by the
              // boundary.
              while (cell_found == false)
                {
                  // define the new point
                  Point<dim> eval_point_iter(
                    eval_point[0] + surf_normal[0] * (nb_step + 1) * step_ratio,
                    eval_point[1] +
                      surf_normal[1] * (nb_step + 1) * step_ratio);
                  /*const auto &cell_iter =
                    GridTools::find_active_cell_around_point(this->dof_handler,
                                                             eval_point_iter);*/
                  // std::cout << "before cell found " << i << std::endl;
                  const auto &cell_iter =
                    find_cell_around_point_with_tree(this->dof_handler,
                                                     eval_point_iter);
                  // std::cout << "cell found " << i<< std::endl;
                  // std::cout << "cell found v index " << cell_vertex_map.first
                  // << std::endl; std::cout << "cell found map " <<
                  // cell_vertex_map.second << std::endl;



                  if (cell_iter->is_artificial() == false)
                    {
                      // const auto
                      // &cell_iter=this->vertices_to_cell[cell_vertex_map.first][cell_vertex_map.second];
                      cell_iter->get_dof_indices(local_dof_indices);
                      // std::cout << "got dof _indices " << std::endl;
                      unsigned int count_small = 0;
                      center_immersed          = particles[p].position;

                      // check if the cell is cut
                      for (unsigned int j = 0; j < local_dof_indices.size();
                           ++j)
                        {
                          // count the number of dof that ar smaller or larger
                          // then the radius of the particles if all the dof are
                          // on one side the cell is not cut by the boundary
                          // meaning we dont have to do anything
                          if ((support_points[local_dof_indices[j]] -
                               center_immersed)
                                .norm() <= particles[p].radius)
                            {
                              ++count_small;
                            }
                        }

                      if (count_small != 0 and
                          count_small != local_dof_indices.size())
                        {
                          cell_found = false;
                        }
                      else
                        {
                          cell_found = true;
                        }

                      // step a bit further away from the boundary.
                      if (cell_found == false)
                        nb_step += 1;
                    }
                  else
                    {
                      break;
                    }
                }


              // when the point is found outside cell that are cut by the
              // boundary we define the 3 point that will be used to create
              // interpolation and  extrapolation of the solution  to evalutate
              // the force on the boundart
              const Point<dim> second_point(
                eval_point[0] + surf_normal[0] * (nb_step + 1) * step_ratio,
                eval_point[1] + surf_normal[1] * (nb_step + 1) * step_ratio);
              const Point<dim> third_point(second_point[0] +
                                             surf_normal[0] * step_ratio,
                                           second_point[1] +
                                             surf_normal[1] * step_ratio);
              const Point<dim> fourth_point(third_point[0] +
                                              surf_normal[0] * step_ratio,
                                            third_point[1] +
                                              surf_normal[1] * step_ratio);



              const auto &cell_2 =
                find_cell_around_point_with_tree(this->dof_handler,
                                                 second_point);



              // Check if the cell is locally owned before doing the evalation.
              // if (cell_vertex_map.first!=vertices_to_cell.size()+1) {
              // const auto
              // &cell_2=this->vertices_to_cell[cell_vertex_map.first][cell_vertex_map.second];
              if (cell_2->is_locally_owned())
                {
                  const auto &cell_3 =
                    find_cell_around_point_with_tree(this->dof_handler,
                                                     third_point);
                  // =
                  // this->vertices_to_cell[cell_vertex_map.first][cell_vertex_map.second];
                  const auto &cell_4 =
                    find_cell_around_point_with_tree(this->dof_handler,
                                                     fourth_point);
                  // const auto &cell_4 =
                  // this->vertices_to_cell[cell_vertex_map.first][cell_vertex_map.second];
                  cell_2->get_dof_indices(local_dof_indices);
                  cell_3->get_dof_indices(local_dof_indices_2);
                  cell_4->get_dof_indices(local_dof_indices_3);
                  // Define the tensors used for the evaluation of the velocity
                  Tensor<1, dim, double> u_1;
                  Tensor<1, dim, double> u_2;
                  Tensor<1, dim, double> u_3;

                  // Define the pressure variables
                  double P1 = 0;
                  double P2 = 0;
                  double P3 = 0;

                  // Define the velocity component of the particle at the
                  // boundary on the reference point we put the reference for
                  // the velocity at the center of the particle this simplifie
                  // the evaluation of the force if the particule is moving
                  u_1[0] = -particles[p].omega[2] * particles[p].radius *
                           sin(i * 2 * PI / (nb_evaluation));
                  u_1[1] = particles[p].omega[2] * particles[p].radius *
                           cos(i * 2 * PI / (nb_evaluation));

                  // Projection of the speed of the boundary on the plan of the
                  // surface used for evaluation
                  double U1 = (surf_vect[0] * u_1[0] + surf_vect[1] * u_1[1]) /
                              surf_vect.norm();


                  // used support function of the cell to define the
                  // interpolation of the velocity
                  Point<dim> second_point_v =
                    immersed_map.transform_real_to_unit_cell(cell_2,
                                                             second_point);
                  Point<dim> third_point_v =
                    immersed_map.transform_real_to_unit_cell(cell_3,
                                                             third_point);
                  Point<dim> fourth_point_v =
                    immersed_map.transform_real_to_unit_cell(cell_4,
                                                             fourth_point);

                  // initialise the component of the velocity
                  u_2[0] = 0;
                  u_2[1] = 0;
                  u_3[0] = 0;
                  u_3[1] = 0;
                  cell_2->get_dof_indices(local_dof_indices);
                  cell_3->get_dof_indices(local_dof_indices_2);
                  cell_4->get_dof_indices(local_dof_indices_3);



                  // define the interpolation of the cell in order to have the
                  // solution at the point previously define
                  for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
                    {
                      auto &present_solution = this->present_solution;
                      const unsigned int component_i =
                        this->fe.system_to_component_index(j).first;
                      if (component_i < dim)
                        {
                          u_2[component_i] +=
                            this->fe.shape_value(j, second_point_v) *
                            present_solution(local_dof_indices[j]);

                          u_3[component_i] +=
                            this->fe.shape_value(j, third_point_v) *
                            present_solution(local_dof_indices_2[j]);
                        }
                      if (component_i == dim)
                        {
                          P1 += this->fe.shape_value(j, second_point_v) *
                                present_solution(local_dof_indices[j]);
                          P2 += this->fe.shape_value(j, third_point_v) *
                                present_solution(local_dof_indices_2[j]);
                          P3 += this->fe.shape_value(j, fourth_point_v) *
                                present_solution(local_dof_indices_3[j]);
                        }
                    }
                  // Evaluate the solution in the reference frame of the
                  // particle
                  u_2[0] = u_2[0] - particles[p].velocity[0];
                  u_2[1] = u_2[1] - particles[p].velocity[1];
                  u_3[0] = u_3[0] - particles[p].velocity[0];
                  u_3[1] = u_3[1] - particles[p].velocity[1];

                  // Project the velocity along the surface
                  double U2 = (surf_vect[0] * u_2[0] + surf_vect[1] * u_2[1]) /
                              surf_vect.norm();
                  double U3 = (surf_vect[0] * u_3[0] + surf_vect[1] * u_3[1]) /
                              surf_vect.norm();

                  // Project the velocity along the normal
                  double U1_r =
                    (surf_normal[0] * u_1[0] + surf_normal[1] * u_1[1]) /
                    surf_normal.norm();
                  double U2_r =
                    (surf_normal[0] * u_2[0] + surf_normal[1] * u_2[1]) /
                    surf_normal.norm();
                  double U3_r =
                    (surf_normal[0] * u_3[0] + surf_normal[1] * u_3[1]) /
                    surf_normal.norm();


                  // Define the 2nd order stencil for the derivative at the
                  // boundary with variable length between the points
                  double du_dn_1 =
                    (U2 / (particles[p].radius +
                           surf_normal.norm() * (nb_step + 1) * step_ratio) -
                     U1 / particles[p].radius) /
                    ((nb_step + 1) * surf_normal.norm() * step_ratio);
                  double du_dn_2 =
                    (U3 / (particles[p].radius +
                           surf_normal.norm() * (nb_step + 2) * step_ratio) -
                     U2 / (particles[p].radius +
                           surf_normal.norm() * (nb_step + 1) * step_ratio)) /
                    (surf_normal.norm() * step_ratio);

                  double du_dr_1 =
                    (U2_r - U1_r) /
                    ((nb_step + 1) * surf_normal.norm() * step_ratio);
                  double du_dr_2 =
                    (U3_r - U2_r) / (surf_normal.norm() * step_ratio);


                  double du_dn = du_dn_1 - (du_dn_2 - du_dn_1) * (nb_step + 1) /
                                             ((nb_step + 1) + 1);
                  double du_dr = du_dr_1 - (du_dr_2 - du_dr_1) * (nb_step + 1) /
                                             ((nb_step + 1) + 1);

                  // Define the 3rd order stencil for the solution at the
                  // boundary of the pressure  with variable length between the
                  // points
                  double P_local = P1 + (nb_step + 1) * (P1 - P2) +
                                   ((nb_step + 2) * (nb_step + 1) / 2) *
                                     ((P1 - P2) - (P2 - P3));


                  // Evaluate the local force on the boundary

                  // First evaluate the viscous force
                  double local_fx_v =
                    ((-mu * du_dr * 2 * surf_normal[0] / surf_normal.norm()) +
                     (-particles[p].radius * mu * du_dn * surf_vect[0] /
                      surf_vect.norm())) *
                    da;
                  double local_fy_v =
                    ((-mu * du_dr * 2 * surf_normal[1] / surf_normal.norm()) +
                     (-particles[p].radius * mu * du_dn * surf_vect[1] /
                      surf_vect.norm())) *
                    da;

                  // Second evaluate the pressure force
                  double local_fx_p_2 =
                    P_local * da * surf_normal[0] / surf_normal.norm();
                  double local_fy_p_2 =
                    P_local * da * surf_normal[1] / surf_normal.norm();


                  // add the local contribution to the global force evaluation
                  fx_v += -local_fx_v;
                  fy_v += -local_fy_v;
                  fx_p_2 += -local_fx_p_2;
                  fy_p_2 += -local_fy_p_2;
                  t_torque += local_fx_v * sin(i * 2 * PI / (nb_evaluation)) *
                                particles[p].radius -
                              local_fy_v * cos(i * 2 * PI / (nb_evaluation)) *
                                particles[p].radius;
                  // nb_eval+=1;
                }
              //}
            }

          // Reduce the solution for each process
          double t_torque_ =
            Utilities::MPI::sum(t_torque, this->mpi_communicator);
          double fx_p_2_ = Utilities::MPI::sum(fx_p_2, this->mpi_communicator);
          double fy_p_2_ = Utilities::MPI::sum(fy_p_2, this->mpi_communicator);
          double fx_v_   = Utilities::MPI::sum(fx_v, this->mpi_communicator);
          double fy_v_   = Utilities::MPI::sum(fy_v, this->mpi_communicator);
          // unsigned int nb_eval_total   = Utilities::MPI::sum(nb_eval,
          // this->mpi_communicator);


          // Present the solution of the force on the boundary of the particle p
          if (this->this_mpi_process == 0)
            {
              if (this->simulation_parameters.forces_parameters.verbosity ==
                  Parameters::Verbosity::verbose)
                {
                  std::cout << "+------------------------------------------+"
                            << std::endl;
                  std::cout << "|  Force  summary particle : " << p
                            << "             |" << std::endl;
                  std::cout << "+------------------------------------------+"
                            << std::endl;

                  std::cout << "particle : " << p
                            << " total_torque :" << t_torque_ << std::endl;
                  std::cout << "particle : " << p
                            << " total_torque :" << t_torque_ << std::endl;
                  std::cout << "fx_P: " << fx_p_2_ << std::endl;
                  std::cout << "fy_P: " << fy_p_2_ << std::endl;
                  std::cout << "fx_v: " << fx_v_ << std::endl;
                  std::cout << "fy_v: " << fy_v_ << std::endl;


                  table_t[p].add_value("particle ID", p);
                  if (this->simulation_parameters.simulation_control.method !=
                      Parameters::SimulationControl::TimeSteppingMethod::steady)
                    table_t[p].add_value(
                      "time", this->simulation_control->get_current_time());
                  table_t[p].add_value("T_z", t_torque_);
                  table_t[p].set_precision(
                    "T_z", this->simulation_parameters.simulation_control.log_precision);



                  table_f[p].add_value("particle ID", p);
                  if (this->simulation_parameters.simulation_control.method !=
                      Parameters::SimulationControl::TimeSteppingMethod::steady)
                    table_f[p].add_value(
                      "time", this->simulation_control->get_current_time());
                  table_f[p].add_value("f_x", fx_p_2_ + fx_v_);
                  table_f[p].add_value("f_y", fy_p_2_ + fy_v_);

                  table_f[p].set_precision(
                    "f_x", this->simulation_parameters.simulation_control.log_precision);
                  table_f[p].set_precision(
                    "f_y", this->simulation_parameters.simulation_control.log_precision);
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


  // same structure as for the 2d case but used 3d variables  so there is 1 more
  // vector on the surface for the evaluation
  if (dim == 3)
    {
      QGauss<dim>   q_formula(this->fe.degree + 1);
      FEValues<dim> fe_values(this->fe, q_formula, update_quadrature_points);

      double         mu = this->simulation_parameters.physical_properties.viscosity;
      MappingQ1<dim> immersed_map;
      std::map<types::global_dof_index, Point<dim>> support_points;
      DoFTools::map_dofs_to_support_points(immersed_map,
                                           this->dof_handler,
                                           support_points);


      std::vector<types::global_dof_index> local_dof_indices(
        this->fe.dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices_2(
        this->fe.dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices_3(
        this->fe.dofs_per_cell);
      unsigned int nb_evaluation =
        this->simulation_parameters.particlesParameters.nb_force_eval;

      // the number of evaluation is round up to the closest square number so
      // there is the same number of evaluation in theta and phi direction
      nb_evaluation = ceil(pow(nb_evaluation, 0.5));

      for (unsigned int p = 0; p < particles.size(); ++p)
        {
          const double center_x = particles[p].position[0];
          const double center_y = particles[p].position[1];
          const double center_z = particles[p].position[2];

          double torque_x = 0;
          double torque_y = 0;
          double torque_z = 0;

          // unsigned int nb_eval =0;

          double fx_v = 0;
          double fy_v = 0;
          double fz_v = 0;

          double fx_p_2 = 0;
          double fy_p_2 = 0;
          double fz_p_2 = 0;

          for (unsigned int i = 0; i < nb_evaluation; ++i)
            {
              for (unsigned int j = 0; j < nb_evaluation; ++j)
                {
                  Tensor<1, dim, double> surf_normal;
                  Tensor<1, dim, double> surf_vect_1;
                  Tensor<1, dim, double> surf_vect_2;

                  double theta  = (i + 0.5) * PI / (nb_evaluation);
                  double phi    = j * 2 * PI / (nb_evaluation);
                  double dtheta = PI / (nb_evaluation);
                  double dphi   = 2 * PI / (nb_evaluation);

                  double da =
                    particles[p].radius * particles[p].radius * dphi *
                    (-cos(theta + dtheta / 2) + cos(theta - dtheta / 2));

                  const Point<dim> eval_point(
                    particles[p].radius * sin(theta) * cos(phi) + center_x,
                    particles[p].radius * sin(theta) * sin(phi) + center_y,
                    particles[p].radius * cos(theta) + center_z);

                  double step_ratio = 0.5;
                  surf_normal[0]    = dr * (sin(theta) * cos(phi));
                  surf_normal[1]    = dr * (sin(theta) * sin(phi));
                  surf_normal[2]    = dr * (cos(theta));

                  surf_vect_1[0] = dr * (sin(theta + PI / 2) * cos(phi));
                  surf_vect_1[1] = dr * (sin(theta + PI / 2) * sin(phi));
                  surf_vect_1[2] = dr * (cos(theta + PI / 2));

                  surf_vect_2 = cross_product_3d(surf_normal, surf_vect_1);
                  surf_vect_2 = dr * surf_vect_2 / surf_vect_2.norm();


                  unsigned int nb_step    = 0;
                  bool         cell_found = false;

                  while (cell_found == false)
                    {
                      Point<dim> eval_point_2(
                        eval_point[0] +
                          surf_normal[0] * (nb_step + 1) * step_ratio,
                        eval_point[1] +
                          surf_normal[1] * (nb_step + 1) * step_ratio,
                        eval_point[2] +
                          surf_normal[2] * (nb_step + 1) * step_ratio);

                      const auto &cell_iter =
                        find_cell_around_point_with_tree(this->dof_handler,
                                                         eval_point_2);

                      if (cell_iter->is_artificial() == false)
                        {
                          cell_iter->get_dof_indices(local_dof_indices);

                          unsigned int count_small = 0;
                          if (dim == 3)
                            {
                              center_immersed(0) = particles[p].position[0];
                              center_immersed(1) = particles[p].position[1];
                              center_immersed(2) = particles[p].position[2];
                            }
                          for (unsigned int j = 0; j < local_dof_indices.size();
                               ++j)
                            {
                              // Count the number of dof that ar smaller or
                              // larger then the radius of the particles if all
                              // the dof are on one side the cell is not cut by
                              // the boundary meaning we dont have to do
                              // anything

                              if ((support_points[local_dof_indices[j]] -
                                   center_immersed)
                                    .norm() <= particles[p].radius)
                                {
                                  ++count_small;
                                }
                            }

                          if (count_small != 0 and
                              count_small != local_dof_indices.size())
                            {
                              cell_found = false;
                            }
                          else
                            {
                              cell_found = true;
                            }


                          if (cell_found == false)
                            nb_step += 1;
                        }
                      else
                        {
                          break;
                        }
                    }

                  const Point<dim> second_point(
                    eval_point[0] + surf_normal[0] * (nb_step + 1) * step_ratio,
                    eval_point[1] + surf_normal[1] * (nb_step + 1) * step_ratio,
                    eval_point[2] +
                      surf_normal[2] * (nb_step + 1) * step_ratio);

                  const Point<dim> third_point(
                    second_point[0] + surf_normal[0] * step_ratio,
                    second_point[1] + surf_normal[1] * step_ratio,
                    second_point[2] + surf_normal[2] * step_ratio);

                  const Point<dim> fourth_point(
                    third_point[0] + surf_normal[0] * step_ratio,
                    third_point[1] + surf_normal[1] * step_ratio,
                    third_point[2] + surf_normal[2] * step_ratio);

                  const auto &cell_2 =
                    find_cell_around_point_with_tree(this->dof_handler,
                                                     second_point);

                  // if (cell_vertex_map.first!=vertices_to_cell.size()+1) {
                  // const auto &cell_2 =
                  // this->vertices_to_cell[cell_vertex_map.first][cell_vertex_map.second];
                  if (cell_2->is_locally_owned())
                    {
                      const auto &cell_3 =
                        find_cell_around_point_with_tree(this->dof_handler,
                                                         third_point);
                      // const auto &cell_3 =
                      // this->vertices_to_cell[cell_vertex_map.first][cell_vertex_map.second];
                      const auto &cell_4 =
                        find_cell_around_point_with_tree(this->dof_handler,
                                                         fourth_point);
                      // const auto &cell_4 =
                      // this->vertices_to_cell[cell_vertex_map.first][cell_vertex_map.second];
                      cell_2->get_dof_indices(local_dof_indices);
                      cell_3->get_dof_indices(local_dof_indices_2);
                      cell_4->get_dof_indices(local_dof_indices_3);

                      // Define the tensor used for the velocity evaluation.
                      Tensor<1, dim, double> u_1;
                      Tensor<1, dim, double> u_2;
                      Tensor<1, dim, double> u_3;
                      double                 P1 = 0;
                      double                 P2 = 0;
                      double                 P3 = 0;

                      // Define the velocity component of the particle at the
                      // boundary on the reference point in 3 d the only
                      // rotation is around the z axis
                      u_1[0] = particles[p].omega[1] * particles[p].radius *
                                 surf_normal[2] / surf_normal.norm() -
                               particles[p].omega[2] * particles[p].radius *
                                 surf_normal[1] / surf_normal.norm();
                      u_1[1] = particles[p].omega[2] * particles[p].radius *
                                 surf_normal[0] / surf_normal.norm() -
                               particles[p].omega[0] * particles[p].radius *
                                 surf_normal[2] / surf_normal.norm();
                      u_1[2] = particles[p].omega[0] * particles[p].radius *
                                 surf_normal[1] / surf_normal.norm() -
                               particles[p].omega[1] * particles[p].radius *
                                 surf_normal[0] / surf_normal.norm();

                      // Projection of the speed of the boundary on the plan of
                      // the surface used for evaluation
                      double U1_1 =
                        (surf_vect_1[0] * u_1[0] + surf_vect_1[1] * u_1[1] +
                         surf_vect_1[2] * u_1[2]) /
                        surf_vect_1.norm();

                      double U1_2 =
                        (surf_vect_2[0] * u_1[0] + surf_vect_2[1] * u_1[1] +
                         surf_vect_2[2] * u_1[2]) /
                        surf_vect_2.norm();


                      // Used support function of the cell to define the
                      // interpolation of the velocity
                      Point<dim> second_point_v =
                        immersed_map.transform_real_to_unit_cell(cell_2,
                                                                 second_point);
                      Point<dim> third_point_v =
                        immersed_map.transform_real_to_unit_cell(cell_3,
                                                                 third_point);
                      Point<dim> fourth_point_v =
                        immersed_map.transform_real_to_unit_cell(cell_4,
                                                                 fourth_point);


                      cell_3->get_dof_indices(local_dof_indices_2);
                      for (unsigned int j = 0; j < local_dof_indices.size();
                           ++j)
                        {
                          const unsigned int component_i =
                            this->fe.system_to_component_index(j).first;
                          auto &present_solution = this->present_solution;
                          if (component_i < dim)
                            {
                              u_2[component_i] +=
                                this->fe.shape_value(j, second_point_v) *
                                present_solution(local_dof_indices[j]);

                              u_3[component_i] +=
                                this->fe.shape_value(j, third_point_v) *
                                present_solution(local_dof_indices_2[j]);
                            }
                          if (component_i == dim)
                            {
                              P1 += this->fe.shape_value(j, second_point_v) *
                                    present_solution(local_dof_indices[j]);
                              P2 += this->fe.shape_value(j, third_point_v) *
                                    present_solution(local_dof_indices_2[j]);
                              P3 += this->fe.shape_value(j, fourth_point_v) *
                                    present_solution(local_dof_indices_3[j]);
                            }
                        }

                      // Evaluate the solution in the reference frame of the
                      // particle
                      u_2[0] = u_2[0] - particles[p].velocity[0];
                      u_2[1] = u_2[1] - particles[p].velocity[1];
                      u_2[2] = u_2[2] - particles[p].velocity[2];
                      u_3[0] = u_3[0] - particles[p].velocity[0];
                      u_3[1] = u_3[1] - particles[p].velocity[1];
                      u_3[2] = u_3[2] - particles[p].velocity[2];

                      double U2_1 =
                        (surf_vect_1[0] * u_2[0] + surf_vect_1[1] * u_2[1] +
                         surf_vect_1[2] * u_2[2]) /
                        surf_vect_1.norm();
                      double U3_1 =
                        (surf_vect_1[0] * u_3[0] + surf_vect_1[1] * u_3[1] +
                         surf_vect_1[2] * u_3[2]) /
                        surf_vect_1.norm();
                      double U2_2 =
                        (surf_vect_2[0] * u_2[0] + surf_vect_2[1] * u_2[1] +
                         surf_vect_2[2] * u_2[2]) /
                        surf_vect_2.norm();
                      double U3_2 =
                        (surf_vect_2[0] * u_3[0] + surf_vect_2[1] * u_3[1] +
                         surf_vect_2[2] * u_3[2]) /
                        surf_vect_2.norm();

                      double U1_r =
                        (surf_normal[0] * u_1[0] + surf_normal[1] * u_1[1] +
                         surf_normal[2] * u_1[2]) /
                        surf_normal.norm();
                      double U2_r =
                        (surf_normal[0] * u_2[0] + surf_normal[1] * u_2[1] +
                         surf_normal[2] * u_2[2]) /
                        surf_normal.norm();
                      double U3_r =
                        (surf_normal[0] * u_3[0] + surf_normal[1] * u_3[1] +
                         surf_normal[2] * u_3[2]) /
                        surf_normal.norm();

                      double du_dn_1_1 =
                        (U2_1 /
                           (particles[p].radius +
                            surf_normal.norm() * (nb_step + 1) * step_ratio) -
                         U1_1 / particles[p].radius) /
                        ((nb_step + 1) * surf_normal.norm() * step_ratio);

                      double du_dn_2_1 =
                        (U3_1 /
                           (particles[p].radius +
                            surf_normal.norm() * (nb_step + 2) * step_ratio) -
                         U2_1 /
                           (particles[p].radius +
                            surf_normal.norm() * (nb_step + 1) * step_ratio)) /
                        (surf_normal.norm() * step_ratio);

                      double du_dn_1_2 =
                        (U2_2 /
                           (particles[p].radius +
                            surf_normal.norm() * (nb_step + 1) * step_ratio) -
                         U1_2 / particles[p].radius) /
                        ((nb_step + 1) * surf_normal.norm() * step_ratio);

                      double du_dn_2_2 =
                        (U3_2 /
                           (particles[p].radius +
                            surf_normal.norm() * (nb_step + 2) * step_ratio) -
                         U2_2 /
                           (particles[p].radius +
                            surf_normal.norm() * (nb_step + 1) * step_ratio)) /
                        (surf_normal.norm() * step_ratio);


                      double du_dr_1 =
                        (U2_r - U1_r) /
                        ((nb_step + 1) * surf_normal.norm() * step_ratio);

                      double du_dr_2 =
                        (U3_r - U2_r) / (surf_normal.norm() * step_ratio);

                      double du_dn_1 = du_dn_1_1 - (du_dn_2_1 - du_dn_1_1) *
                                                     (nb_step + 1) /
                                                     ((nb_step + 1) + 1);

                      double du_dn_2 = du_dn_1_2 - (du_dn_2_2 - du_dn_1_2) *
                                                     (nb_step + 1) /
                                                     ((nb_step + 1) + 1);

                      double du_dr = du_dr_1 - (du_dr_2 - du_dr_1) *
                                                 (nb_step + 1) /
                                                 ((nb_step + 1) + 1);

                      double P_local = P1 + (nb_step + 1) * (P1 - P2) +
                                       ((nb_step + 2) * (nb_step + 1) / 2) *
                                         ((P1 - P2) - (P2 - P3));

                      double local_fx_v =
                        ((-mu * du_dr * 2 * surf_normal[0] /
                          surf_normal.norm()) +
                         (-particles[p].radius * mu * du_dn_1 * surf_vect_1[0] /
                          surf_vect_1.norm()) +
                         (-particles[p].radius * mu * du_dn_2 * surf_vect_2[0] /
                          surf_vect_2.norm())) *
                        da;

                      double local_fy_v =
                        ((-mu * du_dr * 2 * surf_normal[1] /
                          surf_normal.norm()) +
                         (-particles[p].radius * mu * du_dn_1 * surf_vect_1[1] /
                          surf_vect_1.norm()) +
                         (-particles[p].radius * mu * du_dn_2 * surf_vect_2[1] /
                          surf_vect_2.norm())) *
                        da;

                      double local_fz_v =
                        ((-mu * du_dr * 2 * surf_normal[2] /
                          surf_normal.norm()) +
                         (-particles[p].radius * mu * du_dn_1 * surf_vect_1[2] /
                          surf_vect_1.norm()) +
                         (-particles[p].radius * mu * du_dn_2 * surf_vect_2[2] /
                          surf_vect_2.norm())) *
                        da;

                      double local_fx_p_2 =
                        P_local * da * surf_normal[0] / surf_normal.norm();
                      double local_fy_p_2 =
                        P_local * da * surf_normal[1] / surf_normal.norm();
                      double local_fz_p_2 =
                        P_local * da * surf_normal[2] / surf_normal.norm();

                      fx_v += -local_fx_v;
                      fy_v += -local_fy_v;
                      fz_v += -local_fz_v;
                      fx_p_2 += -local_fx_p_2;
                      fy_p_2 += -local_fy_p_2;
                      fz_p_2 += -local_fz_p_2;

                      torque_x += local_fy_v * surf_normal[2] /
                                    surf_normal.norm() * particles[p].radius -
                                  local_fz_v * surf_normal[1] /
                                    surf_normal.norm() * particles[p].radius;
                      torque_y += local_fz_v * surf_normal[0] /
                                    surf_normal.norm() * particles[p].radius -
                                  local_fx_v * surf_normal[2] /
                                    surf_normal.norm() * particles[p].radius;
                      torque_z += local_fx_v * surf_normal[1] /
                                    surf_normal.norm() * particles[p].radius -
                                  local_fy_v * surf_normal[0] /
                                    surf_normal.norm() * particles[p].radius;
                      // nb_eval+=1;
                    }
                  //}
                }
            }
          double t_torque_x =
            Utilities::MPI::sum(torque_x, this->mpi_communicator);
          double t_torque_y =
            Utilities::MPI::sum(torque_y, this->mpi_communicator);
          double t_torque_z =
            Utilities::MPI::sum(torque_z, this->mpi_communicator);
          double fx_p_2_ = Utilities::MPI::sum(fx_p_2, this->mpi_communicator);
          double fy_p_2_ = Utilities::MPI::sum(fy_p_2, this->mpi_communicator);
          double fz_p_2_ = Utilities::MPI::sum(fz_p_2, this->mpi_communicator);
          double fx_v_   = Utilities::MPI::sum(fx_v, this->mpi_communicator);
          double fy_v_   = Utilities::MPI::sum(fy_v, this->mpi_communicator);
          double fz_v_   = Utilities::MPI::sum(fz_v, this->mpi_communicator);
          // unsigned int nb_eval_total   = Utilities::MPI::sum(nb_eval,
          // this->mpi_communicator);

          if (this->this_mpi_process == 0)
            {
              if (this->simulation_parameters.forces_parameters.verbosity ==
                  Parameters::Verbosity::verbose)
                {
                  std::cout << "+------------------------------------------+"
                            << std::endl;
                  std::cout << "|  Force  summary particle : " << p
                            << "             |" << std::endl;
                  std::cout << "+------------------------------------------+"
                            << std::endl;

                  std::cout << "particle : " << p
                            << " total_torque_x :" << t_torque_x << std::endl;
                  std::cout << "particle : " << p
                            << " total_torque_y :" << t_torque_y << std::endl;
                  std::cout << "particle : " << p
                            << " total_torque_z :" << t_torque_z << std::endl;
                  std::cout << "fx_P: " << fx_p_2_ << std::endl;
                  std::cout << "fy_P: " << fy_p_2_ << std::endl;
                  std::cout << "fz_P: " << fz_p_2_ << std::endl;
                  std::cout << "fx_v: " << fx_v_ << std::endl;
                  std::cout << "fy_v: " << fy_v_ << std::endl;
                  std::cout << "fz_v: " << fz_v_ << std::endl;
                  // std::cout << "nb eval" << nb_eval_total << std::endl;
                  table_t[p].add_value("particle ID", p);
                  if (this->simulation_parameters.simulation_control.method !=
                      Parameters::SimulationControl::TimeSteppingMethod::steady)
                    table_t[p].add_value(
                      "time", this->simulation_control->get_current_time());
                  table_t[p].add_value("T_x", t_torque_x);
                  table_t[p].set_precision(
                    "T_x", this->simulation_parameters.simulation_control.log_precision);
                  table_t[p].add_value("T_y", t_torque_x);
                  table_t[p].set_precision(
                    "T_y", this->simulation_parameters.simulation_control.log_precision);
                  table_t[p].add_value("T_z", t_torque_x);
                  table_t[p].set_precision(
                    "T_z", this->simulation_parameters.simulation_control.log_precision);



                  table_f[p].add_value("particle ID", p);
                  if (this->simulation_parameters.simulation_control.method !=
                      Parameters::SimulationControl::TimeSteppingMethod::steady)
                    table_f[p].add_value(
                      "time", this->simulation_control->get_current_time());

                  table_f[p].add_value("f_x", fx_p_2_ + fx_v_);
                  table_f[p].add_value("f_y", fy_p_2_ + fy_v_);

                  table_f[p].set_precision(
                    "f_x", this->simulation_parameters.simulation_control.log_precision);
                  table_f[p].set_precision(
                    "f_y", this->simulation_parameters.simulation_control.log_precision);

                  table_f[p].add_value("f_z", fz_p_2_ + fz_v_);
                  table_f[p].set_precision(
                    "f_z", this->simulation_parameters.simulation_control.log_precision);
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
              this->simulation_parameters.particlesParameters.ib_force_output_file + "." +
              Utilities::int_to_string(p, 2) + ".dat";
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
  auto &present_solution = this->present_solution;
  if (this->simulation_control->is_output_iteration())
    this->write_output_results(present_solution);

  // Calculate error with respect to analytical solution
  if (!firstIter && this->simulation_parameters.analytical_solution->calculate_error())
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

  QGauss<dim>         quadrature_formula(this->number_quadrature_points + 1);
  const MappingQ<dim> mapping(this->velocity_fem_degree,
                              this->simulation_parameters.fem_parameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  const unsigned int dofs_per_cell =
    this->fe.dofs_per_cell; // This gives you dofs per cell
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
  MappingQ1<dim>                                immersed_map;
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(immersed_map,
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
          bool check_error = true;
          for (unsigned int p = 0; p < particles.size(); ++p)
            {
              unsigned int count_small = 0;
              center_immersed(0)       = particles[p].position[0];
              center_immersed(1)       = particles[p].position[1];
              if (dim == 3)
                {
                  center_immersed(2) = particles[p].position[2];
                }
              for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
                {
                  // Count the number of dof that ar smaller or larger then the
                  // radius of the particles if all the dof are on one side the
                  // cell is not cut by the boundary meaning we dont have to do
                  // anything
                  if ((support_points[local_dof_indices[j]] - center_immersed)
                        .norm() <= particles[p].radius)
                    {
                      ++count_small;
                    }
                }
              if (count_small != 0 and count_small != local_dof_indices.size())
                {
                  check_error = false;
                }
            }

          if (check_error)
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
GLSSharpNavierStokesSolver<dim>::sharp_edge()
{
  // This function defines a Immersed Boundary based on the sharp edge method on
  // a hyper_shere of dim=2 or dim=3

  TimerOutput::Scope t(this->computing_timer, "assemble_sharp");
  using numbers::PI;
  Point<dim>                                                  center_immersed;
  Point<dim>                                                  pressure_bridge;
  std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
    active_neighbors_set;
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
    active_neighbors_2;


  // Define a map to all dof and it's support point
  MappingQ1<dim>                                immersed_map;
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(immersed_map,
                                       this->dof_handler,
                                       support_points);

  // Initalize fe value object in order to do calculation with it later
  QGauss<dim>        q_formula(this->number_quadrature_points);
  FEValues<dim>      fe_values(this->fe,
                          q_formula,
                          update_quadrature_points | update_JxW_values);
  const unsigned int dofs_per_cell   = this->fe.dofs_per_cell;
  unsigned int       vertex_per_cell = GeometryInfo<dim>::vertices_per_cell;


  unsigned int n_q_points = q_formula.size();

  // Define multiple local_dof_indices one for the cell iterator one for the
  // cell with the second point for the sharp edge stancil and one for
  // manipulation on the neighbors cell.
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_2(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_3(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_4(dofs_per_cell);
  std::set<unsigned int>               clear_line;

  // Define minimal cell length
  double dr = GridTools::minimal_cell_diameter(*this->triangulation) / sqrt(2);

  // Define cell iterator
  const auto &cell_iterator = this->dof_handler.active_cell_iterators();

  // Loop on all the cell to define if the sharp edge cut them
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_locally_owned())
        {
          double sum_line = 0;
          fe_values.reinit(cell);
          cell->get_dof_indices(local_dof_indices);
          std::vector<int> set_pressure_cell;
          set_pressure_cell.resize(particles.size());

          // Define the order of magnitude for the stencil.
          for (unsigned int qf = 0; qf < n_q_points; ++qf)
            sum_line += fe_values.JxW(qf);

          // Loop over all particle  to see if one of them is cutting this cell
          for (unsigned int p = 0; p < particles.size(); ++p)
            {
              unsigned int count_small = 0;
              center_immersed          = particles[p].position;
              pressure_bridge =
                particles[p].position - particles[p].pressure_location;

              for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
                {
                  // Count the number of dof that are smaller or larger then the
                  // radius of the particles if all the dof are on one side the
                  // cell is not cut by the boundary meaning we dont have to do
                  // anything
                  if ((support_points[local_dof_indices[j]] - center_immersed)
                        .norm() <= particles[p].radius)
                    ++count_small;
                }

              // Impose the pressure inside the particle if the inside of the
              // particle is solved

              bool cell_found = false;
              try
                {
                  // Define the cell and check if the point is inside of the
                  // cell
                  const Point<dim, double> p_cell =
                    immersed_map.transform_real_to_unit_cell(cell,
                                                             pressure_bridge);
                  const double dist_2 =
                    GeometryInfo<dim>::distance_to_unit_cell(p_cell);

                  if (dist_2 == 0)
                    {
                      // If the point is in this cell then the dist is equal
                      // to 0 and we have found our cell
                      cell_found = true;
                    }
                }
              // may cause error if the point is not in cell
              catch (
                const typename MappingQGeneric<dim>::ExcTransformationFailed &)
                {}

              if (cell_found)
                {
                  // clear the line in the matrix
                  unsigned int inside_index = local_dof_indices[dim];
                  for (unsigned int vi = 0; vi < vertex_per_cell; ++vi)
                    {
                      unsigned int v_index = cell->vertex_index(vi);
                      active_neighbors_set = this->vertices_to_cell[v_index];
                      for (unsigned int m = 0; m < active_neighbors_set.size();
                           m++)
                        {
                          const auto &cell_3 = active_neighbors_set[m];
                          cell_3->get_dof_indices(local_dof_indices_3);
                          for (unsigned int o = 0;
                               o < local_dof_indices_3.size();
                               ++o)
                            {
                              if (std::find(local_dof_indices_3.begin(),
                                            local_dof_indices_3.end(),
                                            inside_index) !=
                                  local_dof_indices_3.end())
                                {
                                  for (unsigned int o = 0;
                                       o < local_dof_indices_3.size();
                                       ++o)
                                    {
                                      this->system_matrix.set(
                                        inside_index,
                                        local_dof_indices_3[o],
                                        0);
                                    }
                                }
                            }
                        }
                    }

                  // this->system_matrix.clear_row(inside_index);
                  // set new equation for the first pressure dof of the
                  // cell. this is the new reference pressure inside a
                  // particle
                  this->system_matrix.set(inside_index,
                                          local_dof_indices[dim],
                                          sum_line);
                  auto &system_rhs = this->system_rhs;
                  system_rhs(inside_index) =
                    0 - this->local_evaluation_point(inside_index) * sum_line;
                }



              // If the cell is cut by the IB the count is not 0 or the
              // number of total dof in a cell

              if (count_small != 0 && count_small != local_dof_indices.size())
                {
                  // If we are here the cell is cut by the immersed boundary
                  // loops on the dof that reprensant the velocity  component
                  // and pressure separately
                  for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
                    {
                      const unsigned int component_i =
                        this->fe.system_to_component_index(i).first;

                      if (component_i < dim)
                        {
                          // We are working on the velocity of th
                          // loops on the dof that are for vx or vy separately
                          // loops on all the dof of the the cell that represent
                          // a specific component
                          // define which dof is going to be redefine
                          unsigned int global_index_overwrite =
                            local_dof_indices[i];

                          // Define the distance vector between the
                          // immersed boundary and the dof support point
                          // for each dof
                          Tensor<1, dim, double> vect_dist =
                            (support_points[local_dof_indices[i]] -
                             center_immersed -
                             particles[p].radius *
                               (support_points[local_dof_indices[i]] -
                                center_immersed) /
                               (support_points[local_dof_indices[i]] -
                                center_immersed)
                                 .norm());

                          // Define the length ratio that represent the
                          // zone used for the stencil. The length is
                          // defined as the length between the dof and
                          // the IB
                          unsigned int length_ratio    = 8;
                          double       length_fraction = 1. / length_ratio;
                          double       tp_ratio        = 1. / 2.;
                          double       fp_ratio        = 3. / 4.;

                          if (this->simulation_parameters.particlesParameters.order == 3)
                            {
                              tp_ratio = 1. / 3.;
                              fp_ratio = 2. / 3.;
                            }



                          // Define the other points for the stencil
                          // (IB point, original dof and the other
                          // points) this goes up to a 5 point stencil.
                          Point<dim, double> first_point(
                            support_points[local_dof_indices[i]] - vect_dist);

                          Point<dim, double> second_point(
                            support_points[local_dof_indices[i]] +
                            vect_dist * length_fraction);

                          Point<dim, double> third_point(
                            support_points[local_dof_indices[i]] +
                            vect_dist * length_fraction * tp_ratio);

                          Point<dim, double> fourth_point(
                            support_points[local_dof_indices[i]] +
                            vect_dist * length_fraction * fp_ratio);

                          Point<dim, double> fifth_point(
                            support_points[local_dof_indices[i]] +
                            vect_dist * length_fraction * 1 / 4);

                          double dof_2;
                          double sp_2;

                          double dof_3;
                          double sp_3;
                          double tp_3;

                          double dof_4;
                          double sp_4;
                          double tp_4;
                          double fp_4;

                          double dof_5;
                          double fp2_5;
                          double tp_5;
                          double fp1_5;
                          double sp_5;

                          // Define the stencil coefficient in function
                          // of the length ratio. This will be
                          // automatically generated in future version
                          if (length_ratio == 4)
                            {
                              dof_2 = 5;
                              sp_2  = -4;

                              dof_3 = 45;
                              sp_3  = 36;
                              tp_3  = -80;

                              dof_4 = 455;
                              sp_4  = -364;
                              tp_4  = -1260;
                              fp_4  = 1170;

                              dof_5 = 4845;
                              fp2_5 = -18240;
                              tp_5  = 25840;
                              fp1_5 = -16320;
                              sp_5  = 3876;
                            }
                          else if (length_ratio == 2)
                            {
                              dof_2 = 3;
                              sp_2  = -2;

                              dof_3 = 15;
                              sp_3  = 10;
                              tp_3  = -24;

                              dof_4 = 84;
                              sp_4  = -56;
                              tp_4  = -216;
                              fp_4  = 189;

                              dof_5 = 495;
                              fp2_5 = -1760;
                              tp_5  = 2376;
                              fp1_5 = -1440;
                              sp_5  = 330;
                            }
                          else if (length_ratio == 8)
                            {
                              dof_2 = 9;
                              sp_2  = -8;

                              dof_3 = 153;
                              sp_3  = 136;
                              tp_3  = -288;

                              dof_4 = 2925;
                              sp_4  = -2600;
                              tp_4  = -8424;
                              fp_4  = 8100;

                              dof_5 = 58905;
                              fp2_5 = -228480;
                              tp_5  = 332640;
                              fp1_5 = -215424;
                              sp_5  = 52360;
                            }


                          // Define the vertex associated with the dof
                          unsigned int cell_found = 0;
                          bool         break_bool = false;
                          for (unsigned int vi = 0; vi < vertex_per_cell; ++vi)
                            {
                              unsigned int v_index = cell->vertex_index(vi);

                              // Get a cell iterator for all the cell
                              // neighbors of that vertex
                              active_neighbors_set =
                                this->vertices_to_cell[v_index];
                              unsigned int n_active_cells =
                                active_neighbors_set.size();

                              // Loops on those cell to find in which of
                              // them the new point for or sharp edge
                              // stencil is
                              for (unsigned int cell_index = 0;
                                   cell_index < n_active_cells;
                                   ++cell_index)
                                {
                                  try
                                    {
                                      // Define the cell and check if
                                      // the point is inside of the cell
                                      const Point<dim, double> p_cell =
                                        immersed_map
                                          .transform_real_to_unit_cell(
                                            active_neighbors_set[cell_index],
                                            second_point);
                                      const double dist_2 = GeometryInfo<
                                        dim>::distance_to_unit_cell(p_cell);

                                      // define the cell and check if
                                      // the point is inside of the cell
                                      if (dist_2 == 0)
                                        {
                                          // If the point is in this
                                          // cell then the dist is equal
                                          // to 0 and we have found our
                                          // cell
                                          cell_found = cell_index;
                                          break_bool = true;
                                          active_neighbors =
                                            active_neighbors_set;
                                          break;
                                        }
                                    }
                                  // may cause error if the point is not
                                  // in cell
                                  catch (const typename MappingQGeneric<
                                         dim>::ExcTransformationFailed &)
                                    {}
                                }
                            }

                          auto &cell_2       = active_neighbors[cell_found];
                          bool  skip_stencil = false;

                          if (break_bool == false)
                            {
                              std::cout << "cell not found around point "
                                        << std::endl;
                              std::cout << "cell index " << cell_found
                                        << std::endl;
                              std::cout << "second point  " << second_point
                                        << std::endl;
                              cell_2 = GridTools::find_active_cell_around_point(
                                this->dof_handler, second_point);
                              cell_2->get_dof_indices(local_dof_indices_2);
                              std::cout
                                << "dof point  "
                                << support_points[global_index_overwrite]
                                << std::endl;
                            }



                          // We have or next cell needed to complete the
                          // stencil

                          // Define the unit cell points for the points
                          // used in the stencil for extrapolation.
                          Point<dim> first_point_v =
                            immersed_map.transform_real_to_unit_cell(
                              cell_2, first_point);
                          Point<dim> second_point_v =
                            immersed_map.transform_real_to_unit_cell(
                              cell_2, second_point);
                          Point<dim> third_point_v =
                            immersed_map.transform_real_to_unit_cell(
                              cell_2, third_point);
                          Point<dim> fourth_point_v =
                            immersed_map.transform_real_to_unit_cell(
                              cell_2, fourth_point);
                          Point<dim> fifth_point_v =
                            immersed_map.transform_real_to_unit_cell(
                              cell_2, fifth_point);

                          cell_2->get_dof_indices(local_dof_indices_2);

                          // Clear the current line of this dof  by
                          // looping on the neighbors cell of this dof
                          // and clear all the associated dof
                          for (unsigned int vi = 0; vi < vertex_per_cell; ++vi)
                            {
                              unsigned int v_index = cell->vertex_index(vi);
                              active_neighbors_set =
                                this->vertices_to_cell[v_index];
                              for (unsigned int m = 0;
                                   m < active_neighbors_set.size();
                                   m++)
                                {
                                  const auto &cell_3 = active_neighbors_set[m];
                                  cell_3->get_dof_indices(local_dof_indices_3);
                                  if (std::find(local_dof_indices_3.begin(),
                                                local_dof_indices_3.end(),
                                                global_index_overwrite) !=
                                      local_dof_indices_3.end())
                                    {
                                      for (unsigned int o = 0;
                                           o < local_dof_indices_3.size();
                                           ++o)
                                        {
                                          this->system_matrix.set(
                                            global_index_overwrite,
                                            local_dof_indices_3[o],
                                            0);
                                        }
                                    }
                                }
                            }

                          // Check if the DOF intersect the IB
                          bool do_rhs = false;
                          if (cell_2 == cell)
                            {
                              skip_stencil = true;
                              this->system_matrix.set(global_index_overwrite,
                                                      global_index_overwrite,
                                                      sum_line);
                              auto &system_rhs = this->system_rhs;
                              system_rhs(global_index_overwrite) = 0;
                              // Tolerence to define a intersection of
                              // the DOF and IB
                              if (vect_dist.norm() <= 1e-12 * dr)
                                {
                                  do_rhs = true;
                                }
                              else
                                {
                                  system_rhs(global_index_overwrite) = 0;
                                }
                            }


                          // Define the variable used for the
                          // extrapolation of the actual solution at the
                          // boundary in order to define the correction
                          double local_interp_sol   = 0;
                          double local_interp_sol_2 = 0;
                          double local_interp_sol_3 = 0;
                          double local_interp_sol_4 = 0;

                          // Define the new matrix entry for this dof
                          if (skip_stencil == false)
                            {
                              // First the dof itself
                              for (unsigned int j = 0;
                                   j < local_dof_indices_2.size();
                                   ++j)
                                {
                                  // First the dof itself
                                  const unsigned int component_j =
                                    this->fe.system_to_component_index(j).first;
                                  if (component_j == component_i)
                                    {
                                      if (global_index_overwrite ==
                                          local_dof_indices_2[j])
                                        {
                                          // Define the solution at each
                                          // point used for the stencil and
                                          // applied the stencil for the
                                          // specfic dof. for stencil with
                                          // order of convergence higher
                                          // then 5 the stencil is define
                                          // trough direct extrapolation of
                                          // the cell
                                          auto &evaluation_point =
                                            this->evaluation_point;

                                          if (this->simulation_parameters.particlesParameters
                                                .order == 1)
                                            {
                                              this->system_matrix.set(
                                                global_index_overwrite,
                                                local_dof_indices_2[j],
                                                sp_2 *
                                                    this->fe.shape_value(
                                                      j, second_point_v) *
                                                    sum_line +
                                                  dof_2 * sum_line);
                                              local_interp_sol +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, second_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                            }

                                          if (this->simulation_parameters.particlesParameters
                                                .order == 2)
                                            {
                                              this->system_matrix.set(
                                                global_index_overwrite,
                                                local_dof_indices_2[j],
                                                sp_3 *
                                                    this->fe.shape_value(
                                                      j, second_point_v) *
                                                    sum_line +
                                                  dof_3 * sum_line +
                                                  tp_3 *
                                                    this->fe.shape_value(
                                                      j, third_point_v) *
                                                    sum_line);

                                              local_interp_sol +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, second_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_2 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, third_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                            }
                                          if (this->simulation_parameters.particlesParameters
                                                .order == 3)
                                            {
                                              this->system_matrix.set(
                                                global_index_overwrite,
                                                local_dof_indices_2[j],
                                                sp_4 *
                                                    this->fe.shape_value(
                                                      j, second_point_v) *
                                                    sum_line +
                                                  dof_4 * sum_line +
                                                  tp_4 *
                                                    this->fe.shape_value(
                                                      j, third_point_v) *
                                                    sum_line +
                                                  fp_4 *
                                                    this->fe.shape_value(
                                                      j, fourth_point_v) *
                                                    sum_line);

                                              local_interp_sol +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, second_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_2 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, third_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_3 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, fourth_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                            }
                                          if (this->simulation_parameters.particlesParameters
                                                .order > 4)
                                            {
                                              this->system_matrix.set(
                                                global_index_overwrite,
                                                local_dof_indices_2[j],
                                                this->fe.shape_value(
                                                  j, first_point_v) *
                                                  sum_line);
                                              local_interp_sol +=
                                                this->fe.shape_value(
                                                  j, first_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                            }

                                          if (this->simulation_parameters.particlesParameters
                                                .order == 4)
                                            {
                                              this->system_matrix.set(
                                                global_index_overwrite,
                                                local_dof_indices_2[j],
                                                dof_5 * sum_line +
                                                  sp_5 *
                                                    this->fe.shape_value(
                                                      j, second_point_v) *
                                                    sum_line +
                                                  tp_5 *
                                                    this->fe.shape_value(
                                                      j, third_point_v) *
                                                    sum_line +
                                                  fp1_5 *
                                                    this->fe.shape_value(
                                                      j, fourth_point_v) *
                                                    sum_line +
                                                  fp2_5 *
                                                    this->fe.shape_value(
                                                      j, fifth_point_v) *
                                                    sum_line);


                                              local_interp_sol +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, second_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_2 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, third_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_3 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, fourth_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_4 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, fifth_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                            }
                                        }
                                      // Then the third point trough
                                      // interpolation from the dof of the
                                      // cell in which the third point is
                                      else
                                        {
                                          auto &evaluation_point =
                                            this->evaluation_point;
                                          if (this->simulation_parameters.particlesParameters
                                                .order == 1)
                                            {
                                              this->system_matrix.set(
                                                global_index_overwrite,
                                                local_dof_indices_2[j],
                                                sp_2 *
                                                  this->fe.shape_value(
                                                    j, second_point_v) *
                                                  sum_line);
                                              local_interp_sol +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, second_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                            }

                                          if (this->simulation_parameters.particlesParameters
                                                .order == 2)
                                            {
                                              this->system_matrix.set(
                                                global_index_overwrite,
                                                local_dof_indices_2[j],
                                                sp_3 *
                                                    this->fe.shape_value(
                                                      j, second_point_v) *
                                                    sum_line +
                                                  tp_3 *
                                                    this->fe.shape_value(
                                                      j, third_point_v) *
                                                    sum_line);

                                              local_interp_sol +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, second_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_2 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, third_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                            }
                                          if (this->simulation_parameters.particlesParameters
                                                .order == 3)
                                            {
                                              this->system_matrix.set(
                                                global_index_overwrite,
                                                local_dof_indices_2[j],
                                                sp_4 *
                                                    this->fe.shape_value(
                                                      j, second_point_v) *
                                                    sum_line +
                                                  tp_4 *
                                                    this->fe.shape_value(
                                                      j, third_point_v) *
                                                    sum_line +
                                                  fp_4 *
                                                    this->fe.shape_value(
                                                      j, fourth_point_v) *
                                                    sum_line);

                                              local_interp_sol +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, second_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_2 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, third_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_3 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, fourth_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                            }
                                          if (this->simulation_parameters.particlesParameters
                                                .order > 4)
                                            {
                                              this->system_matrix.set(
                                                global_index_overwrite,
                                                local_dof_indices_2[j],
                                                this->fe.shape_value(
                                                  j, first_point_v) *
                                                  sum_line);
                                              local_interp_sol +=
                                                this->fe.shape_value(
                                                  j, first_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                            }
                                          if (this->simulation_parameters.particlesParameters
                                                .order == 4)
                                            {
                                              this->system_matrix.set(
                                                global_index_overwrite,
                                                local_dof_indices_2[j],
                                                sp_5 *
                                                    this->fe.shape_value(
                                                      j, second_point_v) *
                                                    sum_line +
                                                  tp_5 *
                                                    this->fe.shape_value(
                                                      j, third_point_v) *
                                                    sum_line +
                                                  fp1_5 *
                                                    this->fe.shape_value(
                                                      j, fourth_point_v) *
                                                    sum_line +
                                                  fp2_5 *
                                                    this->fe.shape_value(
                                                      j, fifth_point_v) *
                                                    sum_line);

                                              local_interp_sol +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, second_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_2 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, third_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_3 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, fourth_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                              local_interp_sol_4 +=
                                                1 *
                                                this->fe.shape_value(
                                                  j, fifth_point_v) *
                                                sum_line *
                                                evaluation_point(
                                                  local_dof_indices_2[j]);
                                            }
                                        }
                                    }
                                }
                            }



                          // Define the RHS of the stencil used for the
                          // IB
                          if (skip_stencil == false || do_rhs)
                            {
                              // Different boundary condition depending
                              // if the dof is vx ,vy or vz and if the
                              // problem we solve is 2d or 3d.
                              if (component_i == 0)
                                {
                                  double vx      = 0;
                                  double rhs_add = 0;
                                  if (dim == 2)
                                    {
                                      vx = -particles[p].omega[2] *
                                             particles[p].radius *
                                             ((support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed) /
                                              (support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed)
                                                .norm())[1] +
                                           particles[p].velocity[0];
                                    }
                                  if (dim == 3)
                                    {
                                      vx = particles[p].omega[1] *
                                             ((support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed) /
                                              (support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed)
                                                .norm())[2] *
                                             particles[p].radius -
                                           particles[p].omega[2] *
                                             ((support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed) /
                                              (support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed)
                                                .norm())[1] *
                                             particles[p].radius +
                                           particles[p].velocity[0];
                                    }

                                  auto &evaluation_point =
                                    this->evaluation_point;
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      1)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_2 -
                                                local_interp_sol * sp_2;
                                    }
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      2)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_3 -
                                                local_interp_sol * sp_3 -
                                                local_interp_sol_2 * tp_3;
                                    }
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      3)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_4 -
                                                local_interp_sol * sp_4 -
                                                local_interp_sol_2 * tp_4 -
                                                local_interp_sol_3 * fp_4;
                                    }
                                  if (this->simulation_parameters.particlesParameters.order >
                                      4)
                                    rhs_add = -local_interp_sol;

                                  if (this->simulation_parameters.particlesParameters.order ==
                                      4)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_5 -
                                                local_interp_sol * sp_5 -
                                                local_interp_sol_2 * tp_5 -
                                                local_interp_sol_3 * fp1_5 -
                                                local_interp_sol_4 * fp2_5;
                                    }

                                  auto &system_rhs = this->system_rhs;
                                  system_rhs(global_index_overwrite) =
                                    vx * sum_line + rhs_add;
                                  if (do_rhs)
                                    system_rhs(global_index_overwrite) =
                                      vx * sum_line -
                                      evaluation_point(global_index_overwrite) *
                                        sum_line;
                                }
                              else if (component_i == 1)
                                {
                                  double vy      = 0;
                                  double rhs_add = 0;
                                  if (dim == 2)
                                    {
                                      vy = particles[p].omega[2] *
                                             particles[p].radius *
                                             ((support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed) /
                                              (support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed)
                                                .norm())[0] +
                                           particles[p].velocity[1];
                                    }
                                  if (dim == 3)
                                    {
                                      vy = particles[p].omega[2] *
                                             ((support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed) /
                                              (support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed)
                                                .norm())[0] *
                                             particles[p].radius -
                                           particles[p].omega[0] *
                                             ((support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed) /
                                              (support_points
                                                 [local_dof_indices[i]] -
                                               center_immersed)
                                                .norm())[2] *
                                             particles[p].radius +
                                           particles[p].velocity[2];
                                    }

                                  auto &evaluation_point =
                                    this->evaluation_point;
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      1)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_2 -
                                                local_interp_sol * sp_2;
                                    }
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      2)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_3 -
                                                local_interp_sol * sp_3 -
                                                local_interp_sol_2 * tp_3;
                                    }
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      3)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_4 -
                                                local_interp_sol * sp_4 -
                                                local_interp_sol_2 * tp_4 -
                                                local_interp_sol_3 * fp_4;
                                    }
                                  if (this->simulation_parameters.particlesParameters.order >
                                      4)
                                    rhs_add = -local_interp_sol;
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      4)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_5 -
                                                local_interp_sol * sp_5 -
                                                local_interp_sol_2 * tp_5 -
                                                local_interp_sol_3 * fp1_5 -
                                                local_interp_sol_4 * fp2_5;
                                    }

                                  auto &system_rhs = this->system_rhs;
                                  system_rhs(global_index_overwrite) =
                                    vy * sum_line + rhs_add;
                                  if (do_rhs)
                                    system_rhs(global_index_overwrite) =
                                      vy * sum_line -
                                      evaluation_point(global_index_overwrite) *
                                        sum_line;
                                }
                              else if (component_i == 2 && dim == 3)
                                {
                                  double vz =
                                    particles[p].omega[0] *
                                      ((support_points[local_dof_indices[i]] -
                                        center_immersed) /
                                       (support_points[local_dof_indices[i]] -
                                        center_immersed)
                                         .norm())[1] *
                                      particles[p].radius -
                                    particles[p].omega[1] *
                                      ((support_points[local_dof_indices[i]] -
                                        center_immersed) /
                                       (support_points[local_dof_indices[i]] -
                                        center_immersed)
                                         .norm())[0] *
                                      particles[p].radius +
                                    particles[p].velocity[2];

                                  double rhs_add = 0;
                                  auto & evaluation_point =
                                    this->evaluation_point;
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      1)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_2 -
                                                local_interp_sol * sp_2;
                                    }
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      2)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_3 -
                                                local_interp_sol * sp_3 -
                                                local_interp_sol_2 * tp_3;
                                    }
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      3)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_4 -
                                                local_interp_sol * sp_4 -
                                                local_interp_sol_2 * tp_4 -
                                                local_interp_sol_3 * fp_4;
                                    }
                                  if (this->simulation_parameters.particlesParameters.order >
                                      4)
                                    rhs_add = -local_interp_sol;
                                  if (this->simulation_parameters.particlesParameters.order ==
                                      4)
                                    {
                                      rhs_add = -evaluation_point(
                                                  global_index_overwrite) *
                                                  sum_line * dof_5 -
                                                local_interp_sol * sp_5 -
                                                local_interp_sol_2 * tp_5 -
                                                local_interp_sol_3 * fp1_5 -
                                                local_interp_sol_4 * fp2_5;
                                    }

                                  auto &system_rhs = this->system_rhs;
                                  system_rhs(global_index_overwrite) =
                                    vz * sum_line + rhs_add;
                                  if (do_rhs)
                                    system_rhs(global_index_overwrite) =
                                      vz * sum_line -
                                      evaluation_point(global_index_overwrite) *
                                        sum_line;
                                }
                            }
                        }

                      if (component_i == dim)
                        {
                          // Applied equation on dof that have no equation
                          // define for them. those DOF become Dummy dof. This
                          // is usefull for high order cell or when a dof is
                          // only element of cell that are cut.
                          unsigned int vertex_per_cell =
                            GeometryInfo<dim>::vertices_per_cell;

                          unsigned int global_index_overwrite =
                            local_dof_indices[i];
                          bool dummy_dof = true;
                          for (unsigned int vi = 0; vi < vertex_per_cell; ++vi)
                            {
                              unsigned int v_index = cell->vertex_index(vi);
                              active_neighbors_set =
                                this->vertices_to_cell[v_index];
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
                                          unsigned int count_small_1 = 0;
                                          for (unsigned int q = 0;
                                               q < local_dof_indices_3.size();
                                               ++q)
                                            {
                                              // Count the number of dof
                                              // that are smaller or
                                              // larger then the radius
                                              // of the particles if all
                                              // the dof are on one side
                                              // the cell is not cut by
                                              // the boundary meaning we
                                              // dont have to do
                                              // anything
                                              if ((support_points
                                                     [local_dof_indices_3[q]] -
                                                   center_immersed)
                                                    .norm() <=
                                                  particles[p].radius)
                                                {
                                                  ++count_small_1;
                                                }
                                            }

                                          if (count_small_1 == 0 or
                                              count_small_1 ==
                                                local_dof_indices_3.size())
                                            dummy_dof = false;
                                        }
                                    }
                                }
                            }

                          if (dummy_dof)
                            {
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
template <bool                                              assemble_matrix,
          Parameters::SimulationControl::TimeSteppingMethod scheme,
          Parameters::VelocitySource::VelocitySourceType    velocity_source>
void
GLSSharpNavierStokesSolver<dim>::assembleGLS()
{
  auto &system_rhs = this->system_rhs;
  MPI_Barrier(this->mpi_communicator);
  if (assemble_matrix)
    this->system_matrix = 0;
  system_rhs = 0;
  // erase_inertia();
  double         viscosity_ = this->simulation_parameters.physical_properties.viscosity;
  Function<dim> *l_forcing_function = this->forcing_function;

  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(this->velocity_fem_degree,
                              this->simulation_parameters.fem_parameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients |
                            update_hessians);
  const unsigned int  dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int  n_q_points    = quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  FullMatrix<double>               local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>                   local_rhs(dofs_per_cell);
  std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<Tensor<1, dim>>          present_velocity_values(n_q_points);
  std::vector<Tensor<2, dim>>          present_velocity_gradients(n_q_points);
  std::vector<double>                  present_pressure_values(n_q_points);
  std::vector<Tensor<1, dim>>          present_pressure_gradients(n_q_points);
  std::vector<Tensor<1, dim>>          present_velocity_laplacians(n_q_points);
  std::vector<Tensor<2, dim>>          present_velocity_hess(n_q_points);

  Tensor<1, dim> force;

  std::vector<double>         div_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<Tensor<3, dim>> hess_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> laplacian_phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi_p(dofs_per_cell);

  // Values at previous time step for transient schemes
  std::vector<Tensor<1, dim>> p1_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> p2_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> p3_velocity_values(n_q_points);

  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  // support point
  MappingQ1<dim>                                immersed_map;
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(immersed_map,
                                       this->dof_handler,
                                       support_points);

  Point<dim> center_immersed;
  Point<dim> pressure_bridge;


  // Time steps and inverse time steps which is used for numerous calculations
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Vector for the BDF coefficients
  // The coefficients are stored in the following fashion :
  // 0 - n+1
  // 1 - n
  // 2 - n-1
  // 3 - n-2

  Vector<double> bdf_coefs;

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf1)
    bdf_coefs = bdf_coefficients(1, time_steps_vector);

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    bdf_coefs = bdf_coefficients(2, time_steps_vector);

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    bdf_coefs = bdf_coefficients(3, time_steps_vector);


  // Matrix of coefficients for the SDIRK methods
  // The lines store the information required for each step
  // Column 0 always refer to outcome of the step that is being calculated
  // Column 1 always refer to step n
  // Column 2+ refer to intermediary steps
  FullMatrix<double> sdirk_coefs;
  if (is_sdirk2(scheme))
    sdirk_coefs = sdirk_coefficients(2, dt);

  if (is_sdirk3(scheme))
    sdirk_coefs = sdirk_coefficients(3, dt);



  // Element size
  double h;



  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
    {
      bool assemble_bool = true;
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int k = 0; k < particles.size(); ++k)
            {
              unsigned int count_small = 0;
              if (dim == 2)
                {
                  center_immersed(0) = particles[k].position[0];
                  center_immersed(1) = particles[k].position[1];
                  // define arbitrary point on the boundary where the pressure
                  // will be link between the 2 domain
                  pressure_bridge(0) = particles[k].position[0] -
                                       particles[k].pressure_location[0];
                  pressure_bridge(1) = particles[k].position[1] -
                                       particles[k].pressure_location[1];
                }
              else if (dim == 3)
                {
                  center_immersed(0) = particles[k].position[0];
                  center_immersed(1) = particles[k].position[1];
                  center_immersed(2) = particles[k].position[2];
                  // define arbitrary point on the boundary where the pressure
                  // will be link between the 2 domain
                  pressure_bridge(0) = particles[k].position[0] -
                                       particles[k].pressure_location[0];
                  pressure_bridge(1) = particles[k].position[1] -
                                       particles[k].pressure_location[1];
                  pressure_bridge(2) = particles[k].position[2] -
                                       particles[k].pressure_location[2];
                }

              for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
                {
                  // count the number of dof that are smaller or larger then the
                  // radius of the particles if all the dof are on one side the
                  // cell is not cut by the boundary meaning we dont have to do
                  // anything
                  if ((support_points[local_dof_indices[j]] - center_immersed)
                        .norm() <= particles[k].radius)
                    {
                      ++count_small;
                    }
                }

              if (count_small != 0 and count_small != local_dof_indices.size())
                {
                  assemble_bool = false;
                  break;
                }
            }


          if (assemble_bool == true)
            {
              fe_values.reinit(cell);

              if (dim == 2)
                h = std::sqrt(4. * cell->measure() / M_PI) /
                    this->velocity_fem_degree;
              else if (dim == 3)
                h = pow(6 * cell->measure() / M_PI, 1. / 3.) /
                    this->velocity_fem_degree;

              local_matrix = 0;
              local_rhs    = 0;

              // Gather velocity (values, gradient and laplacian)
              auto &evaluation_point = this->evaluation_point;
              fe_values[velocities].get_function_values(
                evaluation_point, present_velocity_values);
              fe_values[velocities].get_function_gradients(
                evaluation_point, present_velocity_gradients);
              fe_values[velocities].get_function_laplacians(
                evaluation_point, present_velocity_laplacians);


              // Gather pressure (values, gradient)
              fe_values[pressure].get_function_values(evaluation_point,
                                                      present_pressure_values);
              fe_values[pressure].get_function_gradients(
                evaluation_point, present_pressure_gradients);


              // Calculate forcing term if there is a forcing function
              if (l_forcing_function)
                l_forcing_function->vector_value_list(
                  fe_values.get_quadrature_points(), rhs_force);

              // Gather the previous time steps depending on the number of
              // stages of the time integration scheme
              if (scheme !=
                  Parameters::SimulationControl::TimeSteppingMethod::steady)
                fe_values[velocities].get_function_values(this->solution_m1,
                                                          p1_velocity_values);

              if (time_stepping_method_has_two_stages(scheme))
                fe_values[velocities].get_function_values(this->solution_m2,
                                                          p2_velocity_values);

              if (time_stepping_method_has_three_stages(scheme))
                fe_values[velocities].get_function_values(this->solution_m3,
                                                          p3_velocity_values);

              // Loop over the quadrature points
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  // Calculation of the magnitude of the velocity for the
                  // stabilization parameter
                  const double u_mag =
                    std::max(present_velocity_values[q].norm(),
                             1e-12 * GLS_u_scale);

                  // Calculation of the GLS stabilization parameter. The
                  // stabilization parameter used is different if the simulation
                  // is steady or unsteady. In the unsteady case it includes the
                  // value of the time-step
                  double tau;
                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::steady)
                    tau =
                      1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                                     9 * std::pow(4 * viscosity_ / (h * h), 2));
                  else
                    tau =
                      1. /
                      std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                                9 * std::pow(4 * viscosity_ / (h * h), 2));

                  if (PSPG == false)
                    tau = 0;
                  // Gather the shape functions, their gradient and their
                  // laplacian for the velocity and the pressure
                  for (unsigned int k = 0; k < dofs_per_cell; ++k)
                    {
                      div_phi_u[k]  = fe_values[velocities].divergence(k, q);
                      grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                      phi_u[k]      = fe_values[velocities].value(k, q);
                      hess_phi_u[k] = fe_values[velocities].hessian(k, q);
                      phi_p[k]      = fe_values[pressure].value(k, q);
                      grad_phi_p[k] = fe_values[pressure].gradient(k, q);

                      for (int d = 0; d < dim; ++d)
                        laplacian_phi_u[k][d] = trace(hess_phi_u[k][d]);
                    }

                  // Establish the force vector
                  for (int i = 0; i < dim; ++i)
                    {
                      const unsigned int component_i =
                        this->fe.system_to_component_index(i).first;
                      force[i] = rhs_force[q](component_i);
                    }

                  // Calculate the divergence of the velocity
                  double present_velocity_divergence =
                    trace(present_velocity_gradients[q]);

                  // Calculate the strong residual for GLS stabilization
                  auto strong_residual =
                    present_velocity_gradients[q] * present_velocity_values[q] +
                    present_pressure_gradients[q] -
                    viscosity_ * present_velocity_laplacians[q] - force;

                  /* Adjust the strong residual in cases where the scheme is
                   transient.
                   The BDF schemes require values at previous time steps which
                   are stored in the p1, p2 and p3 vectors. The SDIRK scheme
                   require the values at the different stages, which are also
                   stored in the same arrays.
                   */

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf1)
                    strong_residual +=
                      bdf_coefs[0] * present_velocity_values[q] +
                      bdf_coefs[1] * p1_velocity_values[q];

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                    strong_residual +=
                      bdf_coefs[0] * present_velocity_values[q] +
                      bdf_coefs[1] * p1_velocity_values[q] +
                      bdf_coefs[2] * p2_velocity_values[q];

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                    strong_residual +=
                      bdf_coefs[0] * present_velocity_values[q] +
                      bdf_coefs[1] * p1_velocity_values[q] +
                      bdf_coefs[2] * p2_velocity_values[q] +
                      bdf_coefs[3] * p3_velocity_values[q];


                  if (is_sdirk_step1(scheme))
                    strong_residual +=
                      sdirk_coefs[0][0] * present_velocity_values[q] +
                      sdirk_coefs[0][1] * p1_velocity_values[q];

                  if (is_sdirk_step2(scheme))
                    {
                      strong_residual +=
                        sdirk_coefs[1][0] * present_velocity_values[q] +
                        sdirk_coefs[1][1] * p1_velocity_values[q] +
                        sdirk_coefs[1][2] * p2_velocity_values[q];
                    }

                  if (is_sdirk_step3(scheme))
                    {
                      strong_residual +=
                        sdirk_coefs[2][0] * present_velocity_values[q] +
                        sdirk_coefs[2][1] * p1_velocity_values[q] +
                        sdirk_coefs[2][2] * p2_velocity_values[q] +
                        sdirk_coefs[2][3] * p3_velocity_values[q];
                    }

                  // Matrix assembly
                  if (assemble_matrix)
                    {
                      // We loop over the column first to prevent recalculation
                      // of the strong jacobian in the inner loop
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          auto strong_jac =
                            (present_velocity_gradients[q] * phi_u[j] +
                             grad_phi_u[j] * present_velocity_values[q] +
                             grad_phi_p[j] - viscosity_ * laplacian_phi_u[j]);

                          if (is_bdf(scheme))
                            strong_jac += phi_u[j] * bdf_coefs[0];
                          if (is_sdirk(scheme))
                            strong_jac += phi_u[j] * sdirk_coefs[0][0];

                          for (unsigned int i = 0; i < dofs_per_cell; ++i)
                            {
                              local_matrix(i, j) +=
                                ( // Momentum terms
                                  viscosity_ * scalar_product(grad_phi_u[j],
                                                              grad_phi_u[i]) +
                                  present_velocity_gradients[q] * phi_u[j] *
                                    phi_u[i] +
                                  grad_phi_u[j] * present_velocity_values[q] *
                                    phi_u[i] -
                                  div_phi_u[i] * phi_p[j] +
                                  // Continuity
                                  phi_p[i] * div_phi_u[j]) *
                                fe_values.JxW(q);

                              // Mass matrix
                              if (is_bdf(scheme))
                                local_matrix(i, j) += phi_u[j] * phi_u[i] *
                                                      bdf_coefs[0] *
                                                      fe_values.JxW(q);

                              if (is_sdirk(scheme))
                                local_matrix(i, j) += phi_u[j] * phi_u[i] *
                                                      sdirk_coefs[0][0] *
                                                      fe_values.JxW(q);


                              local_matrix(i, j) += tau * strong_jac *
                                                    grad_phi_p[i] *
                                                    fe_values.JxW(q);

                              if (PSPG)
                                {
                                  // PSPG TAU term is currently disabled because
                                  // it does not alter the matrix sufficiently
                                  local_matrix(i, j) +=
                                    -tau * tau * tau * 4 / h / h *
                                    (present_velocity_values[q] * phi_u[j]) *
                                    strong_residual * grad_phi_p[i] *
                                    fe_values.JxW(q);
                                }

                              // Jacobian is currently incomplete
                              if (SUPG)
                                {
                                  local_matrix(i, j) +=
                                    tau *
                                    (strong_jac * (grad_phi_u[i] *
                                                   present_velocity_values[q]) +
                                     strong_residual *
                                       (grad_phi_u[i] * phi_u[j])) *
                                    fe_values.JxW(q);

                                  // SUPG TAU term is currently disabled because
                                  // it does not alter the matrix sufficiently
                                  local_matrix(i, j) +=
                                    -strong_residual *
                                    (grad_phi_u[i] *
                                     present_velocity_values[q]) *
                                    tau * tau * tau * 4 / h / h *
                                    (present_velocity_values[q] * phi_u[j]) *
                                    fe_values.JxW(q);
                                }
                            }
                        }
                    }

                  // Assembly of the right-hand side
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      // Navier-Stokes Residual
                      local_rhs(i) +=
                        ( // Momentum
                          -viscosity_ *
                            scalar_product(present_velocity_gradients[q],
                                           grad_phi_u[i]) -
                          present_velocity_gradients[q] *
                            present_velocity_values[q] * phi_u[i] +
                          present_pressure_values[q] * div_phi_u[i] +
                          force * phi_u[i] -
                          // Continuity
                          present_velocity_divergence * phi_p[i]) *
                        fe_values.JxW(q);

                      // Residual associated with BDF schemes
                      if (scheme == Parameters::SimulationControl::
                                      TimeSteppingMethod::bdf1)
                        local_rhs(i) -=
                          bdf_coefs[0] *
                          (present_velocity_values[q] - p1_velocity_values[q]) *
                          phi_u[i] * fe_values.JxW(q);

                      if (scheme == Parameters::SimulationControl::
                                      TimeSteppingMethod::bdf2)
                        local_rhs(i) -=
                          (bdf_coefs[0] *
                             (present_velocity_values[q] * phi_u[i]) +
                           bdf_coefs[1] * (p1_velocity_values[q] * phi_u[i]) +
                           bdf_coefs[2] * (p2_velocity_values[q] * phi_u[i])) *
                          fe_values.JxW(q);

                      if (scheme == Parameters::SimulationControl::
                                      TimeSteppingMethod::bdf3)
                        local_rhs(i) -=
                          (bdf_coefs[0] *
                             (present_velocity_values[q] * phi_u[i]) +
                           bdf_coefs[1] * (p1_velocity_values[q] * phi_u[i]) +
                           bdf_coefs[2] * (p2_velocity_values[q] * phi_u[i]) +
                           bdf_coefs[3] * (p3_velocity_values[q] * phi_u[i])) *
                          fe_values.JxW(q);

                      // Residuals associated with SDIRK schemes
                      if (is_sdirk_step1(scheme))
                        local_rhs(i) -=
                          (sdirk_coefs[0][0] *
                             (present_velocity_values[q] * phi_u[i]) +
                           sdirk_coefs[0][1] *
                             (p1_velocity_values[q] * phi_u[i])) *
                          fe_values.JxW(q);

                      if (is_sdirk_step2(scheme))
                        {
                          local_rhs(i) -=
                            (sdirk_coefs[1][0] *
                               (present_velocity_values[q] * phi_u[i]) +
                             sdirk_coefs[1][1] *
                               (p1_velocity_values[q] * phi_u[i]) +
                             sdirk_coefs[1][2] *
                               (p2_velocity_values[q] * phi_u[i])) *
                            fe_values.JxW(q);
                        }

                      if (is_sdirk_step3(scheme))
                        {
                          local_rhs(i) -=
                            (sdirk_coefs[2][0] *
                               (present_velocity_values[q] * phi_u[i]) +
                             sdirk_coefs[2][1] *
                               (p1_velocity_values[q] * phi_u[i]) +
                             sdirk_coefs[2][2] *
                               (p2_velocity_values[q] * phi_u[i]) +
                             sdirk_coefs[2][3] *
                               (p3_velocity_values[q] * phi_u[i])) *
                            fe_values.JxW(q);
                        }

                      // PSPG GLS term
                      if (PSPG)
                        local_rhs(i) += -tau *
                                        (strong_residual * grad_phi_p[i]) *
                                        fe_values.JxW(q);

                      // SUPG GLS term
                      if (SUPG)
                        {
                          local_rhs(i) +=
                            -tau *
                            (strong_residual *
                             (grad_phi_u[i] * present_velocity_values[q])) *
                            fe_values.JxW(q);
                        }
                    }
                }

              cell->get_dof_indices(local_dof_indices);

              // The non-linear solver assumes that the nonzero constraints have
              // already been applied to the solution
              const AffineConstraints<double> &constraints_used =
                this->zero_constraints;
              // initial_step ? nonzero_constraints : zero_constraints;
              if (assemble_matrix)
                {
                  constraints_used.distribute_local_to_global(
                    local_matrix,
                    local_rhs,
                    local_dof_indices,
                    this->system_matrix,
                    system_rhs);
                }
              else
                {
                  constraints_used.distribute_local_to_global(local_rhs,
                                                              local_dof_indices,
                                                              system_rhs);
                }
            }
          else
            {
              // could assemble someting in the cells tahat are cut  have to
              // code it here
            }
        }
    }


  if (assemble_matrix)
    this->system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_system");

  if (this->simulation_parameters.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::none)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
    }

  else if (this->simulation_parameters.velocitySource.type ==
           Parameters::VelocitySource::VelocitySourceType::srf)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
    }
  vertices_cell_mapping();
  sharp_edge();
}
template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_rhs");

  if (this->simulation_parameters.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::none)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
    }
  if (this->simulation_parameters.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::srf)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
    }
  vertices_cell_mapping();
  sharp_edge();
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::solve()
{
  read_mesh_and_manifolds(this->triangulation,
                          this->simulation_parameters.mesh,
                          this->simulation_parameters.manifolds_parameters,
                          this->simulation_parameters.restart_parameters.restart,
                          this->simulation_parameters.boundary_conditions);

  define_particles();
  this->setup_dofs();

  // To change once refinement is split into two function
  double temp_refine = this->simulation_parameters.mesh_adaptation.refinement_fraction;
  double temp_coarse = this->simulation_parameters.mesh_adaptation.coarsening_fraction;
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


  this->set_initial_condition(this->simulation_parameters.initial_condition->type,
                              this->simulation_parameters.restart_parameters.restart);

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      if (this->simulation_control->is_at_start())
        this->first_iteration();
      else
        {
          refine_ib();
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          this->iterate();
        }

      this->postprocess_fd(false);

      this->finish_time_step();

      if (this->simulation_parameters.particlesParameters.calculate_force_ib)
        force_on_ib();
      write_force_ib();
      MPI_Barrier(this->mpi_communicator);
    }

  if (this->simulation_parameters.particlesParameters.calculate_force_ib)


    this->finish_simulation();
}


// Pre-compile the 2D and 3D versopm solver to ensure that the library is
// valid before we actually compile the final solver
template class GLSSharpNavierStokesSolver<2>;
template class GLSSharpNavierStokesSolver<3>;
