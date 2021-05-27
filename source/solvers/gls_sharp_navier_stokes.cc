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
    // Find all the cells around each vertices
    TimerOutput::Scope t(this->computing_timer, "vertices_to_cell_map");
    vertices_to_cell.clear();
    const auto &cell_iterator = this->dof_handler.active_cell_iterators();

    // // Loop on all the cells and find their vertices to fill the map of sets of cells around each vertex
    for (const auto &cell : cell_iterator)
    {
        if (cell->is_locally_owned() || cell->is_ghost())
        {
            const unsigned int vertices_per_cell =
                    GeometryInfo<dim>::vertices_per_cell;
            for (unsigned int i = 0; i < vertices_per_cell; i++)
            {
                // First obtain vertex index
                unsigned int v_index = cell->vertex_index(i);

                // Insert the cell into the set of cell around that vertex.
                vertices_to_cell[v_index].insert(cell);
            }
        }
    }
}

template <int dim>
typename DoFHandler<dim>::active_cell_iterator
GLSSharpNavierStokesSolver<dim>::find_cell_around_point_with_neighbors(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                                       Point<dim>             point)
{
    //Find the cell around a point based on an initial cell.

    //Find the cells around the initial cell ( cells that share a vertex with the original cell).
    std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors_set = find_cells_around_cell(cell);
    //Loop over that group of cells
    for (unsigned int i = 0; i < active_neighbors_set.size(); ++i){
        bool inside_cell=point_inside_cell(active_neighbors_set[i], point);
            if(inside_cell) {
                return active_neighbors_set[i];
                }


    }
    // The cell is not found near the initial cell so we use the cell tree algorithm instead (much slower).
    std::cout << "Cell not found around " << point << std::endl;
    return find_cell_around_point_with_tree(this->dof_handler,point);

}

template <int dim>
bool
GLSSharpNavierStokesSolver<dim>::point_inside_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,Point<dim>             point)
{
    try
    {
        const Point<dim, double> p_cell =
                this->mapping->transform_real_to_unit_cell(cell, point);
        const double dist = GeometryInfo<dim>::distance_to_unit_cell(p_cell);
        // if the cell contains the point, the distance is equal to 0
        if (dist == 0)
        {
            // The cell is found so we return it and exit the function.

            return true;
        }
    }
    catch (const typename MappingQGeneric<dim>::ExcTransformationFailed &)
    {}
    return false;
}

template <int dim>
std::vector<typename DoFHandler<dim>::active_cell_iterator>
GLSSharpNavierStokesSolver<dim>::find_cells_around_cell(const typename DoFHandler<dim>::active_cell_iterator &cell)
{
    // Find all the cells that share a vertex with a reference cell including the initial cell.
    std::set<typename DoFHandler<dim>::active_cell_iterator> neighbors_cells;
    // Loop over the vertices of the initial cell and find all the cells around each vertex and add them to the set of cells around the reference cell.
    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; i++)
    {
        unsigned int v_index=cell->vertex_index(i);
        neighbors_cells.insert(this->vertices_to_cell[v_index].begin(),this->vertices_to_cell[v_index].end());
    }
    // Transform the set into a vector.
    std::vector<typename DoFHandler<dim>::active_cell_iterator> cells_sharing_vertices(neighbors_cells.begin(),neighbors_cells.end());
    return cells_sharing_vertices;
}


template <int dim>
void
GLSSharpNavierStokesSolver<dim>::clear_line_in_matrix(const typename DoFHandler<dim>::active_cell_iterator &cell, unsigned int dof_index){
    // Clear a line in the matrix based on an initial cell that contains the DOF.
    // This function ensure that if the dof is a ghost all the entry of the matrix will be erased.


    // if the dof is locally owned we can use the default function of deal ii for matrix (most of the time this is ok)
    if (this->locally_owned_dofs.is_element(dof_index)){
        this->system_matrix.clear_row(dof_index);
    }
    else{
        // This DOf is special, it's at a frontier between two processors. Only one of the processors owns it. If we have reached this point, we are not on that processor.
        // In this case  the clear_row function doesn't clear the entire row it only clear DOFs that are owned by the row. This is an issue.
        // To fix that we force all DOFs contain in ghost cells to be put to 0.
        std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors_set = find_cells_around_cell(cell);

        const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
        std::vector<types::global_dof_index> local_dof_indices_iter(dofs_per_cell);
        // Loop over the neighbours cells and erase the entry of the matrix that could be linked to neighbours ghost cells.
        for (unsigned int i = 0; i < active_neighbors_set.size(); ++i) {
            const auto &cell_3 = active_neighbors_set[i];
            cell_3->get_dof_indices(local_dof_indices_iter);
            if (std::find(local_dof_indices_iter.begin(), local_dof_indices_iter.end(), dof_index) !=
            local_dof_indices_iter.end()) {
                for (unsigned int k = 0; k < local_dof_indices_iter.size(); ++k) {
                    this->system_matrix.set(dof_index, local_dof_indices_iter[k], 0);
                }
            }
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
  TimerOutput::Scope t(this->computing_timer, "force_evaluation");
  // Calculate the torque and force on an immersed boundary
  // The boundary is a circle in 2D or a sphere in 3D

  std::vector<typename DoFHandler<dim>::active_cell_iterator>
                                                              active_neighbors_set;
  std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;

  const double min_cell_diameter =
    GridTools::minimal_cell_diameter(*this->triangulation);

  double dr  = (min_cell_diameter) / std::sqrt(2);
  double rho = this->simulation_parameters.particlesParameters.density;
  // Define stuff for later use
  using numbers::PI;
  Point<dim> center_immersed;

  if (dim == 2)
    {
      // Define general stuff useful for the evaluation of force with stencil
      QGauss<dim>   q_formula(this->fe->degree + 1);
      FEValues<dim> fe_values(*this->fe, q_formula, update_quadrature_points);

      double mu = this->simulation_parameters.physical_properties.viscosity;

      std::map<types::global_dof_index, Point<dim>> support_points;
      DoFTools::map_dofs_to_support_points(*this->mapping,
                                           this->dof_handler,
                                           support_points);


      std::vector<types::global_dof_index> local_dof_indices(
        this->fe->dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices_2(
        this->fe->dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices_3(
        this->fe->dofs_per_cell);
      const unsigned int nb_evaluation =
        this->simulation_parameters.particlesParameters.nb_force_eval;

      // Loop on all the particles
      for (unsigned int p = 0; p < particles.size(); ++p)
        {
          // Define the center
          center_immersed = particles[p].position;

          // Initialize the output variable for this particle

          double t_torque = 0;
          // unsigned int nb_eval=0;
          double fx_v   = 0;
          double fy_v   = 0;
          double fx_p_2 = 0;
          double fy_p_2 = 0;


          // loop on all the evaluation points

          for (unsigned int i = 0; i < nb_evaluation; ++i)
            {
              // Define the normal to the surface evaluated and the vector that
              // is along the surface.
              Tensor<1, dim, double> surf_normal;
              Tensor<1, dim, double> surf_vect;


              surf_normal[0] = dr * cos(i * 2 * PI / (nb_evaluation));
              surf_normal[1] = dr * sin(i * 2 * PI / (nb_evaluation));
              surf_vect[0]   = -surf_normal[1];
              surf_vect[1]   = surf_normal[0];
              double da      = 2 * PI * particles[p].radius / (nb_evaluation);
              // Define the reference point for the surface evaluated.
              // the step ratio is the proportion constant used to multiply the
              // size of the smallest cell (dr) the step made are then:
              // step_ratio * dr
              double           step_ratio = 0.5;
              const Point<dim> eval_point(
                particles[p].radius * cos(i * 2 * PI / (nb_evaluation)) +
                  center_immersed(0),
                particles[p].radius * sin(i * 2 * PI / (nb_evaluation)) +
                  center_immersed(1));

              // Step in the normal direction to the surface until the point
              // used for the IB stencil is not in a cell that is cut by the
              // boundary.
              unsigned int     nb_step    = 0;
              bool             cell_found = false;
              const Point<dim> eval_point_2(
                eval_point[0] + surf_normal[0] * (nb_step + 1) * step_ratio,
                eval_point[1] + surf_normal[1] * (nb_step + 1) * step_ratio);

              // Step in the normal direction to the surface until the point
              // used for the IB stencil is not in a cell that is cut by the
              // boundary.
              while (cell_found == false)
                {
                  // define the new point
                  Point<dim> eval_point_iter(
                    eval_point[0] + surf_normal[0] * (nb_step + 1) * step_ratio,
                    eval_point[1] +
                      surf_normal[1] * (nb_step + 1) * step_ratio);

                  const auto &cell_iter =
                    find_cell_around_point_with_tree(this->dof_handler,
                                                     eval_point_iter);



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
                          // Count the number of dofs that are smaller or larger
                          // then the radius of the particles if all the dofs
                          // are on one side the cell is not cut by the boundary
                          // meaning we don’t have to do anything
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

                      // Step a bit further away from the boundary.
                      if (cell_found == false)
                        nb_step += 1;
                    }
                  else
                    {
                      break;
                    }
                }


              // When the point is found outside the cells that are cut by the
              // boundary we define the 3 points that will be used to create
              // interpolation and  extrapolation of the solution  to evaluate
              // the force on the boundary
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



              // Check if the cell is locally owned before doing the evaluation.
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
                  // the velocity at the center of the particle this simplifies
                  // the evaluation of the force if the particle is moving
                  u_1[0] = -particles[p].omega[2] * particles[p].radius *
                           sin(i * 2 * PI / (nb_evaluation));
                  u_1[1] = particles[p].omega[2] * particles[p].radius *
                           cos(i * 2 * PI / (nb_evaluation));

                  // Projection of the speed of the boundary on the plan of the
                  // surface used for evaluation
                  double U1 = (surf_vect[0] * u_1[0] + surf_vect[1] * u_1[1]) /
                              surf_vect.norm();


                  // Used support function of the cell to define the
                  // interpolation of the velocity
                  Point<dim> second_point_v =
                          this->mapping->transform_real_to_unit_cell(cell_2,
                                                             second_point);
                  Point<dim> third_point_v =
                          this->mapping->transform_real_to_unit_cell(cell_3,
                                                             third_point);
                  Point<dim> fourth_point_v =
                          this->mapping->transform_real_to_unit_cell(cell_4,
                                                             fourth_point);

                  // Initialize the component of the velocity
                  u_2[0] = 0;
                  u_2[1] = 0;
                  u_3[0] = 0;
                  u_3[1] = 0;
                  cell_2->get_dof_indices(local_dof_indices);
                  cell_3->get_dof_indices(local_dof_indices_2);
                  cell_4->get_dof_indices(local_dof_indices_3);



                  // Define the interpolation of the cell in order to have the
                  // solution at the point previously defined
                  for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
                    {
                      auto &present_solution = this->present_solution;
                      const unsigned int component_i =
                        this->fe->system_to_component_index(j).first;
                      if (component_i < dim)
                        {
                          u_2[component_i] +=
                            this->fe->shape_value(j, second_point_v) *
                            present_solution(local_dof_indices[j]);

                          u_3[component_i] +=
                            this->fe->shape_value(j, third_point_v) *
                            present_solution(local_dof_indices_2[j]);
                        }
                      if (component_i == dim)
                        {
                          P1 += this->fe->shape_value(j, second_point_v) *
                                present_solution(local_dof_indices[j]);
                          P2 += this->fe->shape_value(j, third_point_v) *
                                present_solution(local_dof_indices_2[j]);
                          P3 += this->fe->shape_value(j, fourth_point_v) *
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
          particles[p].forces[0]  = fx_p_2_ + fx_v_;
          particles[p].forces[1]  = fy_p_2_ + fy_v_;
          particles[p].torques[2] = t_torque_;
        }
    }


  // Same structure as for the 2d case but used 3d variables  so there is 1 more
  // vector on the surface for the evaluation
  if (dim == 3)
    {
      QGauss<dim>   q_formula(this->fe->degree + 1);
      FEValues<dim> fe_values(*this->fe, q_formula, update_quadrature_points);

      double mu = this->simulation_parameters.physical_properties.viscosity;

      std::map<types::global_dof_index, Point<dim>> support_points;
      DoFTools::map_dofs_to_support_points(*this->mapping,
                                           this->dof_handler,
                                           support_points);


      std::vector<types::global_dof_index> local_dof_indices(
        this->fe->dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices_2(
        this->fe->dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices_3(
        this->fe->dofs_per_cell);
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
                              // Count the number of dofs that are smaller or
                              // larger than the radius of the particles if all
                              // the dofs are on one side the cell is not cut by
                              // the boundary meaning we don’t have to do
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
                      // boundary on the reference point in 3d the only
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
                              this->mapping->transform_real_to_unit_cell(cell_2,
                                                                 second_point);
                      Point<dim> third_point_v =
                              this->mapping->transform_real_to_unit_cell(cell_3,
                                                                 third_point);
                      Point<dim> fourth_point_v =
                              this->mapping->transform_real_to_unit_cell(cell_4,
                                                                 fourth_point);


                      cell_3->get_dof_indices(local_dof_indices_2);
                      for (unsigned int j = 0; j < local_dof_indices.size();
                           ++j)
                        {
                          const unsigned int component_i =
                            this->fe->system_to_component_index(j).first;
                          auto &present_solution = this->present_solution;
                          if (component_i < dim)
                            {
                              u_2[component_i] +=
                                this->fe->shape_value(j, second_point_v) *
                                present_solution(local_dof_indices[j]);

                              u_3[component_i] +=
                                this->fe->shape_value(j, third_point_v) *
                                present_solution(local_dof_indices_2[j]);
                            }
                          if (component_i == dim)
                            {
                              P1 += this->fe->shape_value(j, second_point_v) *
                                    present_solution(local_dof_indices[j]);
                              P2 += this->fe->shape_value(j, third_point_v) *
                                    present_solution(local_dof_indices_2[j]);
                              P3 += this->fe->shape_value(j, fourth_point_v) *
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
            Utilities::MPI::sum(torque_x, this->mpi_communicator) * rho;
          double t_torque_y =
            Utilities::MPI::sum(torque_y, this->mpi_communicator) * rho;
          double t_torque_z =
            Utilities::MPI::sum(torque_z, this->mpi_communicator) * rho;
          double fx_p_2_ =
            Utilities::MPI::sum(fx_p_2, this->mpi_communicator) * rho;
          double fy_p_2_ =
            Utilities::MPI::sum(fy_p_2, this->mpi_communicator) * rho;
          double fz_p_2_ =
            Utilities::MPI::sum(fz_p_2, this->mpi_communicator) * rho;
          double fx_v_ =
            Utilities::MPI::sum(fx_v, this->mpi_communicator) * rho;
          double fy_v_ =
            Utilities::MPI::sum(fy_v, this->mpi_communicator) * rho;
          double fz_v_ =
            Utilities::MPI::sum(fz_v, this->mpi_communicator) * rho;
          // unsigned int nb_eval_total   = Utilities::MPI::sum(nb_eval,
          // this->mpi_communicator);
          particles[p].forces[0]  = fx_p_2_ + fx_v_;
          particles[p].forces[1]  = fy_p_2_ + fy_v_;
          particles[p].forces[2]  = fz_p_2_ + fz_v_;
          particles[p].torques[0] = t_torque_x;
          particles[p].torques[1] = t_torque_y;
          particles[p].torques[2] = t_torque_z;
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
  auto &present_solution = this->present_solution;
  if (this->simulation_control->is_output_iteration())
    this->write_output_results(present_solution);

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
          std::tie(cell_is_cut,std::ignore,local_dof_indices)=cell_cut(cell,local_dof_indices,support_points);

          if (cell_is_cut==false)
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
std::tuple<bool,unsigned int,std::vector<types::global_dof_index> >
GLSSharpNavierStokesSolver<dim>::cell_cut(const typename DoFHandler<dim>::active_cell_iterator &cell, std::vector<types::global_dof_index> &local_dof_indices ,std::map<types::global_dof_index, Point<dim>> &support_points)
{
// Check if a cell is cut and if it's rerun the particle by which it's cut and the local DOFs index.
// The check is done by counting the number of DOFs that is on either side of the boundary define by a particle.

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
        if (nb_dof_inside != 0 && nb_dof_inside != local_dof_indices.size()){
            //Some of the DOFs are inside the boundary, some are outside.
            // This mean that the cell is cut so we return that information and the index of the particle that cut the cell
            // as well as the container containing local DOF of the cell.
            return {true,p,local_dof_indices};
        }
    }
    return {false,0,local_dof_indices};
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
  const unsigned int dofs_per_cell   = this->fe->dofs_per_cell;

  int order=this->simulation_parameters
          .particlesParameters
          .order;



  IBStencils<dim> stencil;
  std::vector<double> ib_coef=stencil.coefficients(order);

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
  double dt=time_steps_vector[0];
  if(Parameters::SimulationControl::TimeSteppingMethod::steady== this->simulation_parameters.simulation_control.method)
      dt=1;


  // impose pressure reference in each of the particle
  for (unsigned int p = 0; p < particles.size(); ++p)
    {
      const auto &cell =find_cell_around_point_with_tree(this->dof_handler,
                                                 particles[p].pressure_location+particles[p].position);

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

          sum_line=sum_line/dt;
          // Clear the line in the matrix
          unsigned int inside_index = local_dof_indices[dim];
          clear_line_in_matrix(cell, inside_index);
          // this->system_matrix.clear_row(inside_index)
          // is not reliable on edge case

          // Set the new equation for the first pressure dofs of the
          // cell. this is the new reference pressure inside a
          // particle

            this->system_matrix.set(inside_index, inside_index, sum_line);
            this->system_rhs(inside_index) = 0 - this->evaluation_point(inside_index) * sum_line;
        }
    }


  // Loop on all the cell to define if the sharp edge cut them
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_locally_owned())
        {
          double sum_line = 0;
          fe_values.reinit(cell);

          // Define the order of magnitude for the stencil.
          for (unsigned int qf = 0; qf < n_q_points; ++qf)
            sum_line += fe_values.JxW(qf);

          sum_line=sum_line/dt;

          // Check if the cell is cut or not by the IB
          auto [cell_is_cut,p,local_dof_indices]=cell_cut(cell,local_dof_indices_2,support_points);

          if (cell_is_cut)
          {
              // If we are here, the cell is cut by the IB.
              // Loops on the dof that represents the velocity  component
              // and pressure separately
              for (unsigned int i = 0; i < local_dof_indices.size(); ++i) {
                  const unsigned int component_i = this->fe->system_to_component_index(i).first;
                  if (component_i < dim) {
                      // We are working on the velocity of the cell cut
                      // loops on the dof that are for vx or vy separately
                      // loops on all the dof of the cell that represent
                      // a specific component

                      // define which dof is going to be redefined
                      unsigned int global_index_overwrite =
                              local_dof_indices[i];

                      // Clear the current line of this dof
                      clear_line_in_matrix(cell, global_index_overwrite);

                      // Define the points for the IB stencil, based on the order and the particle position as well as the DOF position.
                      // Depending on the order, the output variable "point" change definition. In the case of stencil orders 1 to 4 the variable point returns the position of the DOF directly.
                      // In the case of high order stencil, it returns the position of the point that is on the IB.
                      // The variable "interpolation points" return the points used to define the cell used for the stencil definition and
                      // the locations of the points use in the stencil calculation.

                      auto[point, interpolation_points]=stencil.points(order, particles[p],
                                                                       support_points[local_dof_indices[i]]);

                      // Find the cell used for the stencil definition.
                      auto cell_2 = find_cell_around_point_with_neighbors(cell, interpolation_points[
                              stencil.nb_points(order) - 1]);
                      cell_2->get_dof_indices(local_dof_indices_2);

                      bool skip_stencil = false;

                      // Check if the DOF intersect the IB
                      bool dof_on_ib = false;

                      // Check if this dof is a dummy dof or directly on IB and
                      // Check if the point used to define the cell used for the definition of the stencil ("cell_2") is on
                      // a face between the cell that is cut ("cell") and the "cell_2".
                      bool point_in_cell=point_inside_cell(cell,interpolation_points[
                              stencil.nb_points(order) - 1]);
                      if (cell_2 == cell || point_in_cell) {
                          // Give the DOF an approximated value. This help
                          // with pressure shock when the DOF passe from part of the boundary
                          // to the fluid.
                          this->system_matrix.set(global_index_overwrite,
                                                  global_index_overwrite,
                                                  sum_line);
                          skip_stencil = true;

                          // Tolerence to define a intersection of
                          // the DOF and IB
                          if (abs((support_points[local_dof_indices[i]] - particles[p].position).norm() -
                                  particles[p].radius) <= 1e-12 * dr) {
                              dof_on_ib = true;
                          }
                      }
                      // Define the variable used for the
                      // extrapolation of the actual solution at the
                      // boundaries in order to define the correction

                      // Define the unit cell points for the points
                      // used in the stencil.
                      std::vector<Point<dim>> unite_cell_interpolation_points(ib_coef.size());
                      unite_cell_interpolation_points[0] = this->mapping->transform_real_to_unit_cell(cell_2, point);
                      for (unsigned int j = 1;
                           j < ib_coef.size();
                           ++j) {
                          unite_cell_interpolation_points[j] =
                                  this->mapping->transform_real_to_unit_cell(cell_2, interpolation_points[j - 1]);
                      }

                      std::vector<double> local_interp_sol(ib_coef.size());

                      // Define the new matrix entry for this dof
                      if (skip_stencil == false) {
                          for (unsigned int j = 0;
                               j < local_dof_indices_2.size();
                               ++j) {
                              const unsigned int component_j =
                                      this->fe->system_to_component_index(j)
                                              .first;
                              if (component_j == component_i) {

                                  //  Define the solution at each point used for the stencil and applied the stencil for the specific DOF.
                                  //  For stencils of order 4 or higher, the stencil is defined
                                  //  through direct extrapolation of the cell. This can only be done when using a structured mesh
                                  //  as this required a mapping of a point outside of a cell.

                                  // Define the local matrix entries of this DOF based on its contribution of each of the points used in
                                  // the stencil definition and the coefficient associated with this point.
                                  // This loop defined the current solution at the boundary using the same stencil. This is needed to define the residual.
                                  double local_matrix_entry = 0;
                                  for (unsigned int k = 0;
                                  k < ib_coef.size();++k) {
                                      local_matrix_entry += this->fe->shape_value(
                                              j, unite_cell_interpolation_points[k]) * ib_coef[k];
                                      local_interp_sol[k] += this->fe->shape_value(
                                              j, unite_cell_interpolation_points[k]) * this->evaluation_point(local_dof_indices_2[j]);
                                      }
                                  // update the matrix.
                                  this->system_matrix.set(
                                          global_index_overwrite,
                                          local_dof_indices_2[j], local_matrix_entry *
                                          sum_line);

                              }
                          }
                      }

                      // Define the RHS of the equation.

                      double rhs_add = 0;
                      // Different boundary conditions depending
                      // on the component index of the DOF and
                      // the dimension.
                      double v_ib = stencil.ib_velocity(particles[p], support_points[local_dof_indices[i]],
                                                           component_i);

                      for (unsigned int k = 0;k < ib_coef.size();++k) {
                          rhs_add += -local_interp_sol[k] * ib_coef[k] * sum_line;
                      }
                      this->system_rhs(global_index_overwrite) =
                                  v_ib * sum_line + rhs_add;

                      if (dof_on_ib)
                          // Dof is on the immersed boundary
                          this->system_rhs(global_index_overwrite) =
                                  v_ib * sum_line -
                                  this->evaluation_point(global_index_overwrite) *
                                  sum_line;

                      if (skip_stencil && dof_on_ib==false)
                              // Impose the value of the dummy dof. This help
                              // with pressure variation when the IB is
                              // moving.
                              this->system_rhs(global_index_overwrite) =
                                      sum_line * v_ib - this->evaluation_point(
                                              global_index_overwrite) *
                                                        sum_line;
                  }


                  if (component_i == dim)
                        {
                          // Applied equation on dof that have no equation
                          // defined for them. those DOF become Dummy dof. This
                          // is usefull for high order cells or when a dof is
                          // only element of cells that are cut.
                          unsigned int global_index_overwrite =
                            local_dof_indices[i];
                          bool dummy_dof = true;

                          // To check if the pressure dof is a dummy. first check if the matrix entry is close to 0.
                          if(abs(this->system_matrix.el(global_index_overwrite,
                                                     global_index_overwrite))<=1e-16 * dr) {
                              // If the matrix entry on the diagonal of this DOF is close to zero, check if all the cells close are cut.
                              // If it's the case, the DOF is a dummy DOF.
                              active_neighbors_set = find_cells_around_cell(cell);
                              for (unsigned int m = 0; m < active_neighbors_set.size(); m++) {
                                  const auto &cell_3 = active_neighbors_set[m];
                                  cell_3->get_dof_indices(local_dof_indices_3);
                                  for (unsigned int o = 0;
                                       o < local_dof_indices_3.size();
                                       ++o) {
                                      if (global_index_overwrite ==
                                          local_dof_indices_3[o]) {
                                          // cell_3 contain the same dof
                                          // check if this cell is cut if
                                          // it's not cut this dof must not
                                          // be overwritten
                                          bool cell_is_cut;
                                          std::tie(cell_is_cut, std::ignore, std::ignore) = cell_cut(cell_3,
                                                                                                     local_dof_indices_3,
                                                                                                     support_points);

                                          if (cell_is_cut == false) {
                                              dummy_dof = false;
                                              break;
                                          }
                                      }
                                  }
                                  if (dummy_dof == false)
                                      break;
                              }

                              if (dummy_dof) {
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
  double viscosity = this->simulation_parameters.physical_properties.viscosity;
  Function<dim> *l_forcing_function = this->forcing_function;

  QGauss<dim>        quadrature_formula(this->number_quadrature_points);
  FEValues<dim>      fe_values(*this->mapping,
                          *this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients |
                            update_hessians);
  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  FullMatrix<double>               local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>                   local_rhs(dofs_per_cell);
  std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_iter(dofs_per_cell);
  std::vector<Tensor<1, dim>>          present_velocity_values(n_q_points);
  std::vector<Tensor<2, dim>>          present_velocity_gradients(n_q_points);
  std::vector<double>                  present_pressure_values(n_q_points);
  std::vector<Tensor<1, dim>>          present_pressure_gradients(n_q_points);
  std::vector<Tensor<1, dim>>          present_velocity_laplacians(n_q_points);
  std::vector<Tensor<2, dim>>          present_velocity_hess(n_q_points);

  Tensor<1, dim> force;
  Tensor<1, dim> beta_force = this->beta;

  // Velocity dependent source term
  //----------------------------------
  // Angular velocity of the rotating frame. This is always a 3D vector even
  // in 2D.
  Tensor<1, dim> omega_vector;

  double omega_z  = this->simulation_parameters.velocitySource.omega_z;
  omega_vector[0] = this->simulation_parameters.velocitySource.omega_x;
  omega_vector[1] = this->simulation_parameters.velocitySource.omega_y;
  if (dim == 3)
      omega_vector[2] = this->simulation_parameters.velocitySource.omega_z;

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

  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(*this->mapping,
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
  auto & evaluation_point = this->evaluation_point;
  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);


          bool cell_is_cut;
          std::tie(cell_is_cut,std::ignore,local_dof_indices)=cell_cut(cell,local_dof_indices_iter,support_points);


          if (cell_is_cut == false)
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
              fe_values[velocities].get_function_values(evaluation_point,
                                                        present_velocity_values);
              fe_values[velocities].get_function_gradients(
                      evaluation_point, present_velocity_gradients);
              fe_values[velocities].get_function_laplacians(
                      evaluation_point, present_velocity_laplacians);

              // Gather pressure (values, gradient)
              fe_values[pressure].get_function_values(evaluation_point,
                                                      present_pressure_values);
              fe_values[pressure].get_function_gradients(
                      evaluation_point, present_pressure_gradients);

              std::vector<Point<dim>> quadrature_points =
                      fe_values.get_quadrature_points();

              // Calculate forcing term if there is a forcing function
              if (l_forcing_function)
                  l_forcing_function->vector_value_list(quadrature_points, rhs_force);

              // Gather the previous time steps depending on the number of
              // stages of the time integration scheme
              if (scheme !=
                  Parameters::SimulationControl::TimeSteppingMethod::steady)
                  fe_values[velocities].get_function_values(
                          this->previous_solutions[0], p1_velocity_values);

              if (time_stepping_method_has_two_stages(scheme))
                  fe_values[velocities].get_function_values(this->solution_stages[0],
                                                            p2_velocity_values);

              if (time_stepping_method_has_three_stages(scheme))
                  fe_values[velocities].get_function_values(this->solution_stages[1],
                                                            p3_velocity_values);

              if (time_stepping_method_uses_two_previous_solutions(scheme))
                  fe_values[velocities].get_function_values(
                          this->previous_solutions[1], p2_velocity_values);

              if (time_stepping_method_uses_three_previous_solutions(scheme))
                  fe_values[velocities].get_function_values(
                          this->previous_solutions[2], p3_velocity_values);

              // Loop over the quadrature points
              for (unsigned int q = 0; q < n_q_points; ++q)
              {
                  // Gather into local variables the relevant fields
                  const Tensor<1, dim> velocity = present_velocity_values[q];
                  const Tensor<2, dim> velocity_gradient =
                          present_velocity_gradients[q];
                  const double present_velocity_divergence =
                          trace(velocity_gradient);
                  const Tensor<1, dim> p1_velocity = p1_velocity_values[q];
                  const Tensor<1, dim> p2_velocity = p2_velocity_values[q];
                  const Tensor<1, dim> p3_velocity = p3_velocity_values[q];
                  const double current_pressure    = present_pressure_values[q];



                  // Calculation of the magnitude of the velocity for the
                  // stabilization parameter
                  const double u_mag =
                          std::max(velocity.norm(), 1e-12 * GLS_u_scale);

                  // Store JxW in local variable for faster access;
                  const double JxW = fe_values.JxW(q);

                  // Calculation of the GLS stabilization parameter. The
                  // stabilization parameter used is different if the simulation
                  // is steady or unsteady. In the unsteady case it includes the
                  // value of the time-step
                  const double tau =
                          is_steady(scheme) ?
                          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                                         9 * std::pow(4 * viscosity / (h * h), 2)) :
                          1. /
                          std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                                    9 * std::pow(4 * viscosity / (h * h), 2));

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
                              this->fe->system_to_component_index(i).first;
                      force[i] = rhs_force[q](component_i);
                  }
                  // Correct force to include the dynamic forcing term for flow
                  // control
                  force = force + beta_force;

                  // Calculate the strong residual for GLS stabilization
                  auto strong_residual =
                          velocity_gradient * velocity + present_pressure_gradients[q] -
                          viscosity * present_velocity_laplacians[q] - force;

                  if (velocity_source ==
                      Parameters::VelocitySource::VelocitySourceType::srf)
                  {
                      if (dim == 2)
                      {
                          strong_residual +=
                                  2 * omega_z * (-1.) * cross_product_2d(velocity);
                          auto centrifugal =
                                  omega_z * (-1.) *
                                  cross_product_2d(
                                          omega_z * (-1.) *
                                          cross_product_2d(quadrature_points[q]));
                          strong_residual += centrifugal;
                      }
                      else // dim == 3
                      {
                          strong_residual +=
                                  2 * cross_product_3d(omega_vector, velocity);
                          strong_residual += cross_product_3d(
                                  omega_vector,
                                  cross_product_3d(omega_vector, quadrature_points[q]));
                      }
                  }

                  /* Adjust the strong residual in cases where the scheme is
                   transient.
                   The BDF schemes require values at previous time steps which
                   are stored in the p1, p2 and p3 vectors. The SDIRK scheme
                   require the values at the different stages, which are also
                   stored in the same arrays.
                   */

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
                      scheme == Parameters::SimulationControl::TimeSteppingMethod::
                      steady_bdf)
                      strong_residual += bdf_coefs[0] * velocity +
                                         bdf_coefs[1] * p1_velocity_values[q];

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                      strong_residual += bdf_coefs[0] * velocity +
                                         bdf_coefs[1] * p1_velocity +
                                         bdf_coefs[2] * p2_velocity;

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                      strong_residual +=
                              bdf_coefs[0] * velocity + bdf_coefs[1] * p1_velocity +
                              bdf_coefs[2] * p2_velocity + bdf_coefs[3] * p3_velocity;


                  if (is_sdirk_step1(scheme))
                      strong_residual += sdirk_coefs[0][0] * velocity +
                                         sdirk_coefs[0][1] * p1_velocity;

                  if (is_sdirk_step2(scheme))
                  {
                      strong_residual += sdirk_coefs[1][0] * velocity +
                                         sdirk_coefs[1][1] * p1_velocity +
                                         sdirk_coefs[1][2] * p2_velocity;
                  }

                  if (is_sdirk_step3(scheme))
                  {
                      strong_residual += sdirk_coefs[2][0] * velocity +
                                         sdirk_coefs[2][1] * p1_velocity +
                                         sdirk_coefs[2][2] * p2_velocity +
                                         sdirk_coefs[2][3] * p3_velocity;
                  }

                  // Matrix assembly
                  if (assemble_matrix)
                  {
                      // We loop over the column first to prevent recalculation
                      // of the strong jacobian in the inner loop
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      {
                          const auto phi_u_j      = phi_u[j];
                          const auto grad_phi_u_j = grad_phi_u[j];
                          const auto phi_p_j      = phi_p[j];
                          const auto grad_phi_p_j = grad_phi_p[j];



                          auto strong_jac =
                                  (velocity_gradient * phi_u_j + grad_phi_u_j * velocity +
                                   grad_phi_p_j - viscosity * laplacian_phi_u[j]);

                          if (is_bdf(scheme))
                              strong_jac += phi_u_j * bdf_coefs[0];
                          if (is_sdirk(scheme))
                              strong_jac += phi_u_j * sdirk_coefs[0][0];

                          if (velocity_source ==
                              Parameters::VelocitySource::VelocitySourceType::srf)
                          {
                              if (dim == 2)
                                  strong_jac +=
                                          2 * omega_z * (-1.) * cross_product_2d(phi_u_j);
                              else if (dim == 3)
                                  strong_jac +=
                                          2 * cross_product_3d(omega_vector, phi_u_j);
                          }

                          for (unsigned int i = 0; i < dofs_per_cell; ++i)
                          {
                              const auto phi_u_i      = phi_u[i];
                              const auto grad_phi_u_i = grad_phi_u[i];
                              const auto phi_p_i      = phi_p[i];
                              const auto grad_phi_p_i = grad_phi_p[i];


                              local_matrix(i, j) +=
                                      (
                                              // Momentum terms
                                              viscosity *
                                              scalar_product(grad_phi_u_j, grad_phi_u_i) +
                                              velocity_gradient * phi_u_j * phi_u_i +
                                              grad_phi_u_j * velocity * phi_u_i -
                                              div_phi_u[i] * phi_p_j +
                                              // Continuity
                                              phi_p_i * div_phi_u[j]) *
                                      JxW;

                              // Mass matrix
                              if (is_bdf(scheme))
                                  local_matrix(i, j) +=
                                          phi_u_j * phi_u_i * bdf_coefs[0] * JxW;

                              if (is_sdirk(scheme))
                                  local_matrix(i, j) +=
                                          phi_u_j * phi_u_i * sdirk_coefs[0][0] * JxW;

                              // PSPG GLS term
                              local_matrix(i, j) +=
                                      tau * (strong_jac * grad_phi_p_i) * JxW;

                              if (velocity_source == Parameters::VelocitySource::
                              VelocitySourceType::srf)
                              {
                                  if (dim == 2)
                                      local_matrix(i, j) +=
                                              2 * omega_z * (-1.) *
                                              cross_product_2d(phi_u_j) * phi_u_i * JxW;

                                  else if (dim == 3)
                                      local_matrix(i, j) +=
                                              2 * cross_product_3d(omega_vector, phi_u_j) *
                                              phi_u_i * JxW;
                              }


                              // PSPG TAU term is currently disabled because it
                              // does not alter the matrix sufficiently
                              // local_matrix(i, j) +=
                              //  -tau * tau * tau * 4 / h / h *
                              //  (velocity *phi_u_j) *
                              //  strong_residual * grad_phi_p_i *
                              //  fe_values.JxW(q);

                              // Jacobian is currently incomplete
                              if (SUPG)
                              {
                                  local_matrix(i, j) +=
                                          tau *
                                          (strong_jac * (grad_phi_u_i * velocity) +
                                           strong_residual * (grad_phi_u_i * phi_u_j)) *
                                          JxW;

                                  // SUPG TAU term is currently disabled because
                                  // it does not alter the matrix sufficiently
                                  // local_matrix(i, j)
                                  // +=
                                  //   -strong_residual
                                  //   * (grad_phi_u_i
                                  //   *
                                  //   velocity)
                                  //   * tau * tau *
                                  //   tau * 4 / h / h
                                  //   *
                                  //   (velocity
                                  //   *phi_u_j) *
                                  //   fe_values.JxW(q);
                              }
                          }
                      }
                  }

                  // Assembly of the right-hand side
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                      const auto phi_u_i      = phi_u[i];
                      const auto grad_phi_u_i = grad_phi_u[i];
                      const auto phi_p_i      = phi_p[i];
                      const auto grad_phi_p_i = grad_phi_p[i];
                      const auto div_phi_u_i  = div_phi_u[i];


                      // Navier-Stokes Residual
                      local_rhs(i) +=
                              (
                                      // Momentum
                                      -viscosity *
                                      scalar_product(velocity_gradient, grad_phi_u_i) -
                                      velocity_gradient * velocity * phi_u_i +
                                      current_pressure * div_phi_u_i + force * phi_u_i -
                                      // Continuity
                                      present_velocity_divergence * phi_p_i) *
                              JxW;

                      // Residual associated with BDF schemes
                      if (scheme == Parameters::SimulationControl::
                      TimeSteppingMethod::bdf1 ||
                          scheme == Parameters::SimulationControl::
                          TimeSteppingMethod::steady_bdf)
                          local_rhs(i) -=
                                  bdf_coefs[0] * (velocity - p1_velocity) * phi_u_i * JxW;

                      if (scheme ==
                          Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                          local_rhs(i) -= (bdf_coefs[0] * (velocity * phi_u_i) +
                                           bdf_coefs[1] * (p1_velocity * phi_u_i) +
                                           bdf_coefs[2] * (p2_velocity * phi_u_i)) *
                                          JxW;

                      if (scheme ==
                          Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                          local_rhs(i) -= (bdf_coefs[0] * (velocity * phi_u_i) +
                                           bdf_coefs[1] * (p1_velocity * phi_u_i) +
                                           bdf_coefs[2] * (p2_velocity * phi_u_i) +
                                           bdf_coefs[3] * (p3_velocity * phi_u_i)) *
                                          JxW;

                      // Residuals associated with SDIRK schemes
                      if (is_sdirk_step1(scheme))
                          local_rhs(i) -=
                                  (sdirk_coefs[0][0] * (velocity * phi_u_i) +
                                   sdirk_coefs[0][1] * (p1_velocity * phi_u_i)) *
                                  JxW;

                      if (is_sdirk_step2(scheme))
                      {
                          local_rhs(i) -=
                                  (sdirk_coefs[1][0] * (velocity * phi_u_i) +
                                   sdirk_coefs[1][1] * (p1_velocity * phi_u_i) +
                                   sdirk_coefs[1][2] *
                                   (p2_velocity_values[q] * phi_u_i)) *
                                  JxW;
                      }

                      if (is_sdirk_step3(scheme))
                      {
                          local_rhs(i) -=
                                  (sdirk_coefs[2][0] * (velocity * phi_u_i) +
                                   sdirk_coefs[2][1] * (p1_velocity * phi_u_i) +
                                   sdirk_coefs[2][2] * (p2_velocity * phi_u_i) +
                                   sdirk_coefs[2][3] * (p3_velocity * phi_u_i)) *
                                  JxW;
                      }

                      if (velocity_source ==
                          Parameters::VelocitySource::VelocitySourceType::srf)
                      {
                          if (dim == 2)
                          {
                              local_rhs(i) += -2 * omega_z * (-1.) *
                                              cross_product_2d(velocity) * phi_u_i *
                                              JxW;
                              auto centrifugal =
                                      omega_z * (-1.) *
                                      cross_product_2d(
                                              omega_z * (-1.) *
                                              cross_product_2d(quadrature_points[q]));
                              local_rhs(i) += -centrifugal * phi_u_i * JxW;
                          }
                          else if (dim == 3)
                          {
                              local_rhs(i) +=
                                      -2 * cross_product_3d(omega_vector, velocity) *
                                      phi_u_i * JxW;
                              local_rhs(i) +=
                                      -cross_product_3d(
                                              omega_vector,
                                              cross_product_3d(omega_vector,
                                                               quadrature_points[q])) *
                                      phi_u_i * JxW;
                          }
                      }

                      // PSPG GLS term
                      local_rhs(i) += -tau * (strong_residual * grad_phi_p_i) * JxW;

                      // SUPG GLS term
                      if (SUPG)
                      {
                          local_rhs(i) +=
                                  -tau * (strong_residual * (grad_phi_u_i * velocity)) *
                                  JxW;
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
                  constraints_used.distribute_local_to_global(local_matrix,
                                                              local_rhs,
                                                              local_dof_indices,
                                                              this->system_matrix,
                                                              this->system_rhs);
              }
              else
              {
                  constraints_used.distribute_local_to_global(local_rhs,
                                                              local_dof_indices,
                                                              this->system_rhs);
              }
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

  if (this->simulation_parameters.particlesParameters.integrate_motion)
    {
      force_on_ib();
      integrate_particles();
    }
  if (this->simulation_parameters.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::none)
    {
      TimerOutput::Scope t(this->computing_timer, "assemble_system");
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
      TimerOutput::Scope t(this->computing_timer, "assemble_system");
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

  sharp_edge();
}
template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{


  if (this->simulation_parameters.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::none)
    {
        TimerOutput::Scope t(this->computing_timer, "assemble_rhs");
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
      TimerOutput::Scope t(this->computing_timer, "assemble_rhs");
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
  sharp_edge();
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
  vertices_cell_mapping();
  this->simulation_parameters.mesh_adaptation.refinement_fraction = temp_refine;
  this->simulation_parameters.mesh_adaptation.coarsening_fraction = temp_coarse;


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
        this->first_iteration();
      else
        {
          refine_ib();
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          vertices_cell_mapping();
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
