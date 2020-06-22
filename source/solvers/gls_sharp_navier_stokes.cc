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
#include "core/manifolds.h"
#include "core/sdirk.h"
#include "core/time_integration_utilities.h"

// Constructor for class GLSNavierStokesSolver
template <int dim>
GLSSharpNavierStokesSolver<dim>::GLSSharpNavierStokesSolver(
  NavierStokesSolverParameters<dim> &p_nsparam,
  const unsigned int                 p_degreeVelocity,
  const unsigned int                 p_degreePressure)
  : NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>(
      p_nsparam,
      p_degreeVelocity,
      p_degreePressure)
{}

template <int dim>
GLSSharpNavierStokesSolver<dim>::~GLSSharpNavierStokesSolver()
{
  this->dof_handler.clear();
}



template <int dim>
void
GLSSharpNavierStokesSolver<dim>::set_solution_vector(double value)
{
  this->present_solution = value;
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::setup_dofs()
{
  TimerOutput::Scope t(this->computing_timer, "setup_dofs");

  // Clear the preconditioner before the matrix they are associated with is
  // cleared
  amg_preconditioner.reset();
  ilu_preconditioner.reset();

  // Now reset system matrix
  system_matrix.clear();

  this->dof_handler.distribute_dofs(this->fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          this->locally_relevant_dofs);

  const MappingQ<dim>        mapping(this->degreeVelocity_,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValuesExtractors::Vector velocities(0);

  // Non-zero constraints
  {
    this->nonzero_constraints.clear();

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->nonzero_constraints);
    for (unsigned int i_bc = 0; i_bc < this->nsparam.boundary_conditions.size;
         ++i_bc)
      {
        if (this->nsparam.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::noslip)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              ZeroFunction<dim>(dim + 1),
              this->nonzero_constraints,
              this->fe.component_mask(velocities));
          }
        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->nsparam.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              this->nonzero_constraints);
          }
        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::function)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              NavierStokesFunctionDefined<dim>(
                &this->nsparam.boundary_conditions.bcFunctions[i_bc].u,
                &this->nsparam.boundary_conditions.bcFunctions[i_bc].v,
                &this->nsparam.boundary_conditions.bcFunctions[i_bc].w),
              this->nonzero_constraints,
              this->fe.component_mask(velocities));
          }

        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints<DoFHandler<dim>>(
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              this->nsparam.boundary_conditions.periodic_id[i_bc],
              this->nsparam.boundary_conditions.periodic_direction[i_bc],
              this->nonzero_constraints);
          }
      }
  }
  this->nonzero_constraints.close();

  {
    this->zero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->zero_constraints);

    for (unsigned int i_bc = 0; i_bc < this->nsparam.boundary_conditions.size;
         ++i_bc)
      {
        if (this->nsparam.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->nsparam.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              this->zero_constraints);
          }
        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints<DoFHandler<dim>>(
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              this->nsparam.boundary_conditions.periodic_id[i_bc],
              this->nsparam.boundary_conditions.periodic_direction[i_bc],
              this->zero_constraints);
          }
        else // if(nsparam.boundaryConditions.boundaries[i_bc].type==Parameters::noslip
          // || Parameters::function)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              ZeroFunction<dim>(dim + 1),
              this->zero_constraints,
              this->fe.component_mask(velocities));
          }
      }
  }
  this->zero_constraints.close();

  this->present_solution.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                this->mpi_communicator);
  this->solution_m1.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs,
                           this->mpi_communicator);
  this->solution_m2.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs,
                           this->mpi_communicator);
  this->solution_m3.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs,
                           this->mpi_communicator);

  this->newton_update.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->system_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->local_evaluation_point.reinit(this->locally_owned_dofs,
                                      this->mpi_communicator);

  DynamicSparsityPattern dsp(this->locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  this->nonzero_constraints,
                                  false);
#if !(DEAL_II_VERSION_GTE(9, 2, 0))
  SparsityTools::distribute_sparsity_pattern(
    dsp,
    this->dof_handler.compute_n_locally_owned_dofs_per_processor(),
    this->mpi_communicator,
    this->locally_relevant_dofs);
  system_matrix.reinit(this->locally_owned_dofs,
                       this->locally_owned_dofs,
                       dsp,
                       this->mpi_communicator);
#else
  SparsityTools::distribute_sparsity_pattern(
    dsp,
    this->dof_handler.locally_owned_dofs(),
    this->mpi_communicator,
    this->locally_relevant_dofs);
  system_matrix.reinit(this->locally_owned_dofs,
                       this->locally_owned_dofs,
                       dsp,
                       this->mpi_communicator);

#endif

  double global_volume = GridTools::volume(*this->triangulation);

  this->pcout << "   Number of active cells:       "
              << this->triangulation->n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;
  this->pcout << "   Volume of triangulation:      " << global_volume
              << std::endl;
}

template <int dim>
void GLSSharpNavierStokesSolver<dim>::vertices_cell_mapping()
{

    //map the vertex index to the cell that include that vertex used later in which cell a point falls in
    //vertices_to_cell is a vector of vectof of dof handler active cell iterator each element i of the vector is a vector of all the cell in contact with the vertex i

    vertices_to_cell.clear();
    vertices_to_cell.resize(this->dof_handler.n_dofs()/(dim+1));
    const auto &cell_iterator=this->dof_handler.active_cell_iterators();
    //loop on all the cell and
    for (const auto &cell : cell_iterator) {
        if (cell->is_locally_owned() | cell->is_ghost()) {
            unsigned int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
            for (unsigned int i = 0; i < vertices_per_cell; i++) {
                //add this cell as neighbors for all it's vertex
                unsigned int v_index = cell->vertex_index(i);
                std::vector<typename DoFHandler<dim>::active_cell_iterator> adjacent = vertices_to_cell[v_index];
                //can only add the cell if it's a set and not a vector
                std::set<typename DoFHandler<dim>::active_cell_iterator> adjacent_2(adjacent.begin(), adjacent.end());
                adjacent_2.insert(cell);
                //convert back the set to a vector and add it in the vertices_to_cell;
                std::vector<typename DoFHandler<dim>::active_cell_iterator> adjacent_3(adjacent_2.begin(),
                                                                                       adjacent_2.end());
                vertices_to_cell[v_index] = adjacent_3;
            }
        }
    }
}

template <int dim>
void GLSSharpNavierStokesSolver<dim>::define_particles() {

    //define position and velocity of particles
    particles.resize(this->nsparam.particlesParameters.nb);
    // define position of particles
    //x y z
    if (dim ==2) {
        for (unsigned int i=0 ; i< this->nsparam.particlesParameters.nb;++i) {
            particles[i].resize(3 * dim);
            //x y
            particles[i][0] = this->nsparam.particlesParameters.particles[i][0];
            particles[i][1] = this->nsparam.particlesParameters.particles[i][1];
            //Vx Vy
            particles[i][2] = this->nsparam.particlesParameters.particles[i][3];
            particles[i][3] = this->nsparam.particlesParameters.particles[i][4];
            //omega
            particles[i][4] = this->nsparam.particlesParameters.particles[i][8];
            //radius
            particles[i][5] = this->nsparam.particlesParameters.particles[i][9];;
        }
    }

    if (dim ==3) {
        for (unsigned int i=0 ; i< this->nsparam.particlesParameters.nb;++i) {
            particles[i].resize(3 * dim+1);
            //x y
            particles[i][0] = this->nsparam.particlesParameters.particles[i][0];
            particles[i][1] = this->nsparam.particlesParameters.particles[i][1];
            particles[i][2] = this->nsparam.particlesParameters.particles[i][2];
            //Vx Vy
            particles[i][3] = this->nsparam.particlesParameters.particles[i][3];
            particles[i][4] = this->nsparam.particlesParameters.particles[i][4];
            particles[i][5] = this->nsparam.particlesParameters.particles[i][5];
            //omega
            particles[i][6] = this->nsparam.particlesParameters.particles[i][6];
            particles[i][7] = this->nsparam.particlesParameters.particles[i][7];
            particles[i][8] = this->nsparam.particlesParameters.particles[i][8];
            //radius
            particles[i][9] = this->nsparam.particlesParameters.particles[i][9];;
        }
    }
}

template <int dim>
void GLSSharpNavierStokesSolver<dim>::refine_ib() {

    Point<dim> center_immersed;
    MappingQ1<dim> immersed_map;
    std::map< types::global_dof_index, Point< dim >>  	support_points;
    DoFTools::map_dofs_to_support_points(immersed_map,this->dof_handler,support_points);
    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    //define cell iterator
    const auto &cell_iterator=this->dof_handler.active_cell_iterators();
    for (const auto &cell : cell_iterator) {
        if (cell->is_locally_owned()) {
            cell->get_dof_indices(local_dof_indices);
            for (unsigned int p = 0; p < particles.size(); ++p) {
                unsigned int count_small = 0;
                if (dim == 2) {
                    center_immersed(0) = particles[p][0];
                    center_immersed(1) = particles[p][1];

                } else if (dim == 3) {
                    center_immersed(0) = particles[p][0];
                    center_immersed(1) = particles[p][1];
                    center_immersed(2) = particles[p][2];

                }
                for (unsigned int j = 0; j < local_dof_indices.size(); ++j) {
                    //count the number of dof that ar smaller or larger then the radius of the particles
                    //if all the dof are on one side the cell is not cut by the boundary meaning we dont have to do anything
                    if ((support_points[local_dof_indices[j]] - center_immersed).norm() <=
                        particles[p][particles[p].size() - 1] * this->nsparam.particlesParameters.outside_radius &
                        (support_points[local_dof_indices[j]] - center_immersed).norm() >=
                        particles[p][particles[p].size() - 1] * this->nsparam.particlesParameters.inside_radius) {
                        ++count_small;
                    }
                }
                if (count_small > 0) {
                    cell->set_refine_flag();
                }
            }
        }
    }
}


template <int dim>
void GLSSharpNavierStokesSolver<dim>::force_on_ib() {

    // calculate the torque and force on a immersed boundary
    Tensor<1, dim, double> force_vect;
    double dr = (GridTools::minimal_cell_diameter(*this->triangulation) *
                 GridTools::minimal_cell_diameter(*this->triangulation)) / sqrt(2 *
                                                                                (GridTools::minimal_cell_diameter(
                                                                                        *this->triangulation) *
                                                                                 GridTools::minimal_cell_diameter(
                                                                                         *this->triangulation)));
    //std::cout << "dr " << dr  << std::endl;

    using numbers::PI;

    if (dim==2) {
        // loop on all particule
        for (unsigned int p = 0; p < particles.size(); ++p) {

            // define stuff usefull for the evaluation
            const double center_x = particles[p][0];
            const double center_y = particles[p][1];


            QGauss<dim> q_formula(this->fe.degree + 1);
            FEValues<dim> fe_values(this->fe, q_formula, update_quadrature_points);

            double mu = this->nsparam.physical_properties.viscosity;

            MappingQ1<dim> immersed_map;
            std::vector<types::global_dof_index> local_dof_indices(this->fe.dofs_per_cell);
            std::vector<types::global_dof_index> local_dof_indices_2(this->fe.dofs_per_cell);
            std::vector<types::global_dof_index> local_dof_indices_3(this->fe.dofs_per_cell);
            unsigned int nb_evaluation = this->nsparam.particlesParameters.nb_force_eval;
            double t_torque = 0;

            double fx_v = 0;
            double fy_v = 0;
            double fx_p_2 = 0;
            double fy_p_2 = 0;

            // loop on all the evaluation point
            for (unsigned int i = 0; i < nb_evaluation; ++i) {
                // define the normal to the surface evaluated.
                Tensor<1, dim, double> surf_normal;
                Tensor<1, dim, double> surf_vect;
                double step_ratio=2;
                surf_normal[0]=dr * cos(i * 2 * PI / (nb_evaluation));
                surf_normal[1]=dr * sin(i * 2 * PI / (nb_evaluation));
                surf_vect[0]=-surf_normal[1];
                surf_vect[1]=surf_normal[0];
                // define the reference point for the surface evaluated.
                const Point<dim> eval_point(particles[p][5] * cos(i * 2 * PI / (nb_evaluation)) + center_x,
                                            particles[p][5] * sin(i * 2 * PI / (nb_evaluation)) + center_y);

                // define the point used for the evaluation
                // points normal to the surface surface

                // step in the normal direction of the surface until we find a cell that is not the cell where the immersed boundary is found.
                const auto &cell = GridTools::find_active_cell_around_point(this->dof_handler, eval_point);
                Point<dim> eval_point_2(eval_point[0] + surf_normal[0]*step_ratio,
                                              eval_point[1] + surf_normal[1]*step_ratio);
                cell->get_dof_indices(local_dof_indices_3);

                unsigned int nb_step=0;
                bool cell_found=false;
                try {
                    //define the cell and check if the point is inside of the cell
                    Point<dim, double> p_cell = immersed_map.transform_real_to_unit_cell(
                            cell, eval_point_2);
                    double dist_2 = GeometryInfo<dim>::distance_to_unit_cell(p_cell);

                    if (dist_2 != 0) {
                        //if the point is in this cell then the dist is equal to 0 and we have found our cell
                        cell_found = true;
                    }
                }
                    // may cause error if the point is not in cell
                catch (typename MappingQGeneric<dim>::ExcTransformationFailed) {
                }


                while (cell_found==false){
                    Point<dim> eval_point_2(eval_point[0] + surf_normal[0]*(nb_step+1)*step_ratio,
                                                  eval_point[1] + surf_normal[1]*(nb_step+1)*step_ratio);
                    try {
                        //define the cell and check if the point is inside of the cell
                        const Point<dim, double> p_cell = immersed_map.transform_real_to_unit_cell(
                                cell, eval_point_2);
                        const double dist_2 = GeometryInfo<dim>::distance_to_unit_cell(p_cell);

                        if (dist_2 <= 0.000000001*dr) {
                            //if the point is in this cell then the dist is equal to 0 and we have found our cell
                            cell_found = true;
                        }
                    }
                        // may cause error if the point is not in cell
                    catch (typename MappingQGeneric<dim>::ExcTransformationFailed) {
                    }
                    if (cell_found==false)
                        nb_step+=1;

                }

                Point<dim> eval_point_4(eval_point[0] + surf_normal[0]*(nb_step+1)*step_ratio,
                                        eval_point[1] + surf_normal[1]*(nb_step+1)*step_ratio);
                const Point<dim> eval_point_3(eval_point_4[0] + surf_normal[0]*step_ratio,
                                              eval_point_4[1] + surf_normal[1]*step_ratio);
                //std::cout << "eval point " << i <<" "<< eval_point_4<< " nb_step "<< nb_step << std::endl;
                const auto &cell_2 = GridTools::find_active_cell_around_point(this->dof_handler, eval_point_4);
                const auto &cell_3 = GridTools::find_active_cell_around_point(this->dof_handler, eval_point_3);

                cell_3->get_dof_indices(local_dof_indices_2);
                cell_2->get_dof_indices(local_dof_indices);
                if (cell_2->is_locally_owned()) {

                   // define the tensor used for the velocity evaluation.
                    Tensor<1, dim, double> u_1;
                    Tensor<1, dim, double> u_2;
                    Tensor<1, dim, double> u_3;

                    // define the velocity component of the particle at the boundary on the reference point
                    u_1[0] = -particles[p][4]*particles[p][5]*sin(i * 2 * PI / (nb_evaluation))+particles[p][2];
                    u_1[1] = particles[p][4]*particles[p][5]*cos(i * 2 * PI / (nb_evaluation))+particles[p][3];
                    // projection of the speed of the boundary on the plan of the surface used for evalaution
                    double U1 =(surf_vect[0]*u_1[0]+surf_vect[1]*u_1[1])/surf_vect.norm();


                    // used support function of the cell to define the interpolation of the velocity
                    Point<dim> second_point_v = immersed_map.transform_real_to_unit_cell(cell_2, eval_point_4);
                    Point<dim> third_point_v = immersed_map.transform_real_to_unit_cell(cell_3, eval_point_3);


                    cell_3->get_dof_indices(local_dof_indices_2);
                    unsigned  int j=0;
                    while(j < this->fe.dofs_per_cell){

                        u_2[0] += this->fe.shape_value(j, second_point_v) * this->present_solution(local_dof_indices[j]);
                        u_2[1] += this->fe.shape_value(j + 1, second_point_v) * this->present_solution(local_dof_indices[j + 1]);
                        u_3[0] += this->fe.shape_value(j, third_point_v) * this->present_solution(local_dof_indices_2[j]);
                        u_3[1] += this->fe.shape_value(j + 1, third_point_v) * this->present_solution(local_dof_indices_2[j + 1]);

                        if (j < (dim + 1) *pow(1+this->nsparam.fem_parameters.pressureOrder,dim)) {
                            j = j + dim + 1;
                        } else {
                            j = j + dim;
                        }
                    }

                    double U2 = (surf_vect[0]*u_2[0]+surf_vect[1]*u_2[1])/surf_vect.norm();
                    double U3 = (surf_vect[0]*u_3[0]+surf_vect[1]*u_3[1])/surf_vect.norm();
                    double du_dn_1=(U2/(particles[p][5]+surf_normal.norm()*(nb_step+1)*step_ratio)-U1/particles[p][5])/((nb_step+1)*surf_normal.norm()*step_ratio);
                    double du_dn_2=(U3/(particles[p][5]+surf_normal.norm()*(nb_step+2)*step_ratio)-U2/(particles[p][5]+surf_normal.norm()*(nb_step+1)*step_ratio))/(surf_normal.norm()*step_ratio);
                    //std::cout << "U1 " << U1 <<" U2 "<<U2 <<" U3 "<< U3 << std::endl;



                    double du_dn = du_dn_1 - (du_dn_2-du_dn_1)*(nb_step+1)/((nb_step+1)+1) ;
                    //std::cout << "du/dn " << du_dn  << " du/dr " << du_dr << std::endl;

                    double local_fx_v=(du_dn* particles[p][5]) * mu  * 2 * PI * particles[p][5] / (nb_evaluation - 1) *
                                      sin(i * 2 * PI / (nb_evaluation));
                    double local_fy_v=-(du_dn * particles[p][5]) * mu  * 2 * PI * particles[p][5] / (nb_evaluation - 1) *
                                      cos(i * 2 * PI / (nb_evaluation));

                    fx_v += local_fx_v;
                    fy_v += local_fy_v;
                    t_torque += local_fx_v * sin(i * 2 * PI / (nb_evaluation) )* particles[p][5]-local_fy_v*cos(i * 2 * PI / (nb_evaluation))* particles[p][5];
                }
            }

            double t_torque_ =Utilities::MPI::sum(t_torque, this->mpi_communicator);



            fx_p_2 = 0;
            fy_p_2 = 0;
            for (unsigned int i = 0; i < nb_evaluation; ++i) {
                const Point<dim> eval_point(particles[p][5] * cos(i * 2 * PI / (nb_evaluation)) + center_x,
                                            particles[p][5] * sin(i * 2 * PI / (nb_evaluation)) + center_y);
                const Point<dim> eval_point_2(eval_point[0] + 1 * dr * cos(i * 2 * PI / (nb_evaluation)),
                                              eval_point[1] + 1 * dr * sin(i * 2 * PI / (nb_evaluation)));
                const Point<dim> eval_point_3(eval_point[0] + 2 * dr * cos(i * 2 * PI / (nb_evaluation)),
                                              eval_point[1] + 2 * dr * sin(i * 2 * PI / (nb_evaluation)));
                const Point<dim> eval_point_4(eval_point[0] + 3 * dr * cos(i * 2 * PI / (nb_evaluation)),
                                              eval_point[1] + 3 * dr * sin(i * 2 * PI / (nb_evaluation)));
                const auto &cell = GridTools::find_active_cell_around_point(this->dof_handler, eval_point_2);
                if (cell->is_locally_owned()) {
                    const auto &cell2 = GridTools::find_active_cell_around_point(this->dof_handler, eval_point_3);
                    const auto &cell3 = GridTools::find_active_cell_around_point(this->dof_handler, eval_point_4);

                    Point<dim> second_point_v = immersed_map.transform_real_to_unit_cell(cell, eval_point_2);
                    Point<dim> second_point_v_2 = immersed_map.transform_real_to_unit_cell(cell2, eval_point_3);
                    Point<dim> second_point_v_3 = immersed_map.transform_real_to_unit_cell(cell3, eval_point_4);
                    cell->get_dof_indices(local_dof_indices);

                    cell2->get_dof_indices(local_dof_indices_2);
                    cell3->get_dof_indices(local_dof_indices_3);
                    double P_1 = 0;
                    double P_2 = 0;
                    double P_3 = 0;
                    unsigned int j=dim;
                    while(j<(dim + 1) *pow(1+this->nsparam.fem_parameters.pressureOrder,dim)){
                        P_1 += this->fe.shape_value(j, second_point_v) * this->present_solution(local_dof_indices[j]);
                        P_2 += this->fe.shape_value(j, second_point_v_2) *
                               this->present_solution(local_dof_indices_2[j]);
                        P_3 += this->fe.shape_value(j, second_point_v_3) *
                               this->present_solution(local_dof_indices_3[j]);
                        if (j < (dim + 1) *pow(1+this->nsparam.fem_parameters.pressureOrder,dim)) {
                            j = j + dim + 1;
                        } else {
                            j = j + dim;
                        }
                    }
                    double P2 = P_1 + (P_1 - P_2) + ((P_1 - P_2) - (P_2 - P_3));
                    fx_p_2 += P2 * -cos(i * 2 * PI / (nb_evaluation)) * 2 * PI * particles[p][5] / (nb_evaluation - 1);
                    fy_p_2 += P2 * -sin(i * 2 * PI / (nb_evaluation)) * 2 * PI * particles[p][5] / (nb_evaluation - 1);
                }
            }
            double fx_p_2_ =Utilities::MPI::sum(fx_p_2, this->mpi_communicator);
            double fy_p_2_ =Utilities::MPI::sum(fy_p_2, this->mpi_communicator);
            double fx_v_ =Utilities::MPI::sum(fx_v, this->mpi_communicator);
            double fy_v_ =Utilities::MPI::sum(fy_v, this->mpi_communicator);





            if  (this->this_mpi_process == 0){

                force_vect[0]=fx_p_2_+fx_v_;
                force_vect[1]=fy_p_2_+fy_v_;
                if (this->nsparam.forces_parameters.verbosity == Parameters::Verbosity::verbose &&
                    this->this_mpi_process == 0)
                {
                    std::cout <<"particle : "<< p << " total_torque :" << t_torque_ << std::endl;
                    std::cout <<"particle : "<< p << " total_torque :" << t_torque_ << std::endl;
                    std::cout << "fx_P: " << fx_p_2_ << std::endl;
                    std::cout << "fy_P: " << fy_p_2_ << std::endl;
                    std::cout << "fx_v: " << fx_v_ << std::endl;
                    std::cout << "fy_v: " << fy_v_ << std::endl;
                    table_f.add_value("particle ID", p);
                    table_f.add_value("f_x", fx_p_2_+fx_v_);
                    table_f.add_value("f_y", fy_p_2_+fy_v_);

                    table_f.set_precision("f_x",
                                          this->nsparam.forces_parameters.display_precision);
                    table_f.set_precision("f_y",
                                          this->nsparam.forces_parameters.display_precision);
                    if (dim == 3)
                    {
                        table_f.add_value("f_z", fy_p_2_+fy_v_);
                        table_f.set_precision("f_z",
                                              this->nsparam.forces_parameters.display_precision);
                    }
                    std::cout << "+------------------------------------------+" << std::endl;
                    std::cout << "|  Force  summary                          |" << std::endl;
                    std::cout << "+------------------------------------------+" << std::endl;


                }

            }
        }

        table_f.write_text(std::cout);
    }
}




template <int dim>
void GLSSharpNavierStokesSolver<dim>::finish_time_step_particles() {

    if (this->nsparam.analytical_solution->calculate_error()){
        // Update the time of the exact solution to the actual time
        this->exact_solution->set_time(this->simulationControl->get_current_time());
        const double error = this->calculate_L2_error_particles();
        if (this->nsparam.simulation_control.method ==
            Parameters::SimulationControl::TimeSteppingMethod::steady){
            this->error_table.clear_current_row();
            this->error_table.add_value(
                    "cells", this->triangulation->n_global_active_cells());
            this->error_table.add_value("error_velocity", error);
            this->error_table.add_value("error_pressure",  0);

        }
        else{
            this->error_table.clear_current_row();
            this->error_table.add_value("time",
                                        this->simulationControl->get_current_time());
            this->error_table.add_value("error", error);
        }
        if (this->nsparam.analytical_solution->verbosity ==
            Parameters::Verbosity::verbose){
            this->pcout << "L2 error : " << error << std::endl;
        }
    }
}


template <int dim>
double GLSSharpNavierStokesSolver<dim>::calculate_L2_error_particles() {

    TimerOutput::Scope t(this->computing_timer, "error");

    QGauss<dim>         quadrature_formula(this->degreeQuadrature_ + 1);
    const MappingQ<dim> mapping(this->degreeVelocity_,
                                this->nsparam.fem_parameters.qmapping_all);
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
    std::vector<double> div_phi_u(dofs_per_cell);
    std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);

    Function<dim> *l_exact_solution = this->exact_solution;

    Point<dim> center_immersed;
    MappingQ1<dim> immersed_map;
    std::map< types::global_dof_index, Point< dim >>  	support_points;
    DoFTools::map_dofs_to_support_points(immersed_map,this->dof_handler,support_points);

    double l2errorU = 0.;
    double div = 0.;

    // loop over elements
    typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
            .begin_active(),
            endc = this->dof_handler.end();
    for (; cell != endc; ++cell)
    {

        if (cell->is_locally_owned())
        {
            cell->get_dof_indices(local_dof_indices);
            bool check_error=true;
            for (unsigned int p = 0; p < particles.size(); ++p) {
                unsigned int count_small = 0;
                if (dim == 2) {
                    center_immersed(0) = particles[p][0];
                    center_immersed(1) = particles[p][1];

                } else if (dim == 3) {
                    center_immersed(0) = particles[p][0];
                    center_immersed(1) = particles[p][1];
                    center_immersed(2) = particles[p][2];

                }
                for (unsigned int j = 0; j < local_dof_indices.size(); ++j) {
                    //count the number of dof that ar smaller or larger then the radius of the particles
                    //if all the dof are on one side the cell is not cut by the boundary meaning we dont have to do anything
                    if ((support_points[local_dof_indices[j]] - center_immersed).norm() <=particles[p][particles[p].size() - 1] ) {
                        ++count_small;
                    }
                }
                if(count_small != 0 and count_small != local_dof_indices.size()){
                    check_error=false;
                }
            }

            if(check_error) {

                fe_values.reinit(cell);
                fe_values[velocities].get_function_values(this->present_solution,
                                                          local_velocity_values);
                fe_values[pressure].get_function_values(this->present_solution,
                                                        local_pressure_values);
                fe_values[velocities].get_function_gradients(
                        this->evaluation_point, present_velocity_gradients);


                // Retrieve the effective "connectivity matrix" for this element
                cell->get_dof_indices(local_dof_indices);

                // Get the exact solution at all gauss points
                l_exact_solution->vector_value_list(fe_values.get_quadrature_points(),
                                                    q_exactSol);
                for (unsigned int q = 0; q < n_q_points; q++) {


                    double present_velocity_divergence =
                            trace(present_velocity_gradients[q]);
                    div+=present_velocity_divergence* fe_values.JxW(q);


                    // Find the values of x and u_h (the finite element solution) at
                    // the quadrature points
                    double ux_sim = local_velocity_values[q][0];
                    double ux_exact = q_exactSol[q][0];

                    double uy_sim = local_velocity_values[q][1];
                    double uy_exact = q_exactSol[q][1];
                    l2errorU +=
                            (ux_sim - ux_exact) * (ux_sim - ux_exact) * fe_values.JxW(q);
                    l2errorU +=
                            (uy_sim - uy_exact) * (uy_sim - uy_exact) * fe_values.JxW(q);

                    if (dim == 3) {
                        double uz_sim = local_velocity_values[q][2];
                        double uz_exact = q_exactSol[q][2];
                        l2errorU += (uz_sim - uz_exact) * (uz_sim - uz_exact) *
                                    fe_values.JxW(q);
                    }
                }
            }
        }
    }
    l2errorU = Utilities::MPI::sum(l2errorU, this->mpi_communicator);
    div = Utilities::MPI::sum(div, this->mpi_communicator);

    this->pcout << "div: " << div << std::endl;

    return std::sqrt(l2errorU);

}


template <int dim>
void GLSSharpNavierStokesSolver<dim>::sharp_edge(const bool initial_step) {
    //This function define a immersed boundary base on the sharp edge method on a hyper_shere of dim 2 or 3

    //define stuff  in a later version the center of the hyper_sphere would be defined by a particle handler and the boundary condition associeted with it also.

    using numbers::PI;
    Point<dim> center_immersed;
    Point<dim> pressure_bridge;
    std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;
    std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors_set;
    std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors_2;
    //Vector <int> dof_done;
    //dof_done.reinit(this->dof_handler.n_dofs());
    //define a map to all dof and it's support point
    MappingQ1<dim> immersed_map;
    std::map< types::global_dof_index, Point< dim >>  	support_points;
    DoFTools::map_dofs_to_support_points(immersed_map,this->dof_handler,support_points);

    // initalise fe value object in order to do calculation with it later
    QGauss<dim> q_formula(this->degreeQuadrature_);
    FEValues<dim> fe_values(this->fe, q_formula,update_quadrature_points|update_JxW_values);
    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    unsigned int vertex_per_cell = GeometryInfo<dim>::vertices_per_cell;


    unsigned int n_q_points  = q_formula.size();
    // define multiple local_dof_indices one for the cell iterator one for the cell with the second point for
    // the sharp edge stancil and one for manipulation on the neighbors cell.
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices_2(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices_3(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices_4(dofs_per_cell);
    // define minimal cell length

    double dr = (GridTools::minimal_cell_diameter(*this->triangulation) *
                 GridTools::minimal_cell_diameter(*this->triangulation)) / sqrt(2 *
                                                                                (GridTools::minimal_cell_diameter(
                                                                                        *this->triangulation) *
                                                                                 GridTools::minimal_cell_diameter(
                                                                                         *this->triangulation)));


    //define cell iterator
    const auto &cell_iterator=this->dof_handler.active_cell_iterators();


    //loop on all the cell to define if the sharp edge cut them
    for (const auto &cell : cell_iterator) {
        if (cell->is_locally_owned()) {
            double sum_line=0 ;
            fe_values.reinit(cell);
            cell->get_dof_indices(local_dof_indices);
            std::vector<int> set_pressure_cell;
            set_pressure_cell.resize(particles.size());
            // define the order of magnetude for the stencil.
            for(unsigned int qf =0 ; qf<n_q_points  ; ++qf) {
                sum_line += fe_values.JxW(qf);
            }
            // loop over all particle  to see if one of them is cutting this cell
            for (unsigned int p = 0; p < particles.size(); ++p) {
                unsigned int count_small = 0;
                if (dim == 2) {
                    center_immersed(0) = particles[p][0];
                    center_immersed(1) = particles[p][1];
                    // define arbitrary point on the boundary where the pressure will be link between the 2 domain
                    pressure_bridge(0) = particles[p][0]-this->nsparam.particlesParameters.pressure_offset[p][0];
                    pressure_bridge(1) = particles[p][1]-this->nsparam.particlesParameters.pressure_offset[p][1];
                }
                else if (dim == 3) {
                    center_immersed(0) = particles[p][0];
                    center_immersed(1) = particles[p][1];
                    center_immersed(2) = particles[p][2];
                    // define arbitrary point on the boundary where the pressure will be link between the 2 domain
                    pressure_bridge(0) = particles[p][0]-this->nsparam.particlesParameters.pressure_offset[p][0];
                    pressure_bridge(1) = particles[p][1]-this->nsparam.particlesParameters.pressure_offset[p][1];
                    pressure_bridge(2) = particles[p][2]-this->nsparam.particlesParameters.pressure_offset[p][2];

                }

                for (unsigned int j = 0; j < local_dof_indices.size(); ++j) {
                    //count the number of dof that ar smaller or larger then the radius of the particles
                    //if all the dof are on one side the cell is not cut by the boundary meaning we dont have to do anything
                    if ((support_points[local_dof_indices[j]] - center_immersed).norm() <= particles[p][particles[p].size()-1]) {
                        ++count_small;
                    }
                }

                // impose the pressure inside the particle if the inside of the particle is solved
                if(this->nsparam.particlesParameters.assemble_inside & this->nsparam.particlesParameters.P_assemble==Parameters::Particle_Assemble_type::NS ) {
                    bool cell_found = false;
                    try {
                        //define the cell and check if the point is inside of the cell
                        const Point<dim, double> p_cell = immersed_map.transform_real_to_unit_cell(
                                cell, pressure_bridge);
                        const double dist_2 = GeometryInfo<dim>::distance_to_unit_cell(p_cell);

                        if (dist_2 == 0) {
                            //if the point is in this cell then the dist is equal to 0 and we have found our cell
                            cell_found = true;
                        }
                    }
                        // may cause error if the point is not in cell
                    catch (typename MappingQGeneric<dim>::ExcTransformationFailed) {
                    }

                    if (cell_found) {

                        // cleart the line in the matrix
                        unsigned int inside_index = local_dof_indices[dim];
                        this->system_matrix.clear_row(inside_index);
                        for (unsigned int vi = 0; vi < vertex_per_cell; ++vi) {
                            unsigned int v_index = cell->vertex_index(vi);
                            active_neighbors_set = this->vertices_to_cell[v_index];
                            for (unsigned int m = 0; m < active_neighbors_set.size(); m++) {
                                const auto &cell_3 = active_neighbors_set[m];
                                cell_3->get_dof_indices(local_dof_indices_3);
                                for (unsigned int o = 0; o < local_dof_indices_3.size(); ++o) {
                                    this->system_matrix.set(inside_index, local_dof_indices_3[o],
                                                            0);
                                }
                            }
                        }
                        // set new equation for the first pressure dof of the cell. this is the new reference pressure inside a particle
                        system_matrix.set(inside_index, local_dof_indices[dim], sum_line);

                        if (initial_step)
                            this->system_rhs(inside_index) = 0-this->local_evaluation_point(inside_index)*sum_line;
                        else
                            this->system_rhs(inside_index) = 0;

                    }

                }


                //if the cell is cut by the IB the count wont equal 0 or the number of total dof in a cell

                if (count_small != 0 and count_small != local_dof_indices.size()) {

                    //if we are here the cell is cut by the immersed boundary
                    //loops on the dof that reprensant the velocity  component and pressure separately
                    for (unsigned int k = 0; k < dim + 1; ++k) {
                        if (k < dim) {
                            //we are working on the velocity of th
                            //loops on the dof that are for vx or vy separately
                            unsigned int l = k;
                            // loops on all the dof of the the cell that represent a specific component
                            while (l < local_dof_indices.size()) {
                                    if (true) {
                                        // define which dof is going to be redefine
                                        unsigned int global_index_overrigth = local_dof_indices[l];
                                        //define the distance vector between the immersed boundary and the dof support point for each dof
                                        Tensor<1, dim, double> vect_dist = (support_points[local_dof_indices[l]] -
                                                                            center_immersed -
                                                                            particles[p][particles[p].size() - 1] *
                                                                            (support_points[local_dof_indices[l]] -
                                                                             center_immersed) /
                                                                            (support_points[local_dof_indices[l]] -
                                                                             center_immersed).norm());

                                        //define the other points for the stencil ( IB point, original dof and the other points)
                                        // this goes up to a 5 point stencil.

                                        // define the length ratio that represent the zone used for the stencil. the length is define as the length between the dof and the IB
                                        unsigned int length_ratio = 8;

                                        double length_fraction = 1. / length_ratio;

                                        Point<dim, double> second_point(
                                                support_points[local_dof_indices[l]] + vect_dist * length_fraction);
                                        Point<dim, double> third_point(
                                                support_points[local_dof_indices[l]] +
                                                vect_dist * length_fraction * 1 / 2);

                                        Point<dim, double> fourth_point(
                                                support_points[local_dof_indices[l]] +
                                                vect_dist * length_fraction * 3 / 4);

                                        Point<dim, double> fifth_point(
                                                support_points[local_dof_indices[l]] +
                                                vect_dist * length_fraction * 1 / 4);

                                        double dof_2;
                                        double sp_2;

                                        double dof_3;
                                        double sp_3;
                                        double tp_3;

                                        double dof_5;
                                        double fp2_5;
                                        double tp_5;
                                        double fp1_5;
                                        double sp_5;

                                        // define the stencil coefficient in function of the length ratio. this will be automaticly generated in futur version

                                        if (length_ratio == 4) {
                                            dof_2 = 5;
                                            sp_2 = -4;

                                            dof_3 = 45;
                                            sp_3 = 36;
                                            tp_3 = -80;

                                            dof_5 = 4845;
                                            fp2_5 = -18240;
                                            tp_5 = 25840;
                                            fp1_5 = -16320;
                                            sp_5 = 3876;

                                        } else if (length_ratio == 2) {
                                            dof_2 = 3;
                                            sp_2 = -2;

                                            dof_3 = 15;
                                            sp_3 = 10;
                                            tp_3 = -24;

                                            dof_5 = 495;
                                            fp2_5 = -1760;
                                            tp_5 = 2376;
                                            fp1_5 = -1440;
                                            sp_5 = 330;
                                        }
                                        else if (length_ratio == 8) {
                                            dof_2 = 9;
                                            sp_2 = -8;

                                            dof_3 = 153;
                                            sp_3 = 136;
                                            tp_3 = -288;

                                            dof_5 = 58905;
                                            fp2_5 = -228480;
                                            tp_5 = 332640;
                                            fp1_5 = -215424;
                                            sp_5 = 52360;
                                        }





                                        //define the vertex associated with the dof
                                        unsigned int cell_found = 0;
                                        bool break_bool = false;
                                        for (unsigned int vi = 0; vi < vertex_per_cell; ++vi) {
                                            unsigned int v_index = cell->vertex_index(vi);
                                            //get a cell iterator for all the cell neighbors of that vertex
                                            active_neighbors_set = this->vertices_to_cell[v_index];
                                            unsigned int n_active_cells = active_neighbors_set.size();

                                            //loops on those cell to find in which of them the new point for or sharp edge stencil is
                                            for (unsigned int cell_index = 0;
                                                 cell_index < n_active_cells; ++cell_index) {
                                                try {
                                                    //define the cell and check if the point is inside of the cell
                                                    const Point<dim, double> p_cell = immersed_map.transform_real_to_unit_cell(
                                                            active_neighbors_set[cell_index], second_point);
                                                    const double dist_2 = GeometryInfo<dim>::distance_to_unit_cell(
                                                            p_cell);

                                                    //define the cell and check if the point is inside of the cell
                                                    if (dist_2 == 0) {
                                                        //if the point is in this cell then the dist is equal to 0 and we have found our cell
                                                        cell_found = cell_index;
                                                        break_bool = true;
                                                        active_neighbors = active_neighbors_set;
                                                        break;
                                                    }
                                                }
                                                    // may cause error if the point is not in cell
                                                catch (typename MappingQGeneric<dim>::ExcTransformationFailed) {
                                                }
                                            }

                                        }

                                        auto &cell_2 = active_neighbors[cell_found];
                                        bool skip_stencil = false;


                                        if (break_bool == false) {
                                            std::cout << "cell not found around point " << std::endl;
                                            std::cout << "cell index " << cell_found << std::endl;
                                            cell_2 = GridTools::find_active_cell_around_point(this->dof_handler,
                                                                                              second_point);
                                            cell_2->get_dof_indices(local_dof_indices_2);
                                            std::cout << "dof point  " << support_points[global_index_overrigth]
                                                      << std::endl;
                                            std::cout << "second point  " << second_point << std::endl;
                                        }


                                        //we have or next cell needed to complete the stencil

                                        //define the unit cell points for the points used in the stencil for extrapolation.
                                        Point<dim> second_point_v = immersed_map.transform_real_to_unit_cell(cell_2,
                                                                                                             second_point);
                                        Point<dim> third_point_v = immersed_map.transform_real_to_unit_cell(cell_2,
                                                                                                            third_point);
                                        Point<dim> fourth_point_v = immersed_map.transform_real_to_unit_cell(cell_2,
                                                                                                             fourth_point);
                                        Point<dim> fifth_point_v = immersed_map.transform_real_to_unit_cell(cell_2,
                                                                                                            fifth_point);

                                        cell_2->get_dof_indices(local_dof_indices_2);

                                        //clear the current line of this dof  by looping on the neighbors cell of this dof and clear all the associated dof

                                        this->system_matrix.clear_row(global_index_overrigth);

                                        for (unsigned int vi = 0; vi < vertex_per_cell; ++vi) {
                                            unsigned int v_index = cell->vertex_index(vi);
                                            active_neighbors_set = this->vertices_to_cell[v_index];
                                            for (unsigned int m = 0; m < active_neighbors_set.size(); m++) {
                                                const auto &cell_3 = active_neighbors_set[m];
                                                cell_3->get_dof_indices(local_dof_indices_3);
                                                for (unsigned int o = 0; o < local_dof_indices_3.size(); ++o) {

                                                    this->system_matrix.set(global_index_overrigth, local_dof_indices_3[o],
                                                                            0);
                                                }
                                            }
                                        }


                                        // check if the DOF intersect the IB
                                        bool do_rhs = false;
                                        if (cell_2 == cell) {
                                            skip_stencil = true;
                                            this->system_matrix.set(global_index_overrigth, global_index_overrigth,
                                                                    sum_line);
                                            this->system_rhs(global_index_overrigth) = 0;
                                            // Tolerence to define a intersection of the DOF and IB
                                            if (vect_dist.norm() <= 0.000000000001*dr) {
                                                do_rhs = true;
                                            } else {
                                                this->system_rhs(global_index_overrigth) = 0;
                                            }
                                        }


                                        // define the variable used for the extrapolation of the actual solution at the boundary in order to define the correction
                                        double local_interp_sol = 0;
                                        double local_interp_sol_2 = 0;
                                        double local_interp_sol_3 = 0;
                                        double local_interp_sol_4 = 0;
                                        double last_local_interp_sol = 0;
                                        double last_local_interp_sol_2 = 0;
                                        double last_local_interp_sol_3 = 0;
                                        double last_local_interp_sol_4 = 0;
                                        //define the new matrix entry for this dof
                                        if (skip_stencil == false) {
                                            // first the dof itself
                                            unsigned int n = k;
                                            while (n < local_dof_indices_2.size()) {
                                                // first the dof itself
                                                if (global_index_overrigth == local_dof_indices_2[n]) {
                                                    // define the solution at each point used for the stencil and applied the stencil for the specfic dof.

                                                    if (this->nsparam.particlesParameters.order == 2) {
                                                        this->system_matrix.set(global_index_overrigth,
                                                                                local_dof_indices_2[n], sp_2 *
                                                                                                        this->fe.shape_value(
                                                                                                                n,
                                                                                                                second_point_v) *
                                                                                                        sum_line +
                                                                                                        dof_2 *
                                                                                                        sum_line);
                                                        local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) * sum_line *
                                                                        this->evaluation_point(local_dof_indices_2[n]);
                                                        last_local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) *
                                                                this->solution_m1(local_dof_indices_2[n]);
                                                    }

                                                    if (this->nsparam.particlesParameters.order == 3) {
                                                        this->system_matrix.set(global_index_overrigth,
                                                                                local_dof_indices_2[n], sp_3 *
                                                                                                        this->fe.shape_value(
                                                                                                                n,
                                                                                                                second_point_v) *
                                                                                                        sum_line +
                                                                                                        dof_3 *
                                                                                                        sum_line +
                                                                                                        tp_3 *
                                                                                                        this->fe.shape_value(
                                                                                                                n,
                                                                                                                third_point_v) *
                                                                                                        sum_line);

                                                        local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) * sum_line *
                                                                        this->evaluation_point(local_dof_indices_2[n]);
                                                        local_interp_sol_2 +=
                                                                1 * this->fe.shape_value(n, third_point_v) * sum_line *
                                                                        this->evaluation_point(local_dof_indices_2[n]);
                                                        last_local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) *
                                                                        this->solution_m1(local_dof_indices_2[n]);
                                                        last_local_interp_sol_2 +=
                                                                1 * this->fe.shape_value(n, third_point_v) *
                                                                        this->solution_m1(local_dof_indices_2[n]);
                                                    }
                                                    if (this->nsparam.particlesParameters.order > 3) {
                                                        this->system_matrix.set(global_index_overrigth,
                                                                                local_dof_indices_2[n],
                                                                                dof_5 * sum_line + sp_5 *
                                                                                                   this->fe.shape_value(
                                                                                                           n,
                                                                                                           second_point_v) *
                                                                                                   sum_line + tp_5 *
                                                                                                              this->fe.shape_value(
                                                                                                                      n,
                                                                                                                      third_point_v) *
                                                                                                              sum_line +
                                                                                fp1_5 * this->fe.shape_value(n,
                                                                                                             fourth_point_v) *
                                                                                sum_line + fp2_5 *
                                                                                           this->fe.shape_value(n,
                                                                                                                fifth_point_v) *
                                                                                           sum_line);


                                                        local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) * sum_line *
                                                                this->evaluation_point(local_dof_indices_2[n]);
                                                        local_interp_sol_2 +=
                                                                1 * this->fe.shape_value(n, third_point_v) * sum_line *
                                                                this->evaluation_point(local_dof_indices_2[n]);
                                                        local_interp_sol_3 +=
                                                                1 * this->fe.shape_value(n, fourth_point_v) * sum_line *
                                                                this->evaluation_point(local_dof_indices_2[n]);
                                                        local_interp_sol_4 +=
                                                                1 * this->fe.shape_value(n, fifth_point_v) * sum_line *
                                                                this->evaluation_point(local_dof_indices_2[n]);
                                                        last_local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) *
                                                                this->solution_m1(local_dof_indices_2[n]);
                                                        last_local_interp_sol_2 +=
                                                                1 * this->fe.shape_value(n, third_point_v) *
                                                                this->solution_m1(local_dof_indices_2[n]);
                                                        last_local_interp_sol_3 +=
                                                                1 * this->fe.shape_value(n, fourth_point_v) *
                                                                this->solution_m1(local_dof_indices_2[n]);
                                                        last_local_interp_sol_4 +=
                                                                1 * this->fe.shape_value(n, fifth_point_v) *
                                                                this->solution_m1(local_dof_indices_2[n]);

                                                    }
                                                }
                                                    // then the third point trough interpolation from the dof of the cell in which the third point is
                                                else {
                                                    if (this->nsparam.particlesParameters.order == 2) {
                                                        this->system_matrix.set(global_index_overrigth,
                                                                                local_dof_indices_2[n], sp_2 *
                                                                                                        this->fe.shape_value(
                                                                                                                n,
                                                                                                                second_point_v) *
                                                                                                        sum_line);
                                                        local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) * sum_line *
                                                                        this->evaluation_point(local_dof_indices_2[n]);
                                                        last_local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) *
                                                                        this->solution_m1(local_dof_indices_2[n]);
                                                    }

                                                    if (this->nsparam.particlesParameters.order == 3) {
                                                        this->system_matrix.set(global_index_overrigth,
                                                                                local_dof_indices_2[n], sp_3 *
                                                                                                        this->fe.shape_value(
                                                                                                                n,
                                                                                                                second_point_v) *
                                                                                                        sum_line +
                                                                                                        tp_3 *
                                                                                                        this->fe.shape_value(
                                                                                                                n,
                                                                                                                third_point_v) *
                                                                                                        sum_line);

                                                        local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) * sum_line *
                                                                        this->evaluation_point(local_dof_indices_2[n]);
                                                        local_interp_sol_2 +=
                                                                1 * this->fe.shape_value(n, third_point_v) * sum_line *
                                                                        this->evaluation_point(local_dof_indices_2[n]);
                                                        last_local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) *
                                                                        this->solution_m1(local_dof_indices_2[n]);
                                                        last_local_interp_sol_2 +=
                                                                1 * this->fe.shape_value(n, third_point_v) *
                                                                        this->solution_m1(local_dof_indices_2[n]);
                                                    }
                                                    if (this->nsparam.particlesParameters.order > 3) {
                                                        this->system_matrix.set(global_index_overrigth,
                                                                                local_dof_indices_2[n], sp_5 *
                                                                                                        this->fe.shape_value(
                                                                                                                n,
                                                                                                                second_point_v) *
                                                                                                        sum_line +
                                                                                                        tp_5 *
                                                                                                        this->fe.shape_value(
                                                                                                                n,
                                                                                                                third_point_v) *
                                                                                                        sum_line +
                                                                                                        fp1_5 *
                                                                                                        this->fe.shape_value(
                                                                                                                n,
                                                                                                                fourth_point_v) *
                                                                                                        sum_line +
                                                                                                        fp2_5 *
                                                                                                        this->fe.shape_value(
                                                                                                                n,
                                                                                                                fifth_point_v) *
                                                                                                        sum_line);

                                                        local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) * sum_line *
                                                                this->evaluation_point(local_dof_indices_2[n]);
                                                        local_interp_sol_2 +=
                                                                1 * this->fe.shape_value(n, third_point_v) * sum_line *
                                                                        this->evaluation_point(local_dof_indices_2[n]);
                                                        local_interp_sol_3 +=
                                                                1 * this->fe.shape_value(n, fourth_point_v) * sum_line *
                                                                        this->evaluation_point(local_dof_indices_2[n]);
                                                        local_interp_sol_4 +=
                                                                1 * this->fe.shape_value(n, fifth_point_v) * sum_line *
                                                                        this->evaluation_point(local_dof_indices_2[n]);
                                                        last_local_interp_sol +=
                                                                1 * this->fe.shape_value(n, second_point_v) *
                                                                        this->solution_m1(local_dof_indices_2[n]);
                                                        last_local_interp_sol_2 +=
                                                                1 * this->fe.shape_value(n, third_point_v) *
                                                                        this->solution_m1(local_dof_indices_2[n]);
                                                        last_local_interp_sol_3 +=
                                                                1 * this->fe.shape_value(n, fourth_point_v) *
                                                                        this->solution_m1(local_dof_indices_2[n]);
                                                        last_local_interp_sol_4 +=
                                                                1 * this->fe.shape_value(n, fifth_point_v) *
                                                                        this->solution_m1(local_dof_indices_2[n]);
                                                    }
                                                }

                                                if (n < (dim + 1) *
                                                        pow(1 + this->nsparam.fem_parameters.pressureOrder, dim)) {
                                                    n = n + dim + 1;
                                                } else {
                                                    n = n + dim;
                                                }
                                            }
                                        }



                                        // define the rhs of the stencil used for the Ib
                                        if (skip_stencil == false or do_rhs) {
                                            // different boundary condition depending if the dof is vx ,vy or vz and if the problem we solve is 2d or 3d.
                                            if (k == 0) {
                                                if (dim == 2) {
                                                    double vx = -particles[p][4] * particles[p][5] *
                                                                ((support_points[local_dof_indices[l]] -
                                                                  center_immersed) /
                                                                 (support_points[local_dof_indices[l]] -
                                                                  center_immersed).norm())[1] + particles[p][2];
                                                    double rhs_add=0;
                                                    double correction=0;
                                                    double last_interp_last_sol=0;

                                                    if (this->nsparam.particlesParameters.order == 2) {
                                                        last_interp_last_sol =
                                                                this->solution_m1(global_index_overrigth) * dof_2 +
                                                                last_local_interp_sol * sp_2;
                                                        rhs_add = -this->evaluation_point(global_index_overrigth) * sum_line *
                                                                  dof_2 - local_interp_sol * sp_2;
                                                    }
                                                    if (this->nsparam.particlesParameters.order == 3) {
                                                        last_interp_last_sol =
                                                                this->solution_m1(global_index_overrigth) * dof_3 +
                                                                last_local_interp_sol * sp_3 +
                                                                last_local_interp_sol_2 * tp_3;
                                                        rhs_add = -this->evaluation_point(global_index_overrigth) * sum_line *
                                                                  dof_3 - local_interp_sol * sp_3 -
                                                                  local_interp_sol_2 * tp_3;
                                                    }
                                                    if (this->nsparam.particlesParameters.order > 3) {
                                                        last_interp_last_sol =
                                                                this->solution_m1(global_index_overrigth) * dof_5 +
                                                                last_local_interp_sol * sp_5 +
                                                                last_local_interp_sol_2 * tp_5 +
                                                                last_local_interp_sol_3 * fp1_5 +
                                                                last_local_interp_sol_4 * fp2_5;
                                                        rhs_add = -this->evaluation_point(global_index_overrigth) * sum_line *
                                                                  dof_5 - local_interp_sol * sp_5 -
                                                                  local_interp_sol_2 * tp_5 -
                                                                  local_interp_sol_3 * fp1_5 -
                                                                  local_interp_sol_4 * fp2_5;
                                                    }
                                                    // the correction is a value that modified the IB stencil when the patricule is in movement  but is always going to be 0 in this version of the code.

                                                    correction = (last_interp_last_sol - vx) * (1 - particles[p][2] *
                                                                                                    this->nsparam.simulation_control.dt /
                                                                                                    dr) *
                                                                 abs(particles[p][2] *
                                                                     ((support_points[local_dof_indices[l]] -
                                                                       center_immersed) /
                                                                      (support_points[local_dof_indices[l]] -
                                                                       center_immersed).norm())[0] + particles[p][3] *
                                                                                                     ((support_points[local_dof_indices[l]] -
                                                                                                       center_immersed) /
                                                                                                      (support_points[local_dof_indices[l]] -
                                                                                                       center_immersed).norm())[1]);
                                                    correction=0;
                                                    if (!this->simulationControl->is_at_start())
                                                        vx = vx + correction;
                                                    this->system_rhs(global_index_overrigth) = vx * sum_line + rhs_add;
                                                    if (do_rhs)
                                                        this->system_rhs(global_index_overrigth) = vx * sum_line -
                                                                                                   this->evaluation_point(global_index_overrigth) *
                                                                                                   sum_line;
                                                }
                                                if (dim == 3) {
                                                    double vx = particles[p][2];
                                                    if (this->nsparam.particlesParameters.order == 2)
                                                        this->system_rhs(global_index_overrigth) = vx * sum_line -
                                                                                                   this->evaluation_point(global_index_overrigth) *
                                                                                                   sum_line * dof_2 -
                                                                                                   local_interp_sol *
                                                                                                   sp_2;
                                                    if (this->nsparam.particlesParameters.order == 3)
                                                        this->system_rhs(global_index_overrigth) = vx * sum_line -
                                                                                                   this->evaluation_point(global_index_overrigth) *
                                                                                                   sum_line * dof_3 -
                                                                                                   local_interp_sol *
                                                                                                   sp_3 -
                                                                                                   local_interp_sol_2 *
                                                                                                   tp_3;
                                                    if (this->nsparam.particlesParameters.order > 3)
                                                        this->system_rhs(global_index_overrigth) = vx * sum_line -
                                                                                                   this->evaluation_point(global_index_overrigth) *
                                                                                                   sum_line * dof_5 -
                                                                                                   local_interp_sol *
                                                                                                   sp_5 -
                                                                                                   local_interp_sol_2 *
                                                                                                   tp_5 -
                                                                                                   local_interp_sol_3 *
                                                                                                   fp1_5 -
                                                                                                   local_interp_sol_4 *
                                                                                                   fp2_5;
                                                    if (do_rhs)
                                                        this->system_rhs(global_index_overrigth) = vx * sum_line -
                                                                                                   this->evaluation_point(global_index_overrigth) *
                                                                                                   sum_line;
                                                }
                                            } else if (k == 1) {
                                                if (dim == 2) {
                                                    double vy = particles[p][4] * particles[p][5] *
                                                                ((support_points[local_dof_indices[l]] -
                                                                  center_immersed) /
                                                                 (support_points[local_dof_indices[l]] -
                                                                  center_immersed).norm())[0] + particles[p][3];
                                                    double rhs_add = 0;
                                                    double correction = 0;
                                                    double last_interp_last_sol = 0;
                                                    last_interp_last_sol =
                                                            this->solution_m1(global_index_overrigth) * dof_2 +
                                                            last_local_interp_sol * sp_2;

                                                    if (this->nsparam.particlesParameters.order == 2) {
                                                        last_interp_last_sol =
                                                                this->solution_m1(global_index_overrigth) * dof_2 +
                                                                last_local_interp_sol * sp_2;
                                                        rhs_add = -this->evaluation_point(global_index_overrigth)* sum_line *
                                                                  dof_2 - local_interp_sol * sp_2;

                                                    }
                                                    if (this->nsparam.particlesParameters.order == 3) {
                                                        last_interp_last_sol =
                                                                this->solution_m1(global_index_overrigth) * dof_3 +
                                                                last_local_interp_sol * sp_3 +
                                                                last_local_interp_sol_2 * tp_3;
                                                        rhs_add = -this->evaluation_point(global_index_overrigth) * sum_line *
                                                                  dof_3 - local_interp_sol * sp_3 -
                                                                  local_interp_sol_2 * tp_3;
                                                    }
                                                    if (this->nsparam.particlesParameters.order > 3) {
                                                        last_interp_last_sol =
                                                                this->solution_m1(global_index_overrigth) * dof_5 +
                                                                last_local_interp_sol * sp_5 +
                                                                last_local_interp_sol_2 * tp_5 +
                                                                last_local_interp_sol_3 * fp1_5 +
                                                                last_local_interp_sol_4 * fp2_5;
                                                        rhs_add = -this->evaluation_point(global_index_overrigth)* sum_line *
                                                                  dof_5 - local_interp_sol * sp_5 -
                                                                  local_interp_sol_2 * tp_5 -
                                                                  local_interp_sol_3 * fp1_5 -
                                                                  local_interp_sol_4 * fp2_5;
                                                    }
                                                    // the correction is a value that modified the IB stencil when the patricule is in movement  but is always going to be 0 in this version of the code.

                                                    correction = (last_interp_last_sol - vy) * (1 - particles[p][2] *
                                                                                                    this->nsparam.simulation_control.dt /
                                                                                                    dr) *
                                                                 abs(particles[p][2] *
                                                                     ((support_points[local_dof_indices[l]] -
                                                                       center_immersed) /
                                                                      (support_points[local_dof_indices[l]] -
                                                                       center_immersed).norm())[0] + particles[p][3] *
                                                                                                     ((support_points[local_dof_indices[l]] -
                                                                                                       center_immersed) /
                                                                                                      (support_points[local_dof_indices[l]] -
                                                                                                       center_immersed).norm())[1]);
                                                    correction=0;
                                                    if (!this->simulationControl->is_at_start())
                                                        vy = vy + correction;
                                                    this->system_rhs(global_index_overrigth) = vy * sum_line + rhs_add;
                                                    if (do_rhs)
                                                        this->system_rhs(global_index_overrigth) = vy * sum_line -
                                                                                                   this->evaluation_point(global_index_overrigth) *
                                                                                                   sum_line;

                                                }
                                                if (dim == 3) {
                                                    double vy = particles[p][3];
                                                    if (this->nsparam.particlesParameters.order == 2)
                                                        this->system_rhs(global_index_overrigth) = vy * sum_line -
                                                                                                   this->evaluation_point(global_index_overrigth) *
                                                                                                   sum_line * dof_2 -
                                                                                                   local_interp_sol *
                                                                                                   sp_2;
                                                    if (this->nsparam.particlesParameters.order == 3)
                                                        this->system_rhs(global_index_overrigth) = vy * sum_line -
                                                                                                   this->evaluation_point(global_index_overrigth)*
                                                                                                   sum_line * dof_3 -
                                                                                                   local_interp_sol *
                                                                                                   sp_3 -
                                                                                                   local_interp_sol_2 *
                                                                                                   tp_3;
                                                    if (this->nsparam.particlesParameters.order > 3)
                                                        this->system_rhs(global_index_overrigth) = vy * sum_line -
                                                                                                   this->evaluation_point(global_index_overrigth)*
                                                                                                   sum_line * dof_5 -
                                                                                                   local_interp_sol *
                                                                                                   sp_5 -
                                                                                                   local_interp_sol_2 *
                                                                                                   tp_5 -
                                                                                                   local_interp_sol_3 *
                                                                                                   fp1_5 -
                                                                                                   local_interp_sol_4 *
                                                                                                   fp2_5;
                                                    if (do_rhs)
                                                        this->system_rhs(global_index_overrigth) = vy * sum_line -
                                                                                                   this->evaluation_point(global_index_overrigth) *
                                                                                                   sum_line;
                                                }
                                            } else if (k == 2 & dim == 3) {
                                                double vz = particles[p][5];
                                                if (this->nsparam.particlesParameters.order == 2)
                                                    this->system_rhs(global_index_overrigth) = vz * sum_line -
                                                                                               this->evaluation_point(global_index_overrigth)*
                                                                                               sum_line * dof_2 -
                                                                                               local_interp_sol * sp_2;
                                                if (this->nsparam.particlesParameters.order == 3)
                                                    this->system_rhs(global_index_overrigth) = vz * sum_line -
                                                                                               this->evaluation_point(global_index_overrigth) *
                                                                                               sum_line * dof_3 -
                                                                                               local_interp_sol * sp_3 -
                                                                                               local_interp_sol_2 *
                                                                                               tp_3;
                                                if (this->nsparam.particlesParameters.order > 3)
                                                    this->system_rhs(global_index_overrigth) = vz * sum_line -
                                                                                               this->evaluation_point(global_index_overrigth) *
                                                                                               sum_line * dof_5 -
                                                                                               local_interp_sol * sp_5 -
                                                                                               local_interp_sol_2 *
                                                                                               tp_5 -
                                                                                               local_interp_sol_3 *
                                                                                               fp1_5 -
                                                                                               local_interp_sol_4 *
                                                                                               fp2_5;
                                                if (do_rhs)
                                                    this->system_rhs(global_index_overrigth) = vz * sum_line -
                                                                                               this->evaluation_point(global_index_overrigth) *
                                                                                               sum_line;
                                            }
                                        }
                                }

                                if (l < (dim + 1) *pow(1+this->nsparam.fem_parameters.pressureOrder,dim)) {
                                    l = l + dim + 1;
                                } else {
                                    l = l + dim;
                                }
                            }
                        }

                       if (k==dim ){
                           // applied equation on dof that have no equation define for them. those DOF become Dummy dof
                           // This is usefull for high order cell or when a dof is only element of cell that are cuts
                            unsigned int vertex_per_cell = GeometryInfo<dim>::vertices_per_cell;
                            unsigned int l=k;
                            while (l < local_dof_indices.size()) {

                                    bool pressure_impose = true;
                                    for (unsigned int vi = 0; vi < vertex_per_cell; ++vi) {
                                        unsigned int v_index = cell->vertex_index(vi);
                                        active_neighbors_set = this->vertices_to_cell[v_index];
                                        for (unsigned int m = 0; m < active_neighbors_set.size(); m++) {
                                            const auto &cell_3 = active_neighbors_set[m];
                                            cell_3->get_dof_indices(local_dof_indices_3);
                                            for (unsigned int o = 0; o < local_dof_indices_3.size(); ++o) {
                                                if (this->system_matrix.el(local_dof_indices[l], local_dof_indices_3[o]) !=
                                                    0 | this->system_rhs(local_dof_indices[l]) != 0 )
                                                    pressure_impose = false;
                                            }
                                        }
                                    }
                                    if (pressure_impose) {

                                        unsigned int global_index_overrigth = local_dof_indices[l];
                                        this->system_matrix.set(global_index_overrigth, global_index_overrigth,
                                                                sum_line);
                                        this->system_rhs(global_index_overrigth) = 0;
                                    }

                                if (l < (dim + 1) * pow(1 + this->nsparam.fem_parameters.pressureOrder, dim)) {
                                    l = l + dim + 1;
                                } else {
                                    l = l + dim;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    this->system_matrix.compress(VectorOperation::insert);
    this->system_rhs.compress(VectorOperation::insert);
    initial_step_bool=false;
}

template <int dim>
template <bool                                              assemble_matrix,
          Parameters::SimulationControl::TimeSteppingMethod scheme,
          Parameters::VelocitySource::VelocitySourceType    velocity_source>
void
GLSSharpNavierStokesSolver<dim>::assembleGLS()
{

    MPI_Barrier(this->mpi_communicator);
    if (assemble_matrix)
        system_matrix = 0;
    this->system_rhs = 0;
    //erase_inertia();
    double         viscosity_ = this->nsparam.physical_properties.viscosity;
    Function<dim> *l_forcing_function = this->forcing_function;

    QGauss<dim>                      quadrature_formula(this->degreeQuadrature_);
    const MappingQ<dim>              mapping(this->degreeVelocity_,
                                             this->nsparam.fem_parameters.qmapping_all);
    FEValues<dim>                    fe_values(mapping,
                                               this->fe,
                                               quadrature_formula,
                                               update_values | update_quadrature_points |
                                               update_JxW_values | update_gradients |
                                               update_hessians);
    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);
    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> local_rhs(dofs_per_cell);
    std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim + 1));
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);
    std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);
    std::vector<double> present_pressure_values(n_q_points);
    std::vector<Tensor<1, dim>> present_pressure_gradients(n_q_points);
    std::vector<Tensor<1, dim>> present_velocity_laplacians(n_q_points);
    std::vector<Tensor<2, dim>> present_velocity_hess(n_q_points);

    Tensor<1, dim> force;

    std::vector<double> div_phi_u(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
    std::vector<Tensor<3, dim>> hess_phi_u(dofs_per_cell);
    std::vector<Tensor<1, dim>> laplacian_phi_u(dofs_per_cell);
    std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
    std::vector<double> phi_p(dofs_per_cell);
    std::vector<Tensor<1, dim>> grad_phi_p(dofs_per_cell);

    // Values at previous time step for transient schemes
    std::vector<Tensor<1, dim>> p1_velocity_values(n_q_points);
    std::vector<Tensor<1, dim>> p2_velocity_values(n_q_points);
    std::vector<Tensor<1, dim>> p3_velocity_values(n_q_points);

    std::vector<double> time_steps_vector =
            this->simulationControl->get_time_steps_vector();
    // support point
    MappingQ1 <dim> immersed_map;
    std::map<types::global_dof_index, Point<dim >> support_points;
    DoFTools::map_dofs_to_support_points(immersed_map, this->dof_handler, support_points);

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
        bool assemble_bool=true;
        if (cell->is_locally_owned()) {

            cell->get_dof_indices(local_dof_indices);
            for (unsigned int k = 0; k < particles.size(); ++k){
                unsigned int count_small = 0;
                if (dim == 2) {
                    center_immersed(0) = particles[k][0];
                    center_immersed(1) = particles[k][1];
                    // define arbitrary point on the boundary where the pressure will be link between the 2 domain
                    pressure_bridge(0) = particles[k][0]-this->nsparam.particlesParameters.pressure_offset[k][0];
                    pressure_bridge(1) = particles[k][1]-this->nsparam.particlesParameters.pressure_offset[k][1];
                }
                else if (dim == 3) {
                    center_immersed(0) = particles[k][0];
                    center_immersed(1) = particles[k][1];
                    center_immersed(2) = particles[k][2];
                    // define arbitrary point on the boundary where the pressure will be link between the 2 domain
                    pressure_bridge(0) = particles[k][0]-this->nsparam.particlesParameters.pressure_offset[k][0];
                    pressure_bridge(1) = particles[k][1]-this->nsparam.particlesParameters.pressure_offset[k][1];
                    pressure_bridge(2) = particles[k][2]-this->nsparam.particlesParameters.pressure_offset[k][2];
                }

                for (unsigned int j = 0; j < local_dof_indices.size(); ++j) {
                    //count the number of dof that are smaller or larger then the radius of the particles
                    //if all the dof are on one side the cell is not cut by the boundary meaning we dont have to do anything
                    if ((support_points[local_dof_indices[j]] - center_immersed).norm() <= particles[k][particles[k].size()-1]) {
                        ++count_small;
                    }
                }
                if(this->nsparam.particlesParameters.assemble_inside and this->nsparam.particlesParameters.P_assemble==Parameters::Particle_Assemble_type::NS){
                    if (count_small != 0 and count_small!= local_dof_indices.size()){
                        assemble_bool = false;
                        break;
                    }
                }
                else {
                    if (count_small != 0) {
                        assemble_bool = false;
                        break;
                    }
                }
            }


            if (assemble_bool==true ) {



                fe_values.reinit(cell);

                if (dim == 2)
                    h = std::sqrt(4. * cell->measure() / M_PI) / this->degreeVelocity_;
                else if (dim == 3)
                    h =
                            pow(6 * cell->measure() / M_PI, 1. / 3.) / this->degreeVelocity_;

                local_matrix = 0;
                local_rhs = 0;

                // Gather velocity (values, gradient and laplacian)
                fe_values[velocities].get_function_values(this->evaluation_point,
                                                          present_velocity_values);
                fe_values[velocities].get_function_gradients(
                        this->evaluation_point, present_velocity_gradients);
                fe_values[velocities].get_function_laplacians(
                        this->evaluation_point, present_velocity_laplacians);


                // Gather pressure (values, gradient)
                fe_values[pressure].get_function_values(this->evaluation_point,
                                                        present_pressure_values);
                fe_values[pressure].get_function_gradients(
                        this->evaluation_point, present_pressure_gradients);


                // Calculate forcing term if there is a forcing function
                if (l_forcing_function)
                    l_forcing_function->vector_value_list(
                            fe_values.get_quadrature_points(), rhs_force);

                // Gather the previous time steps depending on the number of stages
                // of the time integration scheme
                if (scheme !=
                    Parameters::SimulationControl::TimeSteppingMethod::steady )
                    fe_values[velocities].get_function_values(this->solution_m1,
                                                              p1_velocity_values);

                if (time_stepping_method_has_two_stages(scheme) )
                    fe_values[velocities].get_function_values(this->solution_m2,
                                                              p2_velocity_values);

                if (time_stepping_method_has_three_stages(scheme))
                    fe_values[velocities].get_function_values(this->solution_m3,
                                                              p3_velocity_values);

                // Loop over the quadrature points
                for (unsigned int q = 0; q < n_q_points; ++q) {
                    // Calculation of the magnitude of the velocity for the
                    // stabilization parameter
                    const double u_mag = std::max(present_velocity_values[q].norm(),
                                                  1e-12 * GLS_u_scale);

                    // Calculation of the GLS stabilization parameter. The
                    // stabilization parameter used is different if the simulation is
                    // steady or unsteady. In the unsteady case it includes the value
                    // of the time-step
                    double tau;
                    if (scheme ==
                        Parameters::SimulationControl::TimeSteppingMethod::steady)
                        tau = 1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                                             9 * std::pow(4 * viscosity_ / (h * h), 2));
                    else
                        tau = 1. /
                              std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                                        9 * std::pow(4 * viscosity_ / (h * h), 2));

                    if(PSPG==false)
                        tau=0;
                    // Gather the shape functions, their gradient and their laplacian
                    // for the velocity and the pressure
                    for (unsigned int k = 0; k < dofs_per_cell; ++k) {
                        div_phi_u[k] = fe_values[velocities].divergence(k, q);
                        grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                        phi_u[k] = fe_values[velocities].value(k, q);
                        hess_phi_u[k] = fe_values[velocities].hessian(k, q);
                        phi_p[k] = fe_values[pressure].value(k, q);
                        grad_phi_p[k] = fe_values[pressure].gradient(k, q);

                        for (int d = 0; d < dim; ++d)
                            laplacian_phi_u[k][d] = trace(hess_phi_u[k][d]);
                    }

                    // Establish the force vector
                    for (int i = 0; i < dim; ++i) {
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
                     The BDF schemes require values at previous time steps which are
                     stored in the p1, p2 and p3 vectors. The SDIRK scheme require the
                     values at the different stages, which are also stored in the same
                     arrays.
                     */

                    if (scheme ==
                        Parameters::SimulationControl::TimeSteppingMethod::bdf1 )
                        strong_residual += bdf_coefs[0] * present_velocity_values[q] +
                                           bdf_coefs[1] * p1_velocity_values[q];

                    if (scheme ==
                        Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                        strong_residual += bdf_coefs[0] * present_velocity_values[q] +
                                           bdf_coefs[1] * p1_velocity_values[q] +
                                           bdf_coefs[2] * p2_velocity_values[q];

                    if (scheme ==
                        Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                        strong_residual += bdf_coefs[0] * present_velocity_values[q] +
                                           bdf_coefs[1] * p1_velocity_values[q] +
                                           bdf_coefs[2] * p2_velocity_values[q] +
                                           bdf_coefs[3] * p3_velocity_values[q];


                    if (is_sdirk_step1(scheme))
                        strong_residual +=
                                sdirk_coefs[0][0] * present_velocity_values[q] +
                                sdirk_coefs[0][1] * p1_velocity_values[q];

                    if (is_sdirk_step2(scheme) ) {
                        strong_residual +=
                                sdirk_coefs[1][0] * present_velocity_values[q] +
                                sdirk_coefs[1][1] * p1_velocity_values[q] +
                                sdirk_coefs[1][2] * p2_velocity_values[q];
                    }

                    if (is_sdirk_step3(scheme)) {
                        strong_residual +=
                                sdirk_coefs[2][0] * present_velocity_values[q] +
                                sdirk_coefs[2][1] * p1_velocity_values[q] +
                                sdirk_coefs[2][2] * p2_velocity_values[q] +
                                sdirk_coefs[2][3] * p3_velocity_values[q];
                    }

                    // Matrix assembly
                    if (assemble_matrix) {
                        // We loop over the column first to prevent recalculation of
                        // the strong jacobian in the inner loop
                        for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                            auto strong_jac =
                                    (present_velocity_gradients[q] * phi_u[j] +
                                     grad_phi_u[j] * present_velocity_values[q] +
                                     grad_phi_p[j] - viscosity_ * laplacian_phi_u[j]);

                            if (is_bdf(scheme) )
                                strong_jac += phi_u[j] * bdf_coefs[0];
                            if (is_sdirk(scheme) )
                                strong_jac += phi_u[j] * sdirk_coefs[0][0];

                            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                                local_matrix(i, j) +=
                                        (// Momentum terms
                                                viscosity_ *
                                                scalar_product(grad_phi_u[j], grad_phi_u[i]) +
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

                                if (is_sdirk(scheme) )
                                    local_matrix(i, j) += phi_u[j] * phi_u[i] *
                                                          sdirk_coefs[0][0] *
                                                          fe_values.JxW(q);


                                local_matrix(i, j) +=
                                        tau * strong_jac * grad_phi_p[i] * fe_values.JxW(q);

                                if(PSPG) {
                                    // PSPG TAU term is currently disabled because it does
                                    // not alter the matrix sufficiently
                                    local_matrix(i, j) +=
                                            -tau * tau * tau * 4 / h / h *
                                            (present_velocity_values[q] * phi_u[j]) *
                                            strong_residual * grad_phi_p[i] *
                                            fe_values.JxW(q);
                                }

                                // Jacobian is currently incomplete
                                if (SUPG) {
                                    local_matrix(i, j) +=
                                            tau *
                                            (strong_jac * (grad_phi_u[i] *
                                                           present_velocity_values[q]) +
                                             strong_residual * (grad_phi_u[i] * phi_u[j])) *
                                            fe_values.JxW(q);

                                    // SUPG TAU term is currently disabled because it
                                    // does not alter the matrix sufficiently
                                    local_matrix(i, j)
                                            +=
                                            -strong_residual
                                            * (grad_phi_u[i]
                                               *
                                               present_velocity_values[q])
                                            * tau * tau *
                                            tau * 4 / h / h
                                            *
                                            (present_velocity_values[q]
                                             * phi_u[j]) *
                                            fe_values.JxW(q);
                                }
                            }
                        }
                    }

                    // Assembly of the right-hand side
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                        // Navier-Stokes Residual
                        local_rhs(i) +=
                                (       // Momentum
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
                        if (scheme ==
                            Parameters::SimulationControl::TimeSteppingMethod::bdf1 )
                            local_rhs(i) -=
                                    bdf_coefs[0] *
                                    (present_velocity_values[q] - p1_velocity_values[q]) *
                                    phi_u[i] * fe_values.JxW(q);

                        if (scheme ==
                            Parameters::SimulationControl::TimeSteppingMethod::bdf2 )
                            local_rhs(i) -=
                                    (bdf_coefs[0] * (present_velocity_values[q] * phi_u[i]) +
                                     bdf_coefs[1] * (p1_velocity_values[q] * phi_u[i]) +
                                     bdf_coefs[2] * (p2_velocity_values[q] * phi_u[i])) *
                                    fe_values.JxW(q);

                        if (scheme ==
                            Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                            local_rhs(i) -=
                                    (bdf_coefs[0] * (present_velocity_values[q] * phi_u[i]) +
                                     bdf_coefs[1] * (p1_velocity_values[q] * phi_u[i]) +
                                     bdf_coefs[2] * (p2_velocity_values[q] * phi_u[i]) +
                                     bdf_coefs[3] * (p3_velocity_values[q] * phi_u[i])) *
                                    fe_values.JxW(q);

                        // Residuals associated with SDIRK schemes
                        if (is_sdirk_step1(scheme) )
                            local_rhs(i) -=
                                    (sdirk_coefs[0][0] *
                                     (present_velocity_values[q] * phi_u[i]) +
                                     sdirk_coefs[0][1] * (p1_velocity_values[q] * phi_u[i])) *
                                    fe_values.JxW(q);

                        if (is_sdirk_step2(scheme)) {
                            local_rhs(i) -=
                                    (sdirk_coefs[1][0] *
                                     (present_velocity_values[q] * phi_u[i]) +
                                     sdirk_coefs[1][1] *
                                     (p1_velocity_values[q] * phi_u[i]) +
                                     sdirk_coefs[1][2] *
                                     (p2_velocity_values[q] * phi_u[i])) *
                                    fe_values.JxW(q);
                        }

                        if (is_sdirk_step3(scheme)) {
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
                        if(PSPG)
                            local_rhs(i) +=
                                    -tau * (strong_residual * grad_phi_p[i]) * fe_values.JxW(q);

                        // SUPG GLS term
                        if (SUPG) {
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
                if (assemble_matrix) {
                    constraints_used.distribute_local_to_global(local_matrix,
                                                                local_rhs,
                                                                local_dof_indices,
                                                                system_matrix,
                                                                this->system_rhs);
                } else {
                    constraints_used.distribute_local_to_global(local_rhs,
                                                                local_dof_indices,
                                                                this->system_rhs);
                }
            }
            else if(this->nsparam.particlesParameters.P_assemble==Parameters::Particle_Assemble_type::mass ){
                for (unsigned int q = 0; q < n_q_points; ++q) {
                    if (assemble_matrix) {
                        for (unsigned int i = 0; i< dofs_per_cell; ++i) {
                            for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                                local_matrix(i, j) +=(phi_u[i]*phi_u[j]+phi_p[j] * phi_p[i]) *fe_values.JxW(q);
                            }
                            local_rhs(i) =0;

                        }
                    }
                }

                cell->get_dof_indices(local_dof_indices);

                // The non-linear solver assumes that the nonzero constraints have
                // already been applied to the solution
                const AffineConstraints<double> &constraints_used =
                        this->zero_constraints;
                // initial_step ? nonzero_constraints : zero_constraints;
                if (assemble_matrix) {
                    constraints_used.distribute_local_to_global(local_matrix,
                                                                local_rhs,
                                                                local_dof_indices,
                                                                system_matrix,
                                                                this->system_rhs);
                } else {
                    constraints_used.distribute_local_to_global(local_rhs,
                                                                local_dof_indices,
                                                                this->system_rhs);
                }

            }
            else{
                // could assemble someting in the cells tahat are cut  have to code it here
            }
        }
    }


    if (assemble_matrix)
        system_matrix.compress(VectorOperation::add);
    this->system_rhs.compress(VectorOperation::add);

}

/**
 * Set the initial condition using a L2 or a viscous solver
 **/
template <int dim>
void
GLSSharpNavierStokesSolver<dim>::set_initial_condition(
  Parameters::InitialConditionType initial_condition_type,
  bool                             restart)
{
  if (restart)
    {
      this->pcout << "************************" << std::endl;
      this->pcout << "---> Simulation Restart " << std::endl;
      this->pcout << "************************" << std::endl;
      this->read_checkpoint();
    }
  else if (initial_condition_type ==
           Parameters::InitialConditionType::L2projection)
    {
      assemble_L2_projection();
      solve_system_GMRES(true, 1e-15, 1e-15, true);
      this->present_solution = this->newton_update;
      this->finish_time_step();
      this->postprocess(true);
    }
  else if (initial_condition_type == Parameters::InitialConditionType::nodal)
    {
      this->set_nodal_values();
      this->finish_time_step();
      this->postprocess(true);
    }

  else if (initial_condition_type == Parameters::InitialConditionType::viscous)
    {
      this->set_nodal_values();
      double viscosity = this->nsparam.physical_properties.viscosity;
      this->nsparam.physical_properties.viscosity =
        this->nsparam.initial_condition->viscosity;
      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::steady, false, true);
      this->finish_time_step();
      this->postprocess(true);
      this->nsparam.physical_properties.viscosity = viscosity;
    }
  else
    {
      throw std::runtime_error("GLSNS - Initial condition could not be set");
    }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_L2_projection()
{
  system_matrix    = 0;
  this->system_rhs = 0;
  QGauss<dim>                 quadrature_formula(this->degreeQuadrature_);
  const MappingQ<dim>         mapping(this->degreeVelocity_,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValues<dim>               fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int          dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int          n_q_points    = quadrature_formula.size();
  FullMatrix<double>          local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>              local_rhs(dofs_per_cell);
  std::vector<Vector<double>> initial_velocity(n_q_points,
                                               Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const FEValuesExtractors::Vector     velocities(0);
  const FEValuesExtractors::Scalar     pressure(dim);

  Tensor<1, dim> rhs_initial_velocity_pressure;
  double         rhs_initial_pressure;

  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          local_matrix = 0;
          local_rhs    = 0;
          this->nsparam.initial_condition->uvwp.vector_value_list(
            fe_values.get_quadrature_points(), initial_velocity);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_p[k] = fe_values[pressure].value(k, q);
                  phi_u[k] = fe_values[velocities].value(k, q);
                }

              // Establish the rhs tensor operator
              for (int i = 0; i < dim; ++i)
                {
                  const unsigned int component_i =
                    this->fe.system_to_component_index(i).first;
                  rhs_initial_velocity_pressure[i] =
                    initial_velocity[q](component_i);
                }
              rhs_initial_pressure = initial_velocity[q](dim);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (phi_u[j] * phi_u[i]) * fe_values.JxW(q);
                      local_matrix(i, j) +=
                        (phi_p[j] * phi_p[i]) * fe_values.JxW(q);
                    }
                  local_rhs(i) += (phi_u[i] * rhs_initial_velocity_pressure +
                                   phi_p[i] * rhs_initial_pressure) *
                                  fe_values.JxW(q);
                }
            }

          cell->get_dof_indices(local_dof_indices);
          const AffineConstraints<double> &constraints_used =
            this->nonzero_constraints;
          constraints_used.distribute_local_to_global(local_matrix,
                                                      local_rhs,
                                                      local_dof_indices,
                                                      system_matrix,
                                                      this->system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_system");

  if (this->nsparam.velocitySource.type ==
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
    }

  else if (this->nsparam.velocitySource.type ==
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
    }
    vertices_cell_mapping();
    sharp_edge(initial_step_bool);
}
template <int dim>
void
GLSSharpNavierStokesSolver<dim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_rhs");

  if (this->nsparam.velocitySource.type ==
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
    }
  if (this->nsparam.velocitySource.type ==
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
    }
    vertices_cell_mapping();
    sharp_edge(initial_step_bool);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::solve_linear_system(const bool initial_step,
                                                const bool renewed_matrix)
{
  const double absolute_residual = this->nsparam.linear_solver.minimum_residual;
  const double relative_residual =
    this->nsparam.linear_solver.relative_residual;

  if (this->nsparam.linear_solver.solver ==
      Parameters::LinearSolver::SolverType::gmres)
    solve_system_GMRES(initial_step,
                       absolute_residual,
                       relative_residual,
                       renewed_matrix);
  else if (this->nsparam.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::bicgstab)
    solve_system_BiCGStab(initial_step,
                          absolute_residual,
                          relative_residual,
                          renewed_matrix);
  else if (this->nsparam.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::amg)
    solve_system_AMG(initial_step,
                     absolute_residual,
                     relative_residual,
                     renewed_matrix);
  else if (this->nsparam.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::direct)
    solve_system_direct(initial_step,
                        absolute_residual,
                        relative_residual,
                        renewed_matrix);
  else
    throw(std::runtime_error("This solver is not allowed"));
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::setup_ILU()
{
  TimerOutput::Scope t(this->computing_timer, "setup_ILU");

  const double ilu_fill = this->nsparam.linear_solver.ilu_precond_fill;
  const double ilu_atol = this->nsparam.linear_solver.ilu_precond_atol;
  const double ilu_rtol = this->nsparam.linear_solver.ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix, preconditionerOptions);
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::setup_AMG()
{
  TimerOutput::Scope t(this->computing_timer, "setup_AMG");

  std::vector<std::vector<bool>> constant_modes;
  // Constant modes include pressure since everything is in the same matrix
  std::vector<bool> velocity_components(dim + 1, true);
  velocity_components[dim] = true;
  DoFTools::extract_constant_modes(this->dof_handler,
                                   velocity_components,
                                   constant_modes);

  TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
  amg_data.constant_modes = constant_modes;

  const bool elliptic              = false;
  bool       higher_order_elements = false;
  if (this->degreeVelocity_ > 1)
    higher_order_elements = true;
  const unsigned int n_cycles = this->nsparam.linear_solver.amg_n_cycles;
  const bool         w_cycle  = this->nsparam.linear_solver.amg_w_cycles;
  const double       aggregation_threshold =
    this->nsparam.linear_solver.amg_aggregation_threshold;
  const unsigned int smoother_sweeps =
    this->nsparam.linear_solver.amg_smoother_sweeps;
  const unsigned int smoother_overlap =
    this->nsparam.linear_solver.amg_smoother_overlap;
  const bool                                        output_details = false;
  const char *                                      smoother_type  = "ILU";
  const char *                                      coarse_type    = "ILU";
  TrilinosWrappers::PreconditionAMG::AdditionalData preconditionerOptions(
    elliptic,
    higher_order_elements,
    n_cycles,
    w_cycle,
    aggregation_threshold,
    constant_modes,
    smoother_sweeps,
    smoother_overlap,
    output_details,
    smoother_type,
    coarse_type);

  Teuchos::ParameterList              parameter_ml;
  std::unique_ptr<Epetra_MultiVector> distributed_constant_modes;
  preconditionerOptions.set_parameters(parameter_ml,
                                       distributed_constant_modes,
                                       system_matrix);
  const double ilu_fill = this->nsparam.linear_solver.amg_precond_ilu_fill;
  const double ilu_atol = this->nsparam.linear_solver.amg_precond_ilu_atol;
  const double ilu_rtol = this->nsparam.linear_solver.amg_precond_ilu_rtol;
  parameter_ml.set("smoother: ifpack level-of-fill", ilu_fill);
  parameter_ml.set("smoother: ifpack absolute threshold", ilu_atol);
  parameter_ml.set("smoother: ifpack relative threshold", ilu_rtol);

  parameter_ml.set("coarse: ifpack level-of-fill", ilu_fill);
  parameter_ml.set("coarse: ifpack absolute threshold", ilu_atol);
  parameter_ml.set("coarse: ifpack relative threshold", ilu_rtol);
  amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionAMG>();
  amg_preconditioner->initialize(system_matrix, parameter_ml);
}
template <int dim>
void
GLSSharpNavierStokesSolver<dim>::solve_system_direct(const bool   initial_step,
                                                     const double absolute_residual,
                                                     const double relative_residual,
                                                     const bool   renewed_matrix)
{
    const AffineConstraints<double> &constraints_used =
            initial_step ? this->nonzero_constraints : this->zero_constraints;
    const double linear_solver_tolerance =
            std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

    TrilinosWrappers::MPI::Vector completely_distributed_solution(
            this->locally_owned_dofs, this->mpi_communicator);

    SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                                 linear_solver_tolerance,
                                 true,
                                 true);
    TrilinosWrappers::SolverDirect solver(solver_control);

    if (renewed_matrix || !ilu_preconditioner)
        setup_ILU();
    solver.initialize(system_matrix);
    solver.solve(completely_distributed_solution,this->system_rhs);
    constraints_used.distribute(completely_distributed_solution);
    this->newton_update = completely_distributed_solution;
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::solve_system_GMRES(const bool   initial_step,
                                               const double absolute_residual,
                                               const double relative_residual,
                                               const bool   renewed_matrix)
{
  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linear_solver.residual_precision)
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);
  TrilinosWrappers::SolverGMRES solver(solver_control);

  if (renewed_matrix || !ilu_preconditioner)
    setup_ILU();

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix,
                 completely_distributed_solution,
                 this->system_rhs,
                 *ilu_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }
  }
  constraints_used.distribute(completely_distributed_solution);
  this->newton_update = completely_distributed_solution;
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::solve_system_BiCGStab(
  const bool   initial_step,
  const double absolute_residual,
  const double relative_residual,
  const bool   renewed_matrix)
{
  TimerOutput::Scope t(this->computing_timer, "solve");

  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);
  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linear_solver.residual_precision)
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);
  TrilinosWrappers::SolverBicgstab solver(solver_control);

  if (renewed_matrix || !ilu_preconditioner)
    setup_ILU();

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix,
                 completely_distributed_solution,
                 this->system_rhs,
                 *ilu_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }
    constraints_used.distribute(completely_distributed_solution);
    this->newton_update = completely_distributed_solution;
  }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::solve_system_AMG(const bool   initial_step,
                                             const double absolute_residual,
                                             const double relative_residual,
                                             const bool   renewed_matrix)
{
  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;

  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);
  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linear_solver.residual_precision)
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);
  TrilinosWrappers::SolverGMRES solver(solver_control);

  if (renewed_matrix || !amg_preconditioner)
    setup_AMG();

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix,
                 completely_distributed_solution,
                 this->system_rhs,
                 *amg_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }

    constraints_used.distribute(completely_distributed_solution);

    this->newton_update = completely_distributed_solution;
  }
}

template <int dim>
void
GLSSharpNavierStokesSolver<dim>::solve()
{
  read_mesh_and_manifolds(this->triangulation,
                          this->nsparam.mesh,
                          this->nsparam.manifolds_parameters,
                          this->nsparam.boundary_conditions);

  define_particles();

  for (unsigned int i= 0; i< this->nsparam.particlesParameters.initial_refinement;++i) {
      refine_ib();
      NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
      refine_mesh();
  }

  this->setup_dofs();
  this->set_initial_condition(this->nsparam.initial_condition->type,
                              this->nsparam.restart_parameters.restart);
  initial_step_bool=true;

  while (this->simulationControl->integrate())
    {
      this->simulationControl->print_progression(this->pcout);
      if (this->simulationControl->is_at_start())
        this->first_iteration();
      else
        {
          refine_ib();
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          this->iterate();
        }

      this->postprocess(false);
      this->finish_time_step();

      force_on_ib();
      finish_time_step_particles();

      iter_ib+=1;
      initial_step_bool=true;
      MPI_Barrier(this->mpi_communicator);
    }

  this->finish_simulation();
}


// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class GLSSharpNavierStokesSolver<2>;
template class GLSSharpNavierStokesSolver<3>;
