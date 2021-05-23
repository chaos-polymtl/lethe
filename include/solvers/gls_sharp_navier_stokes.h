/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#ifndef LETHE_GLSSHARPNS_H
#define LETHE_GLSSHARPNS_H

#include <core/ib_particle.h>
#include <core/ib_stencils.h>
#include <solvers/gls_navier_stokes.h>

using namespace dealii;

/**
 * A solver class for the Navier-Stokes equation using GLS stabilization and
 * Sharp-Edge immersed boundaries
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Lucka Barbeau, Bruno Blais, 2020
 */

template <int dim>
class GLSSharpNavierStokesSolver : public GLSNavierStokesSolver<dim>
{
public:
  GLSSharpNavierStokesSolver(SimulationParameters<dim> &nsparam);
  ~GLSSharpNavierStokesSolver();

  void
  solve();

private:
  template <bool                                              assemble_matrix,
            Parameters::SimulationControl::TimeSteppingMethod scheme,
            Parameters::VelocitySource::VelocitySourceType    velocity_source>
  void
  assembleGLS();


  // BB - TODO This explanation needs to be made clearer. Adjacent, Adjacent_2
  // and Adjacent_3 needs to be renamed if possible to a clearer notation

  /**
   */

    /**
    * @brief
    * Map the vertex index to the cell that includes that vertex.
    * This map is used to find all the cell close to a specific vertex.
    */
    void
    vertices_cell_mapping();

    /**
    * @brief
    * Defines the particle structure and value based on the parameter file.
    * This structure gives access to the position, velocity, force and other proprieties of each IB particle.
    * All the variables defined for each of the particle are described in the class: IBparticle.
    * see: \include\core\ib_particle.h
    */
    void
    define_particles();

    /**
    * @brief
    * Evaluate the forces applied on each of the IB particle.
    */
    void
    force_on_ib();

    /**
    * @brief
    * Modify the system matrix to impose IB condition using the sharp_edge approach.
    * The detail of this approach are presented in :
     * L. Barbeau, S. Étienne, C. Béguin & B. Blais,
     * «Development of a high-order continuous Galerkin sharp-interface immersed boundary method and
     * its application to incompressible flow problems,» Computers & Fluids, 2020,
     * in press, ref. CAF-D-20-00773
    */
    void
    sharp_edge();

    /**
    * @brief
    * Write in a specifique file for each of the paticles its forces, velocity, position at each time step.
    * LB - TODO refactor the output format of the file
    */
    void
    write_force_ib();

    /**
    * @brief
    * Integrate the particle velocity and position based on the forces and torques
    * and applies the next value to the particle.
    */
    void
    integrate_particles();

    /**
    * @brief
    * Store the solution of the particles dynamics parameters for integration.
    * Defines the table to store the history of each of the particles.
    */
    void
    finish_time_step_particles();

    /**
    * @brief
    * Evaluate the L2 error on the computational domain if an analytical solution is given.
    * This function is slightly different from its standard GLS counterpart as the cells that are cut by an immersed boundary
    * should not be taken into account in the error evaluation. See "computation domain" definition in:
    * L. Barbeau, S. Étienne, C. Béguin & B. Blais,
    * «Development of a high-order continuous Galerkin sharp-interface immersed boundary method and
    * its application to incompressible flow problems,» Computers & Fluids, 2020,
    * in press, ref. CAF-D-20-0077
    */
    double
    calculate_L2_error_particles();




    /**
    * @brief
    * Same function as its standard GLS counterpart but it used the error evaluation that takes into account the particle’s position.
    */
    virtual void
    postprocess_fd(bool firstIter) override;

    /**
    * @brief
    * Allow a refinement around each of the particles.
    * The zone where the cells will be refined is defined by a ring in 2D and a shell in 3D.
    *The outside and inside radius of the ring\shell is defined in relation to the diameter of the particle by the immersed boundaries parameter:
    * "refine mesh inside radius factor" and "refine mesh outside radius factor".
    * These factors multiply the radius of the particle to define the outside and inside radius of the ring\shell.
    */
    void
    refine_ib();

    /**
    * @brief
    * Return a bool to define if a cell is cut by an IB particle and the local DOFs of the cell for later us.
    * If the cell is cut, the function will return the id of the particle that cut it, else it returns 0.
     *
     * @param cell , the cell that we verify whether it is cut or not.
     *
     * @param local_dof_indices, a container for the local dof indices used during the function call.
     *
     * @param support_points, a mapping of support points for the DOFs.
    *
    */
    std::tuple<bool,unsigned int, std::vector<types::global_dof_index>>
    cell_cut(const typename DoFHandler<dim>::active_cell_iterator &cell,std::vector<types::global_dof_index> &local_dof_indices, std::map<types::global_dof_index, Point<dim>> &support_points);

    /**
    * @brief
    * Return the cell around a point based on a initial guess of a closed cell (look in the neighbors of this cell)
     *
     * @param cell , The initial cell. We suspect the point of being in one of the neighbours of this cell.
     *
     * @param point, The point that we want to find the cell that contains it
    */
    typename DoFHandler<dim>::active_cell_iterator
    find_cell_around_point_with_neighbors(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       Point<dim>             point);

    /**
    * @brief
    Return a bool that describes  if a cell contains a specific point
    *
    * @param cell , The initial cell. We want to check if the point is inside.
    *
    * @param point, The point that we check
    */
    bool
    point_inside_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,Point<dim>  point);

    /**
    * @brief
    *Return a vector of cells around a cell including vertex neighbors
     *
     * @param cell , The initial cell. we want to know all the cells that share a vertex with this cell.
    */
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    find_cells_around_cell(const typename DoFHandler<dim>::active_cell_iterator &cell);

    /**
    * @brief
    *Clear all the line of dof even if the dof is not owned but it is ghost
     *
     * @param cell , The initial cell from which we access the DOF,
     *
     * @param dof_index , The dof index for which we want to clear the line in the matrix.
     *
    */
    void
    clear_line_in_matrix(const typename DoFHandler<dim>::active_cell_iterator &cell, unsigned int dof_index);


    // Modified version of assemble_matrix_and_rhs to include the presence of
    // extra steps. For more detail see the same function in the gls_navier_stokes.h solver.
    void
    assemble_matrix_and_rhs(
      const Parameters::SimulationControl::TimeSteppingMethod
        time_stepping_method) override;

    void
    assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                   time_stepping_method) override;

    /**
     * Members
     */
private:
  std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>>
                               vertices_to_cell;
  const bool                   SUPG        = true;
  const bool                   PSPG        = true;
  const double                 GLS_u_scale = 1;
  std::vector<IBParticle<dim>> particles;


  std::vector<TableHandler> table_f;
  std::vector<TableHandler> table_t;
};


#endif
