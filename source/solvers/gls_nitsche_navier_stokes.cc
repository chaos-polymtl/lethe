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
 * Author: Bruno Blais, Carole-Anne Daunais, Valérie Bibeau, Polytechnique
 Montreal, 2020-
 */

#include "solvers/gls_nitsche_navier_stokes.h"

#include <deal.II/particles/data_out.h>

#include "core/bdf.h"
#include "core/grids.h"
#include "core/manifolds.h"
#include "core/sdirk.h"
#include "core/time_integration_utilities.h"
#include <core/utilities.h>

// Constructor for class GLSNitscheNavierStokesSolver
template <int dim, int spacedim>
GLSNitscheNavierStokesSolver<dim, spacedim>::GLSNitscheNavierStokesSolver(
  NavierStokesSolverParameters<spacedim> &p_nsparam,
  const unsigned int                      p_degreeVelocity,
  const unsigned int                      p_degreePressure)
  : GLSNavierStokesSolver<spacedim>(p_nsparam,
                                    p_degreeVelocity,
                                    p_degreePressure)
  , solid(this->nsparam.nitsche, this->triangulation, p_degreeVelocity)
{}

template <int dim, int spacedim>
template <bool assemble_matrix>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_nitsche_restriction()
{
  std::shared_ptr<Particles::ParticleHandler<spacedim>> solid_ph =
    solid.get_solid_particle_handler();

  TimerOutput::Scope t(this->computing_timer, "Assemble Nitsche terms");

  const unsigned int dofs_per_cell = this->fe.dofs_per_cell;

  std::vector<types::global_dof_index> fluid_dof_indices(dofs_per_cell);

  FullMatrix<double>     local_matrix(dofs_per_cell, dofs_per_cell);
  dealii::Vector<double> local_rhs(dofs_per_cell);

  Tensor<1, spacedim> velocity;
  Function<spacedim> *solid_velocity = solid.get_solid_velocity();

  // Penalization terms
  const auto penalty_parameter =
    1.0 / GridTools::minimal_cell_diameter(*this->triangulation);
  double beta = this->nsparam.nitsche->beta;

  // Loop over all local particles
  auto particle = solid_ph->begin();
  while (particle != solid_ph->end())
    {
      local_matrix = 0;
      local_rhs    = 0;


      const auto &cell = particle->get_surrounding_cell(*this->triangulation);
      const auto &dh_cell =
        typename DoFHandler<spacedim>::cell_iterator(*cell, &this->dof_handler);
      dh_cell->get_dof_indices(fluid_dof_indices);

      const auto pic = solid_ph->particles_in_cell(cell);
      Assert(pic.begin() == particle, ExcInternalError());
      for (const auto &p : pic)
        {
          velocity           = 0;
          const auto &ref_q  = p.get_reference_location();
          const auto &real_q = p.get_location();
          const auto &JxW    = p.get_properties()[0];

          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
              const auto comp_k = this->fe.system_to_component_index(k).first;
              if (comp_k < spacedim)
                {
                  // Get the velocity at non-quadrature point (particle in
                  // fluid)
                  velocity[comp_k] +=
                    this->evaluation_point[fluid_dof_indices[k]] *
                    this->fe.shape_value(k, ref_q);
                }
            }
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const auto comp_i = this->fe.system_to_component_index(i).first;
              if (comp_i < spacedim)
                {
                  if (assemble_matrix)
                    {
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          const auto comp_j =
                            this->fe.system_to_component_index(j).first;
                          if (comp_i == comp_j)
                            local_matrix(i, j) +=
                              penalty_parameter * beta *
                              this->fe.shape_value(i, ref_q) *
                              this->fe.shape_value(j, ref_q) * JxW;
                        }
                    }
                  local_rhs(i) += -penalty_parameter * beta * velocity[comp_i] *
                                    this->fe.shape_value(i, ref_q) * JxW +
                                  penalty_parameter * beta *
                                    solid_velocity->value(real_q, comp_i) *
                                    this->fe.shape_value(i, ref_q) * JxW;
                }
            }
        }
      const AffineConstraints<double> &constraints_used =
        this->zero_constraints;
      constraints_used.distribute_local_to_global(local_matrix,
                                                  local_rhs,
                                                  fluid_dof_indices,
                                                  this->system_matrix,
                                                  this->system_rhs);
      particle = pic.end();
    }
  this->system_matrix.compress(VectorOperation::add);
  this->system_rhs.compress(VectorOperation::add);
}

template <int dim, int spacedim>
Tensor<1, spacedim>
GLSNitscheNavierStokesSolver<dim, spacedim>::calculate_forces_on_solid()
{
  std::shared_ptr<Particles::ParticleHandler<spacedim>> solid_ph =
    solid.get_solid_particle_handler();

  const unsigned int dofs_per_cell = this->fe.dofs_per_cell;

  std::vector<types::global_dof_index> fluid_dof_indices(dofs_per_cell);

  Tensor<2, spacedim> velocity_gradient;
  double pressure;
  Tensor<1, spacedim>   normal_vector;
  Tensor<2, spacedim>   fluid_stress;
  Tensor<2, spacedim>   fluid_pressure;
  Tensor<1, spacedim>   force; // to be changed for a vector of tensors when allowing multiple solids
  const double viscosity = this->nsparam.physical_properties.viscosity;

  // Loop over all local particles
  auto particle = solid_ph->begin();
  while (particle != solid_ph->end())
    {
      const auto &cell = particle->get_surrounding_cell(*this->triangulation);
      const auto &dh_cell =
        typename DoFHandler<spacedim>::cell_iterator(*cell, &this->dof_handler);
      dh_cell->get_dof_indices(fluid_dof_indices);

      const auto pic = solid_ph->particles_in_cell(cell);
      Assert(pic.begin() == particle, ExcInternalError());
      for (const auto &p : pic)
        {
          velocity_gradient            = 0;
          pressure                     = 0;
          const auto &ref_q            = p.get_reference_location();
          const auto &JxW              = p.get_properties()[0];
          normal_vector[0] = -p.get_properties()[1];
          normal_vector[1] = -p.get_properties()[2];
          if (spacedim == 3)
          {
            normal_vector[2] = -p.get_properties()[3];
          }

          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
              const auto comp_k = this->fe.system_to_component_index(k).first;
              if (comp_k < spacedim)
                {
                  // Get the velocity at non-quadrature point (particle in
                  // fluid)
                  velocity_gradient[comp_k] +=
                    double(this->evaluation_point[fluid_dof_indices[k]]) *
                    this->fe.shape_grad(k, ref_q);
                }
              else
                {
                  // Get the pressure at non-quadrature point
                  pressure += 
                    this->evaluation_point[fluid_dof_indices[k]] *
                    this->fe.shape_value(k, ref_q);
                }
            }

          for (int d = 0; d < dim; ++d)
            {
              fluid_pressure[d][d] = pressure;
            }
          fluid_stress = viscosity * velocity_gradient + transpose(velocity_gradient) - fluid_pressure;
          force += fluid_stress * normal_vector * JxW;
        }
      
      particle = pic.end();
    }
    return force;
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::postprocess_solid_forces()
{
  TimerOutput::Scope t(this->computing_timer, "Calculate forces on solid");

  std::vector<Tensor<1, spacedim>>   force(1,  this->calculate_forces_on_solid()); // hard coded, has to be changed for when allowing more than 1 solid


  if (this->nsparam.nitsche->verbosity == Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
      {
        std::cout << std::endl;
        const std::vector<unsigned int> solid_indices(1,1); // hard coded, has to be changed for when allowing more than 1 solid

        std::string independent_column_names = "Solid ID";

        std::vector<std::string> dependent_column_names;
        dependent_column_names.push_back("f_x");
        dependent_column_names.push_back("f_y");
        if (spacedim == 3)
          dependent_column_names.push_back("f_z");

        TableHandler table =
          make_table_scalars_tensors(solid_indices,
                                    independent_column_names,
                                    force,
                                    dependent_column_names,
                                    this->nsparam.forces_parameters.display_precision);

        std::cout << "+------------------------------------------+" << std::endl;
        std::cout << "|  Force on solid summary                  |" << std::endl;
        std::cout << "+------------------------------------------+" << std::endl;
        table.write_text(std::cout);
      }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  this->GLSNavierStokesSolver<spacedim>::assemble_matrix_and_rhs(
    time_stepping_method);

  assemble_nitsche_restriction<true>();
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  this->GLSNavierStokesSolver<spacedim>::assemble_rhs(time_stepping_method);

  assemble_nitsche_restriction<false>();
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::solve()
{
  read_mesh_and_manifolds(this->triangulation,
                          this->nsparam.mesh,
                          this->nsparam.manifolds_parameters,
                          this->nsparam.boundary_conditions);

  this->setup_dofs();
  this->set_initial_condition(this->nsparam.initial_condition->type,
                              this->nsparam.restart_parameters.restart);

  while (this->simulationControl->integrate())
    {
      this->simulationControl->print_progression(this->pcout);
      if (this->nsparam.nitsche->enable_particles_motion)
        solid.integrate_velocity(this->simulationControl->get_time_step());
      if (this->simulationControl->is_at_start())
        this->first_iteration();
      else
        {
          this->refine_mesh();
          this->iterate();
        }
      this->postprocess(false);
      if (this->nsparam.nitsche->calculate_force_on_solid)
      {
        postprocess_solid_forces();
      }
      this->finish_time_step();
    }


  this->finish_simulation();
}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class GLSNitscheNavierStokesSolver<2>;
template class GLSNitscheNavierStokesSolver<2, 3>;
template class GLSNitscheNavierStokesSolver<3>;
