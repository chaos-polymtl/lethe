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

#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/particles/data_out.h>

#include <core/solutions_output.h>
#include <core/utilities.h>

#include "core/bdf.h"
#include "core/grids.h"
#include "core/manifolds.h"
#include "core/sdirk.h"
#include "core/time_integration_utilities.h"

// Constructor for class GLSNitscheNavierStokesSolver
template <int dim, int spacedim>
GLSNitscheNavierStokesSolver<dim, spacedim>::GLSNitscheNavierStokesSolver(
  SimulationParameters<spacedim> &p_nsparam)
  : GLSNavierStokesSolver<spacedim>(p_nsparam)
  , solid(this->simulation_parameters.nitsche,
          this->triangulation,
          p_nsparam.fem_parameters.velocity_order)
  , fe_ht(1)
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
  double beta = this->simulation_parameters.nitsche->beta;

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
                  auto &evaluation_point = this->evaluation_point;
                  velocity[comp_k] += evaluation_point[fluid_dof_indices[k]] *
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
      auto &system_rhs = this->system_rhs;
      constraints_used.distribute_local_to_global(local_matrix,
                                                  local_rhs,
                                                  fluid_dof_indices,
                                                  this->system_matrix,
                                                  system_rhs);
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
  double              pressure;
  Tensor<1, spacedim> normal_vector;
  Tensor<2, spacedim> fluid_stress;
  Tensor<2, spacedim> fluid_pressure;
  Tensor<1, spacedim> force; // to be changed for a vector of tensors when
                             // allowing multiple solids
  const double viscosity =
    this->simulation_parameters.physical_properties.viscosity;

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

      // Generate FEField functoin to evaluate values and gradients
      // at the particle location
      auto &evaluation_point = this->evaluation_point;
      Functions::FEFieldFunction<spacedim,
                                 DoFHandler<spacedim>,
                                 TrilinosWrappers::MPI::Vector>
        fe_field(this->dof_handler, evaluation_point);

      fe_field.set_active_cell(dh_cell);

      for (const auto &p : pic)
        {
          velocity_gradient = 0;
          pressure          = 0;
          const auto &q     = p.get_location();
          const auto &JxW   = p.get_properties()[0];
          normal_vector[0]  = -p.get_properties()[1];
          normal_vector[1]  = -p.get_properties()[2];
          if (spacedim == 3)
            {
              normal_vector[2] = -p.get_properties()[3];
            }

          for (int k = 0; k < spacedim; ++k)
            {
              velocity_gradient[k] = fe_field.gradient(q, k);
            }

          pressure = fe_field.value(q, 3);

          for (int d = 0; d < dim; ++d)
            {
              fluid_pressure[d][d] = pressure;
            }

          fluid_stress =
            viscosity * (velocity_gradient + transpose(velocity_gradient)) -
            fluid_pressure;
          force += fluid_stress * normal_vector * JxW;
        }

      particle = pic.end();
    }
  force = Utilities::MPI::sum(force, this->mpi_communicator);
  return force;
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::postprocess_solid_forces()
{
  TimerOutput::Scope t(this->computing_timer, "Calculate forces on solid");

  std::vector<Tensor<1, spacedim>> force(
    1, this->calculate_forces_on_solid()); // hard coded, has to be changed for
                                           // when allowing more than 1 solid


  if (this->simulation_parameters.nitsche->verbosity ==
        Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
    {
      std::cout << std::endl;
      const std::vector<unsigned int> solid_indices(
        1,
        1); // hard coded, has to be changed for when allowing more than 1 solid

      std::string independent_column_names = "Solid ID";

      std::vector<std::string> dependent_column_names;
      dependent_column_names.push_back("f_x");
      dependent_column_names.push_back("f_y");
      if (spacedim == 3)
        dependent_column_names.push_back("f_z");

      TableHandler table = make_table_scalars_tensors(
        solid_indices,
        independent_column_names,
        force,
        dependent_column_names,
        this->simulation_parameters.simulation_control.log_precision);

      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Force on solid summary                  |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      table.write_text(std::cout);
    }

  std::string filename =
    this->simulation_parameters.nitsche->force_output_name + ".dat";
  std::ofstream output(filename.c_str());


  solid_forces_table.add_value("time",
                               this->simulation_control->get_current_time());
  solid_forces_table.add_value("f_x", force[0][0]);
  solid_forces_table.add_value("f_y", force[0][1]);
  if (dim == 3)
    solid_forces_table.add_value("f_z", force[0][2]);
  else
    solid_forces_table.add_value("f_z", 0.);

  // Precision
  solid_forces_table.set_precision(
    "f_x", this->simulation_parameters.forces_parameters.output_precision);
  solid_forces_table.set_precision(
    "f_y", this->simulation_parameters.forces_parameters.output_precision);
  solid_forces_table.set_precision(
    "f_z", this->simulation_parameters.forces_parameters.output_precision);
  solid_forces_table.set_precision(
    "time", this->simulation_parameters.forces_parameters.output_precision);

  solid_forces_table.write_text(output);
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::solve()
{
  read_mesh_and_manifolds(
    this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  this->setup_dofs();
  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);
  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      if (this->simulation_parameters.nitsche->enable_particles_motion)
        {
          if (this->simulation_control->is_at_start())
            {
              solid.initial_setup();
              solid.setup_particles();
              std::shared_ptr<Particles::ParticleHandler<spacedim>> solid_ph =
                solid.get_solid_particle_handler();
              output_solid_particles(solid_ph);
              output_solid_triangulation();
            }
          solid.integrate_velocity(this->simulation_control->get_time_step());
          solid.move_solid_triangulation(
            this->simulation_control->get_time_step());
        }
      if (this->simulation_control->is_at_start())
        this->first_iteration();
      else
        {
          this->refine_mesh();
          this->iterate();
        }

      this->postprocess(false);
      // postprocess_ht();
      if (this->simulation_parameters.nitsche->calculate_force_on_solid &&
          dim == 2 && spacedim == 3)
        {
          postprocess_solid_forces();
        }

      if (this->simulation_control->is_output_iteration())
        {
          std::shared_ptr<Particles::ParticleHandler<spacedim>> solid_ph =
            solid.get_solid_particle_handler();
          output_solid_particles(solid_ph);
          output_solid_triangulation();
        }

      this->finish_time_step();
      if (this->simulation_parameters.multiphysics.heat_transfer)
        {
          //          finish_time_step_ht();
        }
    }
  if (this->simulation_parameters.test.enabled)
    solid.print_particle_positions();
  this->finish_simulation();
  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      //      finish_ht();
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::output_solid_particles(
  std::shared_ptr<Particles::ParticleHandler<spacedim>> particle_handler)
{
  Particles::DataOut<spacedim, spacedim> particles_out;
  particles_out.build_patches(*particle_handler);

  const std::string folder = this->simulation_control->get_output_path();
  const std::string solution_name =
    this->simulation_control->get_output_name() + "_solid_particles";
  const unsigned int iter        = this->simulation_control->get_step_number();
  const double       time        = this->simulation_control->get_current_time();
  const unsigned int group_files = this->simulation_control->get_group_files();

  write_vtu_and_pvd<0, spacedim>(pvdhandler_solid_particles,
                                 *(&particles_out),
                                 folder,
                                 solution_name,
                                 time,
                                 iter,
                                 group_files,
                                 this->mpi_communicator);
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::output_solid_triangulation()
{
  DataOut<dim, DoFHandler<dim, spacedim>> data_out;
  DoFHandler<dim, spacedim> &solid_dh = solid.get_solid_dof_handler();
  data_out.attach_dof_handler(solid_dh);

  data_out.build_patches();

  const std::string folder = this->simulation_control->get_output_path();
  const std::string solution_name =
    this->simulation_control->get_output_name() + "_solid_triangulation";
  const unsigned int iter        = this->simulation_control->get_step_number();
  const double       time        = this->simulation_control->get_current_time();
  const unsigned int group_files = this->simulation_control->get_group_files();

  write_vtu_and_pvd<dim>(pvdhandler_solid_triangulation,
                         data_out,
                         folder,
                         solution_name,
                         time,
                         iter,
                         group_files,
                         this->mpi_communicator);
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

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library
// is valid before we actually compile the solver This greatly helps with
// debugging
template class GLSNitscheNavierStokesSolver<2>;
template class GLSNitscheNavierStokesSolver<2, 3>;
template class GLSNitscheNavierStokesSolver<3>;
