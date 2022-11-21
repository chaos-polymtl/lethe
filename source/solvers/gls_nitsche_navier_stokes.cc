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

#include <core/bdf.h>
#include <core/grids.h>
#include <core/manifolds.h>
#include <core/sdirk.h>
#include <core/solutions_output.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/gls_nitsche_navier_stokes.h>

#include <deal.II/base/multithread_info.h>

#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/particles/data_out.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

// Constructor for class GLSNitscheNavierStokesSolver
template <int dim, int spacedim>
GLSNitscheNavierStokesSolver<dim, spacedim>::GLSNitscheNavierStokesSolver(
  SimulationParameters<spacedim> &p_nsparam)
  : GLSNavierStokesSolver<spacedim>(p_nsparam)
{
  const unsigned int n_solids =
    this->simulation_parameters.nitsche->number_solids;

  for (unsigned int i_solid = 0; i_solid < n_solids; ++i_solid)
    {
      solids.push_back(std::make_shared<SolidBase<dim, spacedim>>(
        this->simulation_parameters.nitsche->nitsche_solids[i_solid],
        this->triangulation,
        this->mapping));
    }

  pvdhandler_solid_particles.resize(n_solids);

  pvdhandler_solid_triangulation.resize(n_solids);

  solid_forces_table.resize(n_solids);
  solid_torques_table.resize(n_solids);
}

template <int dim, int spacedim>
template <bool assemble_matrix>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_nitsche_restriction()
{
  TimerOutput::Scope t(this->computing_timer, "assemble Nitsche restriction");

  Assert(
    !this->simulation_parameters.physical_properties_manager.is_non_newtonian(),
    RequiresConstantViscosity("assemble_nitsche_restriction"));

  // Viscosity for stabilization constant
  const double viscosity =
    this->simulation_parameters.physical_properties_manager.viscosity_scale;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    {
      std::shared_ptr<Particles::ParticleHandler<spacedim>> &solid_ph =
        solids[i_solid]->get_solid_particle_handler();


      const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

      std::vector<types::global_dof_index> fluid_dof_indices(dofs_per_cell);

      FullMatrix<double>     local_matrix(dofs_per_cell, dofs_per_cell);
      dealii::Vector<double> local_rhs(dofs_per_cell);

      Tensor<1, spacedim> velocity;
      Function<spacedim> *solid_velocity =
        solids[i_solid]->get_solid_velocity();

      // Penalization terms
      const double beta =
        this->simulation_parameters.nitsche->nitsche_solids[i_solid]->beta;

      // Loop over all local particles
      auto particle = solid_ph->begin();
      while (particle != solid_ph->end())
        {
          local_matrix = 0;
          local_rhs    = 0;


#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
          const auto &cell =
            particle->get_surrounding_cell(*this->triangulation);
#else
          const auto &cell = particle->get_surrounding_cell();
#endif
          double h_cell = 0;
          if (dim == 2)
            h_cell = std::sqrt(4. * cell->measure() / M_PI) /
                     this->velocity_fem_degree;
          else if (dim == 3)
            h_cell = pow(6 * cell->measure() / M_PI, 1. / 3.) /
                     this->velocity_fem_degree;
          const double penalty_parameter =
            1. / std::pow(h_cell * h_cell, double(dim) / double(spacedim));
          const auto &dh_cell =
            typename DoFHandler<spacedim>::cell_iterator(*cell,
                                                         &this->dof_handler);
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
                  const auto comp_k =
                    this->fe->system_to_component_index(k).first;
                  if (comp_k < spacedim)
                    {
                      // Get the velocity at non-quadrature point (particle in
                      // fluid)
                      auto &evaluation_point = this->evaluation_point;
                      velocity[comp_k] +=
                        evaluation_point[fluid_dof_indices[k]] *
                        this->fe->shape_value(k, ref_q);
                    }
                }

              const double u_mag = velocity.norm();


              const double tau =
                this->simulation_control->get_assembly_method() ==
                    Parameters::SimulationControl::TimeSteppingMethod::steady ?
                  1. / std::sqrt(
                         std::pow(2. * u_mag / h_cell, 2) +
                         9 * std::pow(4 * viscosity / (h_cell * h_cell), 2)) :
                  1. / std::sqrt(
                         std::pow(sdt, 2) + std::pow(2. * u_mag / h_cell, 2) +
                         9 * std::pow(4 * viscosity / (h_cell * h_cell), 2));

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const auto comp_i =
                    this->fe->system_to_component_index(i).first;
                  if (comp_i < spacedim)
                    {
                      if (assemble_matrix)
                        {
                          for (unsigned int j = 0; j < dofs_per_cell; ++j)
                            {
                              const auto comp_j =
                                this->fe->system_to_component_index(j).first;
                              if (comp_i == comp_j)
                                {
                                  local_matrix(i, j) +=
                                    penalty_parameter * beta *
                                    this->fe->shape_value(i, ref_q) *
                                    this->fe->shape_value(j, ref_q) * JxW;

                                  // Partial jacobian of the strong residual
                                  // term so that Nitsche's IB is taken into
                                  // account in the stabilization
                                  local_matrix(i, j) +=
                                    penalty_parameter * beta *
                                    this->fe->shape_value(j, ref_q) *
                                    (this->fe->shape_grad(i, ref_q) *
                                     velocity) *
                                    tau * JxW;
                                }
                            }
                        }
                      // Regular residual
                      local_rhs(i) += -penalty_parameter * beta *
                                        velocity[comp_i] *
                                        this->fe->shape_value(i, ref_q) * JxW +
                                      penalty_parameter * beta *
                                        solid_velocity->value(real_q, comp_i) *
                                        this->fe->shape_value(i, ref_q) * JxW;

                      // Residual due taking into account Nitsche IB in the
                      // strong residual used for GLS stabilization
                      const double strong_residual =
                        -penalty_parameter * beta *
                        (solid_velocity->value(real_q, comp_i) -
                         velocity[comp_i]);

                      const double test_function =
                        this->fe->shape_grad(i, ref_q) * velocity;

                      local_rhs(i) +=
                        -tau * (strong_residual)*test_function * JxW;
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
}

template <>
Tensor<1, 3>
GLSNitscheNavierStokesSolver<2, 3>::calculate_forces_on_solid(
  const unsigned int i_solid)
{
  std::shared_ptr<Particles::ParticleHandler<3>> &solid_ph =
    solids[i_solid]->get_solid_particle_handler();

  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

  std::vector<types::global_dof_index> fluid_dof_indices(dofs_per_cell);

  Tensor<2, 3> velocity_gradient;
  double       pressure;
  Tensor<1, 3> normal_vector;
  Tensor<2, 3> fluid_stress;
  Tensor<2, 3> fluid_pressure;
  Tensor<1, 3> force; // to be changed for a vector of tensors when
  // allowing multiple solids


  Assert(
    !this->simulation_parameters.physical_properties_manager.is_non_newtonian(),
    RequiresConstantViscosity("assemble_nitsche_restriction"));

  const double viscosity =
    this->simulation_parameters.physical_properties_manager.viscosity_scale;

  // Loop over all local particles
  auto particle = solid_ph->begin();
  while (particle != solid_ph->end())
    {
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
      const auto &cell = particle->get_surrounding_cell(*this->triangulation);
#else
      const auto &cell = particle->get_surrounding_cell();
#endif
      const auto &dh_cell =
        typename DoFHandler<3>::cell_iterator(*cell, &this->dof_handler);
      dh_cell->get_dof_indices(fluid_dof_indices);

      const auto pic = solid_ph->particles_in_cell(cell);
      Assert(pic.begin() == particle, ExcInternalError());

      // Generate FEField function to evaluate values and gradients
      // at the particle location
      auto &evaluation_point = this->evaluation_point;

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
      Functions::
        FEFieldFunction<3, DoFHandler<3, 3>, TrilinosWrappers::MPI::Vector>
          fe_field(this->dof_handler, evaluation_point);
#else
      Functions::FEFieldFunction<3, TrilinosWrappers::MPI::Vector> fe_field(
        this->dof_handler, evaluation_point);
#endif


      fe_field.set_active_cell(dh_cell);

      for (const auto &p : pic)
        {
          velocity_gradient = 0;
          pressure          = 0;
          const auto &q     = p.get_location();
          const auto &JxW   = p.get_properties()[0];
          normal_vector[0]  = -p.get_properties()[1];
          normal_vector[1]  = -p.get_properties()[2];
          normal_vector[2]  = -p.get_properties()[3];

          for (int k = 0; k < 3; ++k)
            {
              velocity_gradient[k] = fe_field.gradient(q, k);
            }

          pressure = fe_field.value(q, 3);

          for (int d = 0; d < 2; ++d)
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
Tensor<1, spacedim>
GLSNitscheNavierStokesSolver<dim, spacedim>::calculate_forces_on_solid(
  const unsigned int i_solid)
{
  std::shared_ptr<Particles::ParticleHandler<dim, spacedim>> &solid_ph =
    solids[i_solid]->get_solid_particle_handler();

  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

  std::vector<types::global_dof_index> fluid_dof_indices(dofs_per_cell);

  // Penalization terms
  const double beta =
    this->simulation_parameters.nitsche->nitsche_solids[i_solid]->beta;
  Tensor<1, spacedim> velocity;
  Function<spacedim> *solid_velocity = solids[i_solid]->get_solid_velocity();
  Tensor<1, spacedim> force;
  for (unsigned int i = 0; i < spacedim; ++i)
    force[i] = 0;

  // Loop over all local particles
  auto particle = solid_ph->begin();
  while (particle != solid_ph->end())
    {
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
      const auto &cell = particle->get_surrounding_cell(*this->triangulation);
#else
      const auto &cell = particle->get_surrounding_cell();
#endif
      double h_cell = 0;
      if (dim == 2)
        h_cell =
          std::sqrt(4. * cell->measure() / M_PI) / this->velocity_fem_degree;
      else if (dim == 3)
        h_cell =
          pow(6 * cell->measure() / M_PI, 1. / 3.) / this->velocity_fem_degree;
      const double penalty_parameter =
        1. / std::pow(h_cell * h_cell, double(dim) / double(spacedim));
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
              const auto comp_k = this->fe->system_to_component_index(k).first;
              if (comp_k < spacedim)
                {
                  // Get the velocity at non-quadrature point (particle in
                  // fluid)
                  auto &evaluation_point = this->evaluation_point;
                  velocity[comp_k] += evaluation_point[fluid_dof_indices[k]] *
                                      this->fe->shape_value(k, ref_q);
                }
            }
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const auto comp_i = this->fe->system_to_component_index(i).first;
              if (comp_i < spacedim)
                {
                  force[comp_i] +=
                    penalty_parameter * beta * this->fe->shape_value(i, ref_q) *
                    JxW *
                    (solid_velocity->value(real_q, comp_i) - velocity[comp_i]);
                }
            }
        }
      particle = pic.end();
    }
  force = Utilities::MPI::sum(force, this->mpi_communicator);
  return force;
}

template <int dim, int spacedim>
Tensor<1, 3>
GLSNitscheNavierStokesSolver<dim, spacedim>::calculate_torque_on_solid(
  const unsigned int i_solid)
{
  std::shared_ptr<Particles::ParticleHandler<spacedim>> &solid_ph =
    solids[i_solid]->get_solid_particle_handler();

  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

  std::vector<types::global_dof_index> fluid_dof_indices(dofs_per_cell);

  // Penalization terms
  const double beta =
    this->simulation_parameters.nitsche->nitsche_solids[i_solid]->beta;
  Tensor<1, spacedim> velocity;
  Function<spacedim> *solid_velocity = solids[i_solid]->get_solid_velocity();


  Tensor<1, 3> torque;
  torque = 0;

  // Todo center of rotation should be parameter passed.
  Point<spacedim> center_of_rotation =
    this->simulation_parameters.nitsche->nitsche_solids[i_solid]->cor;

  // Loop over all local particles
  auto particle = solid_ph->begin();
  while (particle != solid_ph->end())
    {
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
      const auto &cell = particle->get_surrounding_cell(*this->triangulation);
#else
      const auto &cell = particle->get_surrounding_cell();
#endif
      double h_cell = 0;
      if (dim == 2)
        h_cell =
          std::sqrt(4. * cell->measure() / M_PI) / this->velocity_fem_degree;
      else if (dim == 3)
        h_cell =
          pow(6 * cell->measure() / M_PI, 1. / 3.) / this->velocity_fem_degree;
      const double penalty_parameter =
        1. / std::pow(h_cell * h_cell, double(dim) / double(spacedim));
      const auto &dh_cell =
        typename DoFHandler<spacedim>::cell_iterator(*cell, &this->dof_handler);
      dh_cell->get_dof_indices(fluid_dof_indices);

      const auto pic = solid_ph->particles_in_cell(cell);
      Assert(pic.begin() == particle, ExcInternalError());
      for (const auto &p : pic)
        {
          Tensor<1, spacedim> force;
          force              = 0;
          velocity           = 0;
          const auto &ref_q  = p.get_reference_location();
          const auto &real_q = p.get_location();
          const auto &JxW    = p.get_properties()[0];

          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
              const auto comp_k = this->fe->system_to_component_index(k).first;
              if (comp_k < spacedim)
                {
                  // Get the velocity at non-quadrature point (particle in
                  // fluid)
                  auto &evaluation_point = this->evaluation_point;
                  velocity[comp_k] += evaluation_point[fluid_dof_indices[k]] *
                                      this->fe->shape_value(k, ref_q);
                }
            }
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const auto comp_i = this->fe->system_to_component_index(i).first;
              if (comp_i < spacedim)
                {
                  force[comp_i] -=
                    penalty_parameter * beta * this->fe->shape_value(i, ref_q) *
                    JxW *
                    (solid_velocity->value(real_q, comp_i) - velocity[comp_i]);
                }
            }
          // Calculate torque on particle location
          auto distance = real_q - center_of_rotation;

          if (dim == 2)
            {
              torque[0] = 0.;
              torque[1] = 0.;
              torque[2] += distance[0] * force[1] - distance[1] * force[0];
            }
          else if (dim == 3)
            {
              torque[0] += distance[1] * force[2] - distance[2] * force[1];
              torque[1] += distance[2] * force[0] - distance[0] * force[2];
              torque[2] += distance[0] * force[1] - distance[1] * force[0];
            }
        }
      particle = pic.end();
    }
  torque = Utilities::MPI::sum(torque, this->mpi_communicator);
  return torque;
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::postprocess_solid_forces(
  const unsigned int i_solid,
  bool               first_solid_forces)
{
  TimerOutput::Scope t(this->computing_timer, "calculate_force_on_solid");

  std::vector<Tensor<1, spacedim>> force;
  std::vector<unsigned int>        solid_indices;

  force.push_back(this->calculate_forces_on_solid(i_solid));
  solid_indices.push_back(i_solid);

  if (this->simulation_parameters.nitsche->verbosity ==
        Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
    {
      std::string              independent_column_names = "Solid ID";
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

      // printed only for the first i_solid for which calculate_forces is true
      // to be improved so that the column_names are not repeated
      if (first_solid_forces)
        {
          std::cout << std::endl;
          std::cout << "+------------------------------------------+"
                    << std::endl;
          std::cout << "|  Force on solid summary                  |"
                    << std::endl;
          std::cout << "+------------------------------------------+"
                    << std::endl;
        }
      table.write_text(std::cout);
    }

  if (this->simulation_control->is_steady())
    {
      solid_forces_table[i_solid].add_value(
        "cells", this->triangulation->n_global_active_cells());
    }
  else
    {
      solid_forces_table[i_solid].add_value(
        "time", this->simulation_control->get_current_time());
      solid_forces_table[i_solid].set_precision(
        "time", this->simulation_parameters.forces_parameters.output_precision);
    }
  solid_forces_table[i_solid].add_value("f_x", force[0][0]);
  solid_forces_table[i_solid].add_value("f_y", force[0][1]);
  if (dim == 3)
    solid_forces_table[i_solid].add_value("f_z", force[0][2]);
  else
    solid_forces_table[i_solid].add_value("f_z", 0.);

  // Precision
  solid_forces_table[i_solid].set_precision(
    "f_x", this->simulation_parameters.forces_parameters.output_precision);
  solid_forces_table[i_solid].set_precision(
    "f_y", this->simulation_parameters.forces_parameters.output_precision);
  solid_forces_table[i_solid].set_precision(
    "f_z", this->simulation_parameters.forces_parameters.output_precision);

  std::string filename_force =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.nitsche->nitsche_solids[i_solid]
      ->force_output_name +
    "_" + Utilities::int_to_string(i_solid, 2) + ".dat";
  std::ofstream output_force(filename_force.c_str());

  solid_forces_table[i_solid].write_text(output_force);
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::postprocess_solid_torques(
  const unsigned int i_solid,
  bool               first_solid_torques)
{
  TimerOutput::Scope t(this->computing_timer, "calculate_torque_on_solid");

  std::vector<Tensor<1, 3>> torque;
  std::vector<unsigned int> solid_indices;

  torque.push_back(this->calculate_torque_on_solid(i_solid));
  solid_indices.push_back(i_solid);

  if (this->simulation_parameters.nitsche->verbosity ==
        Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
    {
      std::string              independent_column_names = "Solid ID";
      std::vector<std::string> dependent_column_names;
      dependent_column_names.push_back("T_x");
      dependent_column_names.push_back("T_y");
      dependent_column_names.push_back("T_z");

      TableHandler table = make_table_scalars_tensors(
        solid_indices,
        independent_column_names,
        torque,
        dependent_column_names,
        this->simulation_parameters.simulation_control.log_precision);

      // printed only for the first i_solid for which calculate_torques is true
      // to be improved so that the column_names are not repeated
      if (first_solid_torques)
        {
          std::cout << std::endl;
          std::cout << "+------------------------------------------+"
                    << std::endl;
          std::cout << "|  Torque on solids summary                |"
                    << std::endl;
          std::cout << "+------------------------------------------+"
                    << std::endl;
        }
      table.write_text(std::cout);
    }

  if (this->simulation_control->is_steady())
    {
      solid_torques_table[i_solid].add_value(
        "cells", this->triangulation->n_global_active_cells());
    }
  else
    {
      solid_torques_table[i_solid].add_value(
        "time", this->simulation_control->get_current_time());
      solid_torques_table[i_solid].set_precision(
        "time", this->simulation_parameters.forces_parameters.output_precision);
    }
  solid_torques_table[i_solid].add_value("T_x", torque[0][0]);
  solid_torques_table[i_solid].add_value("T_y", torque[0][1]);
  solid_torques_table[i_solid].add_value("T_z", torque[0][2]);

  // Precision
  solid_torques_table[i_solid].set_precision(
    "T_x", this->simulation_parameters.forces_parameters.output_precision);
  solid_torques_table[i_solid].set_precision(
    "T_y", this->simulation_parameters.forces_parameters.output_precision);
  solid_torques_table[i_solid].set_precision(
    "T_z", this->simulation_parameters.forces_parameters.output_precision);

  if (this->this_mpi_process == 0)
    {
      std::string filename_torque =
        this->simulation_parameters.simulation_control.output_folder +
        this->simulation_parameters.nitsche->nitsche_solids[i_solid]
          ->torque_output_name +
        "_" + Utilities::int_to_string(i_solid, 2) + ".dat";
      std::ofstream output_torque(filename_torque.c_str());

      solid_torques_table[i_solid].write_text(output_torque);
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::solve()
{
  MultithreadInfo::set_thread_limit(1);

  // Fluid setup
  read_mesh_and_manifolds(
    *this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  this->setup_dofs();

  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);

  // Solid setup
  if (!this->simulation_parameters.restart_parameters.restart)
    {
      TimerOutput::Scope t(this->computing_timer,
                           "Nitsche solid mesh and particles");
      for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
        {
          solids[i_solid]->initial_setup();

          // Output initial configuration, if output_frequency!=0
          if (this->simulation_control->is_output_iteration())
            {
              output_solid_particles(i_solid);
              output_solid_triangulation(i_solid);
            }
        }
    }

  // Time integration
  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);

      {
        TimerOutput::Scope t(this->computing_timer, "Nitsche particles motion");
        for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
          {
            if (this->simulation_parameters.nitsche->nitsche_solids[i_solid]
                  ->enable_particles_motion)
              {
                // Particle and solid displacement is explicit, thus it must go
                // from t to t+dt
                const double time_step =
                  this->simulation_control->get_time_step();
                const double initial_time =
                  this->simulation_control->get_current_time() - time_step;

                solids[i_solid]->update_temperature_time(time_step);
                solids[i_solid]->integrate_velocity(time_step, initial_time);
                solids[i_solid]->move_solid_triangulation(time_step,
                                                          initial_time);
              }
          }
        if constexpr (dim == spacedim)
          {
            // Parse the nitsche solids to the multiphysics interface
            this->multiphysics->set_solid(&solids);
          }
      }
      if (this->simulation_control->is_at_start())
        this->iterate();
      else
        {
          this->refine_mesh();
          this->iterate();
        }

      this->postprocess(false);
      bool first_solid_forces(true);
      bool first_solid_torques(true);

      for (unsigned int i_solid = 0;
           i_solid < this->simulation_parameters.nitsche->number_solids;
           ++i_solid)
        {
          if (this->simulation_parameters.nitsche->nitsche_solids[i_solid]
                ->calculate_force_on_solid)
            {
              postprocess_solid_forces(i_solid, first_solid_forces);
              first_solid_forces = false;
            }
          if (this->simulation_parameters.nitsche->nitsche_solids[i_solid]
                ->calculate_torque_on_solid)
            {
              // for 22 and 33 simulations
              //(to be implemented for 23 simulations)
              if (dim == spacedim)
                {
                  postprocess_solid_torques(i_solid, first_solid_torques);
                  first_solid_torques = false;
                }
            }

          if (this->simulation_control->is_output_iteration())
            {
              output_solid_particles(i_solid);
              output_solid_triangulation(i_solid);
            }
        }
      this->finish_time_step();
    }
  if (this->simulation_parameters.test.enabled)
    {
      for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
        {
          solids[i_solid]->print_particle_positions();
        }
    }
  this->finish_simulation();
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::output_solid_particles(
  const unsigned int i_solid)
{
  std::shared_ptr<Particles::ParticleHandler<spacedim>> &particle_handler =
    solids[i_solid]->get_solid_particle_handler();
  Particles::DataOut<spacedim, spacedim> particles_out;
  particles_out.build_patches(*particle_handler);

  const std::string folder = this->simulation_control->get_output_path();
  const std::string solution_name =
    this->simulation_control->get_output_name() + "_solid_particles_" +
    Utilities::int_to_string(i_solid, 2);
  const unsigned int iter        = this->simulation_control->get_step_number();
  const double       time        = this->simulation_control->get_current_time();
  const unsigned int group_files = this->simulation_control->get_group_files();

  write_vtu_and_pvd<0, spacedim>(pvdhandler_solid_particles[i_solid],
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
GLSNitscheNavierStokesSolver<dim, spacedim>::output_solid_triangulation(
  const unsigned int i_solid)
{
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
  DataOut<dim, DoFHandler<dim, spacedim>> data_out;
#else
  DataOut<dim, spacedim> data_out;
#endif
  DoFHandler<dim, spacedim> &solid_dh =
    solids[i_solid]->get_solid_dof_handler();
  data_out.attach_dof_handler(solid_dh);

  DoFHandler<dim, spacedim> &displacement_dh =
    solids[i_solid]->get_displacement_dof_handler();
  data_out.attach_dof_handler(displacement_dh);

  std::vector<std::string> solution_names(spacedim, "displacement");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      spacedim, DataComponentInterpretation::component_is_part_of_vector);
  TrilinosWrappers::MPI::Vector &displacement_vector =
    solids[i_solid]->get_displacement_vector();

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
  data_out.add_data_vector(
    displacement_vector,
    solution_names,
    DataOut<dim, DoFHandler<dim, spacedim>>::type_dof_data,
    data_component_interpretation);
#else
  data_out.add_data_vector(displacement_vector,
                           solution_names,
                           DataOut<dim, spacedim>::type_dof_data,
                           data_component_interpretation);
#endif



  data_out.build_patches();

  const std::string folder = this->simulation_control->get_output_path();
  const std::string solution_name =
    this->simulation_control->get_output_name() + "_solid_triangulation_" +
    Utilities::int_to_string(i_solid, 2);
  const unsigned int iter        = this->simulation_control->get_step_number();
  const double       time        = this->simulation_control->get_current_time();
  const unsigned int group_files = this->simulation_control->get_group_files();

  write_vtu_and_pvd<dim>(pvdhandler_solid_triangulation[i_solid],
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
GLSNitscheNavierStokesSolver<dim, spacedim>::refine_mesh()
{
  if (this->simulation_control->get_step_number() %
        this->simulation_parameters.mesh_adaptation.frequency ==
      0)
    {
      // If no adaptation is to be carried out, get out of the function
      if (this->simulation_parameters.mesh_adaptation.type ==
          Parameters::MeshAdaptation::Type::none)
        return;

      // Prepare the solid particle handlers for the mesh refinement
      // All type of refinement except none will require that the particle
      // handler be prepared for refinement
      for (const auto &solid : solids)
        solid->get_solid_particle_handler()
          ->prepare_for_coarsening_and_refinement();


      if (this->simulation_parameters.mesh_adaptation.type ==
          Parameters::MeshAdaptation::Type::kelly)
        this->refine_mesh_kelly();

      else if (this->simulation_parameters.mesh_adaptation.type ==
               Parameters::MeshAdaptation::Type::uniform)
        this->refine_mesh_uniform();

      // Unpact them after refinement has occured
      for (const auto &solid : solids)
        solid->get_solid_particle_handler()
          ->unpack_after_coarsening_and_refinement();
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_matrix_and_rhs()
{
  this->GLSNavierStokesSolver<
    spacedim>::assemble_system_matrix_without_preconditioner();

  this->GLSNavierStokesSolver<spacedim>::assemble_system_rhs();

  assemble_nitsche_restriction<true>();
  this->setup_preconditioner();
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_rhs()
{
  this->GLSNavierStokesSolver<spacedim>::assemble_system_rhs();

  assemble_nitsche_restriction<false>();
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::write_checkpoint()
{
  // Prepare solid particles for serialization
  for (auto solid : solids)
    solid->get_solid_particle_handler()->prepare_for_serialization();

  // Call regular checkpointing routine
  this->GLSNavierStokesSolver<spacedim>::write_checkpoint();

  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;

  // Write data about paraview format (pvd)
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      this->simulation_control->save(prefix);
      // Navier-Stokes
      this->pvdhandler.save(prefix);
      // Nitche
      for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
        {
          pvdhandler_solid_particles[i_solid].save(
            prefix + "_solid_particles_" +
            Utilities::int_to_string(i_solid, 2));
          pvdhandler_solid_triangulation[i_solid].save(
            prefix + "_solid_triangulation_" +
            Utilities::int_to_string(i_solid, 2));
        }
    }

  // Write solid base checkpoint
  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    {
      std::string filename =
        this->simulation_parameters.simulation_control.output_folder +
        this->simulation_parameters.restart_parameters.filename + "_solid_" +
        Utilities::int_to_string(i_solid, 2);
      this->solids[i_solid]->write_checkpoint(filename);
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::read_checkpoint()
{
  this->GLSNavierStokesSolver<spacedim>::read_checkpoint();

  TimerOutput::Scope t(this->computing_timer,
                       "Reset Nitsche solid mesh and particles");


  // Reload initial configurations
  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    {
      std::string prefix =
        this->simulation_parameters.simulation_control.output_folder +
        this->simulation_parameters.restart_parameters.filename + "_solid_" +
        Utilities::int_to_string(i_solid, 2);

      solids[i_solid]->read_checkpoint(prefix);
    }

  // Reload particle and solid pvd handlers
  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;
  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    {
      // Load Paraview pvd handler for solid triangulation and particles
      pvdhandler_solid_particles[i_solid].read(
        prefix + "_solid_particles_" + Utilities::int_to_string(i_solid, 2));
      pvdhandler_solid_triangulation[i_solid].read(
        prefix + "_solid_triangulation_" +
        Utilities::int_to_string(i_solid, 2));

      // Refill force and torque table from checkpoint
      if (this->simulation_parameters.nitsche->nitsche_solids[i_solid]
            ->calculate_force_on_solid)
        {
          std::string filename_force =
            this->simulation_parameters.simulation_control.output_folder +
            this->simulation_parameters.nitsche->nitsche_solids[i_solid]
              ->force_output_name +
            "_" + Utilities::int_to_string(i_solid, 2) + ".dat";
          fill_table_from_file(solid_forces_table[i_solid], filename_force);
        }
      if (this->simulation_parameters.nitsche->nitsche_solids[i_solid]
            ->calculate_torque_on_solid)
        {
          std::string filename_torque =
            this->simulation_parameters.simulation_control.output_folder +
            this->simulation_parameters.nitsche->nitsche_solids[i_solid]
              ->torque_output_name +
            "_" + Utilities::int_to_string(i_solid, 2) + ".dat";
          fill_table_from_file(solid_torques_table[i_solid], filename_torque);
        }
    }
}


// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library
// is valid before we actually compile the solver This greatly helps with
// debugging
template class GLSNitscheNavierStokesSolver<2>;
template class GLSNitscheNavierStokesSolver<2, 3>;
template class GLSNitscheNavierStokesSolver<3>;
