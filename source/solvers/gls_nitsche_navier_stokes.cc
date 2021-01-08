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
  NavierStokesSolverParameters<spacedim> &p_nsparam)
  : GLSNavierStokesSolver<spacedim>(p_nsparam)
  , solid(this->nsparam.nitsche,
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
                  auto &evaluation_point = this->get_evaluation_point();
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
      auto &system_rhs = this->get_system_rhs();
      constraints_used.distribute_local_to_global(local_matrix,
                                                  local_rhs,
                                                  fluid_dof_indices,
                                                  this->system_matrix,
                                                  system_rhs);
      particle = pic.end();
    }
  this->system_matrix.compress(VectorOperation::add);
  auto &system_rhs = this->get_system_rhs();
  system_rhs.compress(VectorOperation::add);
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

      // Generate FEField functoin to evaluate values and gradients
      // at the particle location
      auto &evaluation_point = this->get_evaluation_point();
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


  if (this->nsparam.nitsche->verbosity == Parameters::Verbosity::verbose &&
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
        this->nsparam.simulation_control.log_precision);

      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Force on solid summary                  |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      table.write_text(std::cout);
    }

  std::string   filename = this->nsparam.nitsche->force_output_name + ".dat";
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
    "f_x", this->nsparam.forces_parameters.output_precision);
  solid_forces_table.set_precision(
    "f_y", this->nsparam.forces_parameters.output_precision);
  solid_forces_table.set_precision(
    "f_z", this->nsparam.forces_parameters.output_precision);
  solid_forces_table.set_precision(
    "time", this->nsparam.forces_parameters.output_precision);

  solid_forces_table.write_text(output);
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  if (PhysicsSolver<TrilinosWrappers::MPI::Vector>::get_current_physics() ==
      Parameters::Multiphysics::fluid)
    {
      this->GLSNavierStokesSolver<spacedim>::assemble_matrix_and_rhs(
        time_stepping_method);

      assemble_nitsche_restriction<true>();
    }

  if (PhysicsSolver<TrilinosWrappers::MPI::Vector>::get_current_physics() ==
      Parameters::Multiphysics::heat)
    {
      assemble_matrix_and_rhs_ht(time_stepping_method);
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  if (PhysicsSolver<TrilinosWrappers::MPI::Vector>::get_current_physics() ==
      Parameters::Multiphysics::fluid)
    {
      this->GLSNavierStokesSolver<spacedim>::assemble_rhs(time_stepping_method);

      assemble_nitsche_restriction<false>();
    }

  if (PhysicsSolver<TrilinosWrappers::MPI::Vector>::get_current_physics() ==
      Parameters::Multiphysics::heat)
    {
      assemble_matrix_and_rhs_ht(time_stepping_method);
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::solve()
{
  read_mesh_and_manifolds(this->triangulation,
                          this->nsparam.mesh,
                          this->nsparam.manifolds_parameters,
                          this->nsparam.restart_parameters.restart,
                          this->nsparam.boundary_conditions);

  this->setup_dofs();
  this->set_initial_condition(this->nsparam.initial_condition->type,
                              this->nsparam.restart_parameters.restart);
  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      if (this->nsparam.nitsche->enable_particles_motion)
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
      postprocess_ht();
      if (this->nsparam.nitsche->calculate_force_on_solid && dim == 2 &&
          spacedim == 3)
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
      if (this->nsparam.multiphysics.heat_transfer)
        {
          finish_time_step_ht();
        }
    }
  if (this->nsparam.test.enabled)
    solid.print_particle_positions();
  this->finish_simulation();
  if (this->nsparam.multiphysics.heat_transfer)
    {
      finish_ht();
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
GLSNitscheNavierStokesSolver<dim, spacedim>::setup_dofs()
{
  this->setup_dofs_cfd();
  if (this->nsparam.multiphysics.heat_transfer)
    {
      this->setup_dofs_ht();
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::setup_dofs_ht()
{
  // implementation similar to deal.ii step-6
  this->dof_handler_ht.initialize(*(this->triangulation), this->fe_ht);

  locally_owned_dofs_ht = dof_handler_ht.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(dof_handler_ht,
                                          locally_relevant_dofs_ht);

  auto &solution_ht =
    this->get_present_solution(Parameters::Multiphysics::heat);
  solution_ht.reinit(locally_owned_dofs_ht,
                     locally_relevant_dofs_ht,
                     this->mpi_communicator);

  auto &system_rhs_ht = this->get_system_rhs(Parameters::Multiphysics::heat);
  system_rhs_ht.reinit(locally_owned_dofs_ht);

  auto &newton_update_ht =
    this->get_newton_update(Parameters::Multiphysics::heat);
  newton_update_ht.reinit(locally_owned_dofs_ht);

  TrilinosWrappers::MPI::Vector &local_evaluation_point_ht =
    this->get_local_evaluation_point(Parameters::Multiphysics::heat);
  local_evaluation_point_ht.reinit(this->locally_owned_dofs,
                                   this->mpi_communicator);

  // Previous solutions for transient schemes
  solution_ht_m1.reinit(locally_owned_dofs_ht,
                        locally_relevant_dofs_ht,
                        this->mpi_communicator);
  solution_ht_m2.reinit(locally_owned_dofs_ht,
                        locally_relevant_dofs_ht,
                        this->mpi_communicator);
  solution_ht_m3.reinit(locally_owned_dofs_ht,
                        locally_relevant_dofs_ht,
                        this->mpi_communicator);

  // Non-zero constraints
  auto &nonzero_constraints_ht =
    this->get_nonzero_constraints(Parameters::Multiphysics::ID::heat);

  {
    nonzero_constraints_ht.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler_ht,
                                            nonzero_constraints_ht);

    for (unsigned int i_bc = 0;
         i_bc < this->nsparam.boundary_conditions_ht.size;
         ++i_bc)
      {
        // Dirichlet condition : imposed temperature at i_bc
        if (this->nsparam.boundary_conditions_ht.type[i_bc] ==
            BoundaryConditions::BoundaryType::temperature)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler_ht,
              this->nsparam.boundary_conditions_ht.id[i_bc],
              dealii::Functions::ConstantFunction<spacedim>(
                this->nsparam.boundary_conditions_ht.value[i_bc]),
              nonzero_constraints_ht);
          }
      }
  }
  nonzero_constraints_ht.close();

  // Boundary conditions for Newton correction
  {
    zero_constraints_ht.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler_ht,
                                            zero_constraints_ht);

    for (unsigned int i_bc = 0;
         i_bc < this->nsparam.boundary_conditions_ht.size;
         ++i_bc)
      {
        if (this->nsparam.boundary_conditions_ht.type[i_bc] ==
            BoundaryConditions::BoundaryType::temperature)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler_ht,
              this->nsparam.boundary_conditions_ht.id[i_bc],
              Functions::ZeroFunction<spacedim>(),
              zero_constraints_ht);
          }
      }
  }
  zero_constraints_ht.close();

  // Sparse matrices initialization
  DynamicSparsityPattern dsp_ht(this->dof_handler_ht.n_dofs());
  DoFTools::make_sparsity_pattern(this->dof_handler_ht,
                                  dsp_ht,
                                  nonzero_constraints_ht,
                                  /*keep_constrained_dofs = */ true);
  sparsity_pattern_ht.copy_from(dsp_ht);
  system_matrix_ht.reinit(sparsity_pattern_ht);
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::assemble_matrix_and_rhs_ht(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  const double density       = this->nsparam.physical_properties.density;
  const double specific_heat = this->nsparam.physical_properties.specific_heat;
  const double thermal_conductivity =
    this->nsparam.physical_properties.thermal_conductivity;

  auto &solution_ht =
    this->get_present_solution(Parameters::Multiphysics::heat);
  system_matrix_ht    = 0;
  auto &system_rhs_ht = this->get_system_rhs(Parameters::Multiphysics::heat);
  system_rhs_ht       = 0;

  // Vector for the BDF coefficients
  // The coefficients are stored in the following fashion :
  // 0 - n+1
  // 1 - n
  // 2 - n-1
  // 3 - n-2
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  Vector<double> bdf_coefs;

  if (time_stepping_method ==
        Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
      time_stepping_method ==
        Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
    bdf_coefs = bdf_coefficients(1, time_steps_vector);

  if (time_stepping_method ==
      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    bdf_coefs = bdf_coefficients(2, time_steps_vector);

  if (time_stepping_method ==
      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    bdf_coefs = bdf_coefficients(3, time_steps_vector);

  auto &source_term = this->nsparam.sourceTerm->heat_transfer_source;

  const QGauss<spacedim> quadrature_formula(fe_ht.degree + 1);
  FEValues<spacedim>     fe_values_ht(fe_ht,
                                  quadrature_formula,
                                  update_values | update_gradients |
                                    update_quadrature_points |
                                    update_JxW_values);

  auto &evaluation_point =
    this->get_evaluation_point(Parameters::Multiphysics::heat);

  auto &velocity_evaluation_point =
    this->get_evaluation_point(Parameters::Multiphysics::fluid);


  const unsigned int dofs_per_cell_ht = fe_ht.dofs_per_cell;

  FullMatrix<double> cell_matrix(dofs_per_cell_ht, dofs_per_cell_ht);
  Vector<double>     cell_rhs(dofs_per_cell_ht);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell_ht);
  const unsigned int                   n_q_points = quadrature_formula.size();
  std::vector<Tensor<1, spacedim>>     temperature_gradients(n_q_points);
  std::vector<double>                  source_term_values(n_q_points);


  //  const MappingQ<spacedim> mapping(this->velocity_fem_degree,
  //                                   this->nsparam.fem_parameters.qmapping_all);
  FEValues<spacedim> fe_values_flow(this->fe,
                                    quadrature_formula,
                                    update_values | update_quadrature_points |
                                      update_gradients);

  // FaceValues for Robin boundary condition
  QGauss<spacedim - 1>   face_quadrature_formula_ht(this->fe_ht.degree + 1);
  FEFaceValues<spacedim> fe_face_values_ht(fe_ht,
                                           face_quadrature_formula_ht,
                                           update_values |
                                             update_quadrature_points |
                                             update_JxW_values);

  // Velocity values
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  std::vector<Tensor<1, spacedim>> velocity_values(n_q_points);
  std::vector<Tensor<2, spacedim>> velocity_gradient_values(n_q_points);

  // Values for backward Euler scheme
  std::vector<double> present_temperature_values(n_q_points);
  std::vector<double> present_face_temperature_values(
    face_quadrature_formula_ht.size());

  std::vector<double> p1_temperature_values(n_q_points);
  std::vector<double> p2_temperature_values(n_q_points);
  std::vector<double> p3_temperature_values(n_q_points);
  std::vector<double> p4_temperature_values(n_q_points);

  for (const auto &cell : dof_handler_ht.active_cell_iterators())
    {
      cell_matrix = 0;
      cell_rhs    = 0;
      fe_values_ht.reinit(cell);

      fe_values_ht.get_function_gradients(evaluation_point,
                                          temperature_gradients);


      typename DoFHandler<spacedim>::active_cell_iterator velocity_cell(
        &(*this->triangulation),
        cell->level(),
        cell->index(),
        &this->dof_handler);
      fe_values_flow.reinit(velocity_cell);
      fe_values_flow[velocities].get_function_values(velocity_evaluation_point,
                                                     velocity_values);
      fe_values_flow[velocities].get_function_gradients(
        velocity_evaluation_point, velocity_gradient_values);

      // Gather present value
      fe_values_ht.get_function_values(solution_ht, present_temperature_values);

      // Gather the previous time steps for heat transfer depending on
      // the number of stages of the time integration method
      if (time_stepping_method !=
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        fe_values_ht.get_function_values(this->solution_ht_m1,
                                         p1_temperature_values);

      if (time_stepping_method_has_two_stages(time_stepping_method))
        fe_values_ht.get_function_values(this->solution_ht_m2,
                                         p2_temperature_values);

      if (time_stepping_method_has_three_stages(time_stepping_method))
        fe_values_ht.get_function_values(this->solution_ht_m3,
                                         p3_temperature_values);

      source_term.value_list(fe_values_ht.get_quadrature_points(),
                             source_term_values);


      // assembling local matrix and right hand side
      for (const unsigned int q : fe_values_ht.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values_ht.dof_indices())
            {
              for (const unsigned int j : fe_values_ht.dof_indices())
                {
                  // Weak form for : - k * laplacian T + rho * cp * u * grad T -
                  // f - grad(u)*grad(u) =0
                  cell_matrix(i, j) +=
                    (thermal_conductivity * fe_values_ht.shape_grad(i, q) *
                       fe_values_ht.shape_grad(j, q) +
                     density * specific_heat * fe_values_ht.shape_value(i, q) *
                       velocity_values[q] * fe_values_ht.shape_grad(j, q)) *
                    fe_values_ht.JxW(q); // JxW

                  // Mass matrix for transient simulation
                  if (is_bdf(time_stepping_method))
                    cell_matrix(i, j) += density * specific_heat *
                                         fe_values_ht.shape_value(j, q) *
                                         fe_values_ht.shape_value(i, q) *
                                         bdf_coefs[0] * fe_values_ht.JxW(q);
                }

              // rhs for : - k * laplacian T + rho * cp * u * grad T - f -
              // grad(u)*grad(u) = 0
              cell_rhs(i) -=
                (thermal_conductivity * fe_values_ht.shape_grad(i, q) *
                   temperature_gradients[q] +
                 density * specific_heat * fe_values_ht.shape_value(i, q) *
                   velocity_values[q] * temperature_gradients[q] -
                 source_term_values[q] * fe_values_ht.shape_value(i, q) -
                 fe_values_ht.shape_value(i, q) *
                   scalar_product(velocity_gradient_values[q] +
                                    transpose(velocity_gradient_values[q]),
                                  transpose(velocity_gradient_values[q]))) *
                fe_values_ht.JxW(q); // JxW

              // Residual associated with BDF schemes
              if (time_stepping_method ==
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
                  time_stepping_method == Parameters::SimulationControl::
                                            TimeSteppingMethod::steady_bdf)
                cell_rhs(i) -= (bdf_coefs[0] * present_temperature_values[q] +
                                bdf_coefs[1] * p1_temperature_values[q]) *
                               fe_values_ht.shape_value(i, q) *
                               fe_values_ht.JxW(q); // *phi_u[i]*JxW

              if (time_stepping_method ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                cell_rhs(i) -= (bdf_coefs[0] * present_temperature_values[q] +
                                bdf_coefs[1] * p1_temperature_values[q] +
                                bdf_coefs[2] * p2_temperature_values[q]) *
                               fe_values_ht.shape_value(i, q) *
                               fe_values_ht.JxW(q); // *phi_u[i]*JxW

              if (time_stepping_method ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                cell_rhs(i) -= (bdf_coefs[0] * present_temperature_values[q] +
                                bdf_coefs[1] * p1_temperature_values[q] +
                                bdf_coefs[2] * p2_temperature_values[q] +
                                bdf_coefs[3] * p3_temperature_values[q]) *
                               fe_values_ht.shape_value(i, q) *
                               fe_values_ht.JxW(q); // *phi_u[i]*JxW
            }

        } // end loop on quadrature points

      // Robin boundary condition, loop on faces (Newton's cooling law)
      // implementation similar to deal.ii step-7
      for (unsigned int i_bc = 0;
           i_bc < this->nsparam.boundary_conditions_ht.size;
           ++i_bc)
        {
          if (this->nsparam.boundary_conditions_ht.type[i_bc] ==
              BoundaryConditions::BoundaryType::convection)
            {
              if (cell->is_locally_owned())
                {
                  for (unsigned int face = 0;
                       face < GeometryInfo<dim>::faces_per_cell;
                       face++)
                    {
                      if (cell->face(face)->at_boundary() &&
                          (cell->face(face)->boundary_id() ==
                           this->nsparam.boundary_conditions_ht.id[i_bc]))
                        {
                          fe_face_values_ht.reinit(cell, face);
                          fe_face_values_ht.get_function_values(
                            solution_ht, present_face_temperature_values);
                          {
                            for (const unsigned int q :
                                 fe_face_values_ht.quadrature_point_indices())
                              {
                                for (const unsigned int i :
                                     fe_values_ht.dof_indices())
                                  {
                                    for (const unsigned int j :
                                         fe_values_ht.dof_indices())
                                      {
                                        // Weak form modification
                                        cell_matrix(i, j) +=
                                          fe_face_values_ht.shape_value(j, q) *
                                          fe_face_values_ht.shape_value(i, q) *
                                          this->nsparam.boundary_conditions_ht
                                            .value[i_bc] *
                                          fe_face_values_ht.JxW(q);
                                      }
                                    // Residual
                                    cell_rhs(i) -=
                                      fe_face_values_ht.shape_value(i, q) *
                                      this->nsparam.boundary_conditions_ht
                                        .value[i_bc] *
                                      (present_face_temperature_values[q] -
                                       this->nsparam.boundary_conditions_ht
                                         .Tenv[i_bc]) *
                                      fe_face_values_ht.JxW(q);
                                  }
                              }
                          }
                        }
                    }
                }
            }
        } // end loop for Robin condition

      // transfer cell contribution into global objects
      cell->get_dof_indices(local_dof_indices);
      zero_constraints_ht.distribute_local_to_global(cell_matrix,
                                                     cell_rhs,
                                                     local_dof_indices,
                                                     system_matrix_ht,
                                                     system_rhs_ht);
    } // end loop active cell
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::finish_time_step_ht()
{
  if (this->nsparam.simulation_control.method !=
      Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      this->solution_ht_m3 = this->solution_ht_m2;
      this->solution_ht_m2 = this->solution_ht_m1;
      auto &solution_ht =
        this->get_present_solution(Parameters::Multiphysics::heat);
      this->solution_ht_m1 = solution_ht;
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::output_field_hook(
  DataOut<spacedim> &data_out)
{
  if (this->nsparam.multiphysics.heat_transfer)
    {
      auto &present_temperature =
        this->get_present_solution(Parameters::Multiphysics::heat);
      data_out.add_data_vector(dof_handler_ht,
                               present_temperature,
                               "temperature");
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::set_initial_condition(
  Parameters::InitialConditionType initial_condition_type,
  bool                             restart)
{
  this->set_initial_condition_cfd(initial_condition_type, restart);
  if (this->nsparam.multiphysics.heat_transfer)
    {
      set_initial_condition_ht();
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::solve_linear_system(
  const bool initial_step,
  const bool renewed_matrix)
{
  if (PhysicsSolver<TrilinosWrappers::MPI::Vector>::get_current_physics() ==
      Parameters::Multiphysics::fluid)
    {
      this->solve_linear_system_cfd(initial_step, renewed_matrix);
    }

  if (PhysicsSolver<TrilinosWrappers::MPI::Vector>::get_current_physics() ==
      Parameters::Multiphysics::heat)
    {
      solve_linear_system_ht();
    }
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::set_initial_condition_ht()
{
  MappingQ<spacedim> mapping(fe_ht.degree);
  auto &newton_update = this->get_newton_update(Parameters::Multiphysics::heat);
  VectorTools::interpolate(mapping,
                           this->dof_handler_ht,
                           Functions::ZeroFunction<spacedim>(),
                           newton_update);
  auto &nonzero_constraints =
    this->get_nonzero_constraints(Parameters::Multiphysics::heat);
  nonzero_constraints.distribute(newton_update);
  auto &present_solution =
    this->get_present_solution(Parameters::Multiphysics::heat);
  present_solution = newton_update;

  this->finish_time_step();
  if (this->nsparam.multiphysics.heat_transfer)
    {
      finish_time_step_ht();
    }
  this->postprocess(true);
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::solve_linear_system_ht()
{
  auto &system_rhs_ht = this->get_system_rhs(Parameters::Multiphysics::heat);

  const double absolute_residual = this->nsparam.linear_solver.minimum_residual;
  const double relative_residual =
    this->nsparam.linear_solver.relative_residual;

  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs_ht.l2_norm(), absolute_residual);

  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const double ilu_fill = this->nsparam.linear_solver.ilu_precond_fill;
  const double ilu_atol = this->nsparam.linear_solver.ilu_precond_atol;
  const double ilu_rtol = this->nsparam.linear_solver.ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  TrilinosWrappers::PreconditionILU ilu_preconditioner;

  ilu_preconditioner.initialize(system_matrix_ht, preconditionerOptions);

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    locally_owned_dofs_ht, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false, this->nsparam.linear_solver.max_krylov_vectors);


  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);


  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix_ht,
                 completely_distributed_solution,
                 system_rhs_ht,
                 ilu_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }

    zero_constraints_ht.distribute(completely_distributed_solution);


    auto &newton_update =
      this->get_newton_update(Parameters::Multiphysics::heat);
    newton_update = completely_distributed_solution;
  }
}


template <int dim, int spacedim>
double
GLSNitscheNavierStokesSolver<dim, spacedim>::calculate_l2_error_ht(
  const DoFHandler<spacedim> &         dof_handler,
  const TrilinosWrappers::MPI::Vector &evaluation_point,
  const Function<spacedim> &           exact_solution,
  const Parameters::FEM &              fem_parameters,
  const MPI_Comm &                     mpi_communicator)
{
  const FiniteElement<spacedim> &fe = dof_handler.get_fe();


  QGauss<spacedim>         quadrature_formula(fe.degree + 2);
  const MappingQ<spacedim> mapping(fe.degree, fem_parameters.qmapping_all);
  FEValues<spacedim>       fe_values(mapping,
                               fe,
                               quadrature_formula,
                               update_values | update_gradients |
                                 update_quadrature_points | update_JxW_values);



  const unsigned int dofs_per_cell =
    fe.dofs_per_cell; // This gives you dofs per cell
  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points = quadrature_formula.size();

  std::vector<double> q_exact_solution(n_q_points);
  std::vector<double> q_scalar_values(n_q_points);

  double l2error = 0.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(evaluation_point, q_scalar_values);

          // Retrieve the effective "connectivity matrix" for this element
          cell->get_dof_indices(local_dof_indices);

          // Get the exact solution at all gauss points
          exact_solution.value_list(fe_values.get_quadrature_points(),
                                    q_exact_solution);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              double sim   = q_scalar_values[q];
              double exact = q_exact_solution[q];
              l2error += (sim - exact) * (sim - exact) * fe_values.JxW(q);
            }
        }
    }
  l2error = Utilities::MPI::sum(l2error, mpi_communicator);
  return l2error;
}

template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::postprocess_ht()
{
  if (this->nsparam.analytical_solution->calculate_error() == true &&
      this->nsparam.multiphysics.heat_transfer)
    {
      double temperature_error =
        calculate_l2_error_ht(dof_handler_ht,
                              this->get_evaluation_point(
                                Parameters::Multiphysics::heat),
                              this->nsparam.analytical_solution->temperature,
                              this->nsparam.fem_parameters,
                              this->mpi_communicator);

      this->error_table_ht.add_value(
        "cells", this->triangulation->n_global_active_cells());
      this->error_table_ht.add_value("error_temperature", temperature_error);

      if (this->nsparam.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error temperature : " << temperature_error
                      << std::endl;
        }
    }
}


template <int dim, int spacedim>
void
GLSNitscheNavierStokesSolver<dim, spacedim>::finish_ht()
{
  if (this->this_mpi_process == 0 &&
      this->nsparam.analytical_solution->verbosity ==
        Parameters::Verbosity::verbose &&
      this->nsparam.multiphysics.heat_transfer)
    {
      error_table_ht.omit_column_from_convergence_rate_evaluation("cells");
      error_table_ht.evaluate_all_convergence_rates(
        ConvergenceTable::reduction_rate_log2);
      error_table_ht.set_scientific("error_temperature", true);

      error_table_ht.write_text(std::cout);
    }
}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library
// is valid before we actually compile the solver This greatly helps with
// debugging
template class GLSNitscheNavierStokesSolver<2>;
template class GLSNitscheNavierStokesSolver<2, 3>;
template class GLSNitscheNavierStokesSolver<3>;
