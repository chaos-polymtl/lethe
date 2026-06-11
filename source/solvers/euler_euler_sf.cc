#include <core/grids.h>

#include <solvers/euler_euler_sf.h>

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

// using namespace dealii;

template <int dim>
EulerEulerOneWay<dim>::EulerEulerOneWay(
  CFDDEMSimulationParameters<dim> &fluid_parameters,
  const SolidPhaseParameters<dim> &solid_parameters)
  : FluidDynamicsVANS<dim>(fluid_parameters)
  , solid_solver(solid_parameters,
                 this->simulation_control,
                 dynamic_cast<parallel::distributed::Triangulation<dim> &>(
                   *this->triangulation),
                 this->mpi_communicator)
  , euler_void_fraction(this->mpi_communicator, true)
{}


template <int dim>
void
EulerEulerOneWay<dim>::pass_fluid_solution_to_solid()
{
  solid_solver.set_fluid_velocity_field(this->dof_handler,
                                        *this->mapping,
                                        this->present_solution);
}


template <int dim>
void
EulerEulerOneWay<dim>::update_void_fraction_from_solid()
{
  const TrilinosWrappers::MPI::Vector &alpha_s =
    solid_solver.get_solid_volume_fraction();

  euler_void_fraction.set_solid_volume_fraction(alpha_s);
  euler_void_fraction.calculate_alpha_f();

  const TrilinosWrappers::MPI::Vector &alpha_f =
    euler_void_fraction.get_alpha_f();


  this->void_fraction_manager.void_fraction_locally_relevant = alpha_f;
}

template <int dim>
void
EulerEulerOneWay<dim>::assemble_fluid_drag_exchange_rhs()
{
  fluid_drag_rhs = 0;

  QGauss<dim> quadrature_formula(this->number_quadrature_points);

  FEValues<dim> fluid_fe_values(*this->mapping,
                                *this->fe,
                                quadrature_formula,
                                update_values | update_quadrature_points |
                                  update_JxW_values);

  FEValues<dim> solid_fe_values(*this->mapping,
                                solid_solver.get_fe(),
                                quadrature_formula,
                                update_values);

  const FEValuesExtractors::Vector fluid_velocities(0);
  const FEValuesExtractors::Vector solid_velocities(0);
  const FEValuesExtractors::Scalar solid_alpha(dim);

  const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell();
  const unsigned int n_q           = quadrature_formula.size();

  FullMatrix<double> cell_matrix;
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<Tensor<1, dim>> u_f_q(n_q);
  std::vector<Tensor<1, dim>> u_s_q(n_q);
  std::vector<double>         a_s_q(n_q);

  for (const auto &fluid_cell : this->dof_handler.active_cell_iterators())
    {
      if (!fluid_cell->is_locally_owned())
        continue;

      typename DoFHandler<dim>::active_cell_iterator solid_cell(
        &(*this->triangulation),
        fluid_cell->level(),
        fluid_cell->index(),
        &solid_solver.get_dof_handler());

      fluid_fe_values.reinit(fluid_cell);
      solid_fe_values.reinit(solid_cell);

      cell_rhs = 0;

      fluid_fe_values[fluid_velocities].get_function_values(
        this->present_solution, u_f_q);

      solid_fe_values[solid_velocities].get_function_values(
        solid_solver.get_locally_relevant_solution(), u_s_q);

      solid_fe_values[solid_alpha].get_function_values(
        solid_solver.get_locally_relevant_solution(), a_s_q);

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const Tensor<1, dim> drag_q =
            solid_solver.get_beta() * a_s_q[q] * (u_s_q[q] - u_f_q[q]);

          if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0 &&
              std::abs(fluid_fe_values.quadrature_point(q)[0] - 0.526416) <
                1e-3 &&
              std::abs(fluid_fe_values.quadrature_point(q)[1] - 0.223584) <
                1e-3 &&
              std::abs(fluid_fe_values.quadrature_point(q)[2] - 0.223584) <
                1e-3)
            {
              this->pcout << "fluid drag | x_q = "
                          << fluid_fe_values.quadrature_point(q)
                          << " | alpha_s = " << a_s_q[q]
                          << " | u_s = " << u_s_q[q] << " | u_f = " << u_f_q[q]
                          << " | drag_to_fluid = " << drag_q << std::endl;
            }

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const unsigned int comp_i =
                this->fe->system_to_component_index(i).first;

              if (comp_i < dim)
                {
                  const double v_i =
                    fluid_fe_values[fluid_velocities].value(i, q)[comp_i];

                  cell_rhs(i) += v_i * drag_q[comp_i] * fluid_fe_values.JxW(q);
                }
            }
        }

      fluid_cell->get_dof_indices(local_dof_indices);

      this->zero_constraints.distribute_local_to_global(cell_rhs,
                                                        local_dof_indices,
                                                        fluid_drag_rhs);
    }

  fluid_drag_rhs.compress(VectorOperation::add);
}


template <int dim>
void
EulerEulerOneWay<dim>::solve()
{
  read_mesh_and_manifolds(
    *this->triangulation,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh,
    this->cfd_dem_simulation_parameters.cfd_parameters.manifolds_parameters,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
        .restart ||
      this->cfd_dem_simulation_parameters.void_fraction->read_dem == true,
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions);

  for (const auto &cell : this->triangulation->active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary())
        {
          const auto center = face->center();

          if (std::abs(center[0] - 0.0) < 1e-12)
            face->set_boundary_id(1); // inlet
          else if (std::abs(center[0] - 1.0) < 1e-12)
            face->set_boundary_id(2); // outlet
          else
            face->set_boundary_id(0); // walls
        }

  //   for (const auto &cell : this->triangulation->active_cell_iterators())
  //     for (const auto &face : cell->face_iterators())
  //       if (face->at_boundary())
  //         this->pcout << "face center = " << face->center()
  //                     << " boundary id = " << face->boundary_id() <<
  //                     std::endl;

  this->setup_dofs();

  fluid_drag_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);

  solid_solver.setup();

  this->set_initial_condition(
    this->cfd_dem_simulation_parameters.cfd_parameters.initial_condition->type,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
      .restart);

  this->vertices_cell_mapping();

  this->void_fraction_manager.initialize_void_fraction(
    this->simulation_control->get_current_time());

  pass_fluid_solution_to_solid();

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);

      this->update_boundary_conditions();
      this->multiphysics->update_boundary_conditions();

      if (!this->simulation_control->is_at_start())
        {
          NavierStokesBase<dim, GlobalVectorType, IndexSet>::refine_mesh();
          this->vertices_cell_mapping();
        }

      /*
       * Fluid -> solid.
       */
      pass_fluid_solution_to_solid();

      /*
       * Solid solve.
       */
      solid_solver.advance_one_step();

      /*
       * Solid -> fluid.
       */
      update_void_fraction_from_solid();


      assemble_fluid_drag_exchange_rhs();

      this->pcout << "||fluid_drag_rhs|| = " << fluid_drag_rhs.l2_norm()
                  << std::endl;

      this->set_euler_euler_drag_rhs(fluid_drag_rhs);

      /*
       * Fluid solve.
       * This uses the updated void fraction field.
       */
      this->iterate();



      this->postprocess(false);
      this->finish_time_step_fd();

      this->FluidDynamicsVANS<dim>::monitor_mass_conservation();
    }

  solid_solver.finalize();

  this->finish_simulation();
}



template class EulerEulerOneWay<2>;
template class EulerEulerOneWay<3>;
