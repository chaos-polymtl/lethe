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
  solid_solver.set_fluid_solution_field(this->dof_handler,
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

  this->void_fraction_manager.void_fraction_locally_relevant
    .update_ghost_values();
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

      AssertThrow(&(this->dof_handler.get_triangulation()) ==
                    &(solid_solver.get_dof_handler().get_triangulation()),
                  ExcMessage("Fluid and solid DoF handlers do not use "
                             "the same triangulation."));

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

          if constexpr (dim == 3)
            {
              const auto &point = fluid_fe_values.quadrature_point(q);

              if (Utilities::MPI::this_mpi_process(this->mpi_communicator) ==
                    0 &&
                  std::abs(point[0] - 0.526416) < 1e-3 &&
                  std::abs(point[1] - 0.223584) < 1e-3 &&
                  std::abs(point[2] - 0.223584) < 1e-3)
                {
                  this->pcout << "fluid drag | x_q = " << point
                              << " | alpha_s = " << a_s_q[q]
                              << " | u_s = " << u_s_q[q]
                              << " | u_f = " << u_f_q[q]
                              << " | drag_to_fluid = " << drag_q << std::endl;
                }
            }

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const unsigned int comp_i =
                this->fe->system_to_component_index(i).first;

              if (comp_i < dim)
                {
                  const Tensor<1, dim> v_i =
                    fluid_fe_values[fluid_velocities].value(i, q);

                  cell_rhs(i) += (v_i * drag_q) * fluid_fe_values.JxW(q);
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
EulerEulerOneWay<dim>::inspect_mesh_boundaries() const
{
  AssertThrow(this->mapping != nullptr,
              ExcMessage("The fluid mapping is null."));

  AssertThrow(this->fe != nullptr,
              ExcMessage("The fluid finite element is null."));

  const unsigned int face_quadrature_order = std::max(2u, this->fe->degree + 1);

  const QGauss<dim - 1> face_quadrature(face_quadrature_order);

  FEFaceValues<dim> fe_face_values(*this->mapping,
                                   *this->fe,
                                   face_quadrature,
                                   update_JxW_values);

  std::map<types::boundary_id, double> local_boundary_measure;

  std::map<types::boundary_id, unsigned long long> local_boundary_face_count;

  std::set<unsigned int> local_raw_boundary_ids;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
           ++face)
        {
          if (!cell->face(face)->at_boundary())
            continue;

          const types::boundary_id boundary_id =
            cell->face(face)->boundary_id();

          local_raw_boundary_ids.insert(static_cast<unsigned int>(boundary_id));

          fe_face_values.reinit(cell, face);

          double face_measure = 0.0;

          for (unsigned int q = 0; q < face_quadrature.size(); ++q)
            {
              face_measure += fe_face_values.JxW(q);
            }

          local_boundary_measure[boundary_id] += face_measure;

          local_boundary_face_count[boundary_id] += 1;
        }
    }

  const std::vector<unsigned int> local_ids(local_raw_boundary_ids.begin(),
                                            local_raw_boundary_ids.end());

  const auto gathered_ids =
    Utilities::MPI::all_gather(this->mpi_communicator, local_ids);

  std::set<unsigned int> global_raw_boundary_ids;

  for (const auto &rank_ids : gathered_ids)
    {
      global_raw_boundary_ids.insert(rank_ids.begin(), rank_ids.end());
    }

  AssertThrow(!global_raw_boundary_ids.empty(),
              ExcMessage(
                "The mesh does not contain any external boundary faces."));

  this->pcout << "\nImported mesh boundary information:\n";

  for (const unsigned int raw_id : global_raw_boundary_ids)
    {
      const types::boundary_id boundary_id =
        static_cast<types::boundary_id>(raw_id);

      double local_measure = 0.0;

      const auto measure_it = local_boundary_measure.find(boundary_id);

      if (measure_it != local_boundary_measure.end())
        local_measure = measure_it->second;

      unsigned long long local_face_count = 0;

      const auto count_it = local_boundary_face_count.find(boundary_id);

      if (count_it != local_boundary_face_count.end())
        local_face_count = count_it->second;

      const double global_measure =
        Utilities::MPI::sum(local_measure, this->mpi_communicator);

      const unsigned long long global_face_count =
        Utilities::MPI::sum(local_face_count, this->mpi_communicator);

      this->pcout << "  boundary ID " << raw_id
                  << ": measure = " << global_measure
                  << ", active boundary faces = " << global_face_count << '\n';
    }

  this->pcout << "  number of imported boundary IDs = "
              << global_raw_boundary_ids.size() << "\n"
              << std::endl;
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



  this->setup_dofs();

  inspect_mesh_boundaries();

  fluid_drag_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);

  solid_solver.setup();

  this->set_initial_condition(
    this->cfd_dem_simulation_parameters.cfd_parameters.initial_condition->type,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
      .restart);

  this->vertices_cell_mapping();

  this->void_fraction_manager.initialize_void_fraction(
    this->simulation_control->get_current_time());

  // auto &alpha_f =
  //   this->void_fraction_manager.void_fraction_locally_relevant;

  // for (const auto i : alpha_f.locally_owned_elements())
  //   alpha_f[i] = 1.0;

  // alpha_f.compress(VectorOperation::insert);
  // alpha_f.update_ghost_values();

  // update_void_fraction_from_solid();

  pass_fluid_solution_to_solid();

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);

      this->update_boundary_conditions();
      this->multiphysics->update_boundary_conditions();

      // if (!this->simulation_control->is_at_start())
      //   {
      //     NavierStokesBase<dim, GlobalVectorType, IndexSet>::refine_mesh();
      //     this->vertices_cell_mapping();
      //   }

      /*
       * Fluid -> solid.
       */
      pass_fluid_solution_to_solid();

      /*
       * Solid solve.
       */
      // solid_solver.advance_one_step();

      const bool solid_step_completed = solid_solver.advance_one_step();

      AssertThrow(solid_step_completed,
                  ExcMessage("The solid solver did not complete the "
                             "current coupled time step."));

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
