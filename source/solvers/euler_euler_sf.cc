#include <core/grids.h>

#include <solvers/euler_euler_sf.h>

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

// using namespace dealii;

// template <int dim>
// void
// make_euler_euler_grid(parallel::distributed::Triangulation<dim> &tria,
//                       const EulerEulerGridParameters<dim> &grid_parameters)
// {
//   GridGenerator::subdivided_hyper_rectangle(tria,
//                                             grid_parameters.subdivisions,
//                                             grid_parameters.p1,
//                                             grid_parameters.p2);

//   tria.refine_global(grid_parameters.global_refinement);

//   const double tol = 1e-12;

//   const auto axis_from = [](const std::string &d) -> unsigned int {
//     if (d == "x")
//       return 0;
//     if (d == "y")
//       return 1;
//     if (d == "z")
//       return 2;

//     AssertThrow(false, ExcMessage("direction must be x, y, or z"));
//     return 0;
//   };

//   const unsigned int axis1 = axis_from(grid_parameters.direction1);
//   // const unsigned int axis2 = axis_from(grid_parameters.direction2);

//   for (const auto &cell : tria.active_cell_iterators())
//     {
//       for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
//         {
//           if (cell->face(f)->at_boundary())
//             {
//               const Point<dim> fc = cell->face(f)->center();

//               cell->face(f)->set_boundary_id(0);

//               if (std::fabs(fc[axis1] - grid_parameters.p1[axis1]) < tol)
//                 cell->face(f)->set_boundary_id(1);
//               else if (std::fabs(fc[axis1] - grid_parameters.p2[axis1]) <
//               tol)
//                 cell->face(f)->set_boundary_id(2);
//               //   else if (std::fabs(fc[axis2] - grid_parameters.p1[axis2])
//               <
//               //   tol)
//               //     cell->face(f)->set_boundary_id(1);
//               //   else if (std::fabs(fc[axis2] - grid_parameters.p2[axis2])
//               <
//               //   tol)
//               //     cell->face(f)->set_boundary_id(2);
//             }
//         }
//     }
// }

// template <int dim>
// EulerEulerOneWay<dim>::EulerEulerOneWay(
//   const EulerEulerMeshParameters<dim> &mesh_parameters,
//   const SolidPhaseParameters          &solid_parameters,
//   CFDDEMSimulationParameters<dim>     &fluid_parameters,
//   const MPI_Comm                      &mpi_communicator,
//   const bool                           verbose)
//   : mesh_parameters(mesh_parameters)
//   , mpi_communicator(mpi_communicator)
//   , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
//   , verbose(verbose)
//   , solid_parameters(solid_parameters)
//   , fluid_solver(fluid_parameters)
//   , void_fraction_solver(fluid_solver, mpi_communicator, verbose)
// {}



// template <int dim>
// void
// EulerEulerOneWay<dim>::build_and_pass_alpha_f()
// {
//   if (verbose)
//     pcout << "Building alpha_f and passing it to fluid solver" << std::endl;

//   void_fraction_solver.set_solid_volume_fraction(
//     solid_solver->get_solid_volume_fraction());

//   void_fraction_solver.calculate_alpha_f();
//   void_fraction_solver.pass_alpha_f_to_fluid();
// }



// template <int dim>
// void
// EulerEulerOneWay<dim>::solve_solid_one_step()
// {
//   if (verbose)
//     pcout << "Solving solid phase one step" << std::endl;

//   solid_solver->advance_one_step();
// }

// template <int dim>
// void
// EulerEulerOneWay<dim>::solve_fluid_one_step()
// {
//   if (verbose)
//     pcout << "Solving fluid phase one step" << std::endl;

//   fluid_solver.advance_one_step_external();
// }

// template <int dim>
// void
// EulerEulerOneWay<dim>::run()
// {
//   if (verbose)
//     pcout << "Running Euler-Euler coupling" << std::endl;

//   auto &shared_triangulation = fluid_solver.get_triangulation();

//   if (shared_triangulation.n_global_active_cells() == 0)
//     {
//       if (verbose)
//         pcout << "Building Euler-Euler grid" << std::endl;

//       EulerEulerGridParameters<dim> grid_parameters;
//       grid_parameters.p1 = mesh_parameters.p1;
//       grid_parameters.p2 = mesh_parameters.p2;
//       grid_parameters.subdivisions.assign(dim, 1);
//       grid_parameters.subdivisions[0] = mesh_parameters.nx;
//       if (dim > 1)
//         grid_parameters.subdivisions[1] = mesh_parameters.ny;
//       if (dim > 2)
//         grid_parameters.subdivisions[2] = mesh_parameters.nz;
//       grid_parameters.global_refinement = mesh_parameters.global_refinement;
//       grid_parameters.direction1        = mesh_parameters.direction1;
//       grid_parameters.direction2        = mesh_parameters.direction2;

//       make_euler_euler_grid(shared_triangulation, grid_parameters);
//     }

//   solid_solver = std::make_unique<SolidPhaseSolver<dim>>(solid_parameters,
//                                                          shared_triangulation,
//                                                          mpi_communicator);

//   solid_solver->setup();
//   fluid_solver.setup_for_external_stepping();

//   while (!solid_solver->finished())
//     {
//       pcout << "Coupled step = " << solid_solver->get_step_number()
//             << std::endl;

//       void_fraction_solver.set_solid_volume_fraction(
//         solid_solver->get_solid_volume_fraction());

//       void_fraction_solver.calculate_alpha_f();
//       void_fraction_solver.pass_alpha_f_to_fluid();

//       pcout << "Before fluid step: "
//             << shared_triangulation.n_global_active_cells() << std::endl;

//       const bool fluid_ok = fluid_solver.advance_one_step_external();

//       pcout << "After fluid step:  "
//             << shared_triangulation.n_global_active_cells() << std::endl;

//       AssertThrow(
//         fluid_ok,
//         ExcMessage(
//           "Fluid finished before solid. Check time step and end time."));

//       solid_solver->set_fluid_velocity_field(
//         fluid_solver.get_fluid_dof_handler(),
//         fluid_solver.get_fluid_mapping(),
//         fluid_solver.get_fluid_solution());

//       solve_solid_one_step();
//     }

//   solid_solver->finalize();
//   fluid_solver.finish_external_stepping();
// }

// template void
// make_euler_euler_grid<2>(parallel::distributed::Triangulation<2> &tria,
//                          const EulerEulerGridParameters<2> &grid_parameters);

// template void
// make_euler_euler_grid<3>(parallel::distributed::Triangulation<3> &tria,
//                          const EulerEulerGridParameters<3> &grid_parameters);


// template class EulerEulerOneWay<2>;
// template class EulerEulerOneWay<3>;

template <int dim>
EulerEulerOneWay<dim>::EulerEulerOneWay(
  CFDDEMSimulationParameters<dim> &fluid_parameters,
  const SolidPhaseParameters      &solid_parameters)
  : FluidDynamicsVANS<dim>(fluid_parameters)
  , solid_solver(solid_parameters,
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

  /*
   * IMPORTANT:
   * This direct assignment only works if alpha_f has the same DoF layout as
   * this->void_fraction_manager.void_fraction_locally_relevant.
   *
   * If solid alpha_s is still on the solid mixed DoFHandler, this will fail or
   * give wrong results. Then you must interpolate/project alpha_s onto
   * this->void_fraction_manager.dof_handler first.
   */
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
                                update_values |
                                  update_quadrature_points |
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
                std::abs(fluid_fe_values.quadrature_point(q)[0] - 0.526416) < 1e-3 &&
                std::abs(fluid_fe_values.quadrature_point(q)[1] - 0.223584) < 1e-3 &&
                std::abs(fluid_fe_values.quadrature_point(q)[2] - 0.223584) < 1e-3)
              {
                this->pcout << "fluid drag | x_q = "
                            << fluid_fe_values.quadrature_point(q)
                            << " | alpha_s = " << a_s_q[q]
                            << " | u_s = " << u_s_q[q]
                            << " | u_f = " << u_f_q[q]
                            << " | drag_to_fluid = " << drag_q
                            << std::endl;
              }

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const unsigned int comp_i =
                this->fe->system_to_component_index(i).first;

              if (comp_i < dim)
                {
                  const double v_i =
                    fluid_fe_values[fluid_velocities].value(i, q)[comp_i];

                  cell_rhs(i) += v_i * drag_q[comp_i] *
                                 fluid_fe_values.JxW(q);
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
       * The solid solver uses the latest available fluid velocity field
       * when assembling the drag term.
       */
      pass_fluid_solution_to_solid();

      /*
       * Solid solve.
       * This advances alpha_s and u_s by one solid timestep.
       */
      solid_solver.advance_one_step();

      /*
       * Solid -> fluid.
       * Convert alpha_s to alpha_f = 1 - alpha_s and pass it to VANS.
       */
      update_void_fraction_from_solid();


      assemble_fluid_drag_exchange_rhs();

      this->pcout << "||fluid_drag_rhs|| = "
            << fluid_drag_rhs.l2_norm()
            << std::endl;

      this->set_euler_euler_drag_rhs(fluid_drag_rhs);

      /*
       * Fluid solve.
       * This uses the updated void fraction field.
       */
      this->iterate();

      //   if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
      //     {
      //       const Point<dim> probe =
      //         (dim == 2 ? Point<dim>(0.12, 0.5) : Point<dim>(0.12, 0.5,
      //         0.5));

      //       Vector<double> value(dim + 1);

      //       VectorTools::point_value(this->dof_handler,
      //                                this->present_solution,
      //                                probe,
      //                                value);

      //       this->pcout << "Fluid AFTER solve at probe " << probe << " = ";
      //       for (unsigned int d = 0; d < dim; ++d)
      //         this->pcout << value[d] << " ";
      //       this->pcout << std::endl;
      //     }

      this->postprocess(false);
      this->finish_time_step_fd();

      this->FluidDynamicsVANS<dim>::monitor_mass_conservation();
    }

  solid_solver.finalize();

  this->finish_simulation();
}


// Explicit instantiations
template class EulerEulerOneWay<2>;
template class EulerEulerOneWay<3>;
