#include <solvers/euler_euler_sf.h>

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>

using namespace dealii;

template <int dim>
void
make_euler_euler_grid(parallel::distributed::Triangulation<dim> &tria,
                      const EulerEulerGridParameters<dim> &grid_parameters)
{
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            grid_parameters.subdivisions,
                                            grid_parameters.p1,
                                            grid_parameters.p2);

  tria.refine_global(grid_parameters.global_refinement);

  const double tol = 1e-12;

  const auto axis_from = [](const std::string &d) -> unsigned int {
    if (d == "x")
      return 0;
    if (d == "y")
      return 1;
    if (d == "z")
      return 2;

    AssertThrow(false, ExcMessage("direction must be x, y, or z"));
    return 0;
  };

  const unsigned int axis1 = axis_from(grid_parameters.direction1);
  // const unsigned int axis2 = axis_from(grid_parameters.direction2);

  for (const auto &cell : tria.active_cell_iterators())
    {
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        {
          if (cell->face(f)->at_boundary())
            {
              const Point<dim> fc = cell->face(f)->center();

              cell->face(f)->set_boundary_id(0);

              if (std::fabs(fc[axis1] - grid_parameters.p1[axis1]) < tol)
                cell->face(f)->set_boundary_id(1);
              else if (std::fabs(fc[axis1] - grid_parameters.p2[axis1]) < tol)
                cell->face(f)->set_boundary_id(2);
              //   else if (std::fabs(fc[axis2] - grid_parameters.p1[axis2]) <
              //   tol)
              //     cell->face(f)->set_boundary_id(1);
              //   else if (std::fabs(fc[axis2] - grid_parameters.p2[axis2]) <
              //   tol)
              //     cell->face(f)->set_boundary_id(2);
            }
        }
    }
}

template <int dim>
EulerEulerOneWay<dim>::EulerEulerOneWay(
  const EulerEulerMeshParameters<dim> &mesh_parameters,
  const SolidPhaseParameters          &solid_parameters,
  CFDDEMSimulationParameters<dim>     &fluid_parameters,
  const MPI_Comm                      &mpi_communicator,
  const bool                           verbose)
  : mesh_parameters(mesh_parameters)
  , mpi_communicator(mpi_communicator)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  , verbose(verbose)
  , solid_parameters(solid_parameters)
  , fluid_solver(fluid_parameters)
  , void_fraction_solver(fluid_solver, mpi_communicator, verbose)
{}



// template <int dim>
// void
// EulerEulerOneWay<dim>::solve_solid()
// {
//   if (verbose)
//     pcout << "Solving solid phase" << std::endl;

//   solid_solver->run();
// }

// template <int dim>
// void
// EulerEulerOneWay<dim>::solve_fluid()
// {
//   if (verbose)
//     pcout << "Solving fluid phase" << std::endl;

//   fluid_solver.solve();
// }


template <int dim>
void
EulerEulerOneWay<dim>::build_and_pass_alpha_f()
{
  if (verbose)
    pcout << "Building alpha_f and passing it to fluid solver" << std::endl;

  void_fraction_solver.set_solid_volume_fraction(
    solid_solver->get_solid_volume_fraction());

  void_fraction_solver.calculate_alpha_f();
  void_fraction_solver.pass_alpha_f_to_fluid();
}



template <int dim>
void
EulerEulerOneWay<dim>::solve_solid_one_step()
{
  if (verbose)
    pcout << "Solving solid phase one step" << std::endl;

  solid_solver->advance_one_step();
}

template <int dim>
void
EulerEulerOneWay<dim>::solve_fluid_one_step()
{
  if (verbose)
    pcout << "Solving fluid phase one step" << std::endl;

  fluid_solver.advance_one_step_external();
}

template <int dim>
void
EulerEulerOneWay<dim>::run()
{
  if (verbose)
    pcout << "Running Euler-Euler coupling" << std::endl;

  auto &shared_triangulation = fluid_solver.get_triangulation();

  if (shared_triangulation.n_global_active_cells() == 0)
    {
      if (verbose)
        pcout << "Building Euler-Euler grid" << std::endl;

      EulerEulerGridParameters<dim> grid_parameters;
      grid_parameters.p1 = mesh_parameters.p1;
      grid_parameters.p2 = mesh_parameters.p2;
      grid_parameters.subdivisions.assign(dim, 1);
      grid_parameters.subdivisions[0] = mesh_parameters.nx;
      if (dim > 1)
        grid_parameters.subdivisions[1] = mesh_parameters.ny;
      if (dim > 2)
        grid_parameters.subdivisions[2] = mesh_parameters.nz;
      grid_parameters.global_refinement = mesh_parameters.global_refinement;
      grid_parameters.direction1        = mesh_parameters.direction1;
      grid_parameters.direction2        = mesh_parameters.direction2;

      make_euler_euler_grid(shared_triangulation, grid_parameters);
    }

  solid_solver = std::make_unique<SolidPhaseSolver<dim>>(solid_parameters,
                                                         shared_triangulation,
                                                         mpi_communicator);

  solid_solver->setup();
  fluid_solver.setup_for_external_stepping();

  while (!solid_solver->finished())
    {
      pcout << "Coupled step = " << solid_solver->get_step_number()
            << std::endl;

      void_fraction_solver.set_solid_volume_fraction(
        solid_solver->get_solid_volume_fraction());

      void_fraction_solver.calculate_alpha_f();
      void_fraction_solver.pass_alpha_f_to_fluid();

      pcout << "Before fluid step: "
            << shared_triangulation.n_global_active_cells() << std::endl;

      const bool fluid_ok = fluid_solver.advance_one_step_external();

      pcout << "After fluid step:  "
            << shared_triangulation.n_global_active_cells() << std::endl;

      AssertThrow(
        fluid_ok,
        ExcMessage(
          "Fluid finished before solid. Check time step and end time."));

      solid_solver->set_fluid_velocity_field(
        fluid_solver.get_fluid_dof_handler(),
        fluid_solver.get_fluid_mapping(),
        fluid_solver.get_fluid_solution());

      solve_solid_one_step();
    }

  solid_solver->finalize();
  fluid_solver.finish_external_stepping();
}

template void
make_euler_euler_grid<2>(parallel::distributed::Triangulation<2> &tria,
                         const EulerEulerGridParameters<2> &grid_parameters);

template void
make_euler_euler_grid<3>(parallel::distributed::Triangulation<3> &tria,
                         const EulerEulerGridParameters<3> &grid_parameters);


template class EulerEulerOneWay<2>;
template class EulerEulerOneWay<3>;