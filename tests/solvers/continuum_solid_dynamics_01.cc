// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the solid phase (Euler-Euler model) auxiliary physics
 * through the MultiphysicsInterface. It checks that:
 * - the solid phase is reported as active when driven by a VANS fluid solver;
 * - its DoFHandler is shared through the interface and holds dim + 1
 *   components (solid velocity and solid volume fraction);
 * - its output hook provides the solid velocity and solid volume fraction;
 * - its solution survives a mesh refinement;
 * - its solution can be written to and read from a checkpoint;
 * - enabling it with a standard Navier-Stokes solver throws
 *   ContinuumSolidDynamicsWithoutVANSError.
 */

#include <core/multiphysics.h>
#include <core/output_struct.h>
#include <core/parameters.h>
#include <core/simulation_control.h>

#include <solvers/continuum_solid_dynamics.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <../tests/tests.h>

#include <memory>
#include <string>
#include <variant>
#include <vector>

using namespace dealii;

/**
 * @brief Create a distributed triangulation of the unit hyper cube.
 *
 * @param[in] n_refinements Number of global refinements.
 *
 * @return Shared pointer to the triangulation.
 */
template <int dim>
std::shared_ptr<parallel::distributed::Triangulation<dim>>
make_triangulation(const unsigned int n_refinements)
{
  auto tria = std::make_shared<parallel::distributed::Triangulation<dim>>(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening));
  GridGenerator::hyper_cube(*tria, 0., 1.);
  tria->refine_global(n_refinements);
  return tria;
}

template <int dim>
void
test()
{
  SimulationParameters<dim>     solver_parameters;
  ParameterHandler              dummy_handler;
  Parameters::SizeOfSubsections size_of_subsections;
  size_of_subsections.boundary_conditions = 1;
  size_of_subsections.manifolds           = 1;

  solver_parameters.declare(dummy_handler, size_of_subsections);
  solver_parameters.parse(dummy_handler);

  solver_parameters.multiphysics.fluid_dynamics           = true;
  solver_parameters.multiphysics.continuum_solid_dynamics = true;

  std::shared_ptr<SimulationControl> simulation_control =
    std::make_shared<SimulationControlTransient>(
      solver_parameters.simulation_control);

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  const std::string checkpoint_file = "solid_phase_checkpoint.triangulation";

  // Registration, DoFs, output and checkpoint writing under VANS
  {
    auto                       tria = make_triangulation<dim>(2);
    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout,
                                            FluidDynamicsFormulation::VANS);

    deallog << "Active physics under VANS (expected: fluid, solid phase)"
            << std::endl;
    for (const auto &iphys : multiphysics.get_active_physics())
      deallog << int(iphys) << std::endl;

    multiphysics.setup_dofs();
    multiphysics.set_initial_conditions();

    const DoFHandler<dim> &solid_dof_handler =
      multiphysics.get_dof_handler(PhysicsID::continuum_solid_dynamics);
    deallog << "Solid phase components: "
            << solid_dof_handler.get_fe().n_components() << std::endl;
    deallog << "Solid phase DoFs: " << solid_dof_handler.n_dofs() << std::endl;

    // Output hook
    const auto output_structs = multiphysics.gather_output_hook_global_vector();
    deallog << "Number of output structures: " << output_structs.size()
            << std::endl;
    if (const auto *solution_output =
          std::get_if<OutputStructSolution<dim, GlobalVectorType>>(
            &output_structs[0]))
      for (const auto &name : solution_output->solution_names)
        deallog << "Output field: " << name << std::endl;

    // Checkpoint writing
    multiphysics.write_checkpoint();
    tria->save(checkpoint_file);
  }

  // Checkpoint reading on a new triangulation
  {
    auto tria = make_triangulation<dim>(0);
    tria->load(checkpoint_file);
    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout,
                                            FluidDynamicsFormulation::VANS);
    multiphysics.setup_dofs();
    multiphysics.read_checkpoint();
    deallog << "Solid phase DoFs after restart: "
            << multiphysics.get_dof_handler(PhysicsID::continuum_solid_dynamics)
                 .n_dofs()
            << std::endl;
  }

  // Mesh adaptation
  {
    auto                       tria = make_triangulation<dim>(2);
    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout,
                                            FluidDynamicsFormulation::VANS);
    multiphysics.setup_dofs();
    multiphysics.set_initial_conditions();

    for (const auto &cell : tria->active_cell_iterators())
      if (cell->is_locally_owned())
        cell->set_refine_flag();

    multiphysics.prepare_for_mesh_adaptation();
    tria->execute_coarsening_and_refinement();
    multiphysics.setup_dofs();
    multiphysics.post_mesh_adaptation();

    deallog << "Solid phase DoFs after refinement: "
            << multiphysics.get_dof_handler(PhysicsID::continuum_solid_dynamics)
                 .n_dofs()
            << std::endl;
  }

  // Without VANS, enabling the solid phase must throw
  try
    {
      auto                       tria = make_triangulation<dim>(2);
      MultiphysicsInterface<dim> multiphysics(
        solver_parameters,
        tria,
        simulation_control,
        pcout,
        FluidDynamicsFormulation::standard);
      deallog << "Solid phase did NOT throw without VANS (unexpected)"
              << std::endl;
    }
  catch (const ContinuumSolidDynamicsWithoutVANSError &)
    {
      deallog << "Solid phase correctly threw without VANS" << std::endl;
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test<2>();
      test<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
