// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II includes
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function_signed_distance.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/interface_tools.h>
#include <core/parameters.h>
#include <core/simulation_control.h>

#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  /* This test checks the InterfaceTools::SignedDistanceSolver class from within
  an arbitrary structure that replicates the usual physical solver structure,
  denoted by the background prefix. This test performs verification of the
  signed distance computations for a disk in 2D and a sphere in 3D. It checks
  the L2 and Linf norms of the error on the signed distance field in the narrow
  band.
  */
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  ConditionalOStream pcout(std::cout,
                           (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                            0));



  // Discretize the domain as in an arbitrary background solver
  std::shared_ptr<FiniteElement<dim>> background_fe =
    std::make_shared<FE_Q<dim>>(1);
  std::shared_ptr<MappingQ<dim>> background_mapping =
    std::make_shared<MappingQ<dim>>(background_fe->degree);

  unsigned int number_of_mesh_ref = 3;
  unsigned int initial_ref;

  // Initialize the number of meshes for the spatial convergence analysis
  if (dim == 2)
    initial_ref = 8;

  if (dim == 3)
    initial_ref = 3;

  // Intialize spatial convergence table
  ConvergenceTable error_table;

  for (unsigned int mesh_ref = 0; mesh_ref < number_of_mesh_ref; mesh_ref++)
    {
      // Define min and max refinement levels
      unsigned int max_ref = initial_ref + mesh_ref;
      unsigned int min_ref = max_ref - 2;

      // Triangulation (as a shared_ptr to reproduce an arbitrary background
      // solver architecture)
      std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
        background_triangulation =
          std::make_shared<parallel::distributed::Triangulation<dim>>(
            mpi_communicator,
            typename Triangulation<dim>::MeshSmoothing(
              Triangulation<dim>::smoothing_on_refinement |
              Triangulation<dim>::smoothing_on_coarsening));

      // Constraints for the hanging nodes
      AffineConstraints<double> constraints;

      // Definition of the unit square/cube domain
      Point<dim> p_0 = Point<dim>();
      Point<dim> p_1 = Point<dim>();
      for (int n = 0; n < dim; n++)
        {
          p_0[n] = 0.0;
          p_1[n] = 1.0;
        }
      GridGenerator::hyper_rectangle(*background_triangulation, p_0, p_1);
      background_triangulation->refine_global(min_ref);

      // Initialization of the DoFHandler and relevant IndexSets and vectors
      DoFHandler<dim> background_dof_handler =
        DoFHandler(*background_triangulation);
      background_dof_handler.distribute_dofs(*background_fe);

      IndexSet background_locally_owned_dofs =
        background_dof_handler.locally_owned_dofs();
      IndexSet background_locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(background_dof_handler);

      GlobalVectorType background_level_set(background_locally_owned_dofs,
                                            background_locally_relevant_dofs,
                                            mpi_communicator);
      GlobalVectorType background_level_set_owned(background_locally_owned_dofs,
                                                  mpi_communicator);

      // Define the radius, center, and signed ditance function of the
      // disk/sphere
      const double sphere_radius = 0.23;
      Point<dim>   sphere_center = Point<dim>();
      for (int n = 0; n < dim; n++)
        {
          sphere_center[n] = 0.51;
        }
      auto sphere_function =
        Functions::SignedDistance::Sphere<dim>(sphere_center, sphere_radius);

      // Cell side length
      double h = 1.0 / std::pow(2.0, (max_ref));

      // Define relevant quantities: narrow band half-size (band) and max reinit
      // distance up to which the signed distance is computed
      double band                          = 8.0 * h;
      double max_reinitialization_distance = 2 * band;

      // Do max_ref - min_ref + 1 cycle of adapt refinement to have a
      // well-defined disk/sphere
      for (unsigned int i = 0; i < (max_ref - min_ref + 1); ++i)
        {
          // Set the background level-set field of the disk/sphere
          VectorTools::interpolate(*background_mapping,
                                   background_dof_handler,
                                   sphere_function,
                                   background_level_set_owned);

          background_level_set = background_level_set_owned;

          // Define a phase indicator for the adaptive refinement to reproduce
          // what is done in the CLS solver
          GlobalVectorType background_phase_indicator(
            background_locally_owned_dofs,
            background_locally_relevant_dofs,
            mpi_communicator);
          GlobalVectorType background_phase_indicator_owned(
            background_locally_owned_dofs, mpi_communicator);

          for (auto p : background_locally_relevant_dofs)
            {
              double dist = background_level_set(p);
              if (std::abs(dist) > max_reinitialization_distance)
                dist = max_reinitialization_distance * dist / std::abs(dist);

              background_phase_indicator_owned(p) =
                0.5 - 0.5 * tanh(dist / (band));
            }

          background_phase_indicator = background_phase_indicator_owned;

          // Estimate the error for the mesh refinement base on the phase
          // indicator gradient
          Vector<float> estimated_error_per_cell(
            background_triangulation->n_active_cells());

          KellyErrorEstimator<dim>::estimate(
            *background_mapping,
            background_dof_handler,
            QGauss<dim - 1>(background_fe->degree + 1),
            typename std::map<types::boundary_id,
                              const Function<dim, double> *>(),
            background_phase_indicator,
            estimated_error_per_cell,
            ComponentMask(),
            nullptr,
            0,
            background_triangulation->locally_owned_subdomain());

          // Tag the cells to be refined of coarsen
          parallel::distributed::GridRefinement::
            refine_and_coarsen_fixed_fraction(*background_triangulation,
                                              estimated_error_per_cell,
                                              0.99,
                                              0.001);

          if (background_triangulation->n_levels() > max_ref)
            for (auto &cell :
                 background_triangulation->active_cell_iterators_on_level(
                   max_ref))
              cell->clear_refine_flag();

          for (auto &cell :
               background_triangulation->active_cell_iterators_on_level(
                 min_ref))
            cell->clear_coarsen_flag();

          // Adapt the mesh
          background_triangulation->prepare_coarsening_and_refinement();

          background_triangulation->execute_coarsening_and_refinement();

          background_dof_handler.distribute_dofs(*background_fe);

          background_locally_owned_dofs =
            background_dof_handler.locally_owned_dofs();
          background_locally_relevant_dofs =
            DoFTools::extract_locally_relevant_dofs(background_dof_handler);

          background_level_set.reinit(background_locally_owned_dofs,
                                      background_locally_relevant_dofs,
                                      mpi_communicator);
          background_level_set_owned.reinit(background_locally_owned_dofs,
                                            mpi_communicator);

          // Interpolate the signed distance function on the new mesh
          VectorTools::interpolate(*background_mapping,
                                   background_dof_handler,
                                   sphere_function,
                                   background_level_set_owned);

          background_level_set = background_level_set_owned;
        }

      // Apply hanging node constraints
      constraints.clear();
      constraints.reinit(background_dof_handler.locally_owned_dofs(),
                         DoFTools::extract_locally_relevant_dofs(
                           background_dof_handler));
      DoFTools::make_hanging_node_constraints(background_dof_handler,
                                              constraints);

      constraints.close();

      // Instanciate the SignedDistanceSolver as it would as a member of the
      // background solver class
      std::shared_ptr<
        InterfaceTools::SignedDistanceSolver<dim, GlobalVectorType>>
        signed_distance_solver = std::make_shared<
          InterfaceTools::SignedDistanceSolver<dim, GlobalVectorType>>(
          background_triangulation,
          background_fe,
          max_reinitialization_distance,
          0,
          1.0,
          Parameters::Verbosity::quiet);

      signed_distance_solver->setup_dofs();

      auto background_level_set_0 = background_level_set;

      GlobalVectorType background_signed_distance(
        background_locally_owned_dofs,
        background_locally_relevant_dofs,
        mpi_communicator);

      // Initialization the signed_distance_solver
      signed_distance_solver->set_level_set_from_background_mesh(
        background_dof_handler, background_level_set);

      // Solve the signed_distance field.
      signed_distance_solver->solve();

      // Get the signed_distance field from the SignedDistanceSolver
      auto &signed_distance = signed_distance_solver->get_signed_distance();

      // Interpolate the signed_distance solution from the
      // SignedDistanceSolver DoFHandler to the main solver DoFHandler.
      TrilinosWrappers::MPI::Vector tmp_background_signed_distance(
        background_locally_owned_dofs, mpi_communicator);

      FETools::interpolate(signed_distance_solver->dof_handler,
                           signed_distance,
                           background_dof_handler,
                           constraints,
                           tmp_background_signed_distance);

      constraints.close();

      // Exchange
      background_signed_distance = tmp_background_signed_distance;

      background_level_set = background_signed_distance;

      // Compute the volume enclosed by the 0 level-set
      double volume_1, surface_1;
      std::tie(volume_1, surface_1) =
        InterfaceTools::integrate_volume_and_surface(background_dof_handler,
                                                     *background_fe,
                                                     background_level_set,
                                                     0);

      // Compute the analytical volume
      double exact_volume;
      if (dim == 3)
        exact_volume = 4.0 * M_PI * std::pow(sphere_radius, 3) / 3.0;
      if (dim == 2)
        exact_volume = M_PI * std::pow(sphere_radius, 2);

      // Instanciate objects for error norms computations in the narrow band
      QGauss<dim> quadrature(background_fe->degree + 1);

      FEValues<dim>  fe_values(*background_mapping,
                              *background_fe,
                              quadrature,
                              update_values | update_JxW_values |
                                update_quadrature_points);
      Vector<double> cell_dof_values(background_fe->dofs_per_cell);
      Vector<double> cell_dof_values_analytical(background_fe->dofs_per_cell);

      std::vector<double> cell_quadrature_values(fe_values.n_quadrature_points);

      double error_L2   = 0.0;
      double error_Linf = 0.0;
      double volume     = 0.0;

      // Loop on the owned cells to indentify the ones in the narrow band
      for (const auto &cell : background_dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              fe_values.reinit(cell);

              cell->get_dof_values(background_signed_distance,
                                   cell_dof_values.begin(),
                                   cell_dof_values.end());
              cell->get_dof_values(background_level_set_0,
                                   cell_dof_values_analytical.begin(),
                                   cell_dof_values_analytical.end());

              bool in_band = true;

              // Check if the cell is within the narrow band
              for (unsigned int j = 0; j < background_fe->dofs_per_cell; ++j)
                {
                  if (std::abs(cell_dof_values_analytical[j]) >=
                      max_reinitialization_distance)
                    {
                      in_band = false;
                      break;
                    }
                }

              if (!in_band)
                continue;

              // Get the quadrature point values for error norms integration
              fe_values.get_function_values(background_signed_distance,
                                            cell_quadrature_values);

              for (const unsigned int q : fe_values.quadrature_point_indices())
                {
                  double extact_distance =
                    sphere_function.value(fe_values.quadrature_point(q));

                  // Compute error contribution (squared norm for L2)
                  error_L2 += (cell_quadrature_values[q] - extact_distance) *
                              (cell_quadrature_values[q] - extact_distance) *
                              fe_values.JxW(q);

                  volume += fe_values.JxW(q);

                  error_Linf = std::max(std::abs(cell_quadrature_values[q] -
                                                 extact_distance),
                                        error_Linf);
                }
            }
        }

      // Exchange
      error_L2   = Utilities::MPI::sum(error_L2, mpi_communicator);
      volume     = Utilities::MPI::sum(volume, mpi_communicator);
      error_Linf = Utilities::MPI::max(error_Linf, mpi_communicator);

      error_L2 = sqrt(error_L2 / volume);

      // Fill convergence table
      error_table.add_value("ref. level", max_ref);
      error_table.add_value("cells",
                            background_triangulation->n_global_active_cells());
      error_table.add_value("dofs", background_dof_handler.n_dofs());
      error_table.add_value("error_volume", std::abs(volume_1 - exact_volume));
      error_table.add_value("error_distance_L2", error_L2);
      error_table.add_value("error_distance_Linf", error_Linf);
    }

  // Setup printing format of the table
  error_table.set_precision("error_volume", 3);
  error_table.set_scientific("error_volume", true);

  error_table.set_precision("error_distance_L2", 3);
  error_table.set_scientific("error_distance_L2", true);

  error_table.set_precision("error_distance_Linf", 3);
  error_table.set_scientific("error_distance_Linf", true);

  error_table.omit_column_from_convergence_rate_evaluation("ref. level");
  error_table.omit_column_from_convergence_rate_evaluation("cells");
  error_table.omit_column_from_convergence_rate_evaluation("dofs");

  // Compute the rate of convergence
  error_table.evaluate_all_convergence_rates(
    ConvergenceTable::reduction_rate_log2);

  // Print table
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));
  if (this_mpi_process == 0)
    {
      std::cout << dim << "D spatial convergence study" << std::endl;
      error_table.write_text(std::cout);
    }
}

int
main(int argc, char *argv[])
{
  try
    {
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
