// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II includes
#include <deal.II/base/function_signed_distance.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/interface_tools.h>

#include <../tests/tests.h>

void
test()
{
  /* This test checks the computation of the volume enclosed by the level 0.1 of
  a level-set field using the InterfaceTools::compute_cell_wise_volume function.
  The level-set field of interest is the one describing a sphere. The volume is
  computed for 3 mesh refinements and the test checks the error on the volume
  and the convergence rate of the method (formally 2).
  */
  Triangulation<3> triangulation;
  DoFHandler<3>    dof_handler;
  FE_Q<3>          fe(1);
  MappingQ<3>      mapping(1);

  const Point<3> p_0 = Point<3>({0, 0, 0});
  const Point<3> p_1 = Point<3>({1, 1, 1});

  GridGenerator::hyper_rectangle(triangulation, p_0, p_1);
  triangulation.refine_global(3);

  const Point<3> sphere_center = Point<3>({0.5, 0.5, 0.5});
  const double   sphere_radius = 0.25;
  const double   iso_level     = 0.1;

  Vector<double> error_volume(3);

  // Loop for the mesh convergence study
  for (unsigned int n = 0; n < 3; n++)
    {
      dof_handler.reinit(triangulation);
      dof_handler.distribute_dofs(fe);

      dealii::TrilinosWrappers::MPI::Vector signed_distance;

      signed_distance.reinit(dof_handler.locally_owned_dofs(),
                             triangulation.get_communicator());

      // Set the level-set field of the sphere
      VectorTools::interpolate(
        mapping,
        dof_handler,
        Functions::SignedDistance::Sphere<3>(sphere_center, sphere_radius),
        signed_distance);

      // Compute the surface and volume of the sphere with the NonMatching
      // FEValues
      double volume = 0.0;

      // Compute the volume with the cell-wise routine of InterfaceTools
      FEPointEvaluation<1, 3> fe_point_evaluation(
        mapping, fe, update_jacobians | update_JxW_values);

      double volume_cell_wise_sum = 0.0;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              const unsigned int n_dofs_per_cell =
                cell->get_fe().n_dofs_per_cell();
              Vector<double> cell_dof_level_set_values(n_dofs_per_cell);

              cell->get_dof_values(signed_distance,
                                   cell_dof_level_set_values.begin(),
                                   cell_dof_level_set_values.end());

              const double level_set_correction = -iso_level;
              volume += InterfaceTools::compute_cell_wise_volume(
                fe_point_evaluation,
                cell,
                cell_dof_level_set_values,
                level_set_correction,
                cell->get_fe().degree + 1);
            }
        }

      // Analytical volume
      const double analytical_volume =
        4.0 * M_PI * std::pow(sphere_radius + iso_level, 3) / 3.0;

      // Compute and store the volume error
      error_volume[n] = abs(analytical_volume - volume);

      deallog << "The volume error for ref. lev. " << n + 3
              << " is: " << error_volume[n] << std::endl;
      triangulation.refine_global(1);
    }

  // Compute the rate of convergence
  const double convergence_order =
    log(error_volume[2] / error_volume[1]) / log(0.5);

  deallog << "The convergence rate is: " << convergence_order << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test();
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
