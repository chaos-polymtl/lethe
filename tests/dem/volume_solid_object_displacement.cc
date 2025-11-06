// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Displacing volume solid object (dim=3,spacedim=3).
 */

// Deal.II
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle.h>

// Lethe
#include <core/parameters.h>
#include <core/serial_solid.h>
#include <core/solid_base.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim, int spacedim>
void
test()
{
  // RigidSolidObject
  auto param = std::make_shared<Parameters::RigidSolidObject<spacedim>>();

  // Mesh of the solid
  param->solid_mesh.type               = Parameters::Mesh::Type::dealii;
  param->output_bool                   = false;
  param->solid_mesh.grid_type          = "hyper_cube";
  param->solid_mesh.grid_arguments     = "-0.5 : 0.5 : false";
  param->solid_mesh.initial_refinement = 0;
  param->solid_mesh.simplex            = false;
  param->solid_mesh.translation        = Tensor<1, 3>({0.5, 0.5, 0.5});
  param->solid_mesh.rotation_axis      = Tensor<1, 3>({1., 0., 0.});
  param->solid_mesh.rotation_angle     = 0;
  param->center_of_rotation            = Point<3>({0., -0.5, 0.25});

  param->output_bool = false;

  // Functions
  std::vector<double> translational_vector = {1., 0., 0.};
  std::vector<double> angular_vector       = {-0.39269908169,
                                              0.,
                                              0.}; //  - 0.125 * pi

  param->translational_velocity =
    std::make_shared<Functions::ConstantFunction<dim>>(translational_vector);
  param->angular_velocity =
    std::make_shared<Functions::ConstantFunction<dim>>(angular_vector);

  // SerialSolid
  SerialSolid<dim, spacedim> solid(param, 0);

  double t = 0.1, dt = 0.1, t_end = 1.;
  while (t < t_end)
    {
      solid.move_solid_triangulation(dt, t);
      t += dt;
    }

  // Output
  auto displacement_vector = solid.get_displacement_vector();

  for (unsigned int i = 0; i < (displacement_vector.size() / spacedim); ++i)
    {
      deallog << "Vertex " << i << " displacement: ";
      for (int j = 0; j < spacedim; ++j)
        {
          deallog << displacement_vector[i * spacedim + j] << " ";
        }
      deallog << ";" << std::endl;
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test<3, 3>();
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
