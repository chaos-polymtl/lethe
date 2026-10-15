// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test the wall representation used by the particle-wall contact of the
 * sharp immersed boundary DEM. The wall of each boundary id is represented by
 * the plane tangent to the boundary at the point of the boundary closest to the
 * particle.
 *
 * A sphere at rest overlaps the curved boundary of a ball, which has a single
 * boundary id. The contact force must be along the inward normal of the wall
 * at the boundary quadrature point closest to the particle, which is found here
 * by a search over all the boundary faces. The result must not depend on the
 * number of processes.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <boost/serialization/utility.hpp>

// Lethe
#include <core/ib_particle.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <fem-dem/ib_particles_dem.h>

// Tests
#include <../tests/fem-dem/ib_particles_dem_test_utilities.h>

#include <../tests/tests.h>

#include <cfloat>
#include <utility>
#include <vector>

using namespace dealii;

template <int dim>
void
test()
{
  const MPI_Comm mpi_communicator = MPI_COMM_WORLD;
  const double   tolerance        = 1e-10;

  parallel::distributed::Triangulation<dim> triangulation(mpi_communicator);
  GridGenerator::hyper_ball(triangulation);
  triangulation.refine_global((dim == 2) ? 3 : 2);

  const FESystem<dim> fe(FE_Q<dim>(1), dim + 1);
  DoFHandler<dim>     dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  const MappingQ1<dim>  mapping;
  const QGauss<dim - 1> face_quadrature(2);

  auto parameters =
    make_ib_particles_parameters<dim>((dim == 2) ? "0; 0" : "0; 0; 0");

  // The sphere overlaps the boundary of the ball, away from the symmetry
  // planes of the mesh so that the closest boundary point is unique.
  Point<dim> position;
  position[0] = 0.2;
  position[1] = -0.8;
  if constexpr (dim == 3)
    position[2] = 0.1;
  std::vector<IBParticle<dim>> particles{
    make_ib_sphere<dim>(0, 0.25, position)};

  IBParticlesDEM<dim> ib_dem;
  ib_dem.initialize(parameters, nullptr, nullptr, mpi_communicator, particles);
  ib_dem.update_particles_boundary_contact(particles,
                                           dof_handler,
                                           face_quadrature,
                                           mapping);

  std::vector<Tensor<1, 3>> contact_force(particles.size());
  std::vector<Tensor<1, 3>> contact_torque(particles.size());
  ib_dem.calculate_pw_contact_force(1e-5, contact_force, contact_torque);

  // Find the boundary quadrature point closest to the particle, with the
  // inward normal of the wall at this point, over all the processes.
  FEFaceValues<dim>                 fe_face_values(mapping,
                                   fe,
                                   face_quadrature,
                                   update_quadrature_points |
                                     update_normal_vectors);
  std::pair<double, Tensor<1, dim>> closest_point(DBL_MAX, Tensor<1, dim>());
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;
      for (const auto face : cell->face_indices())
        {
          if (!cell->face(face)->at_boundary())
            continue;
          fe_face_values.reinit(cell, face);
          for (const unsigned int q : fe_face_values.quadrature_point_indices())
            {
              const double level_set =
                particles[0].get_levelset(fe_face_values.quadrature_point(q));
              if (level_set < closest_point.first)
                closest_point = {level_set, -fe_face_values.normal_vector(q)};
            }
        }
    }
  const auto all_closest_points =
    Utilities::MPI::all_gather(mpi_communicator, closest_point);
  for (const auto &process_closest_point : all_closest_points)
    if (process_closest_point.first < closest_point.first)
      closest_point = process_closest_point;

  const Tensor<1, 3> inward_normal = tensor_nd_to_3d(closest_point.second);
  const double       force_norm    = contact_force[0].norm();
  const bool         in_contact    = force_norm > 0.;
  const bool         force_is_normal =
    in_contact &&
    (contact_force[0] / force_norm - inward_normal).norm() < tolerance;

  deallog << "Dimension " << dim << std::endl;
  deallog << "Particle is in contact with the wall: "
          << (in_contact ? "true" : "false") << std::endl;
  deallog
    << "Contact force is normal to the wall at the closest boundary point: "
    << (force_is_normal ? "true" : "false") << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      mpi_initlog();
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
