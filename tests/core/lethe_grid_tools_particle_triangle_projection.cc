// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Lethe
#include <core/lethe_grid_tools.h>

#include <../tests/tests.h>

#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

/**
 * @brief Returns the name of a particle-triangle contact type.
 *
 * @param[in] contact_indicator Type of contact between a particle and a
 * triangle.
 *
 * @return Name of the contact type.
 */
std::string
contact_name(
  const LetheGridTools::ParticleTriangleContactIndicator contact_indicator)
{
  switch (contact_indicator)
    {
      case LetheGridTools::ParticleTriangleContactIndicator::face_contact:
        return "face";
      case LetheGridTools::ParticleTriangleContactIndicator::edge_contact:
        return "edge";
      case LetheGridTools::ParticleTriangleContactIndicator::vertex_contact:
        return "vertex";
      case LetheGridTools::ParticleTriangleContactIndicator::no_contact:
        return "none";
    }
  return "unknown";
}

/**
 * @brief Writes the components of a vector with six decimals. Adding zero to
 * each component writes negative zeros as zeros.
 *
 * @param[in] vector Vector to write.
 *
 * @return The components of the vector separated by spaces.
 */
std::string
format_vector(const Tensor<1, 3> &vector)
{
  std::ostringstream stream;
  stream << std::fixed << std::setprecision(6);
  for (unsigned int d = 0; d < 3; ++d)
    stream << (d == 0 ? "" : " ") << vector[d] + 0.;
  return stream.str();
}

void
test()
{
  /* This test checks the projection of particle centers on a triangle by
  LetheGridTools::find_particle_triangle_projection(). The triangle lies in
  the z = 0 plane and its vertices are (0,0,0), (1,0,0) and (0,1,0).
  Particle centers are placed in the face, edge and vertex regions of the
  triangle, on both sides of its plane. Particle centers farther than the
  cutoff distance from the plane must not pass the distance check, while all
  the other ones must be projected on the closest point of the triangle.
  */

  const std::vector<Point<3>> triangle = {Point<3>(0., 0., 0.),
                                          Point<3>(1., 0., 0.),
                                          Point<3>(0., 1., 0.)};

  const auto triangle_data =
    LetheGridTools::prepare_triangle_projection_data(triangle);

  // Largest distance between a particle center and the plane of the triangle
  const double cutoff_distance = 0.1;

  const std::vector<std::pair<std::string, Point<3>>> particles = {
    {"  Face region, above the plane", Point<3>(0.25, 0.25, 0.05)},
    {"  Face region, under the plane", Point<3>(0.25, 0.25, -0.05)},
    {"  Face region, just closer than the cutoff",
     Point<3>(0.25, 0.25, 0.0999)},
    {"  Face region, just farther than the cutoff",
     Point<3>(0.25, 0.25, 0.1001)},
    {"  Face region, farther than the cutoff under the plane",
     Point<3>(0.25, 0.25, -0.15)},
    {"  Edge region of edge 1-2", Point<3>(0.75, 0.75, 0.05)},
    {"  Edge region of edge 0-1", Point<3>(0.5, -0.2, 0.05)},
    {"  Edge region of edge 0-2", Point<3>(-0.2, 0.5, 0.05)},
    {"  Vertex region of vertex 0", Point<3>(-0.2, -0.2, 0.05)},
    {"  Vertex region of vertex 1", Point<3>(1.3, -0.1, 0.05)},
    {"  Vertex region of vertex 2", Point<3>(-0.1, 1.3, 0.05)},
    {"  In the plane, far from the triangle", Point<3>(2., 2., 0.)}};

  for (const auto &[description, particle_position] : particles)
    {
      const auto [pass_distance_check,
                  projection_point,
                  normal_vector,
                  contact_indicator] =
        LetheGridTools::find_particle_triangle_projection(triangle_data,
                                                          particle_position,
                                                          cutoff_distance);

      deallog << description << std::endl;
      deallog << "  Distance check: " << (pass_distance_check ? "pass" : "fail")
              << ", contact type: " << contact_name(contact_indicator)
              << std::endl;
      deallog << "  Projection point: " << format_vector(projection_point)
              << std::endl;
      deallog << "  Normal vector: " << format_vector(normal_vector)
              << std::endl;
    }
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
