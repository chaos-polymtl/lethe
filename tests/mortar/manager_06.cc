// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Generalized mortar arrangement with unequal cell counts on the two
 * sides of a circular interface (rotor: 6 uniform cells, stator: 10 uniform
 * cells). This produces faces that are cut by 2 and by 3 faces on the other
 * side, exercising the arbitrary-N mortar segment arrangement.
 *
 * The test checks:
 *  - the number of mortars covered by each face (1, 2 or 3),
 *  - the global mortar indices of each face,
 *  - conservation: the sum of the quadrature weights over the rotor faces and
 *    over the stator faces both equal 2*pi*radius,
 *  - the matching invariant: the physical quadrature points of a given global
 *    mortar segment are identical whether the segment is reached from the rotor
 *    or the stator side.
 */

#include <core/mortar_coupling_manager.h>

#include <map>

#include "./tests.h"

// Test-only manager that lets us inject arbitrary (possibly non-uniform and
// unequal) rotor/stator breakpoints directly into the arrangement.
template <int dim>
class CustomMortarManagerCircle : public MyMortarManagerCircle<dim>
{
public:
  template <int dim2>
  CustomMortarManagerCircle(const std::vector<unsigned int> &n_subdivisions,
                            const std::vector<double>       &radius,
                            const Quadrature<dim2>          &quadrature,
                            const std::vector<double>       &stage_heights,
                            const std::vector<double>       &rotor_bp,
                            const std::vector<double>       &stator_bp)
    : MyMortarManagerCircle<dim>(n_subdivisions,
                                 radius,
                                 quadrature,
                                 0.0,
                                 stage_heights)
  {
    this->build_arrangement(rotor_bp, stator_bp);
  }
};

int
main()
{
  const unsigned int        dim                 = 2;
  const unsigned int        n_quadrature_points = 3;
  const std::vector<double> radius              = {1.0, 1.0};
  const std::vector<double> stage_heights       = {0.0, 1.0};

  // Rotor: 6 equal cells; stator: 10 equal cells.
  const unsigned int  n_rotor  = 6;
  const unsigned int  n_stator = 10;
  std::vector<double> rotor_bp, stator_bp;
  for (unsigned int k = 0; k < n_rotor; ++k)
    rotor_bp.push_back(k * 2.0 * numbers::PI / n_rotor);
  for (unsigned int k = 0; k < n_stator; ++k)
    stator_bp.push_back(k * 2.0 * numbers::PI / n_stator);

  const std::vector<unsigned int> n_subdivisions = {n_rotor, 1};

  const CustomMortarManagerCircle<dim> manager(n_subdivisions,
                                               radius,
                                               QGauss<dim>(n_quadrature_points),
                                               stage_heights,
                                               rotor_bp,
                                               stator_bp);

  std::cout << "Total number of mortars: " << manager.get_n_total_mortars()
            << std::endl
            << std::endl;

  // Cyclic midpoints of a breakpoint list give the face centers.
  const auto face_centers = [](const std::vector<double> &bp) {
    std::vector<double> centers;
    for (unsigned int i = 0; i < bp.size(); ++i)
      {
        const double lo = bp[i];
        const double hi = (i + 1 < bp.size()) ? bp[i + 1] : bp[0] + 2 * numbers::PI;
        centers.push_back(0.5 * (lo + hi));
      }
    return centers;
  };

  // Collect, per global mortar index, the physical quadrature points seen from
  // each side, to verify the matching invariant.
  std::map<unsigned int, std::vector<Point<dim>>> points_rotor, points_stator;

  const auto process =
    [&](const std::vector<double> &centers, const bool is_inner) {
      double sum_weights = 0.0;

      for (const double angle : centers)
        {
          const auto p = radius_to_point<dim>(radius[0], angle);

          const auto indices = manager.get_mortar_indices(p, is_inner);
          const auto weights = manager.get_weights(p, is_inner);
          const auto points  = manager.get_points(p, is_inner);

          std::cout << (is_inner ? "rotor " : "stator")
                    << " face, n_mortars = "
                    << manager.get_n_mortars(p, is_inner) << ", indices:";
          for (const auto i : indices)
            std::cout << " " << i;
          std::cout << std::endl;

          for (const double w : weights)
            sum_weights += w;

          // Store the physical points segment by segment.
          for (unsigned int s = 0; s < indices.size(); ++s)
            for (unsigned int q = 0; q < n_quadrature_points; ++q)
              {
                const auto pt = points[s * n_quadrature_points + q];
                (is_inner ? points_rotor : points_stator)[indices[s]]
                  .push_back(pt);
              }
        }

      return sum_weights;
    };

  const double sum_rotor  = process(face_centers(rotor_bp), true);
  std::cout << std::endl;
  const double sum_stator = process(face_centers(stator_bp), false);
  std::cout << std::endl;

  printf("Weight sum rotor : %10.6f\n", sum_rotor);
  printf("Weight sum stator: %10.6f\n", sum_stator);
  printf("Expected 2*pi*r  : %10.6f\n", 2.0 * numbers::PI * radius[0]);

  // Matching invariant: each global mortar segment must have the same physical
  // quadrature points whether reached from the rotor or the stator side.
  double max_match_error = 0.0;
  for (const auto &[index, pts] : points_rotor)
    {
      const auto it = points_stator.find(index);
      if (it == points_stator.end())
        continue;
      const auto &pts_s = it->second;
      for (unsigned int q = 0; q < std::min(pts.size(), pts_s.size()); ++q)
        max_match_error =
          std::max(max_match_error, pts[q].distance(pts_s[q]));
    }
  printf("Max rotor/stator matching error: %10.2e\n", max_match_error);
}
