// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Generalized mortar arrangement with a graded (non-uniform) rotor mesh
 * against a uniform stator mesh on a circular interface. The non-uniform cell
 * widths produce faces cut by 1, 2 and 3 faces on the other side, exercising
 * the arbitrary-N mortar segment arrangement with irregular breakpoints.
 *
 * Same checks as manager_06: per-face mortar count and indices, conservation of
 * the quadrature weights, and the rotor/stator matching invariant.
 */

#include <core/mortar_coupling_manager.h>

#include <map>

#include "./tests.h"

// Test-only manager that injects arbitrary rotor/stator breakpoints.
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

  // Graded rotor: cumulative fractions of 2*pi with widening cells. Stator:
  // 8 equal cells.
  const std::vector<double> rotor_fractions =
    {0.0, 0.05, 0.15, 0.30, 0.50, 0.75, 1.0};
  std::vector<double> rotor_bp;
  for (unsigned int i = 0; i + 1 < rotor_fractions.size(); ++i)
    rotor_bp.push_back(rotor_fractions[i] * 2.0 * numbers::PI);

  const unsigned int  n_stator = 8;
  std::vector<double> stator_bp;
  for (unsigned int k = 0; k < n_stator; ++k)
    stator_bp.push_back(k * 2.0 * numbers::PI / n_stator);

  const std::vector<unsigned int> n_subdivisions = {
    static_cast<unsigned int>(rotor_bp.size()), 1};

  const CustomMortarManagerCircle<dim> manager(n_subdivisions,
                                               radius,
                                               QGauss<dim>(n_quadrature_points),
                                               stage_heights,
                                               rotor_bp,
                                               stator_bp);

  std::cout << "Total number of mortars: " << manager.get_n_total_mortars()
            << std::endl
            << std::endl;

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

          for (unsigned int s = 0; s < indices.size(); ++s)
            for (unsigned int q = 0; q < n_quadrature_points; ++q)
              (is_inner ? points_rotor : points_stator)[indices[s]]
                .push_back(points[s * n_quadrature_points + q]);
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

  double max_match_error = 0.0;
  for (const auto &[index, pts] : points_rotor)
    {
      const auto it = points_stator.find(index);
      if (it == points_stator.end())
        continue;
      const auto &pts_s = it->second;
      for (unsigned int q = 0; q < std::min(pts.size(), pts_s.size()); ++q)
        max_match_error = std::max(max_match_error, pts[q].distance(pts_s[q]));
    }
  printf("Max rotor/stator matching error: %10.2e\n", max_match_error);
}
