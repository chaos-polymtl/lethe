// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Per-stage mortar arrangement on a 3D cylindrical interface, where the
 * circumferential (theta) discretization differs from one axial stage to the
 * next, independently on the rotor and stator sides. This exercises the
 * per-stage arrangement and the prefix-sum global indexing that replaces the
 * former uniform theta-by-z extrusion.
 *
 * Two axial stages are used:
 *   - stage 0 (z in [0, 0.5]): rotor 4 cells, stator 6 cells
 *   - stage 1 (z in [0.5, 1]): rotor 3 cells, stator 5 cells
 *
 * The test checks, per stage:
 *  - the number of mortars covered by each face and their global indices,
 *  - conservation: the sum of the quadrature weights over the rotor faces and
 *    over the stator faces both equal 2*pi*radius*dz_stage,
 *  - the matching invariant: the physical quadrature points of a given global
 *    mortar segment are identical whether reached from the rotor or the stator,
 *  - that the global indices are dense and contiguous, growing stage by stage.
 */

#include <core/mortar_coupling_manager.h>

#include <map>

#include "./tests.h"

// Test-only manager that lets us inject arbitrary per-stage rotor/stator
// breakpoints directly into the arrangement.
template <int dim>
class CustomMortarManagerCircle : public MyMortarManagerCircle<dim>
{
public:
  template <int dim2>
  CustomMortarManagerCircle(
    const std::vector<unsigned int>        &n_subdivisions,
    const std::vector<double>              &radius,
    const Quadrature<dim2>                 &quadrature,
    const std::vector<double>              &stage_heights,
    const std::vector<std::vector<double>> &rotor_bp_per_stage,
    const std::vector<std::vector<double>> &stator_bp_per_stage)
    : MyMortarManagerCircle<dim>(n_subdivisions,
                                 radius,
                                 quadrature,
                                 0.0,
                                 stage_heights)
  {
    this->build_arrangement(rotor_bp_per_stage, stator_bp_per_stage);
  }
};

int
main()
{
  const unsigned int        dim                 = 3;
  const unsigned int        n_quadrature_points = 2;
  const std::vector<double> radius              = {1.0, 1.0};
  const std::vector<double> stage_heights       = {0.0, 0.5, 1.0};

  // Uniform breakpoints for n equal cells.
  const auto uniform = [](const unsigned int n) {
    std::vector<double> v;
    for (unsigned int k = 0; k < n; ++k)
      v.push_back(k * 2.0 * numbers::PI / n);
    return v;
  };

  // Different theta discretization per stage, and unequal rotor/stator counts
  // (so faces are cut by 2 or 3 faces on the other side).
  std::vector<std::vector<double>> rotor_bp(2), stator_bp(2);
  rotor_bp[0]  = uniform(4);
  stator_bp[0] = uniform(6);
  rotor_bp[1]  = uniform(3);
  stator_bp[1] = uniform(5);

  const std::vector<unsigned int> n_subdivisions = {4, 2};

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
        const double hi =
          (i + 1 < bp.size()) ? bp[i + 1] : bp[0] + 2 * numbers::PI;
        centers.push_back(0.5 * (lo + hi));
      }
    return centers;
  };

  // Collect, per global mortar index, the physical quadrature points seen from
  // each side, to verify the matching invariant.
  std::map<unsigned int, std::vector<Point<dim>>> points_rotor, points_stator;

  // Track the set of global indices used, to confirm density/contiguity.
  std::set<unsigned int> all_indices;

  const auto process = [&](const std::vector<double> &centers,
                           const bool                 is_inner,
                           const double               z_center) {
    double sum_weights = 0.0;

    for (const double angle : centers)
      {
        auto p     = radius_to_point<dim>(radius[0], angle);
        p[dim - 1] = z_center;

        const auto indices = manager.get_mortar_indices(p, is_inner);
        const auto weights = manager.get_weights(p, is_inner);
        const auto points  = manager.get_points(p, is_inner);

        std::cout << (is_inner ? "rotor " : "stator")
                  << " face, n_mortars = " << manager.get_n_mortars(p, is_inner)
                  << ", indices:";
        for (const auto i : indices)
          {
            std::cout << " " << i;
            all_indices.insert(i);
          }
        std::cout << std::endl;

        for (const double w : weights)
          sum_weights += w;

        // Store the physical points segment by segment.
        for (unsigned int s = 0; s < indices.size(); ++s)
          for (unsigned int q = 0; q < manager.get_n_quadrature_points(); ++q)
            {
              const auto pt =
                points[s * manager.get_n_quadrature_points() + q];
              (is_inner ? points_rotor : points_stator)[indices[s]].push_back(
                pt);
            }
      }

    return sum_weights;
  };

  for (unsigned int stage = 0; stage < 2; ++stage)
    {
      const double z_lo   = stage_heights[stage];
      const double z_hi   = stage_heights[stage + 1];
      const double z_mid  = 0.5 * (z_lo + z_hi);
      const double dz     = z_hi - z_lo;

      std::cout << "=== Stage " << stage << " (z in [" << z_lo << ", " << z_hi
                << "]) ===" << std::endl;

      const double sum_rotor =
        process(face_centers(rotor_bp[stage]), true, z_mid);
      std::cout << std::endl;
      const double sum_stator =
        process(face_centers(stator_bp[stage]), false, z_mid);
      std::cout << std::endl;

      printf("Weight sum rotor   : %10.6f\n", sum_rotor);
      printf("Weight sum stator  : %10.6f\n", sum_stator);
      printf("Expected 2*pi*r*dz : %10.6f\n",
             2.0 * numbers::PI * radius[0] * dz);
      std::cout << std::endl;
    }

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
        max_match_error = std::max(max_match_error, pts[q].distance(pts_s[q]));
    }
  printf("Max rotor/stator matching error: %10.2e\n", max_match_error);

  // Density/contiguity: the union of all global indices must be exactly
  // {0, 1, ..., get_n_total_mortars() - 1}.
  bool dense = all_indices.size() == manager.get_n_total_mortars();
  unsigned int expected = 0;
  for (const auto i : all_indices)
    dense = dense && (i == expected++);
  std::cout << "Global indices dense and contiguous: " << (dense ? "yes" : "no")
            << std::endl;
}
