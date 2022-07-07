/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019- by the Lethe authors
 *
 * This file is part of the Lethe library.
 *
 * The LEthe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019
 */


// This is a template folder for prototype executables to be created when
// developing totally new Lethe functionnalities

#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>

#include "mockup_systems.h"
#include "write_distances.h"


using namespace dealii;

template <int dim>
std::vector<double>
calculate_distance(std::vector<Point<dim>> triangle,
                   std::vector<Point<dim>> particles)
{
  std::vector<double> distances(particles.size());
  const double        radius = 0.25;
  auto               &p_0    = triangle[0];
  auto               &p_1    = triangle[1];
  auto               &p_2    = triangle[2];

  const Tensor<1, dim> e_0 = p_1 - p_0;
  const Tensor<1, dim> e_1 = p_2 - p_0;

  const Tensor<1, dim> normal      = cross_product_3d(e_0, e_1);
  const double         norm_normal = normal.norm();
  const Tensor<1, dim> unit_normal = normal / norm_normal;

  const double a   = e_0.norm_square();
  const double b   = scalar_product(e_0, e_1);
  const double c   = e_1.norm_square();
  const double det = a * c - b * b;


  // Pre-allocation for speed
  Tensor<1, dim> vector_to_plane;
  Point<dim>     pt_in_triangle;

  unsigned int k = 0;
  for (auto &part : particles)
    {
      vector_to_plane         = p_0 - part;
      double distance_squared = scalar_product(vector_to_plane, unit_normal);

      // If the particle is too far from the plane, set distance squared as an
      // arbitrary distance and continue
      if (distance_squared > (radius * radius))
        {
          distances[k] = std::sqrt(distance_squared);
          ++k;
          continue;
        }

      // Otherwise, do the full calculation taken from Eberly 2003
      const double d = scalar_product(e_0, vector_to_plane);
      const double e = scalar_product(e_1, vector_to_plane);

      // Calculate necessary values;
      double s = b * e - c * d;
      double t = b * d - a * e;
      // std::cout << "s " << s << " t " << t << std::endl;

      // const double f = vector_to_plane.norm_square();
      if (s + t <= det)
        {
          if (s < 0)
            {
              if (t < 0)
                {
                  // Region 4
                  if (d < 0)
                    {
                      t = 0;
                      if (-d >= a)
                        s = 1;
                      else
                        s = -d / a;
                    }
                  else
                    {
                      s = 0;
                      if (e >= 0)
                        t = 0;
                      else if (-e >= c)
                        t = 1;
                      else
                        t = e / c;
                    }
                }
              else
                {
                  // Region 3
                  s = 0;
                  if (e >= 0)
                    t = 0;
                  else if (-e >= c)
                    t = 1;
                  else
                    t = -e / c;
                }
            }
          else if (t < 0)
            {
              // Region 5
              t = 0;
              if (d >= 0)
                s = 0;
              else if (-d >= a)
                s = 1;
              else
                s = -d / a;
            }
          else
            {
              // Region 0
              const double inv_det = 1. / det;
              s *= inv_det;
              t *= inv_det;
            }
        }
      else
        {
          if (s < 0)
            {
              // Region 2
              const double tmp0 = b + d;
              const double tmp1 = c + e;
              if (tmp1 > tmp0)
                {
                  const double numer = tmp1 - tmp0;
                  const double denom = a - 2 * b + c;
                  if (numer >= denom)
                    s = 1;
                  else
                    s = numer / denom;

                  t = 1 - s;
                }
              else
                {
                  s = 0;
                  if (tmp1 <= 0)
                    t = 1;
                  else if (e >= 0)
                    t = 0;
                  else
                    t = -e / c;
                }
            }
          else if (t < 0)
            {
              // Region 6
              const double tmp0 = b + e;
              const double tmp1 = a + d;
              if (tmp1 > tmp0)
                {
                  const double numer = tmp1 - tmp0;
                  const double denom = a - 2 * b + c;
                  if (numer >= denom)
                    t = 1;
                  else
                    t = numer / denom;
                  s = 1 - t;
                }
              else
                {
                  t = 0;
                  if (tmp1 <= 0)
                    s = 1;
                  else if (d >= 0)
                    s = 0;
                  else
                    s = -d / a;
                }
            }
          else
            {
              // Region 1
              const double numer = (c + e) - (b + d);
              if (numer <= 0)
                s = 0;
              else
                {
                  const double denom = a - 2 * b + c;
                  if (numer >= denom)
                    s = 1;
                  else
                    s = numer / denom;
                }
              t = 1 - s;
            }
        }

      pt_in_triangle = p_0 + s * e_0 + t * e_1;
      // std::cout << "pt " << part << " s : " << s << " t : " << t
      //           << " pt in triangle " << pt_in_triangle << std::endl;

      distances[k] = pt_in_triangle.distance(part);
      ++k;
    }
  return distances;
}

int
main()
{
  try
    {
      // Testing case 0
      TimerOutput timer(std::cout,
                        TimerOutput::summary,
                        TimerOutput::wall_times);

      {
        auto case_0 = generate_case_0<3>(21, 2);
        timer.enter_subsection("Case 0");
        auto distances = calculate_distance(case_0.first, case_0.second);
        timer.leave_subsection();
        write_distances(case_0.second, distances, "case_0_z.dat");
      }
      {
        auto case_0 = generate_case_0<3>(21, 1);
        timer.enter_subsection("Case 0");
        auto distances = calculate_distance(case_0.first, case_0.second);
        timer.leave_subsection();
        write_distances(case_0.second, distances, "case_0_y.dat");
      }

      {
        unsigned int n_times = 1e5;
        for (unsigned int t = 0; t < n_times; ++t)
          {
            auto case_1 = generate_case_1<3>(1000, t);
            timer.enter_subsection("Case 1");
            auto distances = calculate_distance(case_1.first, case_1.second);
            timer.leave_subsection();
            if (t == n_times - 1)
              write_distances(case_1.second, distances, "random.dat");
          }
      }
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
