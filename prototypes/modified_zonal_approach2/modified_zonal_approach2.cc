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
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include "mockup_systems.h"
#include "write_distances.h"
#include <deal.II/base/timer.h>


using namespace dealii;


template <int dim>
class ParticleTriangleDistance {
public:
    ParticleTriangleDistance<dim>()
    {    }

    std::vector<double> calculate_distance(std::vector<Point<dim>> triangle, std::vector<Point<dim>> particles)
    {
  const double        radius = 0.25;
  auto               &p0    = triangle[0];
  auto               &p1    = triangle[1];
  auto               &p2    = triangle[2];
   std::vector<double> distances(particles.size());

        const Tensor<1, dim> e0 = p1 - p0;
        const Tensor<1, dim> e1 = p2 - p0;

        const Tensor<1, dim> normal      = cross_product_3d(e0, e1);
        const double         norm_normal = normal.norm();
        const Tensor<1, dim> unit_normal = normal / norm_normal;

        const double a   = e0.norm_square();
        const double b   = scalar_product(e0, e1);
        const double c   = e1.norm_square();

       // Pre-allocation for speed
       Tensor<1, dim> vector_to_plane;
       Point<dim>     pt_in_triangle;

       unsigned int k = 0;

    for (auto &part : particles)
                  {
                    vector_to_plane         = p0 - part;
                    double distance_squared = scalar_product(vector_to_plane, unit_normal);

                    // If the particle is too far from the plane, set distance squared as an
                    // arbitrary distance and continue
                    if (distance_squared > (radius * radius))
                      {
                        distances[k] = std::sqrt(distance_squared);
                        ++k;
                        continue;
                      }


                const double d = scalar_product(e0, vector_to_plane);
                const double e = scalar_product(e1, vector_to_plane);

        auto det = (a * c - b * b);
        auto s = b * e - c * d;
        auto t = b * d - a * e;


        if (s+t <= det)
        {
            if (s <0)
            {
                if (t <0)
                  { auto region4_output =  region4(a,b,c,d,e);
                    s = region4_output.first;
                    t =  region4_output.second;
                }
            else
                {
                    auto region3_output = region3(c,e);
                    s = region3_output.first;
                    t = region3_output.second;
                }
            }
            else if (t < 0)
            {
                auto region5_output = region5(a, d);
                s = region5_output.first;
                t = region5_output.second;
            }
            else
            {
                auto region0_output = region0(det, s, t);
                s = region0_output.first;
                t = region0_output.second;
            }
        }
        else
        {
        if (s < 0 )
        {
            auto region_2_output = region2(a,b,c,d,e);
            s = region_2_output.first;
            t = region_2_output.second;
        }
        else if (t < 0)
        {
            auto region6_output = region6(a,b,c,d,e);
            s = region6_output.first;
            t = region6_output.second;

        }
        else
        {
            auto region1_output = region1(a,b,c,d,e);
            s = region1_output.first;
            t = region1_output.second;
        }
    }


       Point<dim> pt_in_triangle = p0 + s * e0 + t * e1;

       distances[k] = pt_in_triangle.distance(part);
       ++k;
    }
 return distances;
}

private:


    std::pair<double, double> region0(const double &det, double &s, double &t)
    {
        const double inv_det = 1. / det;
        s *= inv_det;
        t *= inv_det;

         return std::make_pair(s,t);
    }

    std::pair<double, double> region1(const double &a, const double &b, const double &c, const double &d, const double &e)
    {
        double s = 0.0;
        double t = 0.0;
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

        return std::make_pair(s,t);
    }

    std::pair<double, double> region2(const double &a, const double &b, const double &c, const double &d, const double &e)
    {
        double t = 0.0;
        double s = 0.0;

        double tmpO = b + d;
         double tmpl = c + e;
         if   (tmpl > tmpO) {
             double numer = tmpl - tmpO;
             double denom = a - 2 * b + c;
             s = (numer >= denom? 1 : numer / denom);
             t = 1 - s;
         }
         else
         {
             s = 0.0;
             t = (tmpl <= 0 ? 1 : (e >= 0 ? 0 : -e / c ) );
         }

         return std::make_pair(s,t);
    }

    std::pair<double, double> region3(const double &c, const double &e)
    {
        double s = 0.0;
        double t = 0.0;

        if (e >= 0)
          t = 0;
        else if (-e >= c)
          t = 1;
        else
          t = -e / c;

        return std::make_pair(s,t);
    }


    std::pair<double, double> region4(const double &a, const double &b, const double &c, const double &d, const double &e)
    {
        double t = 0.0;
        double s = 0.0;

        double tmpO = b + d;
         double tmpl = c + e;
         if   (tmpl > tmpO) {
             double numer = tmpl - tmpO;
             double denom = a - 2 * b + c;
             s = (numer >= denom? 1 : numer / denom);
             t = 1 - s;
         }
         else
         {
             s = 0.0;
             t = (tmpl <= 0 ? 1 : (e >= 0 ? 0 : -e / c ) );
         }

         return std::make_pair(s,t);
    }


    std::pair<double, double> region5(const double &a, const double &d)
    {
        double s = 0.0;
        double t = 0.0;

        if (d >= 0)
          s = 0;
        else if (-d >= a)
          s = 1;
        else
          s = -d / a;

        return std::make_pair(s, t);
    }

    std::pair<double, double> region6(const double &a, const double &b, const double &c, const double &d, const double &e)
    {
        double t = 0.0;
        double s = 0.0;

        double tmpO = b + d;
         double tmpl = c + e;
         if   (tmpl > tmpO) {
             double numer = tmpl - tmpO;
             double denom = a - 2 * b + c;
             s = (numer >= denom? 1 : numer / denom);
             t = 1 - s;
         }
         else
         {
             s = 0.0;
             t = (tmpl <= 0 ? 1 : (e >= 0 ? 0 : -e / c ) );
         }

         return std::make_pair(s,t);
    }
};


int
main()
{
  try
    {
        ParticleTriangleDistance<3> distance_calculator;

        // Testing case 0
        TimerOutput timer(std::cout,
                          TimerOutput::summary,
                          TimerOutput::wall_times);

        {
          auto case_0 = generate_case_0<3>(21, 2);
          timer.enter_subsection("Case 0");
          auto distances = distance_calculator.calculate_distance(case_0.first, case_0.second);
          timer.leave_subsection();
          write_distances(case_0.second, distances, "case_0_z.dat");
        }
        {
          auto case_0 = generate_case_0<3>(21, 1);
          timer.enter_subsection("Case 0");
          auto distances = distance_calculator.calculate_distance(case_0.first, case_0.second);
          timer.leave_subsection();
          write_distances(case_0.second, distances, "case_0_y.dat");
        }

        {
          unsigned int n_times = 1e5;
          for (unsigned int t = 0; t < n_times; ++t)
            {
              auto case_1 = generate_case_1<3>(1000, t);
              timer.enter_subsection("Case 1");
              auto distances = distance_calculator.calculate_distance(case_1.first, case_1.second);
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
