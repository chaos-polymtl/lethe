#include <deal.II/base/point.h>

#include <stdlib.h> /* srand, rand */

#include <vector>

using namespace dealii;


template <int dim>
std::pair<std::vector<Point<dim>>, std::vector<Point<dim>>>
generate_case_0(unsigned int n_pts, unsigned int axis)
{
  Point<dim> p_0({1, 1, 1});
  Point<dim> p_1({1, 2, 1});
  Point<dim> p_2({1, 1, 2});

  // Generate the triangle
  std::vector<Point<dim>> triangle;
  triangle.push_back(p_0);
  triangle.push_back(p_1);
  triangle.push_back(p_2);

  // Generate a line of points varying along the z direction
  srand(0);
  std::vector<Point<dim>> pts;
  for (unsigned int p = 0; p < n_pts; ++p)
    {
      // Generate three random doubles
      double perturbation = 1.5 * p / (n_pts - 1) - 1;
      double x0           = 1.05;
      double x1           = 1.5;
      double x2           = 1.5;

      pts.push_back({x0, x1, x2});
      pts[p][axis] += perturbation;
    }

  return std::pair<std::vector<Point<dim>>, std::vector<Point<dim>>>(triangle,
                                                                     pts);
}


template <int dim>
std::pair<std::vector<Point<dim>>, std::vector<Point<dim>>>
generate_case_1(const unsigned int n_pts, const unsigned int seed)
{
  Point<dim> p_0({1, 1, 1});
  Point<dim> p_1({1, 2, 1});
  Point<dim> p_2({1, 1, 2});

  // Generate the triangle
  std::vector<Point<dim>> triangle;
  triangle.push_back(p_0);
  triangle.push_back(p_1);
  triangle.push_back(p_2);

  // Generate n_pts random points between 0.5 and 1.5 in all directions.
  srand(seed);
  std::vector<Point<dim>> pts;
  for (unsigned int p = 0; p < n_pts; ++p)
    {
      // Generate three random doubles
      double x0 = (double)rand() / (double)RAND_MAX + 0.75;
      double x1 = (double)rand() / (double)RAND_MAX + 0.75;
      double x2 = (double)rand() / (double)RAND_MAX + 0.75;
      pts.push_back({x0, x1, x2});
    }

  return std::pair<std::vector<Point<dim>>, std::vector<Point<dim>>>(triangle,
                                                                     pts);
}

template <int dim>
std::pair<std::vector<Point<dim>>, std::vector<Point<dim>>>
generate_case_2(const unsigned int n_pts, const unsigned int seed)
{
  Point<dim> p_0({1.5, 1.0, 0.0});
  Point<dim> p_1({3.5, 1.0, 0.0});
  Point<dim> p_2({1.0, 1.5, 0.0});

  // Generate the triangle
  std::vector<Point<dim>> triangle;
  triangle.push_back(p_0);
  triangle.push_back(p_1);
  triangle.push_back(p_2);

  // Generate n_pts random points between 0.5 and 1.5 in all directions.
  srand(seed);
  std::vector<Point<dim>> pts;
  for (unsigned int p = 0; p < n_pts; ++p)
    {
      // Generate three random doubles
      double x0 = 4.0 * p / n_pts;
      double x1 = -0.2 * x0 + 1.7;
      double x2 = 0.0;
      pts.push_back({x0, x1, x2});
    }

  return std::pair<std::vector<Point<dim>>, std::vector<Point<dim>>>(triangle,
                                                                     pts);
}
