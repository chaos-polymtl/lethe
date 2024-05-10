#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
 
#include <iostream>
#include <fstream>
#include <cmath>
 
 
using namespace dealii;
 
 
void hypercube_grid()
{
  Triangulation<1> triangulation;
  const double point1(0);
  const double point2(1);

  GridGenerator::hyper_cube(triangulation, point1, point2);
  triangulation.refine_global(4);
 
  std::ofstream out("hyper_cube_test.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to hyper_cube_test.svg" << std::endl;
}
 
void hyperrectangle_grid()
{
  Triangulation<1> triangulation;
  const Point<1> point1(0,0);
  const Point<1> point2(1,1);

 
  GridGenerator::hyper_rectangle(triangulation, point1, point2);
  triangulation.refine_global(4);
 
  std::ofstream out("hyper_rectangle.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to hyper_rectangle.svg" << std::endl;
}
 
void generalcell_grid()
{
  Triangulation<1> triangulation;
  const std::vector<Point<1>> vertice<>;

 
  GridGenerator::hyper_rectangle(triangulation, point1, point2);
  triangulation.refine_global(4);
 
  std::ofstream out("hyper_rectangle.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to hyper_rectangle.svg" << std::endl;
}
 
int main()
{
  first_grid();

}