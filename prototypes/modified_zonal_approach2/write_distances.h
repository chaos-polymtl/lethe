#include <deal.II/base/point.h>
#include <deal.II/base/table_handler.h>

#include <iostream>
#include <vector>

using namespace dealii;
template <int dim>
void
write_distances(const std::vector<Point<dim>> pts,
                const std::vector<double>     distances,
                const std::string             filename)
{
  TableHandler table;

  for (unsigned int i = 0; i < distances.size(); ++i)
    {
      table.add_value("p_x", pts[i][0]);
      table.add_value("p_y", pts[i][1]);
      table.add_value("p_z", pts[i][2]);
      table.add_value("distance", distances[i]);
    }

  std::ofstream output_file(filename);
  table.write_text(output_file);
}
