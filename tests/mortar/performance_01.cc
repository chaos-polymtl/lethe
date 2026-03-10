#include <core/mortar_coupling_manager.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>


using namespace dealii;

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const unsigned int dim            = 2;
  const unsigned int n_subdivisions = 20;
  const double       radius         = 1.0;

  std::vector<Point<dim>> points;

  const double h = 2 * numbers::PI * radius / n_subdivisions;

  for (unsigned int i = 0; i < n_subdivisions; ++i)
    {
      const double phi_i = 2 * numbers::PI / n_subdivisions * i;

      points.emplace_back((radius - h) * std::cos(phi_i),
                          (radius - h) * std::sin(phi_i));
      points.emplace_back((radius - 0) * std::cos(phi_i),
                          (radius - 0) * std::sin(phi_i));

      const double phi_o = 2 * numbers::PI / n_subdivisions * i;

      points.emplace_back((radius + 0) * std::cos(phi_o),
                          (radius + 0) * std::sin(phi_o));
      points.emplace_back((radius + h) * std::cos(phi_o),
                          (radius + h) * std::sin(phi_o));
    }

  std::vector<CellData<dim>> cells;

  for (unsigned int i = 0; i < n_subdivisions; ++i)
    {
      CellData<dim> cell_i;
      cell_i.vertices[0] = 4 * (i + 0) + 0;
      cell_i.vertices[1] = 4 * (i + 0) + 1;
      cell_i.vertices[2] = 4 * ((i + 1) % n_subdivisions) + 0;
      cell_i.vertices[3] = 4 * ((i + 1) % n_subdivisions) + 1;
      cells.emplace_back(cell_i);

      CellData<dim> cell_o;
      cell_o.vertices[0] = 4 * (i + 0) + 2;
      cell_o.vertices[1] = 4 * (i + 0) + 3;
      cell_o.vertices[2] = 4 * ((i + 1) % n_subdivisions) + 2;
      cell_o.vertices[3] = 4 * ((i + 1) % n_subdivisions) + 3;
      cells.emplace_back(cell_o);
    }

  SubCellData sub_cell_data;

  parallel::fullydistributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  Triangulation<dim> tria_serial;
  tria_serial.create_triangulation(points, cells, sub_cell_data);

  const unsigned int n_ranks =
    Utilities::MPI::n_mpi_processes(tria.get_mpi_communicator());
  const unsigned int stride = (n_subdivisions + n_ranks - 1) / n_ranks;

  for (const auto &cell : tria_serial.active_cell_iterators())
    {
      const auto index = cell->active_cell_index() / 2;

      cell->set_subdomain_id(index / stride);
    }

  const auto description =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria_serial, tria.get_mpi_communicator());
  tria.create_triangulation(description);

  Vector<double> ranks(tria.n_active_cells());
  ranks = Utilities::MPI::this_mpi_process(tria.get_mpi_communicator());

  DataOut<dim> data_out;
  data_out.attach_triangulation(tria);
  data_out.add_data_vector(ranks, "ranks");
  data_out.build_patches();
  data_out.write_vtu_in_parallel("grid.vtu", tria.get_mpi_communicator());
}