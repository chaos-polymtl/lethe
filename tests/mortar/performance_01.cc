#include <core/mortar_coupling_manager.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/data_out.h>

#include "./tests.h"


using namespace dealii;

double
time_and_reset(std::chrono::time_point<std::chrono::system_clock> &temp)
{
  const double time_mesh =
    (std::chrono::duration_cast<std::chrono::nanoseconds>(
       std::chrono::system_clock::now() - temp)
       .count() /
     1e9);

  temp = std::chrono::system_clock::now();

  return time_mesh;
}

template <int dim>
void
create_mesh(Triangulation<dim> &tria,
            const unsigned int  n_subdivisions,
            const double        radius,
            const double        rotate_pi)
{
  std::vector<Point<dim>> points;

  const double h = 2 * numbers::PI * radius / n_subdivisions;

  for (unsigned int i = 0; i < n_subdivisions; ++i)
    {
      const double phi_i = 2 * numbers::PI / n_subdivisions * i;

      points.emplace_back((radius - h) * std::cos(phi_i),
                          (radius - h) * std::sin(phi_i));
      points.emplace_back((radius - 0) * std::cos(phi_i),
                          (radius - 0) * std::sin(phi_i));

      const double phi_o = 2 * numbers::PI / n_subdivisions * i + rotate_pi;

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

  Triangulation<dim> tria_serial;
  tria_serial.create_triangulation(points, cells, sub_cell_data);

  for (const auto &cell : tria_serial.active_cell_iterators())
    if (cell->active_cell_index() % 2 == 0)
      cell->face(1)->set_boundary_id(1);
    else
      cell->face(0)->set_boundary_id(2);

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
}

// mpirun -np 8 ./tests/mortar/performance_01.debug/performance_01.debug 128 2
// 0.0
int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const unsigned int dim            = 2;
  const unsigned int n_subdivisions = argc >= 2 ? std::stoi(argv[1]) : 20;
  const unsigned int fe_degree      = argc >= 3 ? std::stoi(argv[2]) : 2;
  const double       radius         = 1.0;
  const double       rotate_pi      = argc >= 4 ? std::stod(argv[3]) : 0.0;
  const unsigned int n_components   = 1;
  using Number                      = double;

  FE_Q<dim>      fe(fe_degree);
  QGauss<dim>    quadrature(fe_degree + 1);
  MappingQ1<dim> mapping;

  // Create mesh
  parallel::fullydistributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  create_mesh(tria, n_subdivisions, radius, rotate_pi);

  // Output mesh
  if (true)
    {
      DataOut<dim>   data_out;
      Vector<double> ranks(tria.n_active_cells());
      ranks = Utilities::MPI::this_mpi_process(tria.get_mpi_communicator());
      data_out.attach_triangulation(tria);
      data_out.add_data_vector(ranks, "ranks");
      data_out.build_patches();
      data_out.write_vtu_in_parallel("grid.vtu", tria.get_mpi_communicator());
    }

  // Create DoFHandler
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  const IndexSet            locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  constraints.reinit(dof_handler.locally_owned_dofs(), locally_relevant_dofs);
  constraints.close();


  // Coupling operator: create
  const std::shared_ptr<MortarManagerBase<dim>> mortar_manager =
    std::make_shared<MortarManagerCircle<dim>>(n_subdivisions,
                                               radius,
                                               quadrature,
                                               rotate_pi);

  const std::shared_ptr<CouplingEvaluationBase<dim, Number>>
    coupling_evaluator =
      std::make_shared<CouplingEvaluationSIPG<dim, n_components, Number>>(
        mapping, dof_handler);

  auto time = std::chrono::system_clock::now();

  double time_setup, time_vmult, time_diag = 0;

  const auto coupling_operator =
    std::make_shared<CouplingOperator<dim, Number>>(mapping,
                                                    dof_handler,
                                                    constraints,
                                                    coupling_evaluator,
                                                    mortar_manager,
                                                    1,
                                                    2,
                                                    1.0);

  time_setup = time_and_reset(time);

  // Coupling operator: apply
  if (true)
    {
      LinearAlgebra::distributed::Vector<Number> src(
        dof_handler.locally_owned_dofs(), dof_handler.get_mpi_communicator());
      LinearAlgebra::distributed::Vector<Number> dst(
        dof_handler.locally_owned_dofs(), dof_handler.get_mpi_communicator());

      time_and_reset(time);
      coupling_operator->vmult_add(dst, src);
      time_vmult = time_and_reset(time);
    }


  // Coupling operator: compute diagonal
  if (true)
    {
      LinearAlgebra::distributed::Vector<Number> diagonal(
        dof_handler.locally_owned_dofs(), dof_handler.get_mpi_communicator());

      time_and_reset(time);
      coupling_operator->add_diagonal_entries(diagonal);
      time_diag = time_and_reset(time);
    }

  if (Utilities::MPI::this_mpi_process(tria.get_mpi_communicator()) == 0)
    std::cout << time_setup << " " << time_vmult << " " << time_diag
              << std::endl;
}