// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests QuadraturePointsDataOut, which is used to visualize the QCM
 * filters. One point patch is built per quadrature point of every cell,
 * located at the real-space quadrature point and carrying the filter radius,
 * the filter diameter and the global active index of the cell.
 *
 * The filter radius returned for each cell depends on the cell index, so that
 * the output also verifies that the radius is evaluated per cell and attached
 * to the quadrature points of that cell.
 */

#include <fem-dem/quadrature_points_data_out.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <../tests/tests.h>

#include <sstream>
#include <string>

template <int dim>
void
test()
{
  // Unit cube refined once: 2^dim cells
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0., 1.);
  triangulation.refine_global(1);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  // Default QCM quadrature for a linear void fraction: 2^dim points per cell
  MappingQ<dim> mapping(1);
  QGauss<dim>   quadrature(2);

  QuadraturePointsDataOut<dim> qp_out;
  qp_out.build_patches(dof_handler, mapping, quadrature, [](const auto &cell) {
    return 0.1 * (cell->active_cell_index() + 1);
  });

  // The gnuplot header contains the date and the time at which the file was
  // written. These lines are dropped to keep the output deterministic, along
  // with the blank lines separating the points.
  std::stringstream gnuplot_output;
  qp_out.write_gnuplot(gnuplot_output);

  deallog << "dim = " << dim << std::endl;
  std::string line;
  while (std::getline(gnuplot_output, line))
    {
      const bool is_comment       = line.rfind("#", 0) == 0;
      const bool is_column_header = line.rfind("# <", 0) == 0;
      if (line.empty() || (is_comment && !is_column_header))
        continue;
      deallog << line << std::endl;
    }
}

int
main()
{
  try
    {
      initlog();
      test<2>();
      test<3>();
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
