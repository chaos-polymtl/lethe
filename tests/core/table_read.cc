// SPDX-FileCopyrightText: Copyright (c) 2021, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II includes
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/particles/data_out.h>


// Lethe
#include <core/utilities.h>

#include <deal.II/numerics/data_out.h>

// Tests (with common definitions)
#include <../tests/tests.h>

/*
 * The test read a table and print the results. We make sure the output of these
 * functions always stay the same.
 */

void
test()
{
  std::string table_file_name = "../table_read_input.dat";

  TableHandler table;
  fill_table_from_file(table, table_file_name);

  table.write_tex(deallog.get_file_stream());

  std::map<std::string, std::vector<double>> vectors;
  fill_vectors_from_file(vectors, table_file_name);

  for (const auto &it : vectors)
    {
      deallog << it.first << std::endl;
      for (const auto &j : it.second)
        deallog << j << std::endl;
    }
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test();
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
