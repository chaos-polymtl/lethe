/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

*
* Author:  Lucka Barbeau , Polytechnique Montreal, 2021-
*/

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
  std::string table_file_name = "table_read_input.dat";

  TableHandler table;
  fill_table_from_file(table, table_file_name);

  table.write_tex(deallog.get_file_stream());

  std::map<std::string, std::vector<double>> vectors;
  fill_vectors_from_file(vectors, table_file_name);

  for (std::map<std::string, std::vector<double>>::iterator it =
         vectors.begin();
       it != vectors.end();
       ++it)
    {
      deallog << it->first << std::endl;
      for (unsigned int j = 0; j < it->second.size(); ++j)
        {
          deallog << it->second[j] << std::endl;
        }
    }
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
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
