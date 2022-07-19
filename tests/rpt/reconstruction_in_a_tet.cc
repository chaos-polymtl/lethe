/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 *
 */

/**
 * @brief Check for a few tetrahedron with detector value what is the
 * reconstructed position in the reference space
 */

#include <deal.II/base/point.h>

#include <../tests/tests.h>
#include <rpt/rpt_fem_reconstruction.h>
#include <rpt/parameters_rpt.h>

#include <vector>


void
test()
{
  std::vector<std::vector<double>> vertex_counts;
  vertex_counts.push_back(std::vector<double>({0, 2, 0, 0}));
  vertex_counts.push_back(std::vector<double>({0, 0, 3, 0}));
  vertex_counts.push_back(std::vector<double>({0, 0, 0, 1}));


  std::vector<double> measured_count({0.3, 0.3, 0.3});
  Vector<double>      solution(3);
  Parameters::RPTFEMReconstructionParameters::FEMCostFunction cost_function = Parameters::RPTFEMReconstructionParameters::FEMCostFunction::absolute;

  solution = assemble_matrix_and_rhs<3>(vertex_counts, measured_count, cost_function);


  deallog << "The final solution is : " << solution[0] << " " << solution[1]
          << " " << solution[2] << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
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
