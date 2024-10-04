// SPDX-FileCopyrightText: Copyright (c) 2021 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Assure cost functions for RPT give expected outcomes.
 */

#include <deal.II/base/point.h>

#include <../tests/tests.h>
#include <rpt/rpt.h>
#include <rpt/rpt_calculating_parameters.h>

#include <vector>

void
test()
{
  // Calculated and measured counts
  std::vector<double> calculated_counts(
    {1, 6, 3, 8, 47, 2, 7, 9, 3, 57, 8, 34, 99, 34, 12});
  std::vector<double> measured_counts(
    {2, 7, 2, 7, 44, 1, 9, 7, 2, 55, 9, 33, 88, 32, 8});

  // Objects construction
  RPTCalculatingParameters rpt_parameters;
  RPT<3>                   RPT(rpt_parameters);
  double                   cost_function;

  // Calculate the 3 cost functions
  RPT.rpt_parameters.tuning_param.cost_function_type =
    Parameters::RPTTuningParameters::CostFunctionType::larachi;
  cost_function =
    RPT.calculate_cost_function(calculated_counts, measured_counts);
  deallog << " Larachi cost function = " << cost_function << std::endl;

  RPT.rpt_parameters.tuning_param.cost_function_type =
    Parameters::RPTTuningParameters::CostFunctionType::l1;
  cost_function =
    RPT.calculate_cost_function(calculated_counts, measured_counts);
  deallog << " L1 cost function = " << cost_function << std::endl;

  RPT.rpt_parameters.tuning_param.cost_function_type =
    Parameters::RPTTuningParameters::CostFunctionType::l2;
  cost_function =
    RPT.calculate_cost_function(calculated_counts, measured_counts);
  deallog << " L2 cost function = " << cost_function << std::endl;
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
