// SPDX-FileCopyrightText: Copyright (c) 2019-2020, 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Includes Lethe
#include <core/parameters.h>
#include <core/sdirk_stage_data.h>
#include <core/simulation_control.h>

#include <solvers/navier_stokes_scratch_data.h>
#include <solvers/physical_properties_manager.h>

// Includes Deal.II manquants
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  // We test in dimension 2
  const unsigned int dim = 2;

  // We create a SimulationControl object with the SDIRK22 method
  Parameters::SimulationControl simulation_control_parameters;
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::sdirk22;

  std::shared_ptr<SimulationControl> simulation_control =
    std::make_shared<SimulationControlTransient>(simulation_control_parameters);

  // We create a PhysicalPropertiesManager object with default parameters
  PhysicalPropertiesManager properties_manager;

  // We create the necessary Deal.II elements for the test
  // These are arbitrary
  FESystem<dim>   fe(FE_Q<dim>(1), dim + 1);
  QGauss<dim>     quadrature(2);
  QGauss<dim - 1> face_quadrature(2);
  MappingQ1<dim>  mapping;

  // We create the NavierStokesScratchData object
  NavierStokesScratchData<dim> scratch_data(simulation_control,
                                            properties_manager,
                                            fe,
                                            quadrature,
                                            mapping,
                                            face_quadrature);

  // 12 digits of precision for the output
  deallog << std::setprecision(12) << std::scientific;

  deallog << "Testing SDIRK22 coefficients" << std::endl;

  // Important note : the nomenclature used for the name of the SDIRK methods
  // are sdirkOrderStage sdirk22 means SDIRK with order 2 and 2 stages, sdirk33
  // means SDIRK with order 3 and 3 stages.
  SDIRKTable         table    = scratch_data.sdirk_table;
  const unsigned int n_stages = table.A.m();

  // Data printed at each stage
  for (unsigned int stage_i = 1; stage_i <= n_stages; ++stage_i)
    {
      SDIRKStageData data(table, stage_i);

      deallog << "\nStage " << stage_i << ":" << std::endl;
      deallog << "  a_ij: ";
      for (const auto &a : data.a_ij)
        deallog << std::setw(12) << a << " ";
      deallog << "\n  c_i : " << data.c_i << std::endl;
    }
}

int
main(int argc, char *argv[])
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
