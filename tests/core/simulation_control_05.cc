/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - by the Lethe authors
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
* Author: Bruno Blais, Polytechnique Montreal, 2020-
*/

/**
 * @brief This test checks that the function update_assembly_method work as expected.
 */

// Lethe
#include <core/parameters.h>
#include <core/simulation_control.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  Parameters::SimulationControl simulation_control_parameters;

  simulation_control_parameters.dt                       = 1;
  simulation_control_parameters.adapt                    = false;
  simulation_control_parameters.maxCFL                   = 1;
  simulation_control_parameters.startup_timestep_scaling = 0.4;
  simulation_control_parameters.bdf_startup_method =
    Parameters::SimulationControl::BDFStartupMethods::multiple_step_bdf;
  simulation_control_parameters.timeEnd                = 4;
  simulation_control_parameters.number_mesh_adaptation = 0;
  simulation_control_parameters.output_name            = "test";
  simulation_control_parameters.subdivision            = 7;
  simulation_control_parameters.output_folder          = "canard";
  simulation_control_parameters.output_frequency       = 8;

  {
    simulation_control_parameters.method =
      Parameters::SimulationControl::TimeSteppingMethod::bdf2;
    SimulationControlTransient simulation_control(
      simulation_control_parameters);


    // Constant time-stepping
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - constant output" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_step_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    while (simulation_control.integrate())
      {
        std::string method;
        if (simulation_control.get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady)
          method = "steady";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
          method = "steady_bdf";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf1)
          method = "bdf1";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf2)
          method = "bdf2";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf3)
          method = "bdf3";
        else
          method = "the assembly method is not initialized";

        deallog << "Iteration : " << simulation_control.get_step_number()
                << "    Time : " << simulation_control.get_current_time()
                << " Assembly method : " << method << std::endl;
      }
  }
  {
    simulation_control_parameters.method =
      Parameters::SimulationControl::TimeSteppingMethod::bdf3;
    SimulationControlTransient simulation_control(
      simulation_control_parameters);


    // Constant time-stepping
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - constant output" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_step_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    while (simulation_control.integrate())
      {
        std::string method;
        if (simulation_control.get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady)
          method = "steady";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
          method = "steady_bdf";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf1)
          method = "bdf1";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf2)
          method = "bdf2";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf3)
          method = "bdf3";
        else
          method = "the assembly method is not initialized";

        deallog << "Iteration : " << simulation_control.get_step_number()
                << "    Time : " << simulation_control.get_current_time()
                << " Assembly method : " << method << std::endl;
      }
  }
  simulation_control_parameters.bdf_startup_method =
    Parameters::SimulationControl::BDFStartupMethods::multiple_step_bdf;
  {
    simulation_control_parameters.method =
      Parameters::SimulationControl::TimeSteppingMethod::bdf2;
    SimulationControlTransient simulation_control(
      simulation_control_parameters);


    // Constant time-stepping
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - constant output" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_step_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    while (simulation_control.integrate())
      {
        std::string method;
        if (simulation_control.get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady)
          method = "steady";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
          method = "steady_bdf";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf1)
          method = "bdf1";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf2)
          method = "bdf2";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf3)
          method = "bdf3";
        else
          method = "the assembly method is not initialized";

        deallog << "Iteration : " << simulation_control.get_step_number()
                << "    Time : " << simulation_control.get_current_time()
                << " Assembly method : " << method << std::endl;
      }
  }
  {
    simulation_control_parameters.method =
      Parameters::SimulationControl::TimeSteppingMethod::bdf3;
    SimulationControlTransient simulation_control(
      simulation_control_parameters);


    // Constant time-stepping
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - constant output" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_step_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    while (simulation_control.integrate())
      {
        std::string method;
        if (simulation_control.get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady)
          method = "steady";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
          method = "steady_bdf";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf1)
          method = "bdf1";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf2)
          method = "bdf2";
        else if (simulation_control.get_assembly_method() ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf3)
          method = "bdf3";
        else
          method = "the assembly method is not initialized";

        deallog << "Iteration : " << simulation_control.get_step_number()
                << "    Time : " << simulation_control.get_current_time()
                << " Assembly method : " << method << std::endl;
      }
  }
}

int
main()
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
}
