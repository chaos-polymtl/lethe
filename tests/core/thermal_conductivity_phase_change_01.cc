// SPDX-FileCopyrightText: Copyright (c) 2022 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the thermal_conductivity_phase_change model.
 * During a phase change the physical properties (including conductivity) of the
 * material changes. This interval is defined by [T_solidus,T_liquidus]. This
 * test verifies that the thermal conductivity is calculated correctly. For
 * temperature T<T_solidus, the thermal conductivity should be
 * thermal_conductivity_s. For temperatures T>T_liquidus the thermal
 * conductivity should be thermal_conductivity_l. For temperature in-between, T
 * \in [T_solidus,T_liquidus] the thermal conductivity should be liquid_fraction
 * * thermal_conductivity_l + (1 - liquid_fraction) * thermal_conductivity_s
 */

// Lethe
#include <core/parameters.h>
#include <core/thermal_conductivity_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beggining" << std::endl;

  Parameters::PhaseChange phase_change_params;
  phase_change_params.thermal_conductivity_l = 10;
  phase_change_params.thermal_conductivity_s = 20;
  phase_change_params.T_solidus              = 1;
  phase_change_params.T_liquidus             = 2;


  ThermalConductivityPhaseChange thermal_conductivity_model(
    phase_change_params);


  deallog << "Testing thermal conductivity" << std::endl;

  double T_0 = 0.5;
  double dT  = 0.1;
  double T   = T_0;
  for (int i = 0; i < 20; ++i)
    {
      std::map<field, double> field_values;
      field_values[field::temperature] = T;

      deallog << "T = " << T
              << " k = " << thermal_conductivity_model.value(field_values)
              << std::endl;
      T += dT;
    }


  deallog << "OK" << std::endl;
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
