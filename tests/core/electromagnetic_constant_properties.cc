// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the constant electromagnetic properties : conductivity, permittivity and permeability.
 */

// Lethe
#include <core/electric_conductivity_model.h>
#include <core/electric_permittivity_model.h>
#include <core/magnetic_permeability_model.h>

#include <../tests/tests.h>

void
test()
{
  // field values can remain empty since the constant models do not depend on
  // any properties.

  std::map<field, double> field_values;

  deallog << "Begin test" << std::endl;


  deallog << "Testing Electric Conductivity" << std::endl;
  ConstantElectricConductivity conductivity_model(5);
  deallog << "Electric Conductivity = "
          << conductivity_model.value(field_values) << std::endl;


  deallog << "Testing Electric Permittivity" << std::endl;
  ConstantElectricPermittivity permittivity_model(6);

  deallog << "Electric Permittivity = "
          << permittivity_model.value(field_values) << std::endl;

  deallog << "Testing Magnetic Permeability" << std::endl;
  ConstantMagneticPermeability permeability_model(7);

  deallog << "Magnetic Permeability = "
          << permeability_model.value(field_values) << std::endl;

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
