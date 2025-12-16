// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the constant electromagnetic conductivity, permittivity and permeability properties.
 */

// Lethe
#include <core/electromagnetic_conductivity_model.h>
#include <core/electromagnetic_permeability_model.h>
#include <core/electromagnetic_permittivity_model.h>

#include <../tests/tests.h>

void
test()
{
  // field values can remain empty since the constant models do not depend on
  // any properties.

  std::map<field, double> field_values;

  deallog << "Begin test" << std::endl;


  deallog << "Testing Electromagnetic Conductivity" << std::endl;
  ConstantElectroMagneticConductivity conductivity_model(5);
  deallog << "Electromagnetic Conductivity = "
          << conductivity_model.value(field_values) << std::endl;


  deallog << "Testing Electromagnetic Permittivity" << std::endl;
  ConstantElectroMagneticPermittivity permittivity_model(6);

  deallog << "Electromagnetic Conductivity = "
          << permittivity_model.value(field_values) << std::endl;

  deallog << "Testing Electromagnetic Permeability" << std::endl;
  ConstantElectroMagneticPermeability permeability_model(7);

  deallog << "Electromagnetic Permeability = "
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
