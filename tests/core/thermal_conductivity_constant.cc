/**
 * @brief Tests the constant thermal conductivity model. This model should always return a constant.
 */

// Lethe
#include <core/thermal_conductivity_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beggining" << std::endl;

  ConstantThermalConductivity thermal_conductivity_model(5);

  deallog << "Testing thermal conductivity - k" << std::endl;

  // field values can remain empty since the constant thermal conductivity does
  // not depend on any fields
  std::map<field, double> field_values;

  deallog << " T = 1    , k = "
          << thermal_conductivity_model.value(field_values) << std::endl;
  deallog << " T = 2    , k = "
          << thermal_conductivity_model.value(field_values) << std::endl;

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
