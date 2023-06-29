/**
 * @brief Tests the isothermal ideal gas density model. This model should always return a constant.
 */

// Lethe
#include <core/density_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;

  DensityIsothermalIdealGas density_model(1.2, 287.05, 293.15);

  // Field values must contain pressure
  std::map<field, double> field_values;

  deallog << "Testing isothermal ideal gas density - rho" << std::endl;
  field_values[field::pressure] = 0;
  deallog << " p = 0    , density = " << density_model.value(field_values)
          << std::endl;

  field_values[field::pressure] = 1000;
  deallog << " p = 1000 , density = " << density_model.value(field_values)
          << std::endl;

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
