/**
 * @brief Tests the constant surface tension model. This model should always return a constant.
 */

// Lethe
#include <core/surface_tension_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;


  SurfaceTensionConstant surface_tension_model(72.86);

  deallog << "Testing surface tension" << std::endl;

  // field values can remain empty since the surface tension density does
  // not depend on any fields
  std::map<field, double> field_values;

  deallog << "Test 1, surface tension = "
          << surface_tension_model.value(field_values) << std::endl;
  deallog << "Test 2, surface tension = "
          << surface_tension_model.value(field_values) << std::endl;

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
