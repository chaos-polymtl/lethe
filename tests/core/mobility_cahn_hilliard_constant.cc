/**
 * @brief Tests the constant mobility model. This model should always return a constant.
 */

// Lethe
#include <core/mobility_cahn_hilliard_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;


  MobilityCahnHilliardModelConstant mobility_model(45.91);

  deallog << "Testing constant mobility" << std::endl;

  // field values can remain empty since the mobility does not depend on any
  // fields for the constant model
  std::map<field, double> field_values;

  deallog << "Test 1, mobility = " << mobility_model.value(field_values)
          << std::endl;
  deallog << "Test 2, mobility = " << mobility_model.value(field_values)
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
