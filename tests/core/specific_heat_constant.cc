/**
 * @brief Tests the constant specific heat model. This model should always return a constant.
 */

// Lethe
#include <core/parameters.h>
#include <core/specific_heat_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beggining" << std::endl;


  SpecificHeatConstant specific_heat_model(5);

  deallog << "Testing specific heat" << std::endl;

  deallog << " T = 1    , Cp = " << specific_heat_model.get_specific_heat(1, 1)
          << std::endl;
  deallog << " T = 2    , Cp = " << specific_heat_model.get_specific_heat(2, 2)
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
