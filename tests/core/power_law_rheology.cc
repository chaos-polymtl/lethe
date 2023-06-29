/**
 * @brief Tests the Power law rheology model. This model should always return a constant.
 */

// Lethe
#include <core/rheological_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beggining" << std::endl;

  PowerLaw rheology_model(5, 0.5, 1e-3);


  // Field values must contain shear rate
  std::map<field, double> field_values;

  deallog << "Testing power law viscosity - nu" << std::endl;
  field_values[field::shear_rate] = 1;
  deallog << " gamma = 1 , nu = " << rheology_model.value(field_values)
          << " , dnu/dgamma analytical = "
          << rheology_model.jacobian(field_values, field::shear_rate)
          << " , dnu/dgamma numerical = "
          << rheology_model.numerical_jacobian(field_values, field::shear_rate)
          << std::endl;
  field_values[field::shear_rate] = 2;
  deallog << " gamma = 2 , nu = " << rheology_model.value(field_values)
          << " , dnu/dgamma analytical = "
          << rheology_model.jacobian(field_values, field::shear_rate)
          << " , dnu/dgamma numerical = "
          << rheology_model.numerical_jacobian(field_values, field::shear_rate)

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
