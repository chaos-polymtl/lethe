/**
 * @brief Tests the constant viscosity rheology model. This model should always return a constant.
 */

// Lethe
#include <core/rheological_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;

  Carreau rheology_model(5, 0, 1, 2, 0.5);


  // Field values must contain shear rate
  std::map<field, double> field_values;

  deallog << "Testing Carreau kinematic viscosity - nu" << std::endl;
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
  field_values[field::shear_rate] = 3;
  deallog << " gamma = 3 , nu = " << rheology_model.value(field_values)
          << " , dnu/dgamma analytical = "
          << rheology_model.jacobian(field_values, field::shear_rate)
          << " , dnu/dgamma numerical = "
          << rheology_model.numerical_jacobian(field_values, field::shear_rate)
          << std::endl;

  deallog << "Testing Carreau dynamic viscosity - mu (gamma = 1)" << std::endl;
  const double dummy_double       = 0;
  field_values[field::shear_rate] = 1;
  double density_ref              = 1;
  deallog << " density_ref = 1, nu = " << rheology_model.value(field_values)
          << ", mu = "
          << rheology_model.get_dynamic_viscosity(
               density_ref, field_values.at(field::shear_rate), dummy_double)
          << std::endl;
  density_ref = 2;
  deallog << " density_ref = 2, nu = " << rheology_model.value(field_values)
          << ", mu = "
          << rheology_model.get_dynamic_viscosity(
               density_ref, field_values.at(field::shear_rate), dummy_double)
          << std::endl;
  density_ref = 3;
  deallog << " density_ref = 3, nu = " << rheology_model.value(field_values)
          << ", mu = "
          << rheology_model.get_dynamic_viscosity(
               density_ref, field_values.at(field::shear_rate), dummy_double)
          << std::endl;

  deallog << "Dynamic viscosity vector values (density_ref = 1)" << std::endl;
  std::vector<double>                  shear_rate_magnitude_vector({1, 2, 3});
  std::map<field, std::vector<double>> field_vectors;
  field_vectors[field::shear_rate] = shear_rate_magnitude_vector;
  unsigned int        n_values     = shear_rate_magnitude_vector.size();
  std::vector<double> dynamic_viscosity_values(n_values);
  density_ref = 1;
  rheology_model.get_dynamic_viscosity_vector(density_ref,
                                              field_vectors,
                                              dynamic_viscosity_values);
  deallog << " gamma = {"
          << Utilities::to_string(shear_rate_magnitude_vector[0]);
  for (unsigned int i = 1; i < n_values; ++i)
    {
      deallog << ", " << Utilities::to_string(shear_rate_magnitude_vector[i]);
    }
  deallog << "}" << std::endl;
  deallog << " mu = {" << dynamic_viscosity_values[0];
  for (unsigned int i = 1; i < n_values; ++i)
    {
      deallog << ", " << dynamic_viscosity_values[i];
    }
  deallog << "}" << std::endl;

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
