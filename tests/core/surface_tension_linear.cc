/**
 * @brief Tests the constant surface tension model. This model should always
 * return sigma_0 + dsigma/dT * T.
 */

// Lethe
#include <core/surface_tension_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;

  Parameters::SurfaceTensionParameters surface_tension_parameters;
  surface_tension_parameters.surface_tension_coefficient = 72.86;
  surface_tension_parameters.surface_tension_gradient    = 0.5;

  SurfaceTensionLinear surface_tension_model(surface_tension_parameters);

  deallog
    << "Testing linear surface tension - sigma (sigma_0 = 72.86, dsigma/dT = 0.5)"
    << std::endl;

  // Field values must contain temperature
  std::map<field, double> field_values;

  field_values[field::temperature] = 0;
  deallog << " T = 0, surface tension = "
          << surface_tension_model.value(field_values)
          << ", dsigma/dT analytical = "
          << surface_tension_model.jacobian(field_values, field::temperature)
          << ", dsigma/dT numerical = "
          << surface_tension_model.numerical_jacobian(field_values,
                                                      field::temperature)
          << std::endl;

  field_values[field::temperature] = 300;
  deallog << " T = 300, surface tension = "
          << surface_tension_model.value(field_values)
          << ", dsigma/dT analytical = "
          << surface_tension_model.jacobian(field_values, field::temperature)
          << ", dsigma/dT numerical = "
          << surface_tension_model.numerical_jacobian(field_values,
                                                      field::temperature)
          << std::endl;

  deallog
    << "Surface tension vector values - sigma (sigma_0 = 72.86, dsigma/dT = 0.5)"
    << std::endl;
  std::vector<double>                  temperature_vector({300, 400, 500, 600});
  std::map<field, std::vector<double>> field_vectors;
  field_vectors[field::temperature] = temperature_vector;
  unsigned int        n_values      = temperature_vector.size();
  std::vector<double> surface_tension_values(n_values);
  std::vector<double> surface_tension_gradient_values(n_values);

  surface_tension_model.vector_value(field_vectors, surface_tension_values);
  surface_tension_model.vector_jacobian(field_vectors,
                                        field::temperature,
                                        surface_tension_gradient_values);

  deallog << " T = {" << Utilities::to_string(temperature_vector[0]);
  for (unsigned int i = 1; i < n_values; ++i)
    deallog << ", " << Utilities::to_string(temperature_vector[i]);
  deallog << "}" << std::endl;
  deallog << " sigma = {" << surface_tension_values[0];
  for (unsigned int i = 1; i < n_values; ++i)
    deallog << ", " << surface_tension_values[i];
  deallog << "}" << std::endl;
  deallog << " dsigma/dT = {" << surface_tension_gradient_values[0];
  for (unsigned int i = 1; i < n_values; ++i)
    deallog << ", " << surface_tension_gradient_values[i];
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
