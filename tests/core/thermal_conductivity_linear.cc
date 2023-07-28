/**
 * @brief Tests the linear thermal conductivity model. This model should always return A+B*T.
 * In this case k=5+10*T.
 */

// Lethe
#include <core/thermal_conductivity_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;

  ThermalConductivityLinear thermal_conductivity_model(5, 10);

  deallog << "Testing linear temperature dependent thermal conductivity - k"
          << std::endl;

  // Field values must contain temperature
  std::map<field, double> field_values;
  field_values[field::temperature] = 1;
  deallog << " T = 1    , k = "
          << thermal_conductivity_model.value(field_values) << std::endl;
  field_values[field::temperature] = 2;
  deallog << " T = 2    , k = "
          << thermal_conductivity_model.value(field_values) << std::endl;

  deallog << " Analytical jacobian "
          << thermal_conductivity_model.jacobian(field_values,
                                                 field::temperature)
          << std::endl;

  deallog << " Numerical jacobian "
          << thermal_conductivity_model.numerical_jacobian(field_values,
                                                           field::temperature)
          << std::endl;

  std::vector<double>                  temperature_vector({1, 2, 3, 4, 5});
  std::map<field, std::vector<double>> field_vectors;
  field_vectors[field::temperature] = temperature_vector;
  unsigned int        n_pts         = temperature_vector.size();
  std::vector<double> thermal_conductivities(n_pts);
  std::vector<double> jacobians(n_pts);
  std::vector<double> numerical_jacobians(n_pts);


  thermal_conductivity_model.vector_value(field_vectors,
                                          thermal_conductivities);

  thermal_conductivity_model.vector_jacobian(field_vectors,
                                             field::temperature,
                                             jacobians);
  thermal_conductivity_model.vector_jacobian(field_vectors,
                                             field::temperature,
                                             numerical_jacobians);
  deallog << " Vector values " << std::endl;
  for (unsigned int i = 0; i < n_pts; ++i)
    {
      deallog << thermal_conductivities[i] << " ";
    }
  deallog << std::endl;
  deallog << " Analytical jacobians " << std::endl;
  for (unsigned int i = 0; i < n_pts; ++i)
    {
      deallog << jacobians[i] << " ";
    }
  deallog << std::endl;

  deallog << " Numerical jacobians " << std::endl;
  for (unsigned int i = 0; i < n_pts; ++i)
    {
      deallog << jacobians[i] << " ";
    }
  deallog << std::endl;

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
