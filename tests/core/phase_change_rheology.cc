/**
 * @brief Tests the phase change viscosity rheology model. If T_liquidus > T > T_solidus this model should
 * always return (T-T_solidus)/(T_liquidus-T_solidus) * kinematic_viscosity_l +
 * [1 - (T-T_solidus)/(T_liquidus-T_solidus)] * kinematic_viscosity_s.
 */

// Lethe
#include <core/rheological_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;

  Parameters::PhaseChange phase_change_parameters;

  // initialise necessary parameters
  phase_change_parameters.kinematic_viscosity_l = 1.0;
  phase_change_parameters.kinematic_viscosity_s = 10.0;
  phase_change_parameters.T_liquidus            = 100.5;
  phase_change_parameters.T_solidus             = 100;

  PhaseChangeRheology rheology_model(phase_change_parameters);


  // Field values must contain temperature
  std::map<field, double> field_values;

  deallog << "Testing phase change kinematic viscosity - nu" << std::endl;
  field_values[field::temperature] = 100.1;
  deallog << " T = 100.1, nu = " << rheology_model.value(field_values)
          << ", dnu/dT analytical = "
          << rheology_model.jacobian(field_values, field::shear_rate)
          << ", dnu/dT numerical = "
          << rheology_model.numerical_jacobian(field_values, field::shear_rate)
          << std::endl;
  field_values[field::temperature] = 100.2;
  deallog << " T = 100.2, nu = " << rheology_model.value(field_values)
          << ", dnu/dT analytical = "
          << rheology_model.jacobian(field_values, field::shear_rate)
          << ", dnu/dT numerical = "
          << rheology_model.numerical_jacobian(field_values, field::shear_rate)
          << std::endl;
  field_values[field::temperature] = 100.3;
  deallog << " T = 100.3, nu = " << rheology_model.value(field_values)
          << ", dnu/dT analytical = "
          << rheology_model.jacobian(field_values, field::shear_rate)
          << ", dnu/dT numerical = "
          << rheology_model.numerical_jacobian(field_values, field::shear_rate)
          << std::endl;

  deallog << "Testing phase change dynamic viscosity - mu (T = 100.2)"
          << std::endl;
  field_values[field::temperature] = 100.2;
  double density_ref               = 1;
  deallog << " density_ref = 1, nu = " << rheology_model.value(field_values)
          << ", mu = "
          << rheology_model.get_dynamic_viscosity(density_ref, field_values)
          << std::endl;
  density_ref = 2;
  deallog << " density_ref = 2, nu = " << rheology_model.value(field_values)
          << ", mu = "
          << rheology_model.get_dynamic_viscosity(density_ref, field_values)
          << std::endl;
  density_ref = 3;
  deallog << " density_ref = 3, nu = " << rheology_model.value(field_values)
          << ", mu = "
          << rheology_model.get_dynamic_viscosity(density_ref, field_values)
          << std::endl;

  deallog << "Dynamic viscosity vector values (density_ref = 1)" << std::endl;
  std::vector<double> temperature_vector({100.1, 100.2, 100.3});
  std::map<field, std::vector<double>> field_vectors;
  field_vectors[field::temperature] = temperature_vector;
  unsigned int        n_values      = temperature_vector.size();
  std::vector<double> dynamic_viscosity_values(n_values);
  density_ref = 1;
  rheology_model.get_dynamic_viscosity_vector(density_ref,
                                              field_vectors,
                                              dynamic_viscosity_values);
  deallog << " T = {" << Utilities::to_string(temperature_vector[0]);
  for (unsigned int i = 1; i < n_values; ++i)
    {
      deallog << ", " << Utilities::to_string(temperature_vector[i]);
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
