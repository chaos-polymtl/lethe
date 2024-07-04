/**
 * @brief Tests the Tanh Levelset diffusivity model.
 */

// Lethe
#include <core/tracer_diffusivity_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;

  TanhLevelsetTracerDiffusivity diffusivity_model(2, 0, 0.05);


  // Field values must contain shear rate
  std::map<field, double> field_values;

  deallog << "Testing diffusivity of immersed solids tanh model" << std::endl;
  field_values[field::levelset] = 0;
  deallog << " levelset = 0, diffusivity =  "
          << diffusivity_model.value(field_values)
          << ", dD/dlevelset numerical = "
          << diffusivity_model.numerical_jacobian(field_values, field::levelset)
          << std::endl;
  field_values[field::levelset] = 100;
  deallog << " levelset = 100, diffusivity =  "
          << diffusivity_model.value(field_values)
          << ", dD/dlevelset numerical = "
          << diffusivity_model.numerical_jacobian(field_values, field::levelset)
          << std::endl;
  field_values[field::levelset] = -100;
  deallog << " levelset = -100, diffusivity =  "
          << diffusivity_model.value(field_values)
          << ", dD/dlevelset numerical = "
          << diffusivity_model.numerical_jacobian(field_values, field::levelset)
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
