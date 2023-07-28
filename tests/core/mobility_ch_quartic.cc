/**
 * @brief Tests the quartic mobility model. This model should always return $$ D(1-\phi^2)^2 $$
 */

// Lethe
#include <core/mobility_ch_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;


  MobilityCahnHilliardModelQuartic mobility_model(10);

  deallog << "Testing quartic mobility" << std::endl;

  // field values can remain empty since the surface tension density does
  // not depend on any fields
  std::map<field, double> field_values;
  field_values[field::phase_order_ch] = 0.5;
  deallog << "Test 1, mobility = " << mobility_model.value(field_values)
          << std::endl;
  field_values[field::phase_order_ch] = -0.3;
  deallog << "Test 2, mobility = " << mobility_model.value(field_values)
          << std::endl;

  deallog << "Test 3, Analytical jacobian "
          << mobility_model.jacobian(field_values, field::phase_order_ch)
          << std::endl;


  std::vector<double>                  phase_vector({-0.9, -0.7, 0, 0.4, 0.8});
  std::map<field, std::vector<double>> field_vectors;
  field_vectors[field::phase_order_ch] = phase_vector;
  unsigned int        n_pts            = phase_vector.size();
  std::vector<double> mobilities(n_pts);
  std::vector<double> jacobians(n_pts);


  mobility_model.vector_value(field_vectors, mobilities);

  mobility_model.vector_jacobian(field_vectors,
                                 field::phase_order_ch,
                                 jacobians);

  deallog << " Vector values " << std::endl;
  for (unsigned int i = 0; i < n_pts; ++i)
    {
      deallog << mobilities[i] << " ";
    }
  deallog << std::endl;
  deallog << " Analytical jacobians " << std::endl;
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
