/**
 * @brief Tests the custom parsed functions.
 */

// Lethe
#include <core/parsed_function_custom.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;

  deallog << "Testing first set" << std::endl;
  // Kinetics for a 3 species system of 2A+3B=>1.5C
  std::string variable_names = "A,B,C";
  std::string constants      = "a=2,b=3,c=1.5";
  std::string expression1    = "-a*(A^a*B^b)";
  std::string expression2    = "-b*(A^a*B^b)";
  std::string expression3    = "+c*(A^a*B^b)";
  std::string expression = expression1 + ";" + expression2 + ";" + expression3;
  std::shared_ptr<ParsedFunctionCustom<3>> test_function =
    std::make_shared<ParsedFunctionCustom<3>>();
  test_function->initialize(variable_names, expression, constants);

  Tensor<1, 3> initial_concentrations{};
  initial_concentrations[0] = 1;
  initial_concentrations[1] = 2;
  initial_concentrations[2] = 0;

  Tensor<1, 3> source_term{};
  test_function->vector_value(initial_concentrations, source_term);
  deallog << " Source term = " << source_term << std::endl;
  for (unsigned int i = 0; i < 3; i++)
    deallog << " Source term [" + std::to_string(i) + "] = "
            << test_function->value(initial_concentrations, i) << std::endl;

  Tensor<2, 3> source_term_gradients{};
  test_function->vector_gradient(initial_concentrations, source_term_gradients);
  deallog << " Source term gradients = " << source_term_gradients << std::endl;
  for (unsigned int i = 0; i < 3; i++)
    deallog << " Source term gradient [" + std::to_string(i) + "] = "
            << test_function->gradient(initial_concentrations, i) << std::endl;



  deallog << "Testing second set" << std::endl;
  // Kinetics for a 5 species system of
  // 2A+3B=>1.5C
  // 1C+1D=>2E
  variable_names          = "A,B,C,D,E";
  constants               = "a=2,b=3,c=1.5,d=1,e=2";
  expression1             = "-a*(A^a*B^b)";
  expression2             = "-b*(A^a*B^b)";
  expression3             = "+c*(A^a*B^b)-c*(C^c*D^d)";
  std::string expression4 = "-d*(C^c*D^d)";
  std::string expression5 = "+e*(C^c*D^d)";
  expression = expression1 + ";" + expression2 + ";" + expression3 + ";" +
               expression4 + ";" + expression5;
  std::shared_ptr<ParsedFunctionCustom<5>> test_function_2 =
    std::make_shared<ParsedFunctionCustom<5>>();
  test_function_2->initialize(variable_names, expression, constants);

  Tensor<1, 5> initial_concentrations_2{};
  initial_concentrations_2[0] = 1;
  initial_concentrations_2[1] = 2;
  initial_concentrations_2[2] = 0.1;
  initial_concentrations_2[3] = 2;
  initial_concentrations_2[4] = 0;

  Tensor<1, 5> source_term_2{};
  test_function_2->vector_value(initial_concentrations_2, source_term_2);
  deallog << " Source term = " << source_term_2 << std::endl;
  for (unsigned int i = 0; i < 5; i++)
    deallog << " Source term [" + std::to_string(i) + "] = "
            << test_function_2->value(initial_concentrations_2, i) << std::endl;

  Tensor<2, 5> source_term_gradients_2{};
  test_function_2->vector_gradient(initial_concentrations_2,
                                   source_term_gradients_2);
  deallog << " Source term gradients = " << source_term_gradients_2
          << std::endl;
  for (unsigned int i = 0; i < 5; i++)
    deallog << " Source term gradient [" + std::to_string(i) + "] = "
            << test_function_2->gradient(initial_concentrations_2, i)
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
