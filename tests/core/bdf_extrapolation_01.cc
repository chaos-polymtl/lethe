/**
 * @brief This code tests the coefficients required for BDF integration
 * of order n.
 */

// Lethe
#include <core/bdf.h>

// Tests (with common definitions)
#include <../tests/tests.h>


void
test()
{
  std::vector<double> time_o1({.1, 0});
  std::vector<double> time_o2({.2, .1, 0});
  std::vector<double> time_o2_uneq({.3, .1, 0});
  std::vector<double> time_o3({.3, .2, .1, 0});

  unsigned int                     n_val   = 2;
  unsigned int                     n_times = 3;
  std::vector<std::vector<double>> values;
  values.resize(3);
  for (unsigned int i = 0; i < values.size(); ++i)
    {
      values[i].resize(n_val);
      for (unsigned int v = 0; v < values[i].size(); ++v)
        {
          values[i][v] = 3 - i + 100 * v;
        }
    }

  for (unsigned int v = 0; v < n_val; ++v)
    {
      deallog << "Value set " << v << std::endl;
      deallog << "----------";
      deallog << std::endl;

      for (unsigned int i = 0; i < n_times; ++i)
        {
          deallog << values[i][v] << " ";
        }
      deallog << std::endl;
    }



  std::vector<double> out_values(n_val);

  {
    bdf_extrapolate(time_o1, values, 1, out_values);
    deallog << std::endl;
    deallog << "---------";
    deallog << std::endl;
    deallog << "Order 1: ";


    for (unsigned int v = 0; v < out_values.size(); ++v)
      {
        deallog << out_values[v] << " ";
      }
    deallog << std::endl;
  }

  {
    bdf_extrapolate(time_o2, values, 2, out_values);
    deallog << std::endl;
    deallog << "--------------------------";
    deallog << std::endl;
    deallog << "Order 2 equal time steps: ";


    for (unsigned int v = 0; v < out_values.size(); ++v)
      {
        deallog << out_values[v] << " ";
      }
    deallog << std::endl;
  }

  {
    bdf_extrapolate(time_o2_uneq, values, 2, out_values);
    deallog << std::endl;
    deallog << "-----------------------------";
    deallog << std::endl;
    deallog << "Order 2 variable time steps: ";


    for (unsigned int v = 0; v < out_values.size(); ++v)
      {
        deallog << out_values[v] << " ";
      }
    deallog << std::endl;
  }

  {
    bdf_extrapolate(time_o3, values, 3, out_values);
    deallog << std::endl;
    deallog << "-----------------------------";
    deallog << std::endl;
    deallog << "Order 3 constant time steps: ";


    for (unsigned int v = 0; v < out_values.size(); ++v)
      {
        deallog << out_values[v] << " ";
      }
    deallog << std::endl;
  }
}

int
main(int argc, char *argv[])
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

  return 0;
}
