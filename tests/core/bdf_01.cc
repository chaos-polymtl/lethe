// check the read and write of simulationcontrol

#include "../tests.h"
#include "core/bdf.h"


void
test()
{
  std::vector<double> dt(5, 0.1);
  dt[1] = 0.2;
  dt[2] = 0.3;
  dt[3] = 0.4;
  dt[4] = 0.5;
  deallog << "Time steps ";
  for (unsigned int i = 0; i < dt.size(); ++i)
    {
      deallog << dt[i] << " ";
    }
  deallog << std::endl;

  Vector<double> order1_coefficients = bdf_coefficients(1, dt);
  deallog << "Order 1 : " << order1_coefficients;

  Vector<double> order2_coefficients = bdf_coefficients(2, dt);
  deallog << "Order 2 : " << order2_coefficients;

  Vector<double> order3_coefficients = bdf_coefficients(3, dt);
  deallog << "Order 3 : " << order3_coefficients;
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
