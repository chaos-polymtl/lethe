// check the read and write of simulationcontrol

#include "../tests.h"
#include "bdf.h"

int main(int argc, char *argv[]) {
  try {
    initlog();
    std::vector<double> dt(5, 0.1);
    dt[1] = 0.2;
    dt[2] = 0.3;
    dt[3] = 0.4;
    dt[4] = 0.5;
    deallog << "Time steps ";
    for (unsigned int i = 0; i < dt.size(); ++i) {
      deallog << dt[i] << " ";
    }
    deallog << std::endl;

    Vector<double> order1_coefficients = bdf_coefficients(1, dt);
    if (!approximatelyEqual(order1_coefficients[0], 10., 1e-8))
      throw std::runtime_error("Error order 1 term 0");
    if (!approximatelyEqual(order1_coefficients[1], -10., 1e-8))
      throw std::runtime_error("Error order 1 term 1");
    deallog << "Order 1 : " << order1_coefficients << std::endl;

    Vector<double> order2_coefficients = bdf_coefficients(2, dt);
    if (!approximatelyEqual(order2_coefficients[0], 13.3333333333333, 1e-8))
      throw std::runtime_error("Error order 2 term 0");
    if (!approximatelyEqual(order2_coefficients[1], -15.0000000000000, 1e-8))
      throw std::runtime_error("Error order 2 term 1");
    if (!approximatelyEqual(order2_coefficients[2], 1.66666666666667, 1e-8))
      throw std::runtime_error("Error order 2 term 2");
    deallog << "Order 2 : " << order2_coefficients << std::endl;

    Vector<double> order3_coefficients = bdf_coefficients(3, dt);
    if (!approximatelyEqual(order3_coefficients[0], 15.0000000000000, 1e-8))
      throw std::runtime_error("Error order 3 term 0");
    if (!approximatelyEqual(order3_coefficients[1], -18.0000000000000, 1e-8))
      throw std::runtime_error("Error order 3 term 1");
    if (!approximatelyEqual(order3_coefficients[2], 3.33333333333333, 1e-8))
      throw std::runtime_error("Error order 3 term 2");
    if (!approximatelyEqual(order3_coefficients[3], -0.333333333333333, 1e-8))
      throw std::runtime_error("Error order 3 term 3");
    deallog << "Order 3 : " << order3_coefficients << std::endl;
  } catch (std::exception &exc) {
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
  } catch (...) {
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
