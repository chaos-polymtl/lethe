/**
 * @brief Tests the composite shape representation.
 *
 */

// Lethe
#include <core/shape.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;

  // Solid parameters
  Tensor<1, 3>                         orientation({0., 0., 0.});
  std::shared_ptr<OpenCascadeShape<3>> hex_step =
    std::make_shared<OpenCascadeShape<3>>("hex_extrusion_cad.step",
                                          Point<3>(),
                                          orientation);

  deallog << "Testing value" << std::endl;
  // Testing value of all shape, to confirm proper implementation
  Point<3> p({1.5, 1.8, 1.75});
  deallog << "Hex extrusion value at (1.5; 1.8; 1.75) = " << hex_step->value(p)
          << std::endl;
  deallog << "OK" << std::endl;

  deallog << "Testing translation rotation" << std::endl;
  Point<3>     translation({0.2, 0., 0.3});
  Tensor<1, 3> orientation_2({1., 1., 1.});
  deallog << "Translation =(0.2; 0.; 0.3)" << std::endl;
  deallog << "Orientation =(1.; 1.; 1.)" << std::endl;
  hex_step->set_position(translation);
  hex_step->set_orientation(orientation_2);
  deallog << "Hex extrusion value at (1.5; 1.8; 1.75) = " << hex_step->value(p)
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
