/**
 * @brief Tests the shape thickening feature.
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
  double                    radius          = 0.5;
  double                    layer_thickness = 0.1;
  Point<3>                  position({0., 0., 0.});
  Tensor<1, 3>              orientation({0., 0., 0.});
  std::shared_ptr<Shape<3>> sphere =
    std::make_shared<Sphere<3>>(radius, position, orientation);

  deallog << "Testing value" << std::endl;
  // Testing value of all shape, to confirm proper implementation
  Point<3> p({1., 0.8, 0.75});
  deallog << " Sphere , SD = " << sphere->value(p) << std::endl;
  sphere->set_layer_thickening(layer_thickness);
  deallog << " Sphere , SD = " << sphere->value(p) << std::endl;
  sphere->set_layer_thickening(-layer_thickness);
  deallog << " Sphere , SD = " << sphere->value(p) << std::endl;

  deallog << "OK" << std::endl;

  deallog << "Testing gradient" << std::endl;
  sphere->set_layer_thickening(0);
  Tensor<1, 3> gradient_sphere = sphere->gradient(p);
  deallog << " Gradient for sphere at p[0] = " << gradient_sphere[0]
          << std::endl;
  deallog << " Gradient for sphere at p[1] = " << gradient_sphere[1]
          << std::endl;
  deallog << " Gradient for sphere at p[2] = " << gradient_sphere[2]
          << std::endl;
  sphere->set_layer_thickening(layer_thickness);
  gradient_sphere = sphere->gradient(p);
  deallog << " Gradient for sphere at p[0] = " << gradient_sphere[0]
          << std::endl;
  deallog << " Gradient for sphere at p[1] = " << gradient_sphere[1]
          << std::endl;
  deallog << " Gradient for sphere at p[2] = " << gradient_sphere[2]
          << std::endl;
  sphere->set_layer_thickening(-layer_thickness);
  gradient_sphere = sphere->gradient(p);
  deallog << " Gradient for sphere at p[0] = " << gradient_sphere[0]
          << std::endl;
  deallog << " Gradient for sphere at p[1] = " << gradient_sphere[1]
          << std::endl;
  deallog << " Gradient for sphere at p[2] = " << gradient_sphere[2]
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
