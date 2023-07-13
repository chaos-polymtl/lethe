/**
 * @brief Tests the mesh controller object
 *
 */

// Lethe
#include <core/mesh_controller.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beginning" << std::endl;

  // Create the mesh controller object
  MeshController mesh_controller(10000);
  deallog << "Testing value" << std::endl;
  // Testing value that are return by the controller for various pseudo number
  // of elements.

  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(1000) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(2000) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(4000) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(8000) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(11000) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(12000) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(11500) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(10500) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(10250) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(10000) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(10000) << std::endl;
  deallog << "coarsening coefficient "
          << mesh_controller.calculate_coarsening_factor(10000) << std::endl;
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
