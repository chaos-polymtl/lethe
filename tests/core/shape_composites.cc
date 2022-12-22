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
  double                    radius = 2;
  Tensor<1, 3>              half_lengths({1., 1., 3.});
  Point<3>                  position({0., 0., 0.});
  Tensor<1, 3>              orientation({0., 0., 0.});
  std::shared_ptr<Shape<3>> sphere =
    std::make_shared<Sphere<3>>(radius, Point<3>(), Point<3>());
  std::shared_ptr<Shape<3>> rectangle =
    std::make_shared<Rectangle<3>>(half_lengths, Point<3>(), Point<3>());

  // Creation of a map of shapes
  std::map<unsigned int, std::shared_ptr<Shape<3>>> shapes;
  shapes[0] = sphere;
  shapes[1] = rectangle;

  // Creation of maps of operations
  typedef std::map<
    unsigned int,
    std::tuple<CompositeShape<3>::BooleanOperation, unsigned int, unsigned int>>
                                              op_map;
  typedef CompositeShape<3>::BooleanOperation bool_op;

  op_map operation_union;
  operation_union[2] = std::make_tuple(bool_op::Union, 0, 1);
  op_map operation_difference;
  operation_difference[2] =
    std::make_tuple(bool_op::Difference, 0, 1); // We substract 0 from 1
  op_map operation_intersection;
  operation_intersection[2] = std::make_tuple(bool_op::Intersection, 0, 1);

  // Initialization of composite shapes
  std::shared_ptr<CompositeShape<3>> composite_union =
    std::make_shared<CompositeShape<3>>(shapes,
                                        operation_union,
                                        position,
                                        orientation);
  std::shared_ptr<CompositeShape<3>> composite_difference =
    std::make_shared<CompositeShape<3>>(shapes,
                                        operation_difference,
                                        position,
                                        orientation);
  std::shared_ptr<CompositeShape<3>> composite_intersection =
    std::make_shared<CompositeShape<3>>(shapes,
                                        operation_intersection,
                                        position,
                                        orientation);

  deallog << "Testing value" << std::endl;
  // Testing value of all shape, to confirm proper implementation
  Point<3> p({1., 0.8, 0.75});
  deallog << " Union , SD = " << composite_union->value(p) << std::endl;
  deallog << " Difference , SD = " << composite_difference->value(p)
          << std::endl;
  deallog << " Intersection , SD = " << composite_intersection->value(p)
          << std::endl;
  deallog << "OK" << std::endl;

  deallog << "Testing copy and translation" << std::endl;
  Point<3>                  translation({0.2, 0., 0.3});
  std::shared_ptr<Shape<3>> copy_composite_union =
    composite_union->static_copy();
  copy_composite_union->set_position(translation);
  deallog << " Position of original composite[0] = "
          << composite_union->get_position()[0] << std::endl;
  deallog << " Position of original composite[1] = "
          << composite_union->get_position()[1] << std::endl;
  deallog << " Position of original composite[2] = "
          << composite_union->get_position()[2] << std::endl;
  deallog << " Translation of copy composite" << std::endl;
  deallog << " Position of copy composite[0] = "
          << copy_composite_union->get_position()[0] << std::endl;
  deallog << " Position of copy composite[1] = "
          << copy_composite_union->get_position()[1] << std::endl;
  deallog << " Position of copy composite[2] = "
          << copy_composite_union->get_position()[2] << std::endl;
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
