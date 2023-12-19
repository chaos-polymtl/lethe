#include <core/shape_parsing.h>

template <int dim>
std::shared_ptr<Shape<dim>>
ShapeGenerator::initialize_shape(const std::string   type,
                                 const std::string   shape_arguments_str,
                                 const Point<dim>   &position,
                                 const Tensor<1, 3> &orientation)
{
  std::shared_ptr<Shape<dim>> shape;
  std::vector<double>         shape_arguments;
  if (type == "rbf" || type == "composite" || type == "opencascade")
    {
      shape = initialize_shape_from_file(type,
                                         shape_arguments_str,
                                         position,
                                         orientation);
    }
  else
    {
      std::vector<std::string> shape_arguments_str_list(
        Utilities::split_string_list(shape_arguments_str, ";"));
      shape_arguments = Utilities::string_to_double(shape_arguments_str_list);
      shape           = initialize_shape_from_vector(type,
                                           shape_arguments,
                                           position,
                                           orientation);
    }
  return shape;
}

template <int dim>
std::shared_ptr<Shape<dim>>
ShapeGenerator::initialize_shape_from_vector(
  const std::string         type,
  const std::vector<double> shape_arguments,
  const Point<dim>         &position,
  const Tensor<1, 3>       &orientation)
{
  std::shared_ptr<Shape<dim>> shape;
  if (type == "sphere")
    shape =
      std::make_shared<Sphere<dim>>(shape_arguments[0], position, orientation);
  else if (type == "hyper rectangle")
    {
      Tensor<1, dim> half_lengths;
      for (unsigned int i = 0; i < dim; ++i)
        {
          half_lengths[i] = shape_arguments[i];
        }
      shape = std::make_shared<HyperRectangle<dim>>(half_lengths,
                                                    position,
                                                    orientation);
    }
  else if (type == "ellipsoid")
    {
      Tensor<1, dim> radii;
      for (unsigned int i = 0; i < dim; ++i)
        {
          radii[i] = shape_arguments[i];
        }
      shape = std::make_shared<Ellipsoid<dim>>(radii, position, orientation);
    }
  else if (type == "superquadric")
    {
      if constexpr (dim == 3)
        {
          Tensor<1, dim> half_lengths{};
          Tensor<1, dim> exponents{};
          for (unsigned int d = 0; d < dim; d++)
            {
              half_lengths[d] = shape_arguments[d];
              exponents[d]    = shape_arguments[d + dim];
            }
          double epsilon;
          if (shape_arguments.size() == 7)
            {
              epsilon = shape_arguments[3 + dim];
            }
          else
            {
              epsilon = 1e-8;
            }
          shape = std::make_shared<Superquadric<dim>>(
            half_lengths, exponents, epsilon, position, orientation);
        }
    }
  else if (type == "torus")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<Torus<dim>>(shape_arguments[0],
                                             shape_arguments[1],
                                             position,
                                             orientation);
    }
  else if (type == "cone")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<Cone<dim>>(shape_arguments[0],
                                            shape_arguments[1],
                                            position,
                                            orientation);
    }
  else if (type == "cylinder")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<Cylinder<dim>>(shape_arguments[0],
                                                shape_arguments[1],
                                                position,
                                                orientation);
    }
  else if (type == "cylindrical tube")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<CylindricalTube<dim>>(shape_arguments[0],
                                                       shape_arguments[1],
                                                       shape_arguments[2],
                                                       position,
                                                       orientation);
    }
  else if (type == "cylindrical helix")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<CylindricalHelix<dim>>(shape_arguments[0],
                                                        shape_arguments[1],
                                                        shape_arguments[2],
                                                        shape_arguments[3],
                                                        position,
                                                        orientation);
    }
  else if (type == "cut hollow sphere")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<CutHollowSphere<dim>>(shape_arguments[0],
                                                       shape_arguments[1],
                                                       shape_arguments[2],
                                                       position,
                                                       orientation);
    }
  else if (type == "death star")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<DeathStar<dim>>(shape_arguments[0],
                                                 shape_arguments[1],
                                                 shape_arguments[2],
                                                 position,
                                                 orientation);
    }
  else if (type == "plane")
    {
      shape = std::make_shared<Plane<dim>>(position, orientation);
    }
  else
    StandardExceptions::ExcNotImplemented();
  return shape;
}

template <int dim>
std::shared_ptr<Shape<dim>>
ShapeGenerator::initialize_shape_from_file(const std::string   type,
                                           const std::string   file_name,
                                           const Point<dim>   &position,
                                           const Tensor<1, 3> &orientation)
{
  std::shared_ptr<Shape<dim>> shape;
  std::vector<double>         shape_arguments;
  if (type == "rbf")
    {
      shape = std::make_shared<RBFShape<dim>>(file_name, position, orientation);
    }
  else if (type == "composite")
    {
      // The following lines retrieve information regarding a
      // composite shape.
      std::map<unsigned int, std::shared_ptr<Shape<dim>>> components;
      std::map<unsigned int,
               std::tuple<typename CompositeShape<dim>::BooleanOperation,
                          unsigned int,
                          unsigned int>>
        operations;
      // In the file, we first loop over all component shapes, then
      // we loop over operations
      std::ifstream myfile(file_name);
      // Read file line by line for section names or arguments
      if (myfile)
        {
          std::string              line;
          std::vector<std::string> column_names;
          std::vector<double>      line_of_data;
          bool                     parsing_shapes     = false;
          bool                     parsing_operations = false;
          // comment marker
          std::string marker = "#";
          while (std::getline(myfile, line))
            {
              // Ignore comments in the line
              size_t markerPos = line.find(marker);
              if (markerPos != std::string::npos)
                {
                  line = line.substr(0, markerPos);
                }
              if (line == "")
                continue;
              if (line == "shapes")
                {
                  parsing_shapes     = true;
                  parsing_operations = false;
                }
              else if (line == "operations")
                {
                  parsing_shapes     = false;
                  parsing_operations = true;
                }
              else
                {
                  std::vector<std::string> list_of_words_base =
                    Utilities::split_string_list(line, ";");
                  std::vector<std::string> list_of_words_clean;
                  for (unsigned int j = 0; j < list_of_words_base.size(); ++j)
                    {
                      if (list_of_words_base[j] != "")
                        {
                          list_of_words_clean.push_back(list_of_words_base[j]);
                        }
                    }
                  if (parsing_shapes)
                    {
                      unsigned int identifier    = stoi(list_of_words_clean[0]);
                      std::string  type          = list_of_words_clean[1];
                      std::string  arguments_str = list_of_words_clean[2];
                      std::replace(arguments_str.begin(),
                                   arguments_str.end(),
                                   ':',
                                   ';');
                      std::string position_str    = list_of_words_clean[3];
                      std::string orientation_str = list_of_words_clean[4];

                      std::vector<std::string> position_str_component =
                        Utilities::split_string_list(position_str, ":");
                      std::vector<std::string> orientation_str_component =
                        Utilities::split_string_list(orientation_str, ":");

                      std::vector<double> temp_position_vec =
                        Utilities::string_to_double(position_str_component);
                      std::vector<double> temp_orientation_vec =
                        Utilities::string_to_double(orientation_str_component);

                      Point<dim> temp_position;
                      Point<3>   temp_orientation =
                        Point<3>({temp_orientation_vec[0],
                                  temp_orientation_vec[1],
                                  temp_orientation_vec[2]});
                      temp_position[0] = temp_position_vec[0];
                      temp_position[1] = temp_position_vec[1];
                      if constexpr (dim == 3)
                        temp_position[2] = temp_position_vec[2];

                      std::shared_ptr<Shape<dim>> shape_temp;
                      shape_temp = ShapeGenerator::initialize_shape(
                        type, arguments_str, Point<dim>(), Tensor<1, 3>());
                      shape_temp->set_position(temp_position);
                      shape_temp->set_orientation(temp_orientation);
                      components[identifier] = shape_temp->static_copy();
                    }
                  else if (parsing_operations)
                    {
                      unsigned int identifier    = stoi(list_of_words_clean[0]);
                      std::string  type          = list_of_words_clean[1];
                      std::string  arguments_str = list_of_words_clean[2];
                      std::vector<std::string> arguments_str_component =
                        Utilities::split_string_list(arguments_str, ":");

                      unsigned int first_shape =
                        stoi(arguments_str_component[0]);
                      unsigned int second_shape =
                        stoi(arguments_str_component[1]);
                      if (type == "union")
                        {
                          operations[identifier] = std::make_tuple(
                            CompositeShape<dim>::BooleanOperation::Union,
                            first_shape,
                            second_shape);
                        }
                      else if (type == "difference")
                        {
                          operations[identifier] = std::make_tuple(
                            CompositeShape<dim>::BooleanOperation::Difference,
                            first_shape,
                            second_shape);
                        }
                      else if (type == "intersection")
                        {
                          operations[identifier] = std::make_tuple(
                            CompositeShape<dim>::BooleanOperation::Intersection,
                            first_shape,
                            second_shape);
                        }
                    }
                }
            }
          myfile.close();
          shape = std::make_shared<CompositeShape<dim>>(components,
                                                        operations,
                                                        position,
                                                        orientation);
        }
      else
        throw std::invalid_argument(file_name);
    }
  else if (type == "opencascade")
    {
      shape = std::make_shared<OpenCascadeShape<dim>>(file_name,
                                                      position,
                                                      orientation);
    }
  return shape;
}

template std::shared_ptr<Shape<2>>
ShapeGenerator::initialize_shape(const std::string   type,
                                 const std::string   arguments,
                                 const Point<2>     &position,
                                 const Tensor<1, 3> &orientation);
template std::shared_ptr<Shape<3>>
ShapeGenerator::initialize_shape(const std::string   type,
                                 const std::string   arguments,
                                 const Point<3>     &position,
                                 const Tensor<1, 3> &orientation);
template std::shared_ptr<Shape<2>>
ShapeGenerator::initialize_shape_from_vector(
  const std::string         type,
  const std::vector<double> shape_arguments,
  const Point<2>           &position,
  const Tensor<1, 3>       &orientation);
template std::shared_ptr<Shape<3>>
ShapeGenerator::initialize_shape_from_vector(
  const std::string         type,
  const std::vector<double> shape_arguments,
  const Point<3>           &position,
  const Tensor<1, 3>       &orientation);
template std::shared_ptr<Shape<2>>
ShapeGenerator::initialize_shape_from_file(const std::string   type,
                                           const std::string   file_name,
                                           const Point<2>     &position,
                                           const Tensor<1, 3> &orientation);
template std::shared_ptr<Shape<3>>
ShapeGenerator::initialize_shape_from_file(const std::string   type,
                                           const std::string   file_name,
                                           const Point<3>     &position,
                                           const Tensor<1, 3> &orientation);
