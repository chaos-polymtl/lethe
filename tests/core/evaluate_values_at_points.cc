// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

#include <deal.II/base/function_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <../tests/tests.h>

using namespace dealii;


/**
 * DESCRIPTION:
 * Tests the evaluate_values_at_points() function located in utilities.h
 * using different a linear and cosine solution fields. Both a scalar and vector
 * fields (\f$dim \in \{2, 3\} \f$) are tested here.
 *
 * For the cosine and cosine gradient functions, since we are using linear
 * elements to approximate them, a small error is expected and observed except
 * when the points correspond to a DoF.
 */

/**
 * @brief Generates a simple linear scalar field function:
 * \f$ f = 1 + x + 2 * y + 3 * z \f$
 *
 * @remark If dim = 2, \f$ z = 0\f$.
 *
 * @tparam dim Spatial dimension of the problem.
 */
template <int dim>
class LinearScalarField : public Function<dim>
{
public:
  LinearScalarField()
    : Function<dim>(1)
  {}

  double
  value(const Point<dim> &p,
        const unsigned int /*component*/ = 0) const override
  {
    double v = 1.0 + p[0] + 2.0 * p[1];
    if constexpr (dim == 3)
      v += 3.0 * p[2];
    return v;
  }
};

/**
 * @brief Generates a simple linear vector field function:
 * - component 0: \f$ f_x = 1+x+y+z \f$
 * - component 1: \f$ f_y = 2*(1+x+y+z) \f$
 * - component 2: \f$ f_z = 3*(1+x+y+z) \f$
 *
 * @remark If dim = 2, \f$ z = 0\f$ and \f$ f_z \f$ does not exist.
 *
 * @tparam dim Spatial dimension of the problem.
 */
template <int dim>
class LinearVectorField : public Function<dim>
{
public:
  LinearVectorField()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const override
  {
    double c = static_cast<double>(component + 1);
    double v = 0.0;
    for (unsigned int d = 0; d < dim; ++d)
      v += c * p[d];
    return c * v;
  }

  void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    for (unsigned int d = 0; d < dim; ++d)
      values[d] = value(p, d);
  }
};

/**
 * @brief Initializes evaluation points depending on the spatial dimension of
 * the problem.
 *
 * @tparam dim Spatial dimension of the problem.
 *
 * @param[in, out] evaluation_points Vector of points where the values of the
 * solution field are to be evaluated.
 */
template <int dim>
void
inline initialize_evaluation_points(std::vector<Point<dim>> &evaluation_points)
{
  if constexpr (dim == 2)
    {
      evaluation_points.emplace_back(Point<2>({0.5, 0.5}));
      evaluation_points.emplace_back(Point<2>({0.1, 0.1}));
      evaluation_points.emplace_back(Point<2>({0.7, 0.3}));
      evaluation_points.emplace_back(Point<2>({0.425, 0.635}));
    }
  else if constexpr (dim == 3)
    {
      evaluation_points.emplace_back(Point<3>({0.5, 0.5, 0.5}));
      evaluation_points.emplace_back(Point<3>({0.1, 0.1, 0.1}));
      evaluation_points.emplace_back(Point<3>({0.7, 0.3, 0.4}));
      evaluation_points.emplace_back(Point<3>({0.425, 0.635, 0.815}));
    }
}

/**
 * @brief Prints results on the deallog for scalar fields.
 *
 * @tparam dim Spatial dimension of the problem.
 *
 * @param[in] evaluation_points Vector of points where the values of the
 * solution field are to be evaluated.
 * @param[in] evaluated_scalar_values Evaluated scalar values at evaluation
 * points.
 * @param[in] scalar_function Function returning a scalar field corresponding to
 * the analytical solution.
 */
template <int dim>
void
inline print_results(const std::vector<Point<dim>> &evaluation_points,
              const std::vector<double>     &evaluated_scalar_values,
              const Function<dim>           &scalar_function)
{
  for (unsigned int p = 0; p < evaluated_scalar_values.size(); ++p)
    {
      const auto &point_p = evaluation_points[p];
      const auto &value_p = evaluated_scalar_values[p];

      if constexpr (dim == 2)
        deallog << "Evaluation point: (" << point_p[0] << ", " << point_p[1]
                << ")" << std::endl;

      else if constexpr (dim == 3)
        deallog << "Evaluation point: (" << point_p[0] << ", " << point_p[1]
                << ", " << point_p[2] << ")" << std::endl;

      deallog << "  Evaluated value: " << value_p << std::endl;
      deallog << "  Function value:  " << scalar_function.value(point_p)
              << std::endl;
      deallog << std::endl;
    }
}

/**
 * @brief Prints results on the deallog for vector fields.
 *
 * @tparam dim Spatial dimension of the problem.
 *
 * @param[in] evaluation_points Vector of points where the values of the
 * solution field are to be evaluated.
 * @param[in] evaluated_vector_values Evaluated vector values at evaluation
 * points.
 * @param[in] vector_function Function returning a vector field of dimension
 * @p dim corresponding to the analytical solution.
 */
template <int dim>
void
inline print_results(
  const std::vector<Point<dim>>             &evaluation_points,
  const std::vector<Tensor<1, dim, double>> &evaluated_vector_values,
  const Function<dim>                       &vector_function)
{
  Vector<double> function_vector_value(dim);
  for (unsigned int p = 0; p < evaluated_vector_values.size(); ++p)
    {
      const auto &point_p = evaluation_points[p];
      const auto &value_p = evaluated_vector_values[p];

      if constexpr (dim == 2)
        deallog << "Evaluation point: (" << point_p[0] << ", " << point_p[1]
                << ")" << std::endl;

      else if constexpr (dim == 3)
        deallog << "Evaluation point: (" << point_p[0] << ", " << point_p[1]
                << ", " << point_p[2] << ")" << std::endl;

      deallog << "  Evaluated value: " << value_p << std::endl;
      vector_function.vector_value(point_p, function_vector_value);
      deallog << "  Function value:  " << function_vector_value << std::endl;
      deallog << std::endl;
    }
}

/**
 * @brief Tests the evaluate_values_at_points function in the utilities with a
 * linear (scalar and vector) field within a
 * GridGenerator::subdivided_hyper_cube of dimensions [0,1]x[0,1] in 2D, and
 * [0,1]x[0,1]x[0,1] in 3D.
 *
 * @tparam dim Spatial dimension of the problem.
 *
 * @param[in] n_repetitions Number of cells to generate in each direction of
 * the GridGenerator::subdivided_hyper_cube.
 * @param[in] scalar_function Function returning a scalar field corresponding to
 * the analytical solution.
 * @param[in] vector_function Function returning a vector field of dimension
 * @p dim corresponding to the analytical solution.
 */
template <int dim>
void
test(const unsigned int   n_repetitions,
     const Function<dim> &scalar_function,
     const Function<dim> &vector_function)
{
  const MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> triangulation(mpi_communicator);
  GridGenerator::subdivided_hyper_cube(triangulation, n_repetitions);

  MappingQ<dim> mapping(1);

  DoFHandler<dim> dof_handler(triangulation);
  DoFHandler<dim> dof_handler_vector(triangulation);

  FE_Q<dim>     fe(1);
  FESystem<dim> fe_vector(FE_Q<dim>(1), dim);

  dof_handler.distribute_dofs(fe);
  dof_handler_vector.distribute_dofs(fe_vector);

  // Initialize and fill solution fields with function
  TrilinosWrappers::MPI::Vector solution_field_scalar(
    dof_handler.locally_owned_dofs(), mpi_communicator);
  VectorTools::interpolate(mapping,
                           dof_handler,
                           scalar_function,
                           solution_field_scalar);

  TrilinosWrappers::MPI::Vector solution_field_vector(
    dof_handler_vector.locally_owned_dofs(), mpi_communicator);
  VectorTools::interpolate(mapping,
                           dof_handler_vector,
                           vector_function,
                           solution_field_vector);

  // Define evaluation points
  std::vector<Point<dim>> evaluation_points;
  initialize_evaluation_points(evaluation_points);

  // Initialize solution vectors
  std::vector<double>                 evaluated_scalar_values;
  std::vector<Tensor<1, dim, double>> evaluated_vector_values;

  // Evaluate values at evaluation points
  evaluate_values_at_points(triangulation,
                            mapping,
                            dof_handler,
                            solution_field_scalar,
                            evaluation_points,
                            evaluated_scalar_values);

  evaluate_values_at_points<dim>(triangulation,
                                 mapping,
                                 dof_handler_vector,
                                 solution_field_vector,
                                 evaluation_points,
                                 evaluated_vector_values);

  // Print results
  deallog << "Scalar field results: " << std::endl;
  print_results(evaluation_points, evaluated_scalar_values, scalar_function);
  deallog << std::endl;
  deallog << "Vector field results: " << std::endl;
  print_results(evaluation_points, evaluated_vector_values, vector_function);
  deallog << std::endl;
}


int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      deallog << "Running dim == 2, linear functions" << std::endl;
      test<2>(8, LinearScalarField<2>(), LinearVectorField<2>());
      deallog << "************************************************************"
              << std::endl;
      deallog << "Running dim == 2, cosine functions" << std::endl;
      test<2>(32,
              Functions::CosineFunction<2>(1),
              Functions::CosineGradFunction<2>());
      deallog << "************************************************************"
              << std::endl;
      deallog << "Running dim == 3, linear functions" << std::endl;
      test<3>(8, LinearScalarField<3>(), LinearVectorField<3>());
      deallog << "************************************************************"
              << std::endl;
      deallog << "Running dim == 3, cosine functions" << std::endl;
      test<3>(32,
              Functions::CosineFunction<3>(1),
              Functions::CosineGradFunction<3>());
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
