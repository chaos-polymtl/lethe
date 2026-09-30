// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef anderson_jackson_filter_test_utilities_h
#define anderson_jackson_filter_test_utilities_h

/**
 * @brief Utilities shared by the unit tests of the Anderson-Jackson filter.
 * They build a uniform velocity-pressure field on a uniform mesh of the unit
 * square (cube), and probe the filtered fields. Every probed quantity is
 * reduced over all the processes, so that the output of the tests does not
 * depend on the number of processes.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/grids.h>
#include <core/periodic_boundary.h>
#include <core/vector.h>

#include <fem-dem/anderson_jackson_filter.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace dealii;

/**
 * @brief Format a number in scientific notation.
 *
 * @param[in] value Number to format.
 *
 * @param[in] precision Number of digits after the decimal point.
 *
 * @return Formatted number.
 */
inline std::string
format(const double value, const unsigned int precision = 6)
{
  std::ostringstream stream;
  stream << std::scientific << std::setprecision(precision) << value;
  return stream.str();
}

/**
 * @brief Format a boolean as true or false.
 *
 * @param[in] value Boolean to format.
 *
 * @return Formatted boolean.
 */
inline std::string
format(const bool value)
{
  return value ? "true" : "false";
}

/**
 * @brief Uniform velocity-pressure field on a uniform mesh of the unit square
 * (cube), optionally periodic in every direction. The velocity is (1, 2) in
 * 2D and (1, 2, 3) in 3D, and the pressure is 4.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class UniformFluidProblem
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] n_refinements Number of global refinements of the unit square
   * (cube).
   *
   * @param[in] velocity_degree Degree of the velocity and pressure elements
   * and of the mapping.
   *
   * @param[in] periodic Make the domain periodic in every direction.
   */
  UniformFluidProblem(const unsigned int n_refinements,
                      const unsigned int velocity_degree,
                      const bool         periodic)
    : triangulation(MPI_COMM_WORLD)
    , fe(FE_Q<dim>(velocity_degree), dim, FE_Q<dim>(velocity_degree), 1)
    , mapping(velocity_degree)
  {
    GridGenerator::hyper_cube(triangulation, 0., 1., true);
    if (periodic)
      {
        for (unsigned int d = 0; d < dim; ++d)
          periodic_boundaries[2 * d] = {static_cast<types::boundary_id>(2 * d +
                                                                        1),
                                        d};
        setup_periodic_boundary_conditions(triangulation, periodic_boundaries);
      }
    triangulation.refine_global(n_refinements);

    dof_handler.reinit(triangulation);
    dof_handler.distribute_dofs(fe);

    const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    GlobalVectorType owned_solution(locally_owned_dofs, MPI_COMM_WORLD);
    VectorTools::interpolate(mapping,
                             dof_handler,
                             Functions::ConstantFunction<dim>(uniform_values()),
                             owned_solution);
    solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    solution = owned_solution;
  }

  /**
   * @brief Return the uniform values of the velocity components followed by
   * the pressure.
   *
   * @return Uniform values of the fields.
   */
  static std::vector<double>
  uniform_values()
  {
    if constexpr (dim == 2)
      return {1., 2., 4.};
    else
      return {1., 2., 3., 4.};
  }

  /// Triangulation of the unit square (cube).
  parallel::distributed::Triangulation<dim> triangulation;

  /// Velocity-pressure finite element.
  const FESystem<dim> fe;

  /// DoFHandler of the velocity-pressure field.
  DoFHandler<dim> dof_handler;

  /// Mapping of the triangulation.
  const MappingQ<dim> mapping;

  /// Pairs of periodic boundaries, empty if the domain is not periodic.
  Parameters::PeriodicBoundaries periodic_boundaries;

  /// Uniform velocity-pressure field with its ghost values.
  GlobalVectorType solution;
};

/**
 * @brief Return the value of a filtered field at the filter center located at
 * a given point.
 *
 * @tparam dim Number of spatial dimensions.
 *
 * @param[in] dof_handler DoFHandler of the filtered fields.
 *
 * @param[in] mapping Mapping of the triangulation.
 *
 * @param[in] field Filtered field.
 *
 * @param[in] point Location of the filter center. It must be a support point
 * of the finite element of the filtered fields.
 *
 * @return Value of the field at the filter center.
 */
template <int dim>
double
value_at_filter_center(
  const DoFHandler<dim>                                      &dof_handler,
  const Mapping<dim>                                         &mapping,
  const typename AndersonJacksonFilter<dim>::FieldVectorType &field,
  const Point<dim>                                           &point)
{
  constexpr double location_tolerance = 1e-10;

  const Quadrature<dim> support_points(
    dof_handler.get_fe().get_unit_support_points());
  FEValues<dim>                        fe_values(mapping,
                          dof_handler.get_fe(),
                          support_points,
                          update_quadrature_points);
  std::vector<types::global_dof_index> dof_indices(
    dof_handler.get_fe().n_dofs_per_cell());

  double value = std::numeric_limits<double>::lowest();
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);
      for (unsigned int i = 0; i < dof_indices.size(); ++i)
        if (fe_values.quadrature_point(i).distance(point) < location_tolerance)
          value = field(dof_indices[i]);
    }

  return Utilities::MPI::max(value, MPI_COMM_WORLD);
}

/**
 * @brief Return the maximum errors of the phase-averaged velocity and
 * pressure with respect to the uniform field, over the filter centers where
 * the averages are defined.
 *
 * @tparam dim Number of spatial dimensions.
 *
 * @param[in] filter Filter that was applied to the uniform field.
 *
 * @return Maximum velocity error and maximum pressure error.
 */
template <int dim>
std::pair<double, double>
maximum_phase_average_errors(const AndersonJacksonFilter<dim> &filter)
{
  const std::vector<double> uniform_values =
    UniformFluidProblem<dim>::uniform_values();
  const auto &fields = filter.get_filtered_fields();

  double velocity_error = 0.;
  double pressure_error = 0.;
  for (const types::global_dof_index dof :
       filter.get_dof_handler().locally_owned_dofs())
    if (fields.valid(dof) > 0.5)
      {
        for (unsigned int d = 0; d < dim; ++d)
          velocity_error =
            std::max(velocity_error,
                     std::abs(fields.velocity[d](dof) - uniform_values[d]));
        pressure_error =
          std::max(pressure_error,
                   std::abs(fields.pressure(dof) - uniform_values[dim]));
      }

  return {Utilities::MPI::max(velocity_error, MPI_COMM_WORLD),
          Utilities::MPI::max(pressure_error, MPI_COMM_WORLD)};
}

/**
 * @brief Check that the phase-averaged velocity and pressure are zero at the
 * filter centers where they are not defined.
 *
 * @tparam dim Number of spatial dimensions.
 *
 * @param[in] filter Filter that was applied.
 *
 * @return True if every undefined average is zero, on every process.
 */
template <int dim>
bool
undefined_averages_are_zero(const AndersonJacksonFilter<dim> &filter)
{
  const auto &fields = filter.get_filtered_fields();

  double maximum_value = 0.;
  for (const types::global_dof_index dof :
       filter.get_dof_handler().locally_owned_dofs())
    if (fields.valid(dof) < 0.5)
      {
        for (unsigned int d = 0; d < dim; ++d)
          maximum_value =
            std::max(maximum_value, std::abs(fields.velocity[d](dof)));
        maximum_value = std::max(maximum_value, std::abs(fields.pressure(dof)));
      }

  // The undefined averages are set exactly to zero by the filter.
  return Utilities::MPI::max(maximum_value, MPI_COMM_WORLD) <= 0.;
}

#endif
