// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the Anderson-Jackson filter on an adaptively refined mesh with
 * hanging nodes, in a channel whose walls are described by the boundary
 * conditions of the fluid dynamics.
 *
 * The unit square is periodic in x, with a noslip wall at y = 0 and a wall
 * moving at the velocity (1, 0), imposed by a function boundary condition, at
 * y = 1. The cells below y = 0.25 are refined once, which creates a line of
 * hanging nodes at y = 0.25. The hanging nodes are not filter centers: their
 * values are interpolated from the constraints.
 *
 * Two fields are filtered:
 * - a uniform field, with the kernel truncated by the walls. The test checks
 *   that the value of every filtered field at every hanging node is the
 *   interpolation of the values of the filter centers that constrain it, and
 *   that the uniform field is recovered at every degree of freedom, hanging
 *   nodes included;
 * - the Couette flow u = (y, 0), with the fluid extended beyond the walls at
 *   the velocities returned by make_wall_velocities(). The test prints the
 *   walls and their velocities, and the filtered velocity and the fluid
 *   volume fraction on both walls and at y = 0.625, where the support of the
 *   kernel lies in the uniform part of the mesh without reaching the walls,
 *   so that the filtered velocity is exact.
 *
 * The output is identical for any number of processes.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/types.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/boundary_conditions.h>
#include <core/grids.h>
#include <core/immersed_solid_classifier.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>
#include <core/shape.h>
#include <core/vector.h>

#include <fem-dem/anderson_jackson_filter.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

// Tests (with common definitions)
#include <../tests/fem-dem/anderson_jackson_filter_test_utilities.h>

#include <../tests/tests.h>

namespace
{
  /// Number of spatial dimensions of the test.
  constexpr int dim = 2;

  /**
   * @brief Couette velocity u = (y, 0) and zero pressure.
   */
  class CouetteFlow : public Function<dim>
  {
  public:
    /**
     * @brief Constructor.
     */
    CouetteFlow()
      : Function<dim>(dim + 1)
    {}

    /**
     * @brief Return a component of the velocity-pressure field.
     *
     * @param[in] point Evaluation point.
     *
     * @param[in] component Component of the field.
     *
     * @return Value of the component.
     */
    double
    value(const Point<dim> &point, const unsigned int component) const override
    {
      return (component == 0) ? point[1] : 0.;
    }
  };

  /**
   * @brief Channel periodic in x, refined below y = 0.25, with the boundary
   * conditions of the fluid dynamics and a velocity-pressure field.
   */
  class RefinedChannel
  {
  public:
    /**
     * @brief Constructor.
     *
     * @param[in] field Velocity-pressure field interpolated on the mesh.
     */
    explicit RefinedChannel(const Function<dim> &field)
      : triangulation(MPI_COMM_WORLD)
      , fe(FE_Q<dim>(1), dim, FE_Q<dim>(1), 1)
      , mapping(1)
    {
      // Boundary conditions of the fluid dynamics: periodic in x (boundaries
      // 0 and 1), noslip at y = 0 (boundary 2) and imposed velocity at y = 1
      // (boundary 3).
      ParameterHandler prm;
      boundary_conditions.declare_parameters(prm, 3);
      prm.parse_input_from_string("subsection boundary conditions\n"
                                  "  set number = 3\n"
                                  "  subsection bc 0\n"
                                  "    set id                 = 0\n"
                                  "    set type               = periodic\n"
                                  "    set periodic id        = 1\n"
                                  "    set periodic direction = 0\n"
                                  "  end\n"
                                  "  subsection bc 1\n"
                                  "    set id   = 2\n"
                                  "    set type = noslip\n"
                                  "  end\n"
                                  "  subsection bc 2\n"
                                  "    set id   = 3\n"
                                  "    set type = function\n"
                                  "    subsection u\n"
                                  "      set Function expression = 1\n"
                                  "    end\n"
                                  "    subsection v\n"
                                  "      set Function expression = 0\n"
                                  "    end\n"
                                  "  end\n"
                                  "end\n");
      boundary_conditions.parse_parameters(prm);

      GridGenerator::hyper_cube(triangulation, 0., 1., true);
      setup_periodic_boundary_conditions(
        triangulation, boundary_conditions.periodic_boundaries);
      triangulation.refine_global(3);
      for (const auto &cell : triangulation.active_cell_iterators())
        if (cell->is_locally_owned() &&
            cell->center()[1] < refined_region_height)
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();

      dof_handler.reinit(triangulation);
      dof_handler.distribute_dofs(fe);

      const IndexSet   locally_owned_dofs = dof_handler.locally_owned_dofs();
      GlobalVectorType owned_solution(locally_owned_dofs, MPI_COMM_WORLD);
      VectorTools::interpolate(mapping, dof_handler, field, owned_solution);
      solution.reinit(locally_owned_dofs,
                      DoFTools::extract_locally_relevant_dofs(dof_handler),
                      MPI_COMM_WORLD);
      solution = owned_solution;
    }

    /// Height below which the cells are refined.
    static constexpr double refined_region_height = 0.25;

    /// Triangulation of the channel.
    parallel::distributed::Triangulation<dim> triangulation;

    /// Velocity-pressure finite element.
    const FESystem<dim> fe;

    /// DoFHandler of the velocity-pressure field.
    DoFHandler<dim> dof_handler;

    /// Mapping of the triangulation.
    const MappingQ<dim> mapping;

    /// Boundary conditions of the fluid dynamics.
    BoundaryConditions::NSBoundaryConditions<dim> boundary_conditions;

    /// Velocity-pressure field with its ghost values.
    GlobalVectorType solution;
  };

  /**
   * @brief Return the parameters of the filter used by the test.
   *
   * @param[in] extend_velocity_beyond_walls Fill the part of the kernel beyond
   * the walls with fluid moving at the velocity of the walls.
   *
   * @return Parameters of a gaussian filter of support radius 0.3.
   */
  Parameters::AndersonJacksonFilter
  make_filter_parameters(const bool extend_velocity_beyond_walls)
  {
    Parameters::AndersonJacksonFilter filter_parameters;
    filter_parameters.kernel_type     = Parameters::FilterKernelType::gaussian;
    filter_parameters.filter_width    = 0.1;
    filter_parameters.gaussian_cutoff = 3.;
    filter_parameters.extend_velocity_beyond_walls =
      extend_velocity_beyond_walls;
    return filter_parameters;
  }

  /**
   * @brief Filter a uniform field with the kernel truncated by the walls and
   * check the values of the filtered fields at the hanging nodes.
   */
  void
  test_hanging_nodes()
  {
    deallog << "Uniform field, kernel truncated by the walls" << std::endl;

    const std::vector<double> uniform_values = {1., 2., 4.};
    RefinedChannel channel{Functions::ConstantFunction<dim>(uniform_values)};

    const std::vector<std::shared_ptr<Shape<dim>>> no_solids;
    const ShapesImmersedSolidClassifier<dim>       classifier(no_solids);

    AndersonJacksonFilter<dim> filter(make_filter_parameters(false),
                                      Parameters::Timer(),
                                      MPI_COMM_WORLD);
    filter.apply(channel.dof_handler,
                 channel.mapping,
                 channel.solution,
                 classifier,
                 channel.boundary_conditions.periodic_boundaries);

    // Hanging node constraints of the space of the filtered fields, rebuilt
    // independently of the filter.
    const DoFHandler<dim> &filter_dof_handler = filter.get_dof_handler();
    const IndexSet locally_owned_dofs = filter_dof_handler.locally_owned_dofs();
    AffineConstraints<double> constraints(
      locally_owned_dofs,
      DoFTools::extract_locally_relevant_dofs(filter_dof_handler));
    DoFTools::make_hanging_node_constraints(filter_dof_handler, constraints);
    constraints.close();

    const auto &fields = filter.get_filtered_fields();
    const std::vector<const AndersonJacksonFilter<dim>::FieldVectorType *>
      all_fields = {&fields.kernel_mass,
                    &fields.fluid_volume_fraction,
                    &fields.solid_volume_fraction,
                    &fields.velocity[0],
                    &fields.velocity[1],
                    &fields.pressure};

    // The value at a hanging node must be the interpolation of the values of
    // the filter centers that constrain it.
    types::global_dof_index n_hanging_nodes     = 0;
    double                  interpolation_error = 0.;
    double                  velocity_error      = 0.;
    double                  pressure_error      = 0.;
    for (const types::global_dof_index dof : locally_owned_dofs)
      {
        for (unsigned int d = 0; d < dim; ++d)
          velocity_error =
            std::max(velocity_error,
                     std::abs(fields.velocity[d](dof) - uniform_values[d]));
        pressure_error =
          std::max(pressure_error,
                   std::abs(fields.pressure(dof) - uniform_values[dim]));

        if (!constraints.is_constrained(dof))
          continue;

        ++n_hanging_nodes;
        for (const auto *field : all_fields)
          {
            double interpolated_value = 0.;
            for (const auto &[master_dof, weight] :
                 *constraints.get_constraint_entries(dof))
              interpolated_value += weight * (*field)(master_dof);
            interpolation_error =
              std::max(interpolation_error,
                       std::abs((*field)(dof)-interpolated_value));
          }
      }
    n_hanging_nodes = Utilities::MPI::sum(n_hanging_nodes, MPI_COMM_WORLD);
    interpolation_error =
      Utilities::MPI::max(interpolation_error, MPI_COMM_WORLD);
    velocity_error = Utilities::MPI::max(velocity_error, MPI_COMM_WORLD);
    pressure_error = Utilities::MPI::max(pressure_error, MPI_COMM_WORLD);

    const auto &statistics = filter.get_statistics();
    deallog << "  number of filter centers              : "
            << statistics.n_filter_centers << std::endl;
    deallog << "  number of hanging nodes               : " << n_hanging_nodes
            << std::endl;
    deallog << "  hanging nodes are interpolated        : "
            << format(interpolation_error < 1e-12) << std::endl;
    deallog << "  velocity is recovered                 : "
            << format(velocity_error < 1e-12) << std::endl;
    deallog << "  pressure is recovered                 : "
            << format(pressure_error < 1e-12) << std::endl;
    deallog << "  undefined phase averages              : "
            << statistics.n_undefined_averages << std::endl;
    deallog << "  minimum kernel mass                   : "
            << format(statistics.minimum_kernel_mass) << std::endl;
    deallog << "  maximum kernel mass                   : "
            << format(statistics.maximum_kernel_mass) << std::endl;
    deallog << "  minimum fluid fraction                : "
            << format(statistics.minimum_fluid_volume_fraction) << std::endl;
    deallog << "  maximum fluid fraction                : "
            << format(statistics.maximum_fluid_volume_fraction) << std::endl;
  }

  /**
   * @brief Filter the Couette flow with the fluid extended beyond the walls
   * at the velocities built from the boundary conditions.
   */
  void
  test_wall_velocities()
  {
    deallog << "Couette flow, fluid extended beyond the walls" << std::endl;

    RefinedChannel channel{CouetteFlow()};

    const std::vector<std::shared_ptr<Shape<dim>>> no_solids;
    const ShapesImmersedSolidClassifier<dim>       classifier(no_solids);

    // The periodic boundaries are not walls.
    const AndersonJacksonFilter<dim>::WallVelocities wall_velocities =
      make_wall_velocities(channel.boundary_conditions, 0.);
    deallog << "  number of walls                       : "
            << wall_velocities.size() << std::endl;
    for (const auto &[id, velocity] : wall_velocities)
      {
        const Point<dim> point(0.5, (id == 2) ? 0. : 1.);
        deallog << "  velocity of the wall " << static_cast<unsigned int>(id)
                << "                : (" << format(velocity->value(point, 0))
                << ", " << format(velocity->value(point, 1)) << ")"
                << std::endl;
      }

    AndersonJacksonFilter<dim> filter(make_filter_parameters(true),
                                      Parameters::Timer(),
                                      MPI_COMM_WORLD);
    filter.apply(channel.dof_handler,
                 channel.mapping,
                 channel.solution,
                 classifier,
                 channel.boundary_conditions.periodic_boundaries,
                 wall_velocities);

    const auto &fields   = filter.get_filtered_fields();
    const auto  value_at = [&](const auto &field, const double y) {
      return value_at_filter_center(filter.get_dof_handler(),
                                    channel.mapping,
                                    field,
                                    Point<dim>(0.5, y));
    };

    constexpr double interior_height = 0.625;
    for (const double y : {0., interior_height, 1.})
      deallog << "  y = " << format(y)
              << ", filtered u = " << format(value_at(fields.velocity[0], y))
              << ", fluid fraction = "
              << format(value_at(fields.fluid_volume_fraction, y)) << std::endl;
    deallog << "  exact away from the walls             : "
            << format(std::abs(value_at(fields.velocity[0], interior_height) -
                               interior_height) < 1e-12)
            << std::endl;
    deallog << "  undefined phase averages              : "
            << filter.get_statistics().n_undefined_averages << std::endl;
  }
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      mpi_initlog();

      test_hanging_nodes();
      test_wall_velocities();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}
