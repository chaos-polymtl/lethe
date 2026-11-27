// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the outlet boundary term of the Jacobian of the
 * matrix-free Navier-Stokes operator. The outlet (directional do-nothing)
 * boundary condition adds the term -(v, beta/h^(p+1) (u·n)_ u) to the residual,
 * where (u·n)_ = min(0, u·n). Its Jacobian, -(v, beta/h^(p+1) (u0·n)_ δu), is
 * linearized around the previous Newton iterate u0, whose values on the
 * boundary faces are stored by the operator for every boundary face batch.
 *
 * Two operators that only differ by the type of their boundary conditions,
 * outlet and none, are applied to the same vector around the same
 * linearization point. All the cell terms of the two operators are the same,
 * so the difference of the two results is the contribution of the outlet term.
 * This contribution is compared with a direct evaluation of the term with
 * FEFaceValues.
 *
 * The normal velocity of the linearization point is negative (backflow) on
 * half of the faces of the cube and positive on the other half, and it varies
 * along each face, so that reading the normal velocity of another face batch
 * changes the result. The operator stores the values of the linearization
 * point for the boundary faces only, which MatrixFree numbers after the inner
 * faces. The mesh is uniform: the penalization of a face batch is computed
 * from the size of a single cell of the batch.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/fluid_dynamics_boundary_conditions.h>
#include <core/parameters.h>
#include <core/simulation_control.h>
#include <core/utilities.h>

#include <solvers/fluid_dynamics_matrix_free_operators.h>
#include <solvers/physical_properties_manager.h>

// Tests
#include <../tests/tests.h>

// Standard library
#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>

/// Polynomial degree of the velocity and of the pressure.
constexpr unsigned int fe_degree = 2;

/// Penalization parameter beta of the outlet boundary conditions.
constexpr double outlet_beta = 2.;

/// Tolerance on the relative difference between the outlet contribution of the
/// operator and its direct evaluation. The two only differ by round-off.
constexpr double relative_tolerance = 1e-10;

/**
 * @brief Velocity and pressure of the linearization point. Its normal velocity
 * is negative on the faces x = 1, y = 0 and z = 1 of the unit cube and positive
 * on the faces x = 0, y = 1 and z = 0.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class LinearizationPoint : public Function<dim>
{
public:
  LinearizationPoint()
    : Function<dim>(dim + 1)
  {}

  /**
   * @brief Evaluate one component at a point.
   *
   * @param[in] p Point at which the function is evaluated.
   * @param[in] component Velocity component, or dim for the pressure.
   *
   * @return Value of the component.
   */
  double
  value(const Point<dim> &p, const unsigned int component) const override
  {
    const double x = p[0];
    const double y = p[1];

    if (component == 0)
      return -(1. + y);
    if (component == 1)
      return 0.5 + x;
    if (dim == 3 && component == 2)
      return -(0.5 + y);
    return 0.;
  }
};

/**
 * @brief Vector to which the operators are applied.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class Direction : public Function<dim>
{
public:
  Direction()
    : Function<dim>(dim + 1)
  {}

  /**
   * @brief Evaluate one component at a point.
   *
   * @param[in] p Point at which the function is evaluated.
   * @param[in] component Velocity component, or dim for the pressure.
   *
   * @return Value of the component.
   */
  double
  value(const Point<dim> &p, const unsigned int component) const override
  {
    const double x = p[0];
    const double y = p[1];

    if (component == 0)
      return 1. + x * y;
    if (component == 1)
      return x - y * y;
    if (dim == 3 && component == 2)
      return 0.5 + x * p[2];
    return x + y;
  }
};

/**
 * @brief Apply the Jacobian of a matrix-free operator, linearized around
 * LinearizationPoint, to Direction.
 *
 * @tparam dim Number of spatial dimensions.
 *
 * @param[in] mapping Mapping of the cells.
 * @param[in] dof_handler DoF handler of the velocity and the pressure.
 * @param[in] boundary_conditions Boundary conditions of the operator.
 * @param[in] simulation_control Steady simulation control.
 * @param[in] properties_manager Physical properties of the fluid.
 * @param[out] result Jacobian applied to Direction.
 */
template <int dim>
void
apply_jacobian(
  const Mapping<dim>    &mapping,
  const DoFHandler<dim> &dof_handler,
  const BoundaryConditions::FluidDynamicsBoundaryConditions<dim>
                                                   &boundary_conditions,
  const std::shared_ptr<SimulationControl>         &simulation_control,
  const std::shared_ptr<PhysicalPropertiesManager> &properties_manager,
  LinearAlgebra::distributed::Vector<double>       &result)
{
  AffineConstraints<double> constraints;
  constraints.close();

  Parameters::Stabilization stabilization_parameters;
  stabilization_parameters.stabilization =
    Parameters::Stabilization::NavierStokesStabilization::pspg_supg;

  NavierStokesStabilizedOperator<dim, double> jacobian;
  jacobian.reinit(mapping,
                  dof_handler,
                  constraints,
                  QGauss<dim>(fe_degree + 1),
                  std::shared_ptr<Function<dim>>(),
                  properties_manager,
                  stabilization_parameters,
                  numbers::invalid_unsigned_int,
                  simulation_control,
                  boundary_conditions,
                  false,
                  false,
                  false);

  LinearAlgebra::distributed::Vector<double> linearization_point;
  jacobian.initialize_dof_vector(linearization_point);
  VectorTools::interpolate(mapping,
                           dof_handler,
                           LinearizationPoint<dim>(),
                           linearization_point);
  linearization_point.update_ghost_values();
  jacobian.evaluate_non_linear_term_and_calculate_tau(linearization_point);

  LinearAlgebra::distributed::Vector<double> direction;
  jacobian.initialize_dof_vector(direction);
  VectorTools::interpolate(mapping, dof_handler, Direction<dim>(), direction);

  jacobian.initialize_dof_vector(result);
  jacobian.vmult(result, direction);
}

/**
 * @brief Evaluate the outlet term of the Jacobian directly with FEFaceValues.
 *
 * @tparam dim Number of spatial dimensions.
 *
 * @param[in] mapping Mapping of the cells.
 * @param[in] dof_handler DoF handler of the velocity and the pressure.
 *
 * @return Outlet term of the Jacobian applied to Direction, for every DoF.
 */
template <int dim>
Vector<double>
evaluate_outlet_term(const Mapping<dim>    &mapping,
                     const DoFHandler<dim> &dof_handler)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();

  Vector<double> linearization_point(dof_handler.n_dofs());
  Vector<double> direction(dof_handler.n_dofs());
  VectorTools::interpolate(mapping,
                           dof_handler,
                           LinearizationPoint<dim>(),
                           linearization_point);
  VectorTools::interpolate(mapping, dof_handler, Direction<dim>(), direction);

  const QGauss<dim - 1>            face_quadrature(fe_degree + 1);
  FEFaceValues<dim>                fe_face_values(mapping,
                                   fe,
                                   face_quadrature,
                                   update_values | update_normal_vectors |
                                     update_JxW_values);
  const FEValuesExtractors::Vector velocities(0);

  std::vector<Tensor<1, dim>> linearization_point_values(
    face_quadrature.size());
  std::vector<Tensor<1, dim>>          direction_values(face_quadrature.size());
  std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());

  Vector<double> outlet_term(dof_handler.n_dofs());

  for (const auto &cell : dof_handler.active_cell_iterators())
    for (const auto face : cell->face_indices())
      if (cell->face(face)->at_boundary())
        {
          fe_face_values.reinit(cell, face);
          fe_face_values[velocities].get_function_values(
            linearization_point, linearization_point_values);
          fe_face_values[velocities].get_function_values(direction,
                                                         direction_values);
          cell->get_dof_indices(dof_indices);

          const double cell_size =
            compute_cell_diameter<dim>(cell->measure(), fe_degree);
          const double penalty_parameter =
            outlet_beta / std::pow(cell_size, fe_degree + 1);

          for (const auto q : fe_face_values.quadrature_point_indices())
            {
              const double normal_outflux =
                std::min(0.,
                         linearization_point_values[q] *
                           fe_face_values.normal_vector(q));

              for (const auto i : fe_face_values.dof_indices())
                outlet_term[dof_indices[i]] +=
                  -penalty_parameter * normal_outflux *
                  (direction_values[q] *
                   fe_face_values[velocities].value(i, q)) *
                  fe_face_values.JxW(q);
            }
        }

  return outlet_term;
}

/**
 * @brief Compare the outlet contribution of the matrix-free Jacobian with its
 * direct evaluation.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0., 1., true);
  triangulation.refine_global(dim == 2 ? 2 : 1);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(FESystem<dim>(FE_Q<dim>(fe_degree), dim + 1));

  const MappingQ<dim> mapping(1);

  ParameterHandler prm;
  Parameters::SimulationControl::declare_parameters(prm);
  Parameters::SimulationControl simulation_control_parameters;
  simulation_control_parameters.parse_parameters(prm);
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::steady;
  const std::shared_ptr<SimulationControl> simulation_control =
    std::make_shared<SimulationControlSteady>(simulation_control_parameters);

  Parameters::PhysicalProperties physical_properties;
  physical_properties.fluids.resize(1);
  physical_properties.number_of_fluids = 1;
  physical_properties.fluids[0].rheological_model =
    Parameters::Material::RheologicalModel::newtonian;
  physical_properties.fluids[0].kinematic_viscosity = 0.1;
  physical_properties.fluids[0].density_model =
    Parameters::Material::DensityModel::constant;
  physical_properties.number_of_solids                = 0;
  physical_properties.number_of_material_interactions = 0;
  const auto properties_manager = std::make_shared<PhysicalPropertiesManager>();
  properties_manager->initialize(physical_properties);

  // Every boundary of the colorized cube is an outlet for the first operator
  // and has no boundary condition for the second one.
  BoundaryConditions::FluidDynamicsBoundaryConditions<dim> outlet_conditions;
  BoundaryConditions::FluidDynamicsBoundaryConditions<dim> no_conditions;
  for (types::boundary_id id = 0; std::cmp_less(id, 2 * dim); ++id)
    {
      outlet_conditions.type[id] =
        BoundaryConditions::FluidDynamicsBoundaryType::outlet;
      outlet_conditions.beta[id] = outlet_beta;
      no_conditions.type[id] =
        BoundaryConditions::FluidDynamicsBoundaryType::none;
      no_conditions.beta[id] = 0.;
    }
  outlet_conditions.number_of_boundary_conditions = 2 * dim;
  no_conditions.number_of_boundary_conditions     = 2 * dim;
  outlet_conditions.fix_pressure_constant         = false;
  no_conditions.fix_pressure_constant             = false;

  LinearAlgebra::distributed::Vector<double> with_outlet;
  LinearAlgebra::distributed::Vector<double> without_outlet;
  apply_jacobian(mapping,
                 dof_handler,
                 outlet_conditions,
                 simulation_control,
                 properties_manager,
                 with_outlet);
  apply_jacobian(mapping,
                 dof_handler,
                 no_conditions,
                 simulation_control,
                 properties_manager,
                 without_outlet);

  const Vector<double> outlet_term = evaluate_outlet_term(mapping, dof_handler);

  Vector<double> difference(dof_handler.n_dofs());
  for (types::global_dof_index i = 0; i < dof_handler.n_dofs(); ++i)
    difference[i] = with_outlet(i) - without_outlet(i) - outlet_term[i];

  const double relative_difference =
    difference.l2_norm() / outlet_term.l2_norm();

  deallog << "Dimension " << dim << std::endl;
  deallog << "Norm of the outlet term of the Jacobian: "
          << outlet_term.l2_norm() << std::endl;
  deallog << "Outlet term of the matrix-free Jacobian matches its direct "
             "evaluation: "
          << (relative_difference < relative_tolerance ? "true" : "false")
          << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      initlog();
      test<2>();
      test<3>();
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
