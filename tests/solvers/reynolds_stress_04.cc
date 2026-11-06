// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the time-averaged Reynolds stresses and turbulent
 * kinetic energy at every degree of freedom, for finite elements that carry
 * more than one degree of freedom per component on some mesh objects (Q3, Q4,
 * unequal velocity and pressure degrees, bubble enrichment), with the vector
 * types and the degree of freedom numberings of the three fluid dynamics
 * solvers.
 *
 * Analytical time-dependent velocity and pressure fields, which differ at every
 * support point, are interpolated at every time step. The reference averages,
 * Reynolds stresses and turbulent kinetic energy are computed independently at
 * the support point of every locally owned degree of freedom, by replaying the
 * time-averaging on the analytical values at that point. The maximal relative
 * error of every averaged quantity is printed for every configuration.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_bubbles.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/vector_tools_interpolate.h>

// Lethe
#include <core/parameters.h>
#include <core/vector.h>

#include <solvers/postprocessing_velocities.h>

// Tests
#include <../tests/tests.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

/**
 * @brief Renumbering of the degrees of freedom applied after their
 * distribution.
 */
enum class Renumbering
{
  none,
  cuthill_mckee,
  component_wise
};

/// Number of time steps over which the fields are averaged.
constexpr unsigned int n_time_steps = 6;

/// Constant time step.
constexpr double time_step = 0.1;

/// Relative errors below this tolerance are printed as zero.
constexpr double error_tolerance = 1e-12;

/**
 * @brief Time dependence shared by all the components of the analytical
 * fields.
 *
 * @param[in] time Time at which the factor is evaluated.
 *
 * @return Value of the time factor.
 */
double
time_factor(const double time)
{
  return 1. + time + time * time;
}

/**
 * @brief Spatial shape of the analytical velocity and pressure fields. The
 * components differ from each other at every point, so a mix-up of two
 * components or of two support points changes the averaged quantities.
 *
 * @param[in] p Point at which the shape is evaluated.
 *
 * @param[in] component Component of the field (velocity components first,
 * pressure last).
 *
 * @return Value of the spatial shape.
 */
template <int dim>
double
spatial_shape(const Point<dim> &p, const unsigned int component)
{
  const double x = p[0];
  const double y = p[1];
  double       z = 0.;
  if constexpr (dim == 3)
    z = p[2];

  if (component == dim)
    return x - z;
  if (component == 0)
    return 1. + x + y * y;
  if (component == 1)
    return 2. - y + z * z * z;
  return 3. + x * z;
}

/**
 * @brief Analytical time-dependent velocity and pressure fields.
 */
template <int dim>
class AnalyticalFields : public Function<dim>
{
public:
  AnalyticalFields()
    : Function<dim>(dim + 1)
  {}

  /**
   * @brief Value of one component of the fields at the current time.
   *
   * @param[in] p Point at which the fields are evaluated.
   *
   * @param[in] component Component of the fields.
   *
   * @return Value of the component.
   */
  double
  value(const Point<dim> &p, const unsigned int component) const override
  {
    return spatial_shape(p, component) * time_factor(this->get_time());
  }
};

/**
 * @brief Reference time-averaged quantities at one point.
 */
template <int dim>
struct ReferenceStatistics
{
  /// Averages of the velocity components and of the pressure.
  std::array<double, dim + 1> average;

  /// Reynolds normal stresses <u'u'>, <v'v'> and <w'w'>.
  std::array<double, dim> normal;

  /// Reynolds shear stresses <u'v'>, <v'w'> and <w'u'> (only <u'v'> in 2D).
  std::array<double, 3> shear;

  /// Turbulent kinetic energy.
  double k;
};

/**
 * @brief Replay the time-averaging of AverageVelocities on the analytical
 * fields at one point, with the averaging starting at the first time step.
 *
 * @param[in] p Point at which the averaged quantities are computed.
 *
 * @return Averaged quantities at the end of the last time step.
 */
template <int dim>
ReferenceStatistics<dim>
compute_reference_statistics(const Point<dim> &p)
{
  ReferenceStatistics<dim>    statistics{};
  std::array<double, dim + 1> sum_value{};
  std::array<double, dim>     sum_normal{};
  std::array<double, 3>       sum_shear{};
  double                      total_time = 0.;

  for (unsigned int step = 0; step < n_time_steps; ++step)
    {
      const double time = step * time_step;
      total_time        = time + time_step;

      std::array<double, dim + 1> value;
      for (unsigned int c = 0; c < dim + 1; ++c)
        {
          value[c] = spatial_shape(p, c) * time_factor(time);
          sum_value[c] += time_step * value[c];
          statistics.average[c] = sum_value[c] / total_time;
        }

      std::array<double, dim> fluctuation;
      for (unsigned int d = 0; d < dim; ++d)
        {
          fluctuation[d] = value[d] - statistics.average[d];
          sum_normal[d] += fluctuation[d] * fluctuation[d] * time_step;
        }

      sum_shear[0] += fluctuation[0] * fluctuation[1] * time_step;
      if constexpr (dim == 3)
        {
          sum_shear[1] += fluctuation[1] * fluctuation[2] * time_step;
          sum_shear[2] += fluctuation[2] * fluctuation[0] * time_step;
        }
    }

  statistics.k = 0.;
  for (unsigned int d = 0; d < dim; ++d)
    {
      statistics.normal[d] = sum_normal[d] / total_time;
      statistics.k += 0.5 * statistics.normal[d];
    }
  for (unsigned int s = 0; s < 3; ++s)
    statistics.shear[s] = sum_shear[s] / total_time;

  return statistics;
}

/**
 * @brief Number of Lagrange degrees of freedom of a scalar base element. The
 * bubble functions of FE_Q_Bubbles are numbered after its Lagrange degrees of
 * freedom.
 *
 * @param[in] base Scalar base element.
 *
 * @return Number of Lagrange degrees of freedom of the base element.
 */
template <int dim>
unsigned int
n_lagrange_dofs(const FiniteElement<dim> &base)
{
  // FE_Q_Bubbles(q) has a polynomial degree of q + 1 and (q + 1)^dim Lagrange
  // degrees of freedom.
  if (dynamic_cast<const FE_Q_Bubbles<dim> *>(&base) != nullptr)
    return Utilities::pow(base.degree, dim);
  return base.n_dofs_per_cell();
}

/**
 * @brief Maximal error of one averaged quantity over the locally owned
 * degrees of freedom.
 */
struct QuantityError
{
  std::string name;
  double      max_error     = 0.;
  double      max_reference = 0.;

  /**
   * @brief Account for one degree of freedom.
   *
   * @param[in] computed Value computed by AverageVelocities.
   *
   * @param[in] reference Reference value.
   */
  void
  add(const double computed, const double reference)
  {
    max_error     = std::max(max_error, std::abs(computed - reference));
    max_reference = std::max(max_reference, std::abs(reference));
  }
};

/**
 * @brief Compute the time-averaged quantities with AverageVelocities for one
 * finite element, degree of freedom numbering and vector type, and print the
 * maximal relative error of every averaged quantity.
 *
 * @param[in] fe Finite element with dim velocity components and one pressure
 * component.
 *
 * @param[in] renumbering Renumbering of the degrees of freedom.
 *
 * @param[in] label Description of the configuration.
 */
template <int dim, typename VectorType, typename DofsType>
void
check(const FiniteElement<dim> &fe,
      const Renumbering         renumbering,
      const std::string        &label)
{
  const MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> triangulation(mpi_communicator);
  GridGenerator::hyper_cube(triangulation, 0., 1.);
  triangulation.refine_global(1);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  // Velocity block and pressure block, as in the block solver
  std::vector<unsigned int> block_component(dim + 1, 0);
  block_component[dim] = 1;

  if (renumbering == Renumbering::cuthill_mckee)
    DoFRenumbering::Cuthill_McKee(dof_handler);
  else if (renumbering == Renumbering::component_wise)
    DoFRenumbering::component_wise(dof_handler, block_component);

  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  DofsType owned_dofs;
  DofsType relevant_dofs;
  if constexpr (std::is_same_v<DofsType, IndexSet>)
    {
      owned_dofs    = locally_owned_dofs;
      relevant_dofs = locally_relevant_dofs;
    }
  else
    {
      const std::vector<types::global_dof_index> dofs_per_block =
        DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
      const types::global_dof_index n_velocity_dofs = dofs_per_block[0];
      const types::global_dof_index n_dofs          = dof_handler.n_dofs();
      owned_dofs    = {locally_owned_dofs.get_view(0, n_velocity_dofs),
                       locally_owned_dofs.get_view(n_velocity_dofs, n_dofs)};
      relevant_dofs = {locally_relevant_dofs.get_view(0, n_velocity_dofs),
                       locally_relevant_dofs.get_view(n_velocity_dofs, n_dofs)};
    }

  AverageVelocities<dim, VectorType, DofsType> average_velocities(dof_handler);
  average_velocities.initialize_vectors(owned_dofs,
                                        relevant_dofs,
                                        mpi_communicator);

  Parameters::PostProcessing<dim> post_processing;
  post_processing.calculate_average_velocities        = true;
  post_processing.initial_time_for_average_velocities = 0.;

  const MappingQ<dim>   mapping(1);
  AnalyticalFields<dim> analytical_fields;
  VectorType            solution;
  solution.reinit(owned_dofs, mpi_communicator);

  for (unsigned int step = 0; step < n_time_steps; ++step)
    {
      const double time = step * time_step;
      analytical_fields.set_time(time);
      VectorTools::interpolate(mapping,
                               dof_handler,
                               analytical_fields,
                               solution);
      average_velocities.calculate_average_velocities(solution,
                                                      post_processing,
                                                      time,
                                                      time_step);
    }

  const std::shared_ptr<VectorType> average =
    average_velocities.get_average_velocities();
  const VectorType &normal_stresses =
    average_velocities.get_reynolds_normal_stresses();
  const VectorType &shear_stresses =
    average_velocities.get_reynolds_shear_stresses();

  // Averaged quantities, in the order in which they are printed
  const std::vector<std::string> component_names =
    (dim == 2) ? std::vector<std::string>{"u", "v", "p"} :
                 std::vector<std::string>{"u", "v", "w", "p"};
  std::vector<QuantityError> errors;
  for (unsigned int c = 0; c < dim + 1; ++c)
    errors.push_back({"<" + component_names[c] + ">"});
  const unsigned int normal_index = errors.size();
  for (unsigned int d = 0; d < dim; ++d)
    errors.push_back(
      {"<" + component_names[d] + "'" + component_names[d] + "'>"});
  const unsigned int k_index = errors.size();
  errors.push_back({"k"});
  const unsigned int shear_index = errors.size();
  errors.push_back({"<u'v'>"});
  if constexpr (dim == 3)
    {
      errors.push_back({"<v'w'>"});
      errors.push_back({"<w'u'>"});
    }
  const unsigned int unused_shear_index = errors.size();
  errors.push_back({"unused shear slots"});

  // Compare every locally owned degree of freedom with the reference
  // computed at its support point. The degrees of freedom shared by several
  // cells are compared more than once, with the same result.
  FEValues<dim>                        fe_values(mapping,
                          fe,
                          Quadrature<dim>(fe.get_unit_support_points()),
                          update_quadrature_points);
  std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);

      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        {
          const types::global_dof_index dof = dof_indices[i];
          if (!locally_owned_dofs.is_element(dof))
            continue;

          const unsigned int component = fe.system_to_component_index(i).first;
          const auto [base, index_within_base] = fe.system_to_base_index(i);
          const bool is_bubble =
            index_within_base >= n_lagrange_dofs(fe.base_element(base.first));

          // The interpolation sets the bubble coefficients to zero, so all the
          // averaged quantities vanish at the bubble degrees of freedom.
          ReferenceStatistics<dim> reference{};
          if (!is_bubble)
            reference =
              compute_reference_statistics(fe_values.quadrature_point(i));

          errors[component].add((*average)(dof), reference.average[component]);

          if (component < dim)
            {
              errors[normal_index + component].add(normal_stresses(dof),
                                                   reference.normal[component]);

              // The shear stresses are stored in the slots of the velocity
              // components: <u'v'> in u, <v'w'> in v and <w'u'> in w. The
              // slot of v is unused in 2D.
              if (component == 0 || dim == 3)
                errors[shear_index + component].add(shear_stresses(dof),
                                                    reference.shear[component]);
              else
                errors[unused_shear_index].add(shear_stresses(dof), 0.);
            }
          else
            {
              errors[k_index].add(normal_stresses(dof), reference.k);
              errors[unused_shear_index].add(shear_stresses(dof), 0.);
            }
        }
    }

  deallog << label << std::endl;
  for (const QuantityError &error : errors)
    {
      const double max_error =
        Utilities::MPI::max(error.max_error, mpi_communicator);
      const double max_reference =
        Utilities::MPI::max(error.max_reference, mpi_communicator);
      // Quantities that vanish everywhere are compared in absolute terms
      const double relative_error =
        (max_reference > 0.) ? max_error / max_reference : max_error;
      const std::string padding(20 - error.name.size(), ' ');
      deallog << "  " << error.name << padding << "max relative error: "
              << ((relative_error < error_tolerance) ? 0. : relative_error)
              << std::endl;
    }
  deallog << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      mpi_initlog();

      using LDVector = LinearAlgebra::distributed::Vector<double>;

      // Elements with a single degree of freedom per component on every mesh
      // object
      check<3, GlobalVectorType, IndexSet>(
        FESystem<3>(FE_Q<3>(1), 3, FE_Q<3>(1), 1),
        Renumbering::cuthill_mckee,
        "3D Q1-Q1, Cuthill-McKee, global vector");
      check<3, GlobalVectorType, IndexSet>(
        FESystem<3>(FE_Q<3>(2), 3, FE_Q<3>(2), 1),
        Renumbering::cuthill_mckee,
        "3D Q2-Q2, Cuthill-McKee, global vector");

      // Matrix-based solver: FESystem(fe_u, dim, fe_p, 1), Cuthill-McKee
      check<3, GlobalVectorType, IndexSet>(
        FESystem<3>(FE_Q<3>(3), 3, FE_Q<3>(3), 1),
        Renumbering::none,
        "3D Q3-Q3, no renumbering, global vector");
      check<3, GlobalVectorType, IndexSet>(
        FESystem<3>(FE_Q<3>(3), 3, FE_Q<3>(3), 1),
        Renumbering::cuthill_mckee,
        "3D Q3-Q3, Cuthill-McKee, global vector");
      check<3, GlobalVectorType, IndexSet>(
        FESystem<3>(FE_Q<3>(2), 3, FE_Q<3>(1), 1),
        Renumbering::cuthill_mckee,
        "3D Q2-Q1, Cuthill-McKee, global vector");
      check<3, GlobalVectorType, IndexSet>(
        FESystem<3>(FE_Q<3>(4), 3, FE_Q<3>(4), 1),
        Renumbering::cuthill_mckee,
        "3D Q4-Q4, Cuthill-McKee, global vector");
      // The pressure shares the vertices, the edge midpoints, the face centers
      // and the cell center of the velocity support points
      check<3, GlobalVectorType, IndexSet>(
        FESystem<3>(FE_Q<3>(4), 3, FE_Q<3>(2), 1),
        Renumbering::cuthill_mckee,
        "3D Q4-Q2, Cuthill-McKee, global vector");
      check<2, GlobalVectorType, IndexSet>(
        FESystem<2>(FE_Q<2>(3), 2, FE_Q<2>(3), 1),
        Renumbering::cuthill_mckee,
        "2D Q3-Q3, Cuthill-McKee, global vector");
      check<3, GlobalVectorType, IndexSet>(
        FESystem<3>(FE_Q_Bubbles<3>(1), 3, FE_Q<3>(1), 1),
        Renumbering::cuthill_mckee,
        "3D Q1 with bubbles-Q1, Cuthill-McKee, global vector");

      // Matrix-free solver: FESystem(FE_Q(k), dim + 1), Cuthill-McKee
      check<3, LDVector, IndexSet>(FESystem<3>(FE_Q<3>(3), 4),
                                   Renumbering::cuthill_mckee,
                                   "3D Q3-Q3, Cuthill-McKee, deal.II vector");
      check<3, LDVector, IndexSet>(FESystem<3>(FE_Q<3>(4), 4),
                                   Renumbering::cuthill_mckee,
                                   "3D Q4-Q4, Cuthill-McKee, deal.II vector");
      check<2, LDVector, IndexSet>(FESystem<2>(FE_Q<2>(3), 3),
                                   Renumbering::cuthill_mckee,
                                   "2D Q3-Q3, Cuthill-McKee, deal.II vector");

      // Block solver: FESystem(fe_u, dim, fe_p, 1), component-wise blocks
      check<3, GlobalBlockVectorType, std::vector<IndexSet>>(
        FESystem<3>(FE_Q<3>(1), 3, FE_Q<3>(1), 1),
        Renumbering::component_wise,
        "3D Q1-Q1, component-wise, block vector");
      check<3, GlobalBlockVectorType, std::vector<IndexSet>>(
        FESystem<3>(FE_Q<3>(3), 3, FE_Q<3>(3), 1),
        Renumbering::component_wise,
        "3D Q3-Q3, component-wise, block vector");
      check<3, GlobalBlockVectorType, std::vector<IndexSet>>(
        FESystem<3>(FE_Q<3>(4), 3, FE_Q<3>(4), 1),
        Renumbering::component_wise,
        "3D Q4-Q4, component-wise, block vector");
      check<3, GlobalBlockVectorType, std::vector<IndexSet>>(
        FESystem<3>(FE_Q<3>(2), 3, FE_Q<3>(1), 1),
        Renumbering::component_wise,
        "3D Q2-Q1, component-wise, block vector");
      check<2, GlobalBlockVectorType, std::vector<IndexSet>>(
        FESystem<2>(FE_Q<2>(3), 2, FE_Q<2>(3), 1),
        Renumbering::component_wise,
        "2D Q3-Q3, component-wise, block vector");
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
