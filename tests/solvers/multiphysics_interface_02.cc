// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests that the void fraction can be shared through the
 * multiphysics interface when the fluid dynamics formulation is VANS. A void
 * fraction field, its previous time steps, its DoFHandler and its mapping are
 * registered with PhysicsID::void_fraction, then read back only through the
 * getters of the interface and integrated over the unit square/cube. Since
 * the fields are linear, they are represented exactly by FE_Q(1) and their
 * integrals are known analytically. The test also checks that:
 * - the interface shares the registered vector instead of copying it, so a
 *   new void fraction computed by the solver is seen by the other physics;
 * - the void fraction is refused by the block getters, since it is only
 *   registered as a regular (non-block) field;
 * - the void fraction is refused when the formulation is not VANS.
 */

// Deal.II includes
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools_interpolate.h>

// Lethe
#include <core/multiphysics.h>
#include <core/parameters.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/simulation_parameters.h>

#include <exception>
#include <memory>
#include <vector>

// Tests
#include <../tests/tests.h>

using namespace dealii;

/**
 * @brief Linear void fraction field eps(x) = intercept + slope * x. Its
 * integral over the unit square/cube [0,1]^dim is intercept + slope / 2.
 */
template <int dim>
class LinearVoidFraction : public Function<dim>
{
public:
  LinearVoidFraction(const double intercept, const double slope)
    : Function<dim>(1)
    , intercept(intercept)
    , slope(slope)
  {}

  double
  value(const Point<dim> &p, const unsigned int component = 0) const override;

private:
  const double intercept;
  const double slope;
};

template <int dim>
double
LinearVoidFraction<dim>::value(const Point<dim> &p,
                               const unsigned int /*component*/) const
{
  return intercept + slope * p[0];
}

/**
 * @brief Interpolate a function into a ghosted vector, the same way the
 * ParticleProjector fills its void fraction: interpolate into a locally owned
 * vector, then copy it into the ghosted vector and update the ghost values.
 *
 * @param[in] mapping Mapping used for the interpolation.
 * @param[in] dof_handler DoFHandler of the field.
 * @param[in] function Function to interpolate.
 * @param[in,out] ghosted_vector Ghosted vector that receives the field.
 */
template <int dim>
void
interpolate_in_ghosted_vector(const Mapping<dim>    &mapping,
                              const DoFHandler<dim> &dof_handler,
                              const Function<dim>   &function,
                              GlobalVectorType      &ghosted_vector)
{
  GlobalVectorType locally_owned_vector(dof_handler.locally_owned_dofs(),
                                        dof_handler.get_mpi_communicator());
  VectorTools::interpolate(mapping,
                           dof_handler,
                           function,
                           locally_owned_vector);
  ghosted_vector = locally_owned_vector;
  ghosted_vector.update_ghost_values();
}

/**
 * @brief Integrate a scalar field over the locally owned cells and sum the
 * result over all MPI processes.
 *
 * @param[in] mapping Mapping of the field.
 * @param[in] dof_handler DoFHandler of the field.
 * @param[in] field Ghosted vector of the field.
 *
 * @return Integral of the field over the whole domain.
 */
template <int dim>
double
integrate_field(const Mapping<dim>     &mapping,
                const DoFHandler<dim>  &dof_handler,
                const GlobalVectorType &field)
{
  const QGauss<dim> quadrature(dof_handler.get_fe().degree + 1);
  FEValues<dim>     fe_values(mapping,
                          dof_handler.get_fe(),
                          quadrature,
                          update_values | update_JxW_values);

  std::vector<double> field_values(quadrature.size());
  double              local_integral = 0.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      fe_values.get_function_values(field, field_values);

      for (unsigned int q = 0; q < quadrature.size(); ++q)
        local_integral += field_values[q] * fe_values.JxW(q);
    }

  return Utilities::MPI::sum(local_integral,
                             dof_handler.get_mpi_communicator());
}

template <int dim>
void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  deallog << "  Dimension " << dim << std::endl;

  // Unit square/cube, refined enough to be split between the MPI processes.
  std::shared_ptr<parallel::distributed::Triangulation<dim>> tria =
    std::make_shared<parallel::distributed::Triangulation<dim>>(
      mpi_communicator);
  GridGenerator::hyper_cube(*tria, 0., 1.);

  // Refine the mesh enough to have a non-trivial number of cells per MPI rank
  tria->refine_global(dim == 2 ? 3 : 2);

  // Only the fluid dynamics is enabled, since the other auxiliary physics that
  // depend on the fluid velocity throw under the VANS formulation.
  SimulationParameters<dim>     solver_parameters;
  ParameterHandler              dummy_handler;
  Parameters::SizeOfSubsections size_of_subsections;
  size_of_subsections.boundary_conditions = 1;
  size_of_subsections.manifolds           = 1;

  solver_parameters.declare(dummy_handler, size_of_subsections);
  solver_parameters.parse(dummy_handler);
  solver_parameters.multiphysics.fluid_dynamics = true;

  std::shared_ptr<SimulationControl> simulation_control =
    std::make_shared<SimulationControlTransient>(
      solver_parameters.simulation_control);

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  // Void fraction fields
  const FE_Q<dim>                  fe(1);
  std::shared_ptr<DoFHandler<dim>> dof_handler =
    std::make_shared<DoFHandler<dim>>(*tria);
  dof_handler->distribute_dofs(fe);
  std::shared_ptr<Mapping<dim>> mapping = std::make_shared<MappingQ<dim>>(1);

  const IndexSet locally_owned_dofs = dof_handler->locally_owned_dofs();
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(*dof_handler);

  // Present void fraction: eps(x) = 0.4 + 0.2 x, integral = 0.5
  std::shared_ptr<GlobalVectorType> present_void_fraction =
    std::make_shared<GlobalVectorType>(locally_owned_dofs,
                                       locally_relevant_dofs,
                                       mpi_communicator);
  interpolate_in_ghosted_vector(*mapping,
                                *dof_handler,
                                LinearVoidFraction<dim>(0.4, 0.2),
                                *present_void_fraction);

  // Previous void fractions for n-1 and for n-2, both constant fields of 0.3
  // and 0.2, respectively
  std::shared_ptr<std::vector<GlobalVectorType>> previous_void_fraction =
    std::make_shared<std::vector<GlobalVectorType>>(2);
  const std::vector<double> previous_values = {0.3, 0.2};
  for (unsigned int i = 0; i < previous_void_fraction->size(); ++i)
    {
      (*previous_void_fraction)[i].reinit(locally_owned_dofs,
                                          locally_relevant_dofs,
                                          mpi_communicator);
      interpolate_in_ghosted_vector(*mapping,
                                    *dof_handler,
                                    Functions::ConstantFunction<dim>(
                                      previous_values[i]),
                                    (*previous_void_fraction)[i]);
    }

  {
    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout,
                                            FluidDynamicsFormulation::VANS);

    // Register the void fraction, as done by the VANS solvers in setup_dofs()
    multiphysics.set_dof_handler(PhysicsID::void_fraction, dof_handler);
    multiphysics.set_mapping(PhysicsID::void_fraction, mapping);
    multiphysics.set_solution(PhysicsID::void_fraction, present_void_fraction);
    multiphysics.set_previous_solutions(PhysicsID::void_fraction,
                                        previous_void_fraction);

    deallog << "  Void fraction registered under VANS" << std::endl;

    // From here on, the fields are only accessed through the interface
    const DoFHandler<dim> &interface_dof_handler =
      multiphysics.get_dof_handler(PhysicsID::void_fraction);
    const Mapping<dim> &interface_mapping =
      multiphysics.get_mapping(PhysicsID::void_fraction);

    deallog << "  Integral of the present void fraction (expected 0.5)  : "
            << integrate_field(interface_mapping,
                               interface_dof_handler,
                               multiphysics.get_solution(
                                 PhysicsID::void_fraction))
            << std::endl;

    const std::vector<GlobalVectorType> &interface_previous_void_fraction =
      multiphysics.get_previous_solutions(PhysicsID::void_fraction);
    deallog << "  Integral of the previous void fraction 0 (expected 0.3) : "
            << integrate_field(interface_mapping,
                               interface_dof_handler,
                               interface_previous_void_fraction[0])
            << std::endl;
    deallog << "  Integral of the previous void fraction 1 (expected 0.2) : "
            << integrate_field(interface_mapping,
                               interface_dof_handler,
                               interface_previous_void_fraction[1])
            << std::endl;

    // The interface must store the registered vector itself, not a copy
    const bool same_vector =
      (&multiphysics.get_solution(PhysicsID::void_fraction) ==
       present_void_fraction.get());
    deallog << "  Interface shares the registered vector             : "
            << (same_vector ? "yes" : "no") << std::endl;

    // A new void fraction computed by the solver must be seen through the
    // interface without registering it again: eps(x) = 0.8 - 0.2 x,
    // integral = 0.7
    interpolate_in_ghosted_vector(*mapping,
                                  *dof_handler,
                                  LinearVoidFraction<dim>(0.8, -0.2),
                                  *present_void_fraction);
    deallog << "  Integral of the updated void fraction (expected 0.7)  : "
            << integrate_field(interface_mapping,
                               interface_dof_handler,
                               multiphysics.get_solution(
                                 PhysicsID::void_fraction))
            << std::endl;

    // The void fraction is not a block physics: the block getters must keep
    // refusing it, even under VANS.
    try
      {
        multiphysics.get_block_solution(PhysicsID::void_fraction);
        deallog << "  get_block_solution did NOT throw (unexpected)"
                << std::endl;
      }
    catch (const std::exception &)
      {
        deallog << "  get_block_solution correctly threw under VANS"
                << std::endl;
      }
  }

  // With the standard formulation, the void fraction does not exist and the
  // interface must refuse it.
  {
    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout,
                                            FluidDynamicsFormulation::standard);
    try
      {
        multiphysics.get_solution(PhysicsID::void_fraction);
        deallog << "  get_solution did NOT throw under the standard "
                   "formulation (unexpected)"
                << std::endl;
      }
    catch (const std::exception &)
      {
        deallog << "  get_solution correctly threw under the standard "
                   "formulation"
                << std::endl;
      }
  }
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      MPILogInitAll                    all;

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
