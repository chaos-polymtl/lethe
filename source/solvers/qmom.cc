// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/dimensionality.h>
#include <core/parameters_multiphysics.h>
#include <core/utilities.h>

#include <solvers/initial_conditions.h>
#include <solvers/qmom.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/data_component_interpretation.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string>
#include <utility>

template <int dim>
QMOM<dim>::QMOM(
  MultiphysicsInterface<dim>      *multiphysics_interface,
  const SimulationParameters<dim> &p_simulation_parameters,
  const std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                           &p_triangulation,
  const std::shared_ptr<SimulationControl> &p_simulation_control)
  : AuxiliaryPhysics<dim, GlobalVectorType>(
      p_simulation_parameters.physics_solving_strategy.at(PhysicsID::qmom))
  , n_moments(
      2 * p_simulation_parameters.multiphysics.qmom_parameters.number_of_nodes)
  , n_components(n_scalars + n_moments)
  , multiphysics(multiphysics_interface)
  , computing_timer(p_triangulation->get_mpi_communicator(),
                    this->pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
  , simulation_parameters(p_simulation_parameters)
  , triangulation(p_triangulation)
  , simulation_control(p_simulation_control)
  , dof_handler(std::make_shared<DoFHandler<dim>>(*triangulation))
  , ethanol_extractor(ethanol_component)
  , free_lipids_extractor(free_lipids_component)
{
  AssertThrow(!simulation_parameters.mesh.simplex,
              ExcMessage("The QMOM physics does not support simplex meshes "
                         "yet. Use a quadrilateral or hexahedral mesh."));

  // Relative tolerance used to compare parameter values
  constexpr double parameter_relative_tolerance = 1e-12;
  const auto       values_are_equal = [](const double first_value,
                                   const double second_value) {
    return std::abs(first_value - second_value) <=
           parameter_relative_tolerance *
             std::max({1., std::abs(first_value), std::abs(second_value)});
  };

  // The QMOM phase properties and operating conditions are given in SI units
  // and are not rescaled by the dimensionality subsection
  const Parameters::Dimensionality &dimensionality =
    simulation_parameters.dimensionality;
  const bool dimensionality_is_si = values_are_equal(dimensionality.length, 1.) &&
                                    values_are_equal(dimensionality.mass, 1.) &&
                                    values_are_equal(dimensionality.time, 1.) &&
                                    values_are_equal(dimensionality.temperature,
                                                     1.);
  AssertThrow(dimensionality_is_si,
              ExcMessage(
                "The QMOM physics only supports SI units, but the "
                "dimensionality subsection sets length = " +
                std::to_string(dimensionality.length) +
                ", mass = " + std::to_string(dimensionality.mass) +
                ", time = " + std::to_string(dimensionality.time) +
                " and temperature = " +
                std::to_string(dimensionality.temperature) +
                ". Remove the dimensionality subsection or set these "
                "parameters to 1."));

  // The initial lipid concentration in the organic phase is set both in the
  // QMOM subsection and in the QMOM initial conditions
  const Parameters::QMOM &qmom_parameters =
    simulation_parameters.multiphysics.qmom_parameters;
  const Parameters::QMOMInitialConditions &qmom_initial_conditions =
    simulation_parameters.initial_condition->qmom;
  AssertThrow(
    values_are_equal(qmom_parameters.initial_concentration_organic,
                     qmom_initial_conditions.initial_lipid_concentration),
    ExcMessage(
      "The initial lipid concentration in the organic phase is set to " +
      std::to_string(qmom_parameters.initial_concentration_organic) +
      " by \"initial concentration organic\" in the QMOM subsection, but to " +
      std::to_string(qmom_initial_conditions.initial_lipid_concentration) +
      " by \"initial lipid concentration\" in the QMOM subsection of the "
      "initial conditions. Set both parameters to the same value."));

  // The QMOM initial conditions are not applied yet: all the fields are
  // initialized to zero
  const Parameters::QMOMInitialConditions default_qmom_initial_conditions;
  AssertThrow(
    qmom_initial_conditions.distribution_shape ==
        default_qmom_initial_conditions.distribution_shape &&
      values_are_equal(
        qmom_initial_conditions.initial_lipid_concentration,
        default_qmom_initial_conditions.initial_lipid_concentration),
    ExcMessage(
      "The QMOM initial conditions are not implemented yet: all the QMOM "
      "fields are initialized to zero. Remove the QMOM subsection from the "
      "initial conditions subsection."));

  this->pcout << std::setprecision(simulation_control->get_log_precision())
              << std::scientific;

  // Ethanol and free lipids concentrations followed by the moments. All the
  // components share the same degree, so that they have the same support
  // points.
  fe = std::make_unique<FESystem<dim>>(
    FE_Q<dim>(simulation_parameters.fem_parameters.qmom_degree), n_components);
  mapping         = std::make_shared<MappingQ<dim>>(fe->degree);
  cell_quadrature = std::make_unique<QGauss<dim>>(fe->degree + 1);
  face_quadrature = std::make_unique<QGauss<dim - 1>>(fe->degree + 1);

  moment_extractors.reserve(n_moments);
  for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
    moment_extractors.emplace_back(first_moment_component + i_moment);

  present_solution = std::make_shared<GlobalVectorType>();
  solution_transfer =
    std::make_unique<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);

  previous_solutions = std::make_shared<std::vector<GlobalVectorType>>(
    maximum_number_of_previous_solutions());
  previous_solutions_transfer.reserve(previous_solutions->size());
  for (unsigned int i_previous_solution = 0;
       i_previous_solution < previous_solutions->size();
       ++i_previous_solution)
    previous_solutions_transfer.emplace_back(
      SolutionTransfer<dim, GlobalVectorType>(*dof_handler));
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
QMOM<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;

  std::vector<std::string> solution_names(n_components);
  solution_names[ethanol_component]     = "ethanol_concentration";
  solution_names[free_lipids_component] = "free_lipids_concentration";
  for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
    solution_names[first_moment_component + i_moment] =
      "moment_" + std::to_string(i_moment);

  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation(n_components,
                             DataComponentInterpretation::component_is_scalar);

  solution_output_structs.emplace_back(
    std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
    *dof_handler,
    *present_solution,
    solution_names,
    component_interpretation);
  return solution_output_structs;
}

template <int dim>
void
QMOM<dim>::finish_simulation()
{
  if (simulation_parameters.timer.type == Parameters::Timer::Type::end)
    {
      announce_string(this->pcout, "QMOM");
      this->pcout << std::defaultfloat;
      computing_timer.print_summary();
      this->pcout << std::scientific;
    }
}

template <int dim>
void
QMOM<dim>::percolate_time_vectors()
{
  for (unsigned int i_previous_solution = previous_solutions->size() - 1;
       i_previous_solution > 0;
       --i_previous_solution)
    (*previous_solutions)[i_previous_solution] =
      (*previous_solutions)[i_previous_solution - 1];
  (*previous_solutions)[0] = *present_solution;
}

template <int dim>
void
QMOM<dim>::postprocess([[maybe_unused]] bool first_iteration)
{}

template <int dim>
void
QMOM<dim>::output_per_iteration_timer()
{
  announce_string(this->pcout, "QMOM");
  this->pcout << std::defaultfloat;
  computing_timer.print_summary();
  this->pcout << std::scientific;
  computing_timer.reset();
}

template <int dim>
void
QMOM<dim>::pre_mesh_adaptation()
{
  solution_transfer->prepare_for_coarsening_and_refinement(*present_solution);
  for (unsigned int i_previous_solution = 0;
       i_previous_solution < previous_solutions->size();
       ++i_previous_solution)
    previous_solutions_transfer[i_previous_solution]
      .prepare_for_coarsening_and_refinement(
        (*previous_solutions)[i_previous_solution]);
}

template <int dim>
void
QMOM<dim>::post_mesh_adaptation()
{
  const auto mpi_communicator = triangulation->get_mpi_communicator();

  GlobalVectorType interpolated_solution(locally_owned_dofs, mpi_communicator);
  solution_transfer->interpolate(interpolated_solution);
  nonzero_constraints.distribute(interpolated_solution);
  *present_solution = interpolated_solution;

  for (unsigned int i_previous_solution = 0;
       i_previous_solution < previous_solutions->size();
       ++i_previous_solution)
    {
      GlobalVectorType interpolated_previous_solution(locally_owned_dofs,
                                                      mpi_communicator);
      previous_solutions_transfer[i_previous_solution].interpolate(
        interpolated_previous_solution);
      nonzero_constraints.distribute(interpolated_previous_solution);
      (*previous_solutions)[i_previous_solution] =
        interpolated_previous_solution;
    }
}

template <int dim>
void
QMOM<dim>::compute_error_estimate(
  [[maybe_unused]] const std::pair<const Variable,
                                   Parameters::MultipleAdaptationParameters>
                                         &ivar,
  [[maybe_unused]] dealii::Vector<float> &estimated_error_per_cell)
{
  AssertThrow(false,
              ExcMessage("Mesh adaptation based on QMOM variables is not "
                         "implemented yet. Remove the QMOM variables from the "
                         "mesh adaptation subsection."));
}

template <int dim>
void
QMOM<dim>::write_checkpoint()
{
  std::vector<const GlobalVectorType *> sol_set_transfer;

  solution_transfer =
    std::make_unique<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);

  sol_set_transfer.emplace_back(&(*present_solution));
  for (const auto &previous_solution : *previous_solutions)
    sol_set_transfer.emplace_back(&previous_solution);

  solution_transfer->prepare_for_serialization(sol_set_transfer);
}

template <int dim>
void
QMOM<dim>::read_checkpoint()
{
  const auto mpi_communicator = triangulation->get_mpi_communicator();
  this->pcout << "Reading QMOM checkpoint" << std::endl;

  std::vector<GlobalVectorType *> input_vectors(1 + previous_solutions->size());
  GlobalVectorType distributed_system(locally_owned_dofs, mpi_communicator);
  input_vectors[0] = &distributed_system;

  std::vector<GlobalVectorType> distributed_previous_solutions;
  distributed_previous_solutions.reserve(previous_solutions->size());
  for (unsigned int i_previous_solution = 0;
       i_previous_solution < previous_solutions->size();
       ++i_previous_solution)
    {
      distributed_previous_solutions.emplace_back(
        GlobalVectorType(locally_owned_dofs, mpi_communicator));
      input_vectors[i_previous_solution + 1] =
        &distributed_previous_solutions[i_previous_solution];
    }

  solution_transfer->deserialize(input_vectors);

  *present_solution = distributed_system;
  for (unsigned int i_previous_solution = 0;
       i_previous_solution < previous_solutions->size();
       ++i_previous_solution)
    (*previous_solutions)[i_previous_solution] =
      distributed_previous_solutions[i_previous_solution];
}

template <int dim>
std::vector<OutputStructTableHandler>
QMOM<dim>::gather_tables()
{
  return std::vector<OutputStructTableHandler>();
}

template <int dim>
void
QMOM<dim>::define_constraints()
{
  nonzero_constraints.clear();
  nonzero_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(*dof_handler, nonzero_constraints);
  nonzero_constraints.close();

  zero_constraints.clear();
  zero_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(*dof_handler, zero_constraints);
  zero_constraints.close();
}

template <int dim>
void
QMOM<dim>::setup_dofs()
{
  dof_handler->distribute_dofs(*fe);
  DoFRenumbering::Cuthill_McKee(*dof_handler);

  const auto mpi_communicator = triangulation->get_mpi_communicator();

  locally_owned_dofs    = dof_handler->locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(*dof_handler);

  present_solution->reinit(locally_owned_dofs,
                           locally_relevant_dofs,
                           mpi_communicator);
  for (auto &solution : *previous_solutions)
    solution.reinit(locally_owned_dofs,
                    locally_relevant_dofs,
                    mpi_communicator);

  evaluation_point.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  newton_update.reinit(locally_owned_dofs, mpi_communicator);
  local_evaluation_point.reinit(locally_owned_dofs, mpi_communicator);

  define_constraints();

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(*dof_handler,
                                  dsp,
                                  nonzero_constraints,
                                  /*keep_constrained_dofs = */ true);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  this->pcout << "   Number of QMOM degrees of freedom: "
              << dof_handler->n_dofs() << std::endl;

  // Share the QMOM fields with the other physics
  multiphysics->set_dof_handler(PhysicsID::qmom, dof_handler);
  multiphysics->set_mapping(PhysicsID::qmom, mapping);
  multiphysics->set_solution(PhysicsID::qmom, present_solution);
  multiphysics->set_previous_solutions(PhysicsID::qmom, previous_solutions);
}

template <int dim>
void
QMOM<dim>::set_initial_conditions()
{
  // Zero (no particles, no solute) until QMOM initial conditions are added to
  // the parameters
  *present_solution = 0.;
  percolate_time_vectors();
}

template <int dim>
void
QMOM<dim>::solve_linear_system()
{
  AssertThrow(false,
              ExcMessage("The QMOM linear solver is not implemented yet."));
}

template <int dim>
void
QMOM<dim>::assemble_system_matrix()
{
  AssertThrow(false,
              ExcMessage("The QMOM system matrix assembly is not implemented "
                         "yet."));
}

template <int dim>
void
QMOM<dim>::assemble_system_rhs()
{
  AssertThrow(false,
              ExcMessage("The QMOM right-hand side assembly is not "
                         "implemented yet."));
}

template <int dim>
void
QMOM<dim>::output_newton_update_norms()
{
  this->pcout << "\t  ||dx||_L2 = " << std::setw(6) << newton_update.l2_norm()
              << std::setw(6)
              << "\t  ||dx||_Linfty = " << newton_update.linfty_norm()
              << std::endl;
}

template <int dim>
double
QMOM<dim>::get_residual_rescale_metric() const
{
  return simulation_parameters.linear_solver.at(PhysicsID::qmom)
             .rescale_residual_by_volume ?
           std::sqrt(GridTools::volume(*triangulation, *mapping)) :
           1.;
}

template class QMOM<2>;
template class QMOM<3>;
