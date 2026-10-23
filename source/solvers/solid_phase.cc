// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/utilities.h>

#include <solvers/solid_phase.h>

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
SolidPhase<dim>::SolidPhase(
  MultiphysicsInterface<dim>      *multiphysics_interface,
  const SimulationParameters<dim> &p_simulation_parameters,
  const std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                           &p_triangulation,
  const std::shared_ptr<SimulationControl> &p_simulation_control)
  : AuxiliaryPhysics<dim, GlobalVectorType>(
      p_simulation_parameters.physics_solving_strategy.at(
        PhysicsID::solid_phase))
  , multiphysics(multiphysics_interface)
  , computing_timer(p_triangulation->get_mpi_communicator(),
                    this->pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
  , simulation_parameters(p_simulation_parameters)
  , triangulation(p_triangulation)
  , simulation_control(p_simulation_control)
  , dof_handler(std::make_shared<DoFHandler<dim>>(*triangulation))
  , solid_velocity_extractor(0)
  , solid_volume_fraction_extractor(dim)
{
  AssertThrow(!simulation_parameters.mesh.simplex,
              ExcMessage("The solid phase physics does not support simplex "
                         "meshes yet."));

  this->pcout << std::setprecision(simulation_control->get_log_precision())
              << std::scientific;

  const unsigned int velocity_degree =
    simulation_parameters.fem_parameters.solid_velocity_degree;
  const unsigned int volume_fraction_degree =
    simulation_parameters.fem_parameters.solid_volume_fraction_degree;

  // Solid velocity (dim components) followed by the solid volume fraction
  fe = std::make_shared<FESystem<dim>>(FE_Q<dim>(velocity_degree),
                                       dim,
                                       FE_Q<dim>(volume_fraction_degree),
                                       1);
  // fe->degree is the highest degree of the two fields
  mapping         = std::make_shared<MappingQ<dim>>(fe->degree);
  cell_quadrature = std::make_shared<QGauss<dim>>(fe->degree + 1);
  face_quadrature = std::make_shared<QGauss<dim - 1>>(fe->degree + 1);

  present_solution = std::make_shared<GlobalVectorType>();
  solution_transfer =
    std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);

  previous_solutions = std::make_shared<std::vector<GlobalVectorType>>(
    maximum_number_of_previous_solutions());
  previous_solutions_transfer.reserve(previous_solutions->size());
  for (unsigned int i = 0; i < previous_solutions->size(); ++i)
    previous_solutions_transfer.emplace_back(
      SolutionTransfer<dim, GlobalVectorType>(*dof_handler));
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
SolidPhase<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;

  std::vector<std::string> solution_names(dim, "solid_velocity");
  solution_names.emplace_back("solid_volume_fraction");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  component_interpretation.push_back(
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
SolidPhase<dim>::finish_simulation()
{
  if (simulation_parameters.timer.type == Parameters::Timer::Type::end)
    {
      announce_string(this->pcout, "Solid Phase");
      this->pcout << std::defaultfloat;
      computing_timer.print_summary();
      this->pcout << std::scientific;
    }
}

template <int dim>
void
SolidPhase<dim>::percolate_time_vectors()
{
  for (unsigned int i = previous_solutions->size() - 1; i > 0; --i)
    (*previous_solutions)[i] = (*previous_solutions)[i - 1];
  (*previous_solutions)[0] = *present_solution;
}

template <int dim>
void
SolidPhase<dim>::postprocess([[maybe_unused]] bool first_iteration)
{}

template <int dim>
void
SolidPhase<dim>::output_per_iteration_timer()
{
  announce_string(this->pcout, "Solid Phase");
  this->pcout << std::defaultfloat;
  computing_timer.print_summary();
  this->pcout << std::scientific;
  computing_timer.reset();
}

template <int dim>
void
SolidPhase<dim>::pre_mesh_adaptation()
{
  solution_transfer->prepare_for_coarsening_and_refinement(*present_solution);
  for (unsigned int i = 0; i < previous_solutions->size(); ++i)
    previous_solutions_transfer[i].prepare_for_coarsening_and_refinement(
      (*previous_solutions)[i]);
}

template <int dim>
void
SolidPhase<dim>::post_mesh_adaptation()
{
  const auto mpi_communicator = triangulation->get_mpi_communicator();

  GlobalVectorType tmp(locally_owned_dofs, mpi_communicator);
  solution_transfer->interpolate(tmp);
  nonzero_constraints.distribute(tmp);
  *present_solution = tmp;

  for (unsigned int i = 0; i < previous_solutions->size(); ++i)
    {
      GlobalVectorType tmp_previous_solution(locally_owned_dofs,
                                             mpi_communicator);
      previous_solutions_transfer[i].interpolate(tmp_previous_solution);
      nonzero_constraints.distribute(tmp_previous_solution);
      (*previous_solutions)[i] = tmp_previous_solution;
    }
}

template <int dim>
void
SolidPhase<dim>::compute_error_estimate(
  [[maybe_unused]] const std::pair<const Variable,
                                   Parameters::MultipleAdaptationParameters>
                                         &ivar,
  [[maybe_unused]] dealii::Vector<float> &estimated_error_per_cell)
{
  AssertThrow(false,
              ExcMessage("Mesh adaptation based on solid phase variables is "
                         "not implemented yet."));
}

template <int dim>
void
SolidPhase<dim>::write_checkpoint()
{
  std::vector<const GlobalVectorType *> sol_set_transfer;

  solution_transfer =
    std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);

  sol_set_transfer.emplace_back(&(*present_solution));
  for (const auto &previous_solution : *previous_solutions)
    sol_set_transfer.emplace_back(&previous_solution);

  solution_transfer->prepare_for_serialization(sol_set_transfer);
}

template <int dim>
void
SolidPhase<dim>::read_checkpoint()
{
  const auto mpi_communicator = triangulation->get_mpi_communicator();
  this->pcout << "Reading solid phase checkpoint" << std::endl;

  std::vector<GlobalVectorType *> input_vectors(1 + previous_solutions->size());
  GlobalVectorType distributed_system(locally_owned_dofs, mpi_communicator);
  input_vectors[0] = &distributed_system;

  std::vector<GlobalVectorType> distributed_previous_solutions;
  distributed_previous_solutions.reserve(previous_solutions->size());
  for (unsigned int i = 0; i < previous_solutions->size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        GlobalVectorType(locally_owned_dofs, mpi_communicator));
      input_vectors[i + 1] = &distributed_previous_solutions[i];
    }

  solution_transfer->deserialize(input_vectors);

  *present_solution = distributed_system;
  for (unsigned int i = 0; i < previous_solutions->size(); ++i)
    (*previous_solutions)[i] = distributed_previous_solutions[i];
}

template <int dim>
std::vector<OutputStructTableHandler>
SolidPhase<dim>::gather_tables()
{
  return std::vector<OutputStructTableHandler>();
}

template <int dim>
void
SolidPhase<dim>::define_constraints()
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
SolidPhase<dim>::setup_dofs()
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

  this->pcout << "   Number of solid phase degrees of freedom: "
              << dof_handler->n_dofs() << std::endl;

  // Share the solid phase fields with the other physics (e.g. VANS)
  multiphysics->set_dof_handler(PhysicsID::solid_phase, dof_handler);
  multiphysics->set_mapping(PhysicsID::solid_phase, mapping);
  multiphysics->set_solution(PhysicsID::solid_phase, present_solution);
  multiphysics->set_previous_solutions(PhysicsID::solid_phase,
                                       previous_solutions);
}

template <int dim>
void
SolidPhase<dim>::set_initial_conditions()
{
  // Zero until solid phase initial conditions are added to the parameters
  *present_solution = 0.;
  percolate_time_vectors();
}

template <int dim>
void
SolidPhase<dim>::solve_linear_system()
{
  AssertThrow(false,
              ExcMessage("The solid phase linear solver is not implemented "
                         "yet."));
}

template <int dim>
void
SolidPhase<dim>::assemble_system_matrix()
{
  AssertThrow(false,
              ExcMessage("The solid phase system matrix assembly is not "
                         "implemented yet."));
}

template <int dim>
void
SolidPhase<dim>::assemble_system_rhs()
{
  AssertThrow(false,
              ExcMessage("The solid phase right-hand side assembly is not "
                         "implemented yet."));
}

template <int dim>
void
SolidPhase<dim>::output_newton_update_norms()
{
  this->pcout << "\t  ||dx||_L2 = " << std::setw(6) << newton_update.l2_norm()
              << std::setw(6)
              << "\t  ||dx||_Linfty = " << newton_update.linfty_norm()
              << std::endl;
}

template <int dim>
double
SolidPhase<dim>::get_residual_rescale_metric() const
{
  return simulation_parameters.linear_solver.at(PhysicsID::solid_phase)
             .rescale_residual_by_volume ?
           std::sqrt(GridTools::volume(*triangulation, *mapping)) :
           1.;
}

template class SolidPhase<2>;
template class SolidPhase<3>;
