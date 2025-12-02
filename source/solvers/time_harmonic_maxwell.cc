// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/time_harmonic_maxwell.h>

using VectorType = GlobalVectorType;

template <int dim>
TimeHarmonicMaxwell<dim>::TimeHarmonicMaxwell(
  MultiphysicsInterface<dim>      *multiphysics_interface,
  const SimulationParameters<dim> &p_simulation_parameters,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> p_triangulation,
  std::shared_ptr<SimulationControl> p_simulation_control)
  : AuxiliaryPhysics<dim, GlobalVectorType>()
  , multiphysics(multiphysics_interface)
  , computing_timer(p_triangulation->get_mpi_communicator(),
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , simulation_parameters(p_simulation_parameters)
  , triangulation(p_triangulation)
  , simulation_control(std::move(p_simulation_control))
  , dof_handler_trial_interior(
      std::make_shared<DoFHandler<dim>>(*triangulation))
  , dof_handler_trial_skeleton(
      std::make_shared<DoFHandler<dim>>(*triangulation))
  , dof_handler_test(std::make_shared<DoFHandler<dim>>(*triangulation))
  , extractor_E_real(0)
  , extractor_E_imag(dim)
  , extractor_H_real(2 * dim)
  , extractor_H_imag(3 * dim)
{
  if (simulation_parameters.mesh.simplex)
    {
      // for simplex meshes
      AssertThrow(
        false,
        ExcMessage(
          "TimeHarmonicMaxwell solver not yet implemented for simplex meshes."));
    }
  else
    {
      AssertThrow(dim == 3,
                  ExcMessage(
                    "TimeHarmonicMaxwell only implemented for 3D problems."));

      AssertThrow(
        simulation_parameters.fem_parameters.electromagnetics_trial_order <
          simulation_parameters.fem_parameters.electromagnetics_test_order,
        ExcMessage(
          "The DPG method requires the test space to be of higher order than the trial space."));

      // Usual case, for quad/hex meshes
      fe_trial_interior = std::make_shared<FESystem<dim>>(
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order) ^
          dim,
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order) ^
          dim,
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order) ^
          dim,
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order) ^
          dim);
      fe_trial_skeleton = std::make_shared<FESystem<dim>>(
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order));
      fe_test = std::make_shared<FESystem<dim>>(
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order));
      mapping = std::make_shared<MappingQ<dim>>(fe_trial_interior->degree);
      cell_quadrature = std::make_shared<QGauss<dim>>(fe_test->degree + 1);
      face_quadrature = std::make_shared<QGauss<dim - 1>>(fe_test->degree + 1);
    }

  // Initialize solution shared_ptr
  present_solution = std::make_shared<GlobalVectorType>();

  // Allocate solution transfer
  solution_transfer = std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(
    *dof_handler_trial_interior);
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
TimeHarmonicMaxwell<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;

  // Interior output setup
  std::vector<std::string> solution_interior_names(dim, "E_real");
  for (unsigned int i = 0; i < dim; ++i)
    {
      solution_interior_names.emplace_back("E_imag");
    }
  for (unsigned int i = 0; i < dim; ++i)
    {
      solution_interior_names.emplace_back("H_real");
    }
  for (unsigned int i = 0; i < dim; ++i)
    {
      solution_interior_names.emplace_back("H_imag");
    }

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    solution_interior_data_component_interpretation(
      4 * dim, DataComponentInterpretation::component_is_part_of_vector);

  solution_output_structs.emplace_back(
    std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
    *this->dof_handler_trial_interior,
    *this->present_solution,
    solution_interior_names,
    solution_interior_data_component_interpretation);

  // Skeleton output setup
  // TODO:  it will need its own writer probably because we need to use
  // DataOutFaces object

  return solution_output_structs;
}

template <int dim>
double
TimeHarmonicMaxwell<dim>::calculate_L2_error()
{
  // TODO
  return 0.0;
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::finish_simulation()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::percolate_time_vectors()
{
  // No time-dependent vectors to percolate in time-harmonic Maxwell
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::modify_solution()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::update_boundary_conditions()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::postprocess(bool first_iteration)
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::pre_mesh_adaptation()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::post_mesh_adaptation()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::write_checkpoint()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::read_checkpoint()
{
  // TODO
}

template <int dim>
std::vector<OutputStructTableHandler>
TimeHarmonicMaxwell<dim>::gather_tables()
{
  // TODO
  return std::vector<OutputStructTableHandler>();
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::compute_kelly(
  const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                        &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::compute_energy_norm(
  const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                        &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::setup_dofs()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::setup_preconditioner()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::define_constraints()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::solve_linear_system()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::assemble_system_matrix()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::assemble_system_rhs()
{
  // TODO
}



template class TimeHarmonicMaxwell<2>;
template class TimeHarmonicMaxwell<3>;
