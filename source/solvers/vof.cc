// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof.h>

#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/error_estimator.h>

#include <sys/stat.h>

#include <cmath>

template <int dim>
VolumeOfFluid<dim>::VolumeOfFluid(
  MultiphysicsInterface<dim>      *multiphysics_interface,
  const SimulationParameters<dim> &p_simulation_parameters,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> p_triangulation,
  std::shared_ptr<SimulationControl> p_simulation_control)
  : AuxiliaryPhysics<dim, GlobalVectorType>(
      p_simulation_parameters.non_linear_solver.at(PhysicsID::VOF))
  , multiphysics(multiphysics_interface)
  , computing_timer(p_triangulation->get_mpi_communicator(),
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , simulation_parameters(p_simulation_parameters)
  , triangulation(p_triangulation)
  , simulation_control(std::move(p_simulation_control))
  , dof_handler(*triangulation)
  , sharpening_threshold(simulation_parameters.multiphysics.vof_parameters
                           .regularization_method.sharpening.threshold)
{
  AssertThrow(
    simulation_parameters.physical_properties_manager.get_number_of_fluids() ==
      2,
    InvalidNumberOfFluid(simulation_parameters.physical_properties_manager
                           .get_number_of_fluids()));

  AssertThrow(((simulation_parameters.fem_parameters.VOF_uses_dg &&
                simulation_parameters.multiphysics.vof_parameters
                    .regularization_method.regularization_method_type ==
                  Parameters::RegularizationMethodType::none) ||
               !simulation_parameters.fem_parameters.VOF_uses_dg),
              UnsupportedRegularization());

  if (simulation_parameters.fem_parameters.VOF_uses_dg &&
      this->simulation_parameters.post_processing.calculate_mass_conservation)
    {
      this->pcout
        << "Warning: DG-VOF and the geometric surface and volume computations are not compatible. Only the global volume will be monitored."
        << std::endl;
    }

  if (simulation_parameters.mesh.simplex)
    {
      // for simplex meshes
      fe = std::make_shared<FE_SimplexP<dim>>(
        simulation_parameters.fem_parameters.VOF_order);
      mapping         = std::make_shared<MappingFE<dim>>(*fe);
      cell_quadrature = std::make_shared<QGaussSimplex<dim>>(fe->degree + 1);
      face_quadrature =
        std::make_shared<QGaussSimplex<dim - 1>>(fe->degree + 1);
    }
  else
    {
      // Usual case, for quad/hex meshes
      if (simulation_parameters.fem_parameters.VOF_uses_dg)
        {
          fe = std::make_shared<FE_DGQ<dim>>(
            simulation_parameters.fem_parameters.VOF_order);
        }
      else
        {
          fe = std::make_shared<FE_Q<dim>>(
            simulation_parameters.fem_parameters.VOF_order);
        }
      // Mapping has to be at least Q1, but DGQ0 is allowed
      mapping = std::make_shared<MappingQ<dim>>(std::max(fe->degree, uint(1)));
      cell_quadrature = std::make_shared<QGauss<dim>>(fe->degree + 1);
      face_quadrature = std::make_shared<QGauss<dim - 1>>(fe->degree + 1);
    }

  // Allocate solution transfer
  solution_transfer =
    std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(dof_handler);

  // Set size of previous solutions using BDF schemes information
  previous_solutions.resize(maximum_number_of_previous_solutions());

  // Prepare previous solutions transfer
  previous_solutions_transfer.reserve(previous_solutions.size());
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions_transfer.emplace_back(
        SolutionTransfer<dim, GlobalVectorType>(this->dof_handler));
    }

  // Check the value of interface sharpness
  if (simulation_parameters.multiphysics.vof_parameters.regularization_method
        .sharpening.interface_sharpness < 1.0)
    this->pcout
      << "Warning: interface sharpness values smaller than 1 smooth the interface instead of sharpening it."
      << std::endl
      << "The interface sharpness value should be set between 1 and 2"
      << std::endl;


  // Change the behavior of the timer for situations when you don't want
  // outputs
  if (simulation_parameters.timer.type == Parameters::Timer::Type::none)
    this->computing_timer.disable_output();

  // Initialize the interface object for subequations to solve
  this->vof_subequations_interface =
    std::make_shared<VOFSubequationsInterface<dim>>(this->simulation_parameters,
                                                    this->pcout,
                                                    this->triangulation,
                                                    this->simulation_control);


  if (simulation_parameters.multiphysics.vof_parameters.regularization_method
        .geometric_interface_reinitialization.enable ||
      simulation_parameters.initial_condition
          ->vof_initial_condition_smoothing ==
        Parameters::VOFInitialConditionType::geometric)
    {
      /* For the VOF solver, the interface is defined as the iso-contour
         \f$\phi\f = 0.5$. Hence, for the  SignedDistanceSolver, we set the
         iso-level to 0.5. We also define the inside part of the domain as
         \f$\phi>0.5\f$. Since the SignedDistanceSolver considers the inside
         of the domain as \f$d<0\f$, we scale the phase fraction with a factor
         of -1.*/
      this->signed_distance_solver = std::make_shared<
        InterfaceTools::SignedDistanceSolver<dim, GlobalVectorType>>(
        triangulation,
        fe,
        simulation_parameters.multiphysics.vof_parameters.regularization_method
          .geometric_interface_reinitialization.max_reinitialization_distance,
        0.5,
        -1.0,
        simulation_parameters.multiphysics.vof_parameters.regularization_method
          .verbosity);
      this->signed_distance_transformation =
        SignedDistanceTransformationBase::model_cast(
          simulation_parameters.multiphysics.vof_parameters
            .regularization_method.geometric_interface_reinitialization);
    }
}


template <int dim>
void
VolumeOfFluid<dim>::assemble_rhs()
{
  assemble_system_rhs();
}


template <int dim>
void
VolumeOfFluid<dim>::setup_assemblers()
{
  AssertThrow(
    is_sdirk(this->simulation_control->get_assembly_method()) == false,
    ExcMessage("The SDIRK scheme is not yet supported for this physics"));

  this->assemblers.clear();

  // Time-stepping schemes
  if (is_bdf(this->simulation_control->get_assembly_method()))
    {
      this->assemblers.emplace_back(
        std::make_shared<VOFAssemblerBDF<dim>>(this->simulation_control));
    }

  // If the VOF solver uses DG, a different set of assemblers is used
  if (simulation_parameters.fem_parameters.VOF_uses_dg)
    {
      this->assemblers.emplace_back(
        std::make_shared<VOFAssemblerDGCore<dim>>());
      this->inner_face_assembler = std::make_shared<VOFAssemblerSIPG<dim>>();
    }
  else
    {
      // Core assembler
      this->assemblers.emplace_back(std::make_shared<VOFAssemblerCore<dim>>(
        this->simulation_control,
        this->simulation_parameters.fem_parameters,
        this->simulation_parameters.multiphysics.vof_parameters));

      // DCDD shock-capturing assembler
      if (this->simulation_parameters.stabilization.vof_dcdd_stabilization)
        this->assemblers.emplace_back(
          std::make_shared<VOFAssemblerDCDDStabilization<dim>>(
            this->simulation_control));
    }
}


template <int dim>
void
VolumeOfFluid<dim>::assemble_system_matrix()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");

  this->system_matrix = 0;
  setup_assemblers();

  if (simulation_parameters.fem_parameters.VOF_uses_dg)
    assemble_system_matrix_dg();
  else
    assemble_system_matrix_cg();
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_system_matrix_cg()
{
  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data =
    VOFScratchData<dim>(this->simulation_control,
                        this->simulation_parameters.physical_properties_manager,
                        *this->fe,
                        *this->cell_quadrature,
                        *this->face_quadrature,
                        *this->mapping,
                        dof_handler_fd->get_fe());

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &VolumeOfFluid::assemble_local_system_matrix,
                  &VolumeOfFluid::copy_local_matrix_to_global_matrix,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  this->system_matrix.compress(VectorOperation::add);
}


template <int dim>
void
VolumeOfFluid<dim>::assemble_system_matrix_dg()
{
  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data =
    VOFScratchData<dim>(this->simulation_control,
                        this->simulation_parameters.physical_properties_manager,
                        *this->fe,
                        *this->cell_quadrature,
                        *this->face_quadrature,
                        *this->mapping,
                        dof_handler_fluid->get_fe());

  StabilizedDGMethodsCopyData copy_data(this->fe->n_dofs_per_cell(),
                                        this->cell_quadrature->size());

  // We first wrap the assembly of the matrix within a cell_worker lambda
  // function. This is only done for compatibility reasons with the MeshWorker
  // paradigm and does not have any functional purpose.
  const auto cell_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        VOFScratchData<dim>                                  &scratch_data,
        StabilizedDGMethodsCopyData                          &copy_data) {
      this->assemble_local_system_matrix(cell, scratch_data, copy_data);
    };

  const auto boundary_worker =
    [&]([[maybe_unused]] const typename DoFHandler<dim>::active_cell_iterator
                                                     &cell,
        [[maybe_unused]] const unsigned int          &face_no,
        [[maybe_unused]] VOFScratchData<dim>         &scratch_data,
        [[maybe_unused]] StabilizedDGMethodsCopyData &copy_data) {};

  const auto face_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                   &face_no,
        const unsigned int                                   &sub_face_no,
        const typename DoFHandler<dim>::active_cell_iterator &neigh_cell,
        const unsigned int                                   &neigh_face_no,
        const unsigned int                                   &neigh_sub_face_no,
        VOFScratchData<dim>                                  &scratch_data,
        StabilizedDGMethodsCopyData                          &copy_data) {
      scratch_data.reinit_internal_face(cell,
                                        face_no,
                                        sub_face_no,
                                        neigh_cell,
                                        neigh_face_no,
                                        neigh_sub_face_no,
                                        this->evaluation_point);

      // Pad copy_data memory for the internal faces elementary matrices
      // BB note : Array could be pre-allocated
      copy_data.face_data.emplace_back();
      auto &copy_data_face = copy_data.face_data.back();
      copy_data_face.face_matrix.reinit(scratch_data.n_interface_dofs,
                                        scratch_data.n_interface_dofs);
      copy_data_face.joint_dof_indices =
        scratch_data.fe_interface_values_vof.get_interface_dof_indices();

      // Gather velocity information at the face to advect properly
      // Get the cell that corresponds to the fluid dynamics
      typename DoFHandler<dim>::active_cell_iterator velocity_cell(
        &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

      // Reinit the internal face velocity within the scratch data
      reinit_face_velocity_with_adequate_solution(velocity_cell,
                                                  face_no,
                                                  scratch_data);

      this->inner_face_assembler->assemble_matrix(scratch_data, copy_data);
    };

  const auto copier = [&](const StabilizedDGMethodsCopyData &copy_data) {
    this->copy_local_matrix_to_global_matrix(copy_data);

    const AffineConstraints<double> &constraints_used = this->zero_constraints;

    for (const auto &cdf : copy_data.face_data)
      {
        constraints_used.distribute_local_to_global(cdf.face_matrix,
                                                    cdf.joint_dof_indices,
                                                    system_matrix);
      }
  };

  MeshWorker::mesh_loop(this->dof_handler.begin_active(),
                        this->dof_handler.end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once |
                          MeshWorker::assemble_ghost_faces_both,
                        boundary_worker,
                        face_worker);
}


template <int dim>
void
VolumeOfFluid<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  VOFScratchData<dim>                                  &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(cell, this->evaluation_point, this->previous_solutions);

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*this->triangulation), cell->level(), cell->index(), dof_handler_fd);

  if (multiphysics->fluid_dynamics_is_block())
    {
      // Check if the post processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::FluidDynamicsInitialConditionType::
              average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_time_average_solution(
              PhysicsID::fluid_dynamics),
            *multiphysics->get_block_previous_solutions(
              PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_block_previous_solutions(
              PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
    }
  else
    {
      // Check if the post-processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::FluidDynamicsInitialConditionType::
              average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_time_average_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_previous_solutions(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_previous_solutions(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
    }

  copy_data.reset();

  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_matrix(scratch_data, copy_data);
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
VolumeOfFluid<dim>::copy_local_matrix_to_global_matrix(
  const StabilizedMethodsCopyData &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_matrix,
                                              copy_data.local_dof_indices,
                                              this->system_matrix);
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");

  this->system_rhs = 0;
  setup_assemblers();

  if (simulation_parameters.fem_parameters.VOF_uses_dg)
    assemble_system_rhs_dg();

  else
    assemble_system_rhs_cg();
}


template <int dim>
void
VolumeOfFluid<dim>::assemble_system_rhs_cg()
{
  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data =
    VOFScratchData<dim>(this->simulation_control,
                        this->simulation_parameters.physical_properties_manager,
                        *this->fe,
                        *this->cell_quadrature,
                        *this->face_quadrature,
                        *this->mapping,
                        dof_handler_fd->get_fe());

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &VolumeOfFluid::assemble_local_system_rhs,
                  &VolumeOfFluid::copy_local_rhs_to_global_rhs,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_system_rhs_dg()
{
  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data =
    VOFScratchData<dim>(this->simulation_control,
                        this->simulation_parameters.physical_properties_manager,
                        *this->fe,
                        *this->cell_quadrature,
                        *this->face_quadrature,
                        *this->mapping,
                        dof_handler_fluid->get_fe());

  StabilizedDGMethodsCopyData copy_data(this->fe->n_dofs_per_cell(),
                                        this->cell_quadrature->size());

  const auto cell_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        VOFScratchData<dim>                                  &scratch_data,
        StabilizedDGMethodsCopyData                          &copy_data) {
      this->assemble_local_system_rhs(cell, scratch_data, copy_data);
    };

  const auto boundary_worker =
    [&]([[maybe_unused]] const typename DoFHandler<dim>::active_cell_iterator
                                                     &cell,
        [[maybe_unused]] const unsigned int          &face_no,
        [[maybe_unused]] VOFScratchData<dim>         &scratch_data,
        [[maybe_unused]] StabilizedDGMethodsCopyData &copy_data) {};

  const auto face_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                   &face_no,
        const unsigned int                                   &sub_face_no,
        const typename DoFHandler<dim>::active_cell_iterator &neigh_cell,
        const unsigned int                                   &neigh_face_no,
        const unsigned int                                   &neigh_sub_face_no,
        VOFScratchData<dim>                                  &scratch_data,
        StabilizedDGMethodsCopyData                          &copy_data)

  {
    scratch_data.reinit_internal_face(cell,
                                      face_no,
                                      sub_face_no,
                                      neigh_cell,
                                      neigh_face_no,
                                      neigh_sub_face_no,
                                      this->evaluation_point);

    copy_data.face_data.emplace_back();
    auto &copy_data_face = copy_data.face_data.back();
    copy_data_face.joint_dof_indices =
      scratch_data.fe_interface_values_vof.get_interface_dof_indices();
    copy_data_face.face_rhs.reinit(scratch_data.n_interface_dofs);

    // Gather velocity information at the face to advect properly.
    // First gather the dof handler for the fluid dynamics.
    const DoFHandler<dim> *dof_handler_fluid =
      multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
    // Get the cell that corresponds to the fluid dynamics
    typename DoFHandler<dim>::active_cell_iterator velocity_cell(
      &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

    reinit_face_velocity_with_adequate_solution(velocity_cell,
                                                face_no,
                                                scratch_data);

    this->inner_face_assembler->assemble_rhs(scratch_data, copy_data);
  };


  const auto copier = [&](const StabilizedDGMethodsCopyData &c) {
    this->copy_local_rhs_to_global_rhs(c);
    const AffineConstraints<double> &constraints_used = this->zero_constraints;
    for (const auto &cdf : c.face_data)
      {
        constraints_used.distribute_local_to_global(cdf.face_rhs,
                                                    cdf.joint_dof_indices,
                                                    system_rhs);
      }
  };

  MeshWorker::mesh_loop(this->dof_handler.begin_active(),
                        this->dof_handler.end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once |
                          MeshWorker::assemble_ghost_faces_both,
                        boundary_worker,
                        face_worker);
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  VOFScratchData<dim>                                  &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(cell, this->evaluation_point, this->previous_solutions);

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*this->triangulation), cell->level(), cell->index(), dof_handler_fd);

  if (multiphysics->fluid_dynamics_is_block())
    {
      // Check if the post-processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::FluidDynamicsInitialConditionType::
              average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_time_average_solution(
              PhysicsID::fluid_dynamics),
            *multiphysics->get_block_previous_solutions(
              PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_block_previous_solutions(
              PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
    }
  else
    {
      // Check if the post-processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::FluidDynamicsInitialConditionType::
              average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_time_average_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_previous_solutions(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_previous_solutions(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
    }

  copy_data.reset();

  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_rhs(scratch_data, copy_data);
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
VolumeOfFluid<dim>::copy_local_rhs_to_global_rhs(
  const StabilizedMethodsCopyData &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_rhs,
                                              copy_data.local_dof_indices,
                                              this->system_rhs);
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
VolumeOfFluid<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;
  std::vector<std::string>                         solution_names(1, "phase");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    solution_component_interpretation(
      1, DataComponentInterpretation::component_is_scalar);

  // Phase fraction
  solution_output_structs.emplace_back(
    std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
    this->dof_handler,
    this->present_solution,
    solution_names,
    solution_component_interpretation);

  // Filter phase fraction
  std::vector<std::string> filtered_solution_names(1, "filtered_phase");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    filtered_solution_component_interpretation(
      1, DataComponentInterpretation::component_is_scalar);
  solution_output_structs.emplace_back(
    std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
    this->dof_handler,
    this->filtered_solution,
    filtered_solution_names,
    filtered_solution_component_interpretation);

  auto vof_parameters = this->simulation_parameters.multiphysics.vof_parameters;

  if ((vof_parameters.surface_tension_force.enable &&
       vof_parameters.surface_tension_force.output_vof_auxiliary_fields) ||
      vof_parameters.regularization_method.algebraic_interface_reinitialization
        .enable)
    {
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        projected_phase_fraction_gradient_component_interpretation(
          dim, DataComponentInterpretation::component_is_scalar);
      for (unsigned int i = 0; i < dim; ++i)
        projected_phase_fraction_gradient_component_interpretation[i] =
          DataComponentInterpretation::component_is_part_of_vector;

      std::vector<std::string> solution_names_new(dim,
                                                  "phase_fraction_gradient");

      solution_output_structs.emplace_back(
        std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
        this->vof_subequations_interface->get_dof_handler(
          VOFSubequationsID::phase_gradient_projection),
        this->vof_subequations_interface->get_solution(
          VOFSubequationsID::phase_gradient_projection),
        solution_names_new,
        projected_phase_fraction_gradient_component_interpretation);

      solution_output_structs.emplace_back(
        std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
        this->vof_subequations_interface->get_dof_handler(
          VOFSubequationsID::curvature_projection),
        this->vof_subequations_interface->get_solution(
          VOFSubequationsID::curvature_projection),
        std::vector<std::string>(1, "curvature"),
        std::vector<DataComponentInterpretation::DataComponentInterpretation>(
          1, DataComponentInterpretation::component_is_scalar));
    }

  if (simulation_parameters.multiphysics.vof_parameters.regularization_method
        .geometric_interface_reinitialization.enable)
    {
      if ((simulation_control->get_step_number() %
             simulation_parameters.multiphysics.vof_parameters
               .regularization_method.frequency !=
           0) ||
          (simulation_control->get_step_number() == 0))
        {
          TimerOutput::Scope t(this->computing_timer, "Signed distance output");

          signed_distance_solver->setup_dofs();

          signed_distance_solver->set_level_set_from_background_mesh(
            dof_handler, this->present_solution);

          signed_distance_solver->solve();
        }
      for (auto &output_struct : signed_distance_solver->gather_output_hook())
        solution_output_structs.push_back(output_struct);

      signed_distance_solver->output_interface_reconstruction(
        "interface_reconstruction_" +
          this->simulation_parameters.simulation_control.output_name,
        this->simulation_parameters.simulation_control.output_folder,
        simulation_control->get_current_time(),
        simulation_control->get_step_number());
    }
  return solution_output_structs;
}

template <int dim>
double
VolumeOfFluid<dim>::calculate_L2_error()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->cell_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values);

  const unsigned int  n_q_points = this->cell_quadrature->size();
  std::vector<double> phase_exact_solution(n_q_points);
  std::vector<double> phase_values(n_q_points);

  auto &exact_solution = simulation_parameters.analytical_solution->phase;
  exact_solution.set_time(this->simulation_control->get_current_time());

  double l2error = 0.;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);
          fe_values_vof.get_function_values(this->present_solution,
                                            phase_values);

          // Get the exact solution at all gauss points
          exact_solution.value_list(fe_values_vof.get_quadrature_points(),
                                    phase_exact_solution);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              double sim   = phase_values[q];
              double exact = phase_exact_solution[q];
              l2error += (sim - exact) * (sim - exact) * fe_values_vof.JxW(q);
            }
        }
    }
  l2error = Utilities::MPI::sum(l2error, mpi_communicator);
  return l2error;
}

template <int dim>
template <typename VectorType>
std::pair<Tensor<1, dim>, Tensor<1, dim>>
VolumeOfFluid<dim>::calculate_barycenter(const GlobalVectorType &solution,
                                         const VectorType       &solution_fd)
{
  const MPI_Comm mpi_communicator = this->triangulation->get_mpi_communicator();

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->cell_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values);

  std::shared_ptr<VolumeOfFluidFilterBase> filter =
    VolumeOfFluidFilterBase::model_cast(
      this->simulation_parameters.multiphysics.vof_parameters.phase_filter);

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  FEValues<dim> fe_values_fd(*this->mapping,
                             dof_handler_fd->get_fe(),
                             *this->cell_quadrature,
                             update_values);

  const unsigned int          n_q_points = this->cell_quadrature->size();
  std::vector<double>         phase_values(n_q_points);
  std::vector<Tensor<1, dim>> velocity_values(n_q_points);
  std::vector<Point<dim>>     quadrature_locations(n_q_points);

  const FEValuesExtractors::Vector velocity(0);

  Tensor<1, dim> barycenter_location;
  Tensor<1, dim> barycenter_velocity;
  double         volume = 0;


  std::map<field, std::vector<double>> fields;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);
          quadrature_locations = fe_values_vof.get_quadrature_points();
          fe_values_vof.get_function_values(solution, phase_values);

          // Get fluid dynamics active cell iterator
          typename DoFHandler<dim>::active_cell_iterator cell_fd(
            &(*(this->triangulation)),
            cell->level(),
            cell->index(),
            dof_handler_fd);

          fe_values_fd.reinit(cell_fd);
          fe_values_fd[velocity].get_function_values(solution_fd,
                                                     velocity_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              const double JxW = fe_values_vof.JxW(q);
              const double filtered_phase_value =
                filter->filter_phase(phase_values[q]);

              volume += (filtered_phase_value)*JxW;
              barycenter_location +=
                (filtered_phase_value)*quadrature_locations[q] * JxW;
              barycenter_velocity +=
                (filtered_phase_value)*velocity_values[q] * JxW;
            }
        }
    }

  volume = Utilities::MPI::sum(volume, mpi_communicator);
  barycenter_location =
    Utilities::MPI::sum(barycenter_location, mpi_communicator) / volume;
  barycenter_velocity =
    Utilities::MPI::sum(barycenter_velocity, mpi_communicator) / volume;

  return std::pair<Tensor<1, dim>, Tensor<1, dim>>(barycenter_location,
                                                   barycenter_velocity);
}

template std::pair<Tensor<1, 2>, Tensor<1, 2>>
VolumeOfFluid<2>::calculate_barycenter<GlobalVectorType>(
  const GlobalVectorType &solution,
  const GlobalVectorType &current_solution_fd);


template std::pair<Tensor<1, 3>, Tensor<1, 3>>
VolumeOfFluid<3>::calculate_barycenter<GlobalVectorType>(
  const GlobalVectorType &solution,
  const GlobalVectorType &current_solution_fd);

template std::pair<Tensor<1, 2>, Tensor<1, 2>>
VolumeOfFluid<2>::calculate_barycenter<GlobalBlockVectorType>(
  const GlobalVectorType      &solution,
  const GlobalBlockVectorType &current_solution_fd);


template std::pair<Tensor<1, 3>, Tensor<1, 3>>
VolumeOfFluid<3>::calculate_barycenter<GlobalBlockVectorType>(
  const GlobalVectorType      &solution,
  const GlobalBlockVectorType &current_solution_fd);


template <int dim>
template <typename VectorType>
void
VolumeOfFluid<dim>::calculate_volume_and_mass(
  const GlobalVectorType          &solution,
  const VectorType                &current_solution_fd,
  const Parameters::FluidIndicator monitored_fluid)
{
  const MPI_Comm mpi_communicator = this->triangulation->get_mpi_communicator();

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->cell_quadrature,
                              update_values | update_JxW_values);

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  QGauss<dim> quadrature_formula(dof_handler_fd->get_fe().degree + 1);

  FEValues<dim> fe_values_fd(*this->mapping,
                             dof_handler_fd->get_fe(),
                             quadrature_formula,
                             update_values);

  const FEValuesExtractors::Scalar pressure(dim);

  const unsigned int  n_q_points = this->cell_quadrature->size();
  std::vector<double> phase_values(n_q_points);
  std::vector<double> density_0(n_q_points);
  std::vector<double> density_1(n_q_points);
  std::vector<double> pressure_values(n_q_points);

  this->volume_monitored = 0.;
  this->mass_monitored   = 0.;

  // Physical properties
  const auto density_models =
    this->simulation_parameters.physical_properties_manager
      .get_density_vector();
  std::map<field, std::vector<double>> fields;
  fields.insert(
    std::pair<field, std::vector<double>>(field::pressure, n_q_points));

  for (const auto &cell_vof : this->dof_handler.active_cell_iterators())
    {
      if (cell_vof->is_locally_owned())
        {
          fe_values_vof.reinit(cell_vof);
          fe_values_vof.get_function_values(solution, phase_values);

          if (this->simulation_parameters.physical_properties_manager
                .field_is_required(field::pressure))
            {
              // Get fluid dynamics active cell iterator
              typename DoFHandler<dim>::active_cell_iterator cell_fd(
                &(*(this->triangulation)),
                cell_vof->level(),
                cell_vof->index(),
                dof_handler_fd);

              fe_values_fd.reinit(cell_fd);
              fe_values_fd[pressure].get_function_values(current_solution_fd,
                                                         pressure_values);
              set_field_vector(field::pressure, pressure_values, fields);
            }
          // Calculate physical properties for the cell
          density_models[0]->vector_value(fields, density_0);
          density_models[1]->vector_value(fields, density_1);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              switch (monitored_fluid)
                {
                  case Parameters::FluidIndicator::fluid0:
                    {
                      this->volume_monitored +=
                        fe_values_vof.JxW(q) * (1 - phase_values[q]);
                      this->mass_monitored += fe_values_vof.JxW(q) *
                                              (1 - phase_values[q]) *
                                              density_0[q];
                      break;
                    }
                  case Parameters::FluidIndicator::fluid1:
                    {
                      this->volume_monitored +=
                        fe_values_vof.JxW(q) * phase_values[q];
                      this->mass_monitored +=
                        fe_values_vof.JxW(q) * phase_values[q] * density_1[q];
                      break;
                    }
                  default:
                    throw std::runtime_error(
                      "Unsupported number of fluids (>2)");
                }
            }
        }
    }

  this->volume_monitored =
    Utilities::MPI::sum(this->volume_monitored, mpi_communicator);
  this->mass_monitored =
    Utilities::MPI::sum(this->mass_monitored, mpi_communicator);
}

template <int dim>
template <typename VectorType>
Tensor<1, dim>
VolumeOfFluid<dim>::calculate_momentum(
  const GlobalVectorType          &solution,
  const VectorType                &current_solution_fd,
  const Parameters::FluidIndicator monitored_fluid)
{
  const MPI_Comm mpi_communicator = this->triangulation->get_mpi_communicator();

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->cell_quadrature,
                              update_values | update_JxW_values);

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  FEValues<dim> fe_values_fd(*this->mapping,
                             dof_handler_fd->get_fe(),
                             *this->cell_quadrature,
                             update_values);

  const FEValuesExtractors::Vector velocity(0);
  const FEValuesExtractors::Scalar pressure(dim);

  const unsigned int          n_q_points = this->cell_quadrature->size();
  std::vector<double>         phase_values(n_q_points);
  std::vector<double>         density_0(n_q_points);
  std::vector<double>         density_1(n_q_points);
  std::vector<double>         pressure_values(n_q_points);
  std::vector<Tensor<1, dim>> velocity_values(n_q_points);

  Tensor<1, dim> total_momentum;

  // Physical properties
  const auto density_models =
    this->simulation_parameters.physical_properties_manager
      .get_density_vector();
  std::map<field, std::vector<double>> fields;
  fields.insert(
    std::pair<field, std::vector<double>>(field::pressure, n_q_points));

  for (const auto &cell_vof : this->dof_handler.active_cell_iterators())
    {
      if (cell_vof->is_locally_owned())
        {
          fe_values_vof.reinit(cell_vof);
          fe_values_vof.get_function_values(solution, phase_values);

          // Get fluid dynamics active cell iterator
          typename DoFHandler<dim>::active_cell_iterator cell_fd(
            &(*(this->triangulation)),
            cell_vof->level(),
            cell_vof->index(),
            dof_handler_fd);
          fe_values_fd.reinit(cell_fd);

          fe_values_fd[velocity].get_function_values(current_solution_fd,
                                                     velocity_values);

          if (this->simulation_parameters.physical_properties_manager
                .field_is_required(field::pressure))
            {
              fe_values_fd[pressure].get_function_values(current_solution_fd,
                                                         pressure_values);
              set_field_vector(field::pressure, pressure_values, fields);
            }



          // Calculate physical properties for the cell
          density_models[0]->vector_value(fields, density_0);
          density_models[1]->vector_value(fields, density_1);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              switch (monitored_fluid)
                {
                  case Parameters::FluidIndicator::fluid0:
                    {
                      total_momentum += fe_values_vof.JxW(q) *
                                        (1 - phase_values[q]) *
                                        velocity_values[q] * density_0[q];
                      break;
                    }
                  case Parameters::FluidIndicator::fluid1:
                    {
                      total_momentum += fe_values_vof.JxW(q) * phase_values[q] *
                                        velocity_values[q] * density_1[q];
                      break;
                    }
                  default:
                    throw std::runtime_error(
                      "Unsupported number of fluids (>2)");
                }
            }
        }
    }

  total_momentum = Utilities::MPI::sum(total_momentum, mpi_communicator);

  return total_momentum;
}

template <int dim>
void
VolumeOfFluid<dim>::finish_simulation()
{
  auto         mpi_communicator = this->triangulation->get_mpi_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (this_mpi_process == 0 &&
      simulation_parameters.analytical_solution->verbosity ==
        Parameters::Verbosity::verbose)
    {
      if (this->simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        this->error_table.omit_column_from_convergence_rate_evaluation("cells");
      else
        this->error_table.omit_column_from_convergence_rate_evaluation("time");

      this->error_table.set_scientific("error_phase", true);
      this->error_table.set_precision(
        "error_phase", this->simulation_control->get_log_precision());
      this->error_table.write_text(std::cout);
    }
}

template <int dim>
void
VolumeOfFluid<dim>::percolate_time_vectors()
{
  for (unsigned int i = this->previous_solutions.size() - 1; i > 0; --i)
    {
      this->previous_solutions[i] = this->previous_solutions[i - 1];
    }
  this->previous_solutions[0] = this->present_solution;
}

template <int dim>
void
VolumeOfFluid<dim>::postprocess(bool first_iteration)
{
  auto         mpi_communicator = this->triangulation->get_mpi_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (simulation_parameters.analytical_solution->calculate_error() &&
      !first_iteration)
    {
      double phase_error = calculate_L2_error();

      if (this->simulation_control->is_steady())
        {
          this->error_table.add_value(
            "cells", this->triangulation->n_global_active_cells());
        }
      else
        {
          this->error_table.add_value(
            "time", this->simulation_control->get_current_time());
        }
      this->error_table.add_value("error_phase", phase_error);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error phase : " << phase_error << std::endl;
        }
    }

  if (this->simulation_parameters.post_processing.calculate_mass_conservation)
    {
      const std::vector<Parameters::FluidIndicator> fluid_indicators = {
        Parameters::FluidIndicator::fluid0, Parameters::FluidIndicator::fluid1};

      const unsigned int n_fluids =
        this->simulation_parameters.physical_properties_manager
          .get_number_of_fluids();

      // Set column names according to dim for volume and mass values
      std::string volume_column_name;
      std::string mass_column_name;

      std::string geometric_volume_column_name;
      std::string area_column_name;

      // To display when verbose
      std::vector<std::string> dependent_column_names;
      dependent_column_names.reserve(n_fluids * 2 + 1);
      std::vector<double> volumes_masses_momentum_and_sharpening_threshold;
      volumes_masses_momentum_and_sharpening_threshold.reserve(n_fluids * 2 +
                                                               1);

      if (this_mpi_process == 0)
        {
          // Set conservation monitoring table
          if (this->simulation_control->is_steady())
            {
              this->table_monitoring_vof.add_value(
                "cells", this->triangulation->n_global_active_cells());
            }
          else
            {
              this->table_monitoring_vof.add_value(
                "time", this->simulation_control->get_current_time());
              this->table_monitoring_vof.set_scientific("time", true);
            }
        }

      double geometric_volume_outside, geometric_volume_inside, surface;

      // The mesh classifier in the InterfaceTools::compute_surface_and_volume
      // is not compatible with other FE than FE_Q
      if (!simulation_parameters.fem_parameters.VOF_uses_dg)
        {
          std::tie(geometric_volume_outside, surface) =
            InterfaceTools::compute_surface_and_volume(
              dof_handler, *fe, this->present_solution, 0.5, mpi_communicator);

          // Compute geometric inside volume (fluid 1)
          const double global_volume =
            GridTools::volume(*this->triangulation, *this->mapping);
          geometric_volume_inside = global_volume - geometric_volume_outside;
        }

      for (unsigned int i = 0; i < n_fluids; i++)
        {
          // Calculate volume and mass
          calculate_volume_and_mass(this->present_solution,
                                    *multiphysics->get_solution(
                                      PhysicsID::fluid_dynamics),
                                    fluid_indicators[i]);

          // Calculate momentum
          const Tensor<1, dim> momentum =
            calculate_momentum(this->present_solution,
                               *multiphysics->get_solution(
                                 PhysicsID::fluid_dynamics),
                               fluid_indicators[i]);

          if (this_mpi_process == 0)
            {
              std::string fluid_id("");

              if (fluid_indicators[i] == Parameters::FluidIndicator::fluid1)
                {
                  fluid_id = "fluid_1";
                }
              else if (fluid_indicators[i] ==
                       Parameters::FluidIndicator::fluid0)
                {
                  fluid_id = "fluid_0";
                }

              std::vector<std::string> momentum_names(
                {"momentum-x_" + fluid_id, "momentum-y_" + fluid_id});
              if constexpr (dim == 3)
                momentum_names.emplace_back("momentum-z_" + fluid_id);

              if constexpr (dim == 2)
                {
                  volume_column_name = "surface_" + fluid_id;
                  mass_column_name   = "mass_per_length_" + fluid_id;
                  geometric_volume_column_name =
                    "geometric_surface_" + fluid_id;
                  area_column_name = "length_" + fluid_id;
                }
              else if constexpr (dim == 3)
                {
                  volume_column_name           = "volume_" + fluid_id;
                  mass_column_name             = "mass_" + fluid_id;
                  geometric_volume_column_name = "geometric_volume_" + fluid_id;
                  area_column_name             = "surface_" + fluid_id;
                }

              // Add "surface" or "volume" column
              this->table_monitoring_vof.add_value(volume_column_name,
                                                   this->volume_monitored);
              this->table_monitoring_vof.set_scientific(volume_column_name,
                                                        true);

              if (!simulation_parameters.fem_parameters.VOF_uses_dg)
                {
                  // Add "geometric surface" or "geometric volume" column
                  if (fluid_id == "fluid_1")
                    {
                      this->table_monitoring_vof.add_value(
                        geometric_volume_column_name, geometric_volume_inside);
                      this->table_monitoring_vof.set_scientific(
                        geometric_volume_column_name, true);
                    }
                  else
                    {
                      this->table_monitoring_vof.add_value(
                        geometric_volume_column_name, geometric_volume_outside);
                      this->table_monitoring_vof.set_scientific(
                        geometric_volume_column_name, true);
                    }
                }

              // Add "mass per length" or "mass" column
              this->table_monitoring_vof.add_value(mass_column_name,
                                                   this->mass_monitored);
              this->table_monitoring_vof.set_scientific(mass_column_name, true);

              if (!simulation_parameters.fem_parameters.VOF_uses_dg)
                {
                  this->table_monitoring_vof.add_value(area_column_name,
                                                       surface);
                  this->table_monitoring_vof.set_scientific(area_column_name,
                                                            true);
                }
              // Add "momentum" columns
              for (int d = 0; d < dim; ++d)
                {
                  this->table_monitoring_vof.add_value(momentum_names[d],
                                                       momentum[d]);
                  this->table_monitoring_vof.set_scientific(momentum_names[d],
                                                            true);
                }

              if (this->simulation_parameters.post_processing.verbosity ==
                  Parameters::Verbosity::verbose)
                {
                  dependent_column_names.emplace_back(volume_column_name);

                  if (!simulation_parameters.fem_parameters.VOF_uses_dg)
                    {
                      dependent_column_names.emplace_back(
                        geometric_volume_column_name);
                    }

                  dependent_column_names.emplace_back(mass_column_name);
                  if (!simulation_parameters.fem_parameters.VOF_uses_dg)
                    {
                      dependent_column_names.emplace_back(area_column_name);
                    }
                  volumes_masses_momentum_and_sharpening_threshold.emplace_back(
                    this->volume_monitored);
                  if (!simulation_parameters.fem_parameters.VOF_uses_dg)
                    {
                      if (fluid_id == "fluid_1")
                        volumes_masses_momentum_and_sharpening_threshold
                          .emplace_back(geometric_volume_inside);
                      else
                        volumes_masses_momentum_and_sharpening_threshold
                          .emplace_back(geometric_volume_outside);
                    }
                  volumes_masses_momentum_and_sharpening_threshold.emplace_back(
                    this->mass_monitored);

                  if (!simulation_parameters.fem_parameters.VOF_uses_dg)
                    {
                      volumes_masses_momentum_and_sharpening_threshold
                        .emplace_back(surface);
                    }
                  for (int d = 0; d < dim; ++d)
                    {
                      dependent_column_names.emplace_back(momentum_names[d]);
                      volumes_masses_momentum_and_sharpening_threshold
                        .emplace_back(momentum[d]);
                    }
                }
            }
        }

      if (this_mpi_process == 0)
        {
          // Add sharpening threshold column
          this->table_monitoring_vof.add_value("sharpening_threshold",
                                               this->sharpening_threshold);

          if (this->simulation_control->get_step_number() %
                this->simulation_parameters.post_processing.output_frequency ==
              0)
            {
              // Save table to .dat
              std::string filename =
                this->simulation_parameters.simulation_control.output_folder +
                this->simulation_parameters.post_processing
                  .mass_conservation_output_name +
                ".dat";
              std::ofstream output(filename.c_str());
              this->table_monitoring_vof.write_text(output);
            }

          // Print on terminal
          if (this->simulation_parameters.post_processing.verbosity ==
              Parameters::Verbosity::verbose)
            {
              volumes_masses_momentum_and_sharpening_threshold.emplace_back(
                this->sharpening_threshold);
              dependent_column_names.emplace_back("sharpening_threshold");

              // Dependent variable columns (volumes, masses and sharpening
              // threshold)
              std::vector<std::vector<double>>
                volumes_masses_and_sharpening_thresholds;
              volumes_masses_and_sharpening_thresholds.emplace_back(
                volumes_masses_momentum_and_sharpening_threshold);

              // Time column
              std::string         independent_column_name = "time";
              std::vector<double> time                    = {
                this->simulation_control->get_current_time()};

              TableHandler table = make_table_scalars_vectors(
                time,
                independent_column_name,
                volumes_masses_and_sharpening_thresholds,
                dependent_column_names,
                this->simulation_parameters.simulation_control.log_precision,
                true);

              std::cout << std::endl;



              announce_string(this->pcout, "VOF Mass Conservation");
              table.write_text(std::cout);
            }
        }
    }

  if (this->simulation_parameters.post_processing.calculate_barycenter)
    {
      // Calculate volume and mass (this->mass_monitored)
      std::pair<Tensor<1, dim>, Tensor<1, dim>> position_and_velocity;

      if (multiphysics->fluid_dynamics_is_block())
        {
          // Check if the post processed variable needs to be calculated with
          // the average velocity profile or the fluid solution.
          if (this->simulation_parameters.initial_condition->type ==
                Parameters::FluidDynamicsInitialConditionType::
                  average_velocity_profile &&
              !this->simulation_parameters.multiphysics.fluid_dynamics &&
              simulation_control->get_current_time() >
                this->simulation_parameters.post_processing
                  .initial_time_for_average_velocities)
            {
              position_and_velocity = calculate_barycenter(
                this->present_solution,
                *multiphysics->get_block_time_average_solution(
                  PhysicsID::fluid_dynamics));
            }
          else
            {
              position_and_velocity =
                calculate_barycenter(this->present_solution,
                                     *multiphysics->get_block_solution(
                                       PhysicsID::fluid_dynamics));
            }
        }
      else
        {
          // Check if the post processed variable needs to be calculated with
          // the average velocity profile or the fluid solution.
          if (this->simulation_parameters.initial_condition->type ==
                Parameters::FluidDynamicsInitialConditionType::
                  average_velocity_profile &&
              !this->simulation_parameters.multiphysics.fluid_dynamics &&
              simulation_control->get_current_time() >
                this->simulation_parameters.post_processing
                  .initial_time_for_average_velocities)
            {
              position_and_velocity =
                calculate_barycenter(this->present_solution,
                                     *multiphysics->get_time_average_solution(
                                       PhysicsID::fluid_dynamics));
            }
          else
            {
              position_and_velocity =
                calculate_barycenter(this->present_solution,
                                     *multiphysics->get_solution(
                                       PhysicsID::fluid_dynamics));
            }
        }
      if (this_mpi_process == 0)
        {
          if (simulation_parameters.post_processing.verbosity ==
              Parameters::Verbosity::verbose)
            {
              std::cout << std::endl;
              std::string independent_column_names = "time";

              std::vector<std::string> dependent_column_names;
              dependent_column_names.emplace_back("x_vof");
              dependent_column_names.emplace_back("y_vof");
              if constexpr (dim == 3)
                dependent_column_names.emplace_back("z_vof");
              dependent_column_names.emplace_back("vx_vof");
              dependent_column_names.emplace_back("vy_vof");
              if constexpr (dim == 3)
                dependent_column_names.emplace_back("vz_vof");

              std::vector<Tensor<1, dim>> position_vector{
                position_and_velocity.first};
              std::vector<Tensor<1, dim>> velocity_vector{
                position_and_velocity.second};

              std::vector<std::vector<Tensor<1, dim>>>
                position_and_velocity_vectors{position_vector, velocity_vector};

              std::vector<double> time = {
                this->simulation_control->get_current_time()};

              TableHandler table = make_table_scalars_tensors(
                time,
                independent_column_names,
                position_and_velocity_vectors,
                dependent_column_names,
                this->simulation_parameters.simulation_control.log_precision,
                true);

              announce_string(this->pcout, "VOF Barycenter");
              table.write_text(std::cout);
            }

          this->table_barycenter.add_value(
            "time", simulation_control->get_current_time());
          this->table_barycenter.set_scientific("time", true);

          this->table_barycenter.add_value("x_vof",
                                           position_and_velocity.first[0]);
          this->table_barycenter.set_scientific("x_vof", true);

          this->table_barycenter.add_value("y_vof",
                                           position_and_velocity.first[1]);
          this->table_barycenter.set_scientific("y_vof", true);

          if constexpr (dim == 3)
            {
              this->table_barycenter.add_value("z_vof",
                                               position_and_velocity.first[2]);
              this->table_barycenter.set_scientific("z_vof", true);
            }

          this->table_barycenter.add_value("vx_vof",
                                           position_and_velocity.second[0]);
          this->table_barycenter.set_scientific("vx_vof", true);

          this->table_barycenter.add_value("vy_vof",
                                           position_and_velocity.second[1]);
          this->table_barycenter.set_scientific("vy_vof", true);

          if constexpr (dim == 3)
            {
              this->table_barycenter.add_value("vz_vof",
                                               position_and_velocity.second[2]);
              this->table_barycenter.set_scientific("vz_vof", true);
            }


          if (this->simulation_control->get_step_number() %
                this->simulation_parameters.post_processing.output_frequency ==
              0)
            {
              // Save table to .dat
              std::string filename =
                this->simulation_parameters.simulation_control.output_folder +
                this->simulation_parameters.post_processing
                  .barycenter_output_name +
                ".dat";
              std::ofstream output(filename.c_str());
              this->table_barycenter.write_text(output);
              output.close();
            }
        }
    }

  if (this->simulation_parameters.timer.type ==
      Parameters::Timer::Type::iteration)
    {
      announce_string(this->pcout, "VOF");
      this->computing_timer.print_summary();
      this->computing_timer.reset();
    }
}


template <int dim>
void
VolumeOfFluid<dim>::modify_solution()
{
  {
    TimerOutput::Scope t(this->computing_timer, "Modify solution");

    auto vof_parameters =
      this->simulation_parameters.multiphysics.vof_parameters;
    // Interface sharpening
    if (vof_parameters.regularization_method.sharpening.enable)
      {
        // Interface sharpening is done at a constant frequency
        if (this->simulation_control->get_step_number() %
              this->simulation_parameters.multiphysics.vof_parameters
                .regularization_method.frequency ==
            0)
          {
            handle_interface_sharpening();
          }
      }
  }

  // Apply algebraic interface reinitialization
  if (simulation_parameters.multiphysics.vof_parameters.regularization_method
        .algebraic_interface_reinitialization.enable &&
      (simulation_control->get_step_number() %
         simulation_parameters.multiphysics.vof_parameters.regularization_method
           .frequency ==
       0))
    reinitialize_interface_with_algebraic_method();

  // Apply geometric interface reinitialization
  if (simulation_parameters.multiphysics.vof_parameters.regularization_method
        .geometric_interface_reinitialization.enable &&
      (simulation_control->get_step_number() %
         simulation_parameters.multiphysics.vof_parameters.regularization_method
           .frequency ==
       0))
    reinitialize_interface_with_geometric_method();

  // Apply filter to phase fraction values
  apply_phase_filter(this->present_solution, this->filtered_solution);

  // Solve phase fraction gradient and curvature projections
  if (simulation_parameters.multiphysics.vof_parameters.surface_tension_force
        .enable)
    {
      this->vof_subequations_interface
        ->set_vof_filtered_solution_and_dof_handler(this->filtered_solution,
                                                    this->dof_handler);
      this->vof_subequations_interface->solve_specific_subequation(
        VOFSubequationsID::phase_gradient_projection);
      this->vof_subequations_interface->solve_specific_subequation(
        VOFSubequationsID::curvature_projection);
    }
}

template <int dim>
void
VolumeOfFluid<dim>::handle_interface_sharpening()
{
  if (this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "Sharpening interface at step "
                  << this->simulation_control->get_step_number() << std::endl;
    }
  if (this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.sharpening.type ==
      Parameters::SharpeningType::adaptive)
    {
      if (this->simulation_parameters.multiphysics.vof_parameters
            .regularization_method.verbosity != Parameters::Verbosity::quiet)
        {
          this->pcout << "   Adapting the sharpening threshold" << std::endl;
        }

      this->sharpening_threshold = find_sharpening_threshold();

      if (this->simulation_parameters.multiphysics.vof_parameters
            .regularization_method.verbosity != Parameters::Verbosity::quiet)
        {
          this->pcout << "   ... final sharpening is : "
                      << this->sharpening_threshold << std::endl;
        }
    }
  else
    {
      // Constant sharpening
      this->sharpening_threshold =
        this->simulation_parameters.multiphysics.vof_parameters
          .regularization_method.sharpening.threshold;
    }

  // Sharpen the interface of all solutions (present and previous)
  sharpen_interface(this->present_solution, this->sharpening_threshold, true);
}


template <int dim>
double
VolumeOfFluid<dim>::find_sharpening_threshold()
{
  // Sharpening threshold (st) search range extrema
  double st_min =
    0.5 - this->simulation_parameters.multiphysics.vof_parameters
            .regularization_method.sharpening.threshold_max_deviation;
  double st_max =
    0.5 + this->simulation_parameters.multiphysics.vof_parameters
            .regularization_method.sharpening.threshold_max_deviation;

  // Useful definitions for readability
  const double mass_deviation_tol =
    this->simulation_parameters.multiphysics.vof_parameters
      .regularization_method.sharpening.tolerance *
    this->mass_first_iteration;
  const unsigned int max_iterations =
    this->simulation_parameters.multiphysics.vof_parameters
      .regularization_method.sharpening.max_iterations;

  const Parameters::FluidIndicator monitored_fluid =
    this->simulation_parameters.multiphysics.vof_parameters
      .regularization_method.sharpening.monitored_fluid;

  unsigned int nb_search_ite = 0;
  // Local variable for the tested sharpening_threshold values

  double mass_deviation_min = calculate_mass_deviation(monitored_fluid, st_min);
  double mass_deviation_max = calculate_mass_deviation(monitored_fluid, st_max);
  double mass_deviation_avg = 0.;
  double st_avg             = 0.;


  // Bissection algorithm to calculate an interface sharpening value that would
  // ensure mass conservation of the monitored phase (do-while loop, see
  // condition below)
  do
    {
      nb_search_ite++;
      // Calculate middle point value
      st_avg = (st_min + st_max) / 2.;

      mass_deviation_avg = calculate_mass_deviation(monitored_fluid, st_avg);

      if (this->simulation_parameters.multiphysics.vof_parameters
            .regularization_method.verbosity != Parameters::Verbosity::quiet)
        {
          this->pcout
            << "   ... step " << nb_search_ite
            << " of the search algorithm, min, avg, max mass deviation is : "
            << mass_deviation_min << " " << mass_deviation_avg << " "
            << mass_deviation_max << std::endl;
        }

      if (mass_deviation_avg * mass_deviation_min < 0)
        {
          st_max             = st_avg;
          mass_deviation_max = mass_deviation_avg;
        }
      else
        {
          st_min             = st_avg;
          mass_deviation_min = mass_deviation_avg;
        }
    }
  while (std::abs(mass_deviation_avg) > mass_deviation_tol &&
         nb_search_ite < max_iterations);

  // Take minimum deviation in between the two endpoints of the last
  // interval searched, if out of the do-while loop because max_iterations is
  // reached
  if (std::abs(mass_deviation_avg) > mass_deviation_tol)
    {
      double mass_deviation_endpoint = 0.;
      double st_tested               = 0;
      if (st_min == st_avg)
        st_tested = st_max;
      else if (st_max == st_avg)
        st_tested = st_min;

      mass_deviation_endpoint =
        calculate_mass_deviation(monitored_fluid, st_tested);

      // Retake st_ave value if mass deviation is not lowered at endpoint
      // values
      if (std::abs(mass_deviation_endpoint) > std::abs(mass_deviation_avg))
        {
          st_tested = st_avg;
        }

      // Output message
      if (std::abs(mass_deviation_endpoint) > mass_deviation_tol)
        {
          this->pcout
            << "  WARNING: Maximum number of iterations (" << nb_search_ite
            << ") reached in the " << std::endl
            << "  adaptive sharpening threshold algorithm, remaining error"
            << std::endl
            << "  on mass conservation is: "
            << (this->mass_monitored - this->mass_first_iteration) /
                 this->mass_first_iteration
            << std::endl
            << "  Consider increasing the sharpening threshold range or the "
            << std::endl
            << "  number of iterations to reach the mass conservation tolerance."
            << std::endl;
        }
    }

  // Output message that mass conservation condition is reached
  if (this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "   ... search algorithm took : " << nb_search_ite
                  << " step(s) " << std::endl
                  << "   ... error on mass conservation reached: "
                  << (this->mass_monitored - this->mass_first_iteration) /
                       this->mass_first_iteration
                  << std::endl;
    }

  return st_avg;
}


template <int dim>
double
VolumeOfFluid<dim>::calculate_mass_deviation(
  const Parameters::FluidIndicator monitored_fluid,
  const double                     sharpening_threshold)
{
  // Copy present solution VOF
  evaluation_point = this->present_solution;

  // Sharpen interface using the tested threshold value
  sharpen_interface(evaluation_point, sharpening_threshold, false);

  // Calculate mass of the monitored phase
  calculate_volume_and_mass(evaluation_point,
                            *multiphysics->get_solution(
                              PhysicsID::fluid_dynamics),
                            monitored_fluid);

  // Calculate relative mass deviation
  double mass_deviation = (this->mass_monitored - this->mass_first_iteration) /
                          this->mass_first_iteration;

  return mass_deviation;
}

template <int dim>
void
VolumeOfFluid<dim>::sharpen_interface(GlobalVectorType &solution,
                                      const double      sharpening_threshold,
                                      const bool sharpen_previous_solutions)
{
  // Limit the phase fractions between 0 and 1
  update_solution_and_constraints(solution);
  if (sharpen_previous_solutions)
    {
      for (auto &previous_solution : previous_solutions)
        update_solution_and_constraints(previous_solution);
    }

  // Assemble matrix and solve the system for interface sharpening
  assemble_L2_projection_interface_sharpening(solution, sharpening_threshold);
  solve_interface_sharpening(solution);

  if (sharpen_previous_solutions)
    {
      for (auto &previous_solution : previous_solutions)
        {
          assemble_L2_projection_interface_sharpening(previous_solution,
                                                      sharpening_threshold);
          solve_interface_sharpening(previous_solution);
        }
    }

  // Re limit the phase fractions between 0 and 1 after interface
  // sharpening
  update_solution_and_constraints(solution);
  if (sharpen_previous_solutions)
    {
      for (auto &previous_solution : previous_solutions)
        update_solution_and_constraints(previous_solution);
    }
}


template <int dim>
void
VolumeOfFluid<dim>::smooth_phase_fraction(GlobalVectorType &solution)
{
  assemble_projection_phase_fraction(solution);
  solve_projection_phase_fraction(solution);
}


template <int dim>
void
VolumeOfFluid<dim>::assemble_projection_phase_fraction(
  GlobalVectorType &solution)
{
  // Get fe values of VOF phase fraction
  FEValues<dim> fe_values_phase_fraction(*this->mapping,
                                         *this->fe,
                                         *this->cell_quadrature,
                                         update_values | update_JxW_values |
                                           update_gradients);

  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

  const unsigned int n_q_points = this->cell_quadrature->size();
  FullMatrix<double> local_matrix_phase_fraction(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_phase_fraction(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<double>         phi_phase_fraction(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_phase_fraction_gradient(dofs_per_cell);

  std::vector<double> phase_values(n_q_points);

  // Reinitialize system matrix and rhs for the projected phase fraction
  system_rhs_phase_fraction    = 0;
  system_matrix_phase_fraction = 0;

  double h;
  double cell_measure;

  const double phase_fraction_diffusion_factor =
    this->simulation_parameters.initial_condition
      ->projection_step_diffusion_factor;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_phase_fraction.reinit(cell);

          // Compute cell diameter
          cell_measure = compute_cell_measure_with_JxW(
            fe_values_phase_fraction.get_JxW_values());
          h = compute_cell_diameter<dim>(cell_measure, fe->degree);

          local_matrix_phase_fraction = 0;
          local_rhs_phase_fraction    = 0;

          // Get phase fraction values
          fe_values_phase_fraction.get_function_values(solution, phase_values);

          double color_function = 0.0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              color_function +=
                phase_values[q] * fe_values_phase_fraction.JxW(q);
            }

          color_function /= cell_measure;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_phase_fraction[k] =
                    fe_values_phase_fraction.shape_value(k, q);
                  phi_phase_fraction_gradient[k] =
                    fe_values_phase_fraction.shape_grad(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_phase_fraction(i, j) +=
                        (phi_phase_fraction[j] * phi_phase_fraction[i] +
                         phase_fraction_diffusion_factor * h * h *
                           scalar_product(phi_phase_fraction_gradient[i],
                                          phi_phase_fraction_gradient[j])) *
                        fe_values_phase_fraction.JxW(q);
                    }

                  // rhs
                  local_rhs_phase_fraction(i) +=
                    phi_phase_fraction[i] * color_function *
                    fe_values_phase_fraction.JxW(q);
                }
            }

          cell->get_dof_indices(local_dof_indices);
          this->nonzero_constraints.distribute_local_to_global(
            local_matrix_phase_fraction,
            local_rhs_phase_fraction,
            local_dof_indices,
            system_matrix_phase_fraction,
            system_rhs_phase_fraction);
        }
    }
  system_matrix_phase_fraction.compress(VectorOperation::add);
  system_rhs_phase_fraction.compress(VectorOperation::add);
}


template <int dim>
void
VolumeOfFluid<dim>::solve_projection_phase_fraction(GlobalVectorType &solution)
{
  // Solve the L2 projection system
  const double linear_solver_tolerance = 1e-13;

  GlobalVectorType completely_distributed_phase_fraction_solution(
    this->locally_owned_dofs, triangulation->get_mpi_communicator());

  SolverControl solver_control(
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF).max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverCG solver(solver_control);

  const unsigned int ilu_fill =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .ilu_precond_rtol;

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix_phase_fraction,
                                 preconditionerOptions);
  solver.solve(system_matrix_phase_fraction,
               completely_distributed_phase_fraction_solution,
               system_rhs_phase_fraction,
               *ilu_preconditioner);

  if (this->simulation_parameters.multiphysics.vof_parameters
        .surface_tension_force.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver (phase fraction) took : "
                  << solver_control.last_step() << " steps " << std::endl;
    }

  this->nonzero_constraints.distribute(
    completely_distributed_phase_fraction_solution);
  solution = completely_distributed_phase_fraction_solution;
}


template <int dim>
void
VolumeOfFluid<dim>::pre_mesh_adaptation()
{
  this->solution_transfer->prepare_for_coarsening_and_refinement(
    this->present_solution);

  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      this->previous_solutions_transfer[i]
        .prepare_for_coarsening_and_refinement(this->previous_solutions[i]);
    }
}


template <int dim>
void
VolumeOfFluid<dim>::post_mesh_adaptation()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  // Set up the vectors for the transfer
  GlobalVectorType tmp(this->locally_owned_dofs, mpi_communicator);

  // Interpolate the solution at time and previous time
  this->solution_transfer->interpolate(tmp);

  // Distribute constraints
  this->nonzero_constraints.distribute(tmp);

  // Fix on the new mesh
  this->present_solution = tmp;

  // Transfer previous solutions
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      GlobalVectorType tmp_previous_solution(this->locally_owned_dofs,
                                             mpi_communicator);
      this->previous_solutions_transfer[i].interpolate(tmp_previous_solution);
      this->nonzero_constraints.distribute(tmp_previous_solution);
      this->previous_solutions[i] = tmp_previous_solution;
    }

  // Apply filter to phase fraction
  apply_phase_filter(this->present_solution, this->filtered_solution);
}

template <int dim>
void
VolumeOfFluid<dim>::compute_kelly(
  const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                        &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  if (ivar.first == Variable::phase)
    {
      const FEValuesExtractors::Scalar phase(0);

      KellyErrorEstimator<dim>::estimate(
        *this->mapping,
        this->dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        this->present_solution,
        estimated_error_per_cell,
        this->fe->component_mask(phase));
    }
}

template <int dim>
std::vector<OutputStructTableHandler>
VolumeOfFluid<dim>::gather_tables()
{
  std::vector<OutputStructTableHandler> table_output_structs;

  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder;
  std::string suffix = ".checkpoint";

  if (this->simulation_parameters.analytical_solution->calculate_error())
    table_output_structs.emplace_back(
      this->error_table,
      prefix + this->simulation_parameters.analytical_solution->get_filename() +
        "_VOF" + suffix);

  if (this->simulation_parameters.post_processing.calculate_mass_conservation)
    table_output_structs.emplace_back(
      this->table_monitoring_vof,
      prefix + 
      this->simulation_parameters.post_processing.mass_conservation_output_name +
        suffix);

  if (this->simulation_parameters.post_processing.calculate_barycenter)
    table_output_structs.emplace_back(
      this->table_barycenter,
       prefix +
        this->simulation_parameters.post_processing.barycenter_output_name +
        suffix);     

  return table_output_structs;

}

template <int dim>
void
VolumeOfFluid<dim>::write_checkpoint()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  std::vector<const GlobalVectorType *> sol_set_transfer;

  solution_transfer =
    std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(dof_handler);

  sol_set_transfer.emplace_back(&this->present_solution);
  for (const auto &previous_solution : this->previous_solutions)
    {
      sol_set_transfer.emplace_back(&previous_solution);
    }
  this->solution_transfer->prepare_for_serialization(sol_set_transfer);

  // Serialize all post-processing tables that are currently used with the VOF
  // solver
  const std::vector<OutputStructTableHandler> &table_output_structs =
    this->gather_tables();
  this->serialize_tables_vector(table_output_structs, mpi_communicator);
}

template <int dim>
void
VolumeOfFluid<dim>::read_checkpoint()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  auto previous_solutions_size = this->previous_solutions.size();
  this->pcout << "Reading VOF checkpoint" << std::endl;

  std::vector<GlobalVectorType *> input_vectors(1 + previous_solutions_size);
  GlobalVectorType                distributed_system(this->locally_owned_dofs,
                                      mpi_communicator);
  input_vectors[0] = &distributed_system;


  std::vector<GlobalVectorType> distributed_previous_solutions;
  distributed_previous_solutions.reserve(previous_solutions_size);
  for (unsigned int i = 0; i < previous_solutions_size; ++i)
    {
      distributed_previous_solutions.emplace_back(
        GlobalVectorType(this->locally_owned_dofs, mpi_communicator));
      input_vectors[i + 1] = &distributed_previous_solutions[i];
    }

  this->solution_transfer->deserialize(input_vectors);

  this->present_solution = distributed_system;
  for (unsigned int i = 0; i < previous_solutions_size; ++i)
    {
      this->previous_solutions[i] = distributed_previous_solutions[i];
    }

  // Apply filter to phase fraction
  apply_phase_filter(this->present_solution, this->filtered_solution);

  // Deserialize all post-processing tables that are currently used with the VOF
  // solver
  std::vector<OutputStructTableHandler> table_output_structs =
    this->gather_tables();
  this->deserialize_tables_vector(table_output_structs, mpi_communicator);
  
  if (this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.sharpening.type ==
      Parameters::SharpeningType::adaptive)
    {
      // Calculate volume and mass
      calculate_volume_and_mass(
        this->present_solution,
        *multiphysics->get_solution(PhysicsID::fluid_dynamics),
        this->simulation_parameters.multiphysics.vof_parameters
          .regularization_method.sharpening.monitored_fluid);

      this->mass_first_iteration = this->mass_monitored;
    }
}


template <int dim>
void
VolumeOfFluid<dim>::setup_dofs()
{
  verify_consistency_of_boundary_conditions();

  auto mpi_communicator = triangulation->get_mpi_communicator();

  // Setup DoFs for all active subequations
  this->vof_subequations_interface->setup_dofs();

  this->dof_handler.distribute_dofs(*this->fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();

  this->locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(this->dof_handler);

  this->present_solution.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                mpi_communicator);

  this->filtered_solution.reinit(this->locally_owned_dofs,
                                 this->locally_relevant_dofs,
                                 mpi_communicator);

  this->level_set.reinit(this->locally_owned_dofs,
                         this->locally_relevant_dofs,
                         mpi_communicator);

  // Previous solutions for transient schemes
  for (auto &solution : this->previous_solutions)
    {
      solution.reinit(this->locally_owned_dofs,
                      this->locally_relevant_dofs,
                      mpi_communicator);
    }

  this->system_rhs.reinit(this->locally_owned_dofs, mpi_communicator);

  this->newton_update.reinit(this->locally_owned_dofs, mpi_communicator);

  this->local_evaluation_point.reinit(this->locally_owned_dofs,
                                      mpi_communicator);

  define_non_zero_constraints();

  // Boundary conditions for Newton correction
  define_zero_constraints();

  // Sparse matrices initialization
  DynamicSparsityPattern dsp(this->locally_relevant_dofs);
  if (simulation_parameters.fem_parameters.VOF_uses_dg)
    {
      DoFTools::make_flux_sparsity_pattern(this->dof_handler,
                                           dsp,
                                           nonzero_constraints,
                                           /*keep_constrained_dofs = */ true);
    }
  else
    {
      DoFTools::make_sparsity_pattern(this->dof_handler,
                                      dsp,
                                      nonzero_constraints,
                                      /*keep_constrained_dofs = */ false);
    }

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             this->locally_owned_dofs,
                                             mpi_communicator,
                                             this->locally_relevant_dofs);

  // Initialization of phase fraction matrices for interface sharpening.
  // system_matrix_phase_fraction is used in
  // assemble_L2_projection_interface_sharpening for assembling the system for
  // sharpening the interface, while complete_system_matrix_phase_fraction is
  // used in update_solution_and_constraints to limit the phase fraction
  // values between 0 and 1. According to step-41, to limit the phase fractions
  // we compute the Lagrange multiplier as the residual of the original linear
  // system, given via the variables complete_system_matrix_phase_fraction and
  // complete_system_rhs_phase_fraction
  system_matrix_phase_fraction.reinit(this->locally_owned_dofs,
                                      this->locally_owned_dofs,
                                      dsp,
                                      mpi_communicator);

  complete_system_matrix_phase_fraction.reinit(this->locally_owned_dofs,
                                               this->locally_owned_dofs,
                                               dsp,
                                               mpi_communicator);

  complete_system_rhs_phase_fraction.reinit(this->locally_owned_dofs,
                                            mpi_communicator);

  // In update_solution_and_constraints (which limits the phase fraction
  // between 0 and 1) nodal_phase_fraction_owned copies the solution, then
  // limits it, and finally updates (rewrites) the solution.
  nodal_phase_fraction_owned.reinit(this->locally_owned_dofs, mpi_communicator);

  // Right hand side of the interface sharpening problem (used in
  // assemble_L2_projection_interface_sharpening).
  system_rhs_phase_fraction.reinit(this->locally_owned_dofs, mpi_communicator);

  this->system_matrix.reinit(this->locally_owned_dofs,
                             this->locally_owned_dofs,
                             dsp,
                             mpi_communicator);

  this->pcout << "   Number of VOF degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;

  // Provide the VOF dof_handler and solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::VOF, &this->dof_handler);
  multiphysics->set_solution(PhysicsID::VOF, &this->present_solution);
  multiphysics->set_filtered_solution(PhysicsID::VOF, &this->filtered_solution);
  multiphysics->set_previous_solutions(PhysicsID::VOF,
                                       &this->previous_solutions);

  if (simulation_parameters.multiphysics.vof_parameters.regularization_method
        .geometric_interface_reinitialization.enable)
    {
      signed_distance_solver->setup_dofs();
    }


  mass_matrix_phase_fraction.reinit(this->locally_owned_dofs,
                                    this->locally_owned_dofs,
                                    dsp,
                                    mpi_communicator);

  assemble_mass_matrix(mass_matrix_phase_fraction);
}

template <int dim>
void
VolumeOfFluid<dim>::update_boundary_conditions()
{
  if (!this->simulation_parameters.boundary_conditions_vof.time_dependent)
    return;

  double time = this->simulation_control->get_current_time();
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_vof.type)
    {
      if (type == BoundaryConditions::BoundaryType::vof_dirichlet)
        {
          this->simulation_parameters.boundary_conditions_vof.phase_fraction
            .at(id)
            ->set_time(time);
        }
    }
  define_non_zero_constraints();
}


template <int dim>
void
VolumeOfFluid<dim>::define_zero_constraints()
{
  // Zero constraints
  this->zero_constraints.clear();
  this->zero_constraints.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          this->zero_constraints);

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_vof.type)
    {
      if (type == BoundaryConditions::BoundaryType::vof_dirichlet)
        {
          VectorTools::interpolate_boundary_values(
            this->dof_handler,
            id,
            Functions::ZeroFunction<dim>(),
            this->zero_constraints);
        }
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
            id,
            this->simulation_parameters.boundary_conditions_vof
              .periodic_neighbor_id.at(id),
            this->simulation_parameters.boundary_conditions_vof
              .periodic_direction.at(id),
            this->zero_constraints);
        }
    }
  this->zero_constraints.close();
}


template <int dim>
void
VolumeOfFluid<dim>::define_non_zero_constraints()
{
  {
    nonzero_constraints.clear();
    nonzero_constraints.reinit(this->locally_owned_dofs,
                               this->locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);

    for (auto const &[id, type] :
         this->simulation_parameters.boundary_conditions_vof.type)
      {
        if (type == BoundaryConditions::BoundaryType::vof_dirichlet)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              id,
              *this->simulation_parameters.boundary_conditions_vof
                 .phase_fraction.at(id),
              nonzero_constraints);
          }
        if (type == BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              id,
              this->simulation_parameters.boundary_conditions_vof
                .periodic_neighbor_id.at(id),
              this->simulation_parameters.boundary_conditions_vof
                .periodic_direction.at(id),
              nonzero_constraints);
          }
      }
  }
  nonzero_constraints.close();
}


template <int dim>
void
VolumeOfFluid<dim>::set_initial_conditions()
{
  VectorTools::interpolate(*this->mapping,
                           this->dof_handler,
                           simulation_parameters.initial_condition->VOF,
                           this->newton_update);
  this->nonzero_constraints.distribute(this->newton_update);
  this->present_solution = this->newton_update;


  if (simulation_parameters.initial_condition
        ->vof_initial_condition_smoothing ==
      Parameters::VOFInitialConditionType::diffusive)
    smooth_phase_fraction(this->present_solution);

  else if (simulation_parameters.initial_condition
             ->vof_initial_condition_smoothing ==
           Parameters::VOFInitialConditionType::geometric)
    reinitialize_interface_with_geometric_method();

  apply_phase_filter(this->present_solution, this->filtered_solution);

  if (this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.sharpening.type ==
      Parameters::SharpeningType::adaptive)
    {
      // Calculate volume and mass
      calculate_volume_and_mass(
        this->present_solution,
        *multiphysics->get_solution(PhysicsID::fluid_dynamics),
        this->simulation_parameters.multiphysics.vof_parameters
          .regularization_method.sharpening.monitored_fluid);

      this->mass_first_iteration = this->mass_monitored;
    }

  // Solve initial phase fraction gradient and curvature projections if solution
  // outputs are requested
  if ((simulation_parameters.multiphysics.vof_parameters.surface_tension_force
         .enable &&
       simulation_parameters.multiphysics.vof_parameters.surface_tension_force
         .output_vof_auxiliary_fields) ||
      simulation_parameters.multiphysics.vof_parameters.regularization_method
        .algebraic_interface_reinitialization.enable)
    {
      this->vof_subequations_interface
        ->set_vof_filtered_solution_and_dof_handler(this->filtered_solution,
                                                    this->dof_handler);
      this->vof_subequations_interface->solve_specific_subequation(
        VOFSubequationsID::phase_gradient_projection);
      this->vof_subequations_interface->solve_specific_subequation(
        VOFSubequationsID::curvature_projection);
    }

  // Reset algebraic interface reinitialization output directory;
  // if it does not exist, create it.
  if (simulation_parameters.multiphysics.vof_parameters.regularization_method
        .algebraic_interface_reinitialization.enable &&
      simulation_parameters.multiphysics.vof_parameters.regularization_method
        .algebraic_interface_reinitialization.output_reinitialization_steps)
    {
      auto mpi_communicator = this->triangulation->get_mpi_communicator();
      const std::string folder =
        this->simulation_parameters.simulation_control.output_folder +
        "/algebraic-reinitialization-steps-output/";

      struct stat buffer;

      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          if (stat(folder.c_str(), &buffer) != 0)
            {
              create_output_folder(folder);
            }
          else
            {
              delete_output_folder(folder);
              create_output_folder(folder);
            }
        }
    }

  percolate_time_vectors();
}


template <int dim>
void
VolumeOfFluid<dim>::solve_linear_system(const bool initial_step,
                                        const bool /*renewed_matrix*/)
{
  TimerOutput::Scope t(this->computing_timer, "Solve linear system");

  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;

  const double absolute_residual =
    simulation_parameters.linear_solver.at(PhysicsID::VOF).minimum_residual;
  const double relative_residual =
    simulation_parameters.linear_solver.at(PhysicsID::VOF).relative_residual;

  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::VOF).verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const unsigned int ilu_fill =
    simulation_parameters.linear_solver.at(PhysicsID::VOF).ilu_precond_fill;
  const double ilu_atol =
    simulation_parameters.linear_solver.at(PhysicsID::VOF).ilu_precond_atol;
  const double ilu_rtol =
    simulation_parameters.linear_solver.at(PhysicsID::VOF).ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  TrilinosWrappers::PreconditionILU ilu_preconditioner;

  ilu_preconditioner.initialize(this->system_matrix, preconditionerOptions);

  GlobalVectorType completely_distributed_solution(this->locally_owned_dofs,
                                                   mpi_communicator);

  SolverControl solver_control(
    simulation_parameters.linear_solver.at(PhysicsID::VOF).max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false,
    simulation_parameters.linear_solver.at(PhysicsID::VOF).max_krylov_vectors);

  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);

  solver.solve(this->system_matrix,
               completely_distributed_solution,
               this->system_rhs,
               ilu_preconditioner);

  if (simulation_parameters.linear_solver.at(PhysicsID::VOF).verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() << std::endl;
    }

  // Update constraints and newton vectors
  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
}

// This function is explained in detail in step-41 of deal.II tutorials
template <int dim>
void
VolumeOfFluid<dim>::update_solution_and_constraints(GlobalVectorType &solution)
{
  // This is a penalty parameter for limiting the phase fraction
  // in the range of [0,1]. According to step 41, this parameter depends
  // on the problem itself and needs to be chosen large enough (for example,
  // there is no convergence using the penalty_parameter = 1)
  const double penalty_parameter = 100;

  GlobalVectorType lambda(this->locally_owned_dofs,
                          this->triangulation->get_mpi_communicator());

  nodal_phase_fraction_owned = solution;

  complete_system_matrix_phase_fraction.residual(lambda,
                                                 nodal_phase_fraction_owned,
                                                 system_rhs_phase_fraction);

  this->bounding_constraints.clear();

  std::vector<bool> dof_touched(this->dof_handler.n_dofs(), false);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
               ++v)
            {
              Assert(this->dof_handler.get_fe().dofs_per_cell ==
                       GeometryInfo<dim>::vertices_per_cell,
                     ExcNotImplemented());
              const unsigned int dof_index = cell->vertex_dof_index(v, 0);
              if (this->locally_owned_dofs.is_element(dof_index))
                {
                  const double solution_value =
                    nodal_phase_fraction_owned(dof_index);
                  if (lambda(dof_index) +
                        penalty_parameter *
                          mass_matrix_phase_fraction(dof_index, dof_index) *
                          (solution_value - this->phase_upper_bound) >
                      0)
                    {
                      this->bounding_constraints.add_line(dof_index);
                      this->bounding_constraints.set_inhomogeneity(
                        dof_index, this->phase_upper_bound);
                      nodal_phase_fraction_owned(dof_index) =
                        this->phase_upper_bound;
                      lambda(dof_index) = 0;
                    }
                  else if (lambda(dof_index) +
                             penalty_parameter *
                               mass_matrix_phase_fraction(dof_index,
                                                          dof_index) *
                               (solution_value - this->phase_lower_bound) <
                           0)
                    {
                      this->bounding_constraints.add_line(dof_index);
                      this->bounding_constraints.set_inhomogeneity(
                        dof_index, this->phase_lower_bound);
                      nodal_phase_fraction_owned(dof_index) =
                        this->phase_lower_bound;
                      lambda(dof_index) = 0;
                    }
                }
            }
        }
    }
  solution = nodal_phase_fraction_owned;
  this->bounding_constraints.close();
}


template <int dim>
void
VolumeOfFluid<dim>::assemble_L2_projection_interface_sharpening(
  GlobalVectorType &solution,
  const double      sharpening_threshold)
{
  const double interface_sharpness =
    this->simulation_parameters.multiphysics.vof_parameters
      .regularization_method.sharpening.interface_sharpness;

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->cell_quadrature,
                              update_values | update_JxW_values);

  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int n_q_points    = this->cell_quadrature->size();
  FullMatrix<double> local_matrix_phase_fraction(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_phase_fraction(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_phase(dofs_per_cell);

  std::vector<double> phase_values(n_q_points);

  system_rhs_phase_fraction    = 0;
  system_matrix_phase_fraction = 0;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);

          local_matrix_phase_fraction = 0;
          local_rhs_phase_fraction    = 0;

          fe_values_vof.get_function_values(solution, phase_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              auto phase_value = phase_values[q];
              // phase_value      = std::min(std::max(phase_value, 0.), 1.);

              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_phase[k] = fe_values_vof.shape_value(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_phase_fraction(i, j) +=
                        (phi_phase[j] * phi_phase[i]) * fe_values_vof.JxW(q);
                    }

                  // $$ (if 0 <= \phi <= c)  {\Phi = c ^ (1 - \alpha) * (\phi
                  // ^ \alpha)}$$
                  // $$ (if c <  \phi <= 1)  {\Phi = 1 - (1 - c) ^ (1 -
                  // \alpha)
                  // * (1 - \phi) ^ \alpha}
                  if (phase_value >= 0 && phase_value <= sharpening_threshold)
                    local_rhs_phase_fraction(i) +=
                      std::pow(sharpening_threshold,
                               (1. - interface_sharpness)) *
                      std::pow(phase_value, interface_sharpness) *
                      phi_phase[i] * fe_values_vof.JxW(q);
                  else
                    {
                      local_rhs_phase_fraction(i) +=
                        (1 -
                         std::pow((1. - sharpening_threshold),
                                  (1. - interface_sharpness)) *
                           std::pow((1. - phase_value), interface_sharpness)) *
                        phi_phase[i] * fe_values_vof.JxW(q);
                    }
                }
            }

          cell->get_dof_indices(local_dof_indices);
          this->nonzero_constraints.distribute_local_to_global(
            local_matrix_phase_fraction,
            local_rhs_phase_fraction,
            local_dof_indices,
            system_matrix_phase_fraction,
            system_rhs_phase_fraction);
        }
    }

  system_matrix_phase_fraction.compress(VectorOperation::add);
  system_rhs_phase_fraction.compress(VectorOperation::add);
}

template <int dim>
void
VolumeOfFluid<dim>::solve_interface_sharpening(GlobalVectorType &solution)
{
  // Solve the L2 projection system
  const double linear_solver_tolerance = 1e-15;

  if (this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.verbosity ==
      Parameters::Verbosity::extra_verbose)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  GlobalVectorType completely_distributed_phase_fraction_solution(
    this->locally_owned_dofs, triangulation->get_mpi_communicator());


  SolverControl solver_control(
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF).max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverCG solver(solver_control);

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const unsigned int ilu_fill =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .ilu_precond_rtol;

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix_phase_fraction,
                                 preconditionerOptions);

  solver.solve(system_matrix_phase_fraction,
               completely_distributed_phase_fraction_solution,
               system_rhs_phase_fraction,
               *ilu_preconditioner);

  if (this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.verbosity ==
      Parameters::Verbosity::extra_verbose)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  this->nonzero_constraints.distribute(
    completely_distributed_phase_fraction_solution);
  solution = completely_distributed_phase_fraction_solution;
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_mass_matrix(
  TrilinosWrappers::SparseMatrix &mass_matrix)
{
  QGauss<dim> quadrature_formula(this->cell_quadrature->size());

  FEValues<dim> fe_values_vof(*this->fe,
                              quadrature_formula,
                              update_values | update_JxW_values);


  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);
          cell_matrix = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_matrix(i, i) +=
                (fe_values_vof.shape_value(i, q) *
                 fe_values_vof.shape_value(i, q) * fe_values_vof.JxW(q));
          cell->get_dof_indices(local_dof_indices);
          this->nonzero_constraints.distribute_local_to_global(
            cell_matrix, local_dof_indices, mass_matrix);
        }
    }
}


template <int dim>
void
VolumeOfFluid<dim>::apply_phase_filter(
  const GlobalVectorType &original_solution,
  GlobalVectorType       &filtered_solution)
{
  // Initializations
  auto mpi_communicator = this->triangulation->get_mpi_communicator();
  GlobalVectorType filtered_solution_owned(this->locally_owned_dofs,
                                           mpi_communicator);
  filtered_solution_owned = original_solution;
  filtered_solution.reinit(original_solution);

  // Create filter object
  filter = VolumeOfFluidFilterBase::model_cast(
    this->simulation_parameters.multiphysics.vof_parameters.phase_filter);

  // Apply filter to the solution
  for (auto p : this->locally_owned_dofs)
    {
      filtered_solution_owned[p] =
        filter->filter_phase(filtered_solution_owned[p]);
    }
  filtered_solution = filtered_solution_owned;

  if (this->simulation_parameters.multiphysics.vof_parameters.phase_filter
        .verbosity == Parameters::Verbosity::verbose)
    {
      this->pcout << "Filtered phase values: " << std::endl;
      for (const double filtered_phase : filtered_solution)
        {
          this->pcout << filtered_phase << std::endl;
        }
    }
}

template <int dim>
void
VolumeOfFluid<dim>::reinitialize_interface_with_algebraic_method()
{
  TimerOutput::Scope t(this->computing_timer, "Algebraic reinitialization");

  // Reinitialize previous VOF solution
  // (this is only coherent with BDF1 and BDF2)
  if (this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.frequency > 1)
    {
      auto mpi_communicator = this->triangulation->get_mpi_communicator();

      GlobalVectorType previous_reinitialized_solution_owned(
        this->locally_owned_dofs, mpi_communicator);
      GlobalVectorType previous_filtered_solution(this->locally_owned_dofs,
                                                  this->locally_relevant_dofs,
                                                  mpi_communicator);

      // Apply filter to previous solution
      apply_phase_filter(this->previous_solutions[0],
                         previous_filtered_solution);

      // Set VOF information in the VOF subequations interface
      this->vof_subequations_interface
        ->set_vof_filtered_solution_and_dof_handler(previous_filtered_solution,
                                                    this->dof_handler);
      this->vof_subequations_interface->set_vof_solution(
        this->previous_solutions[0]);

      // Solve phase gradient and curvature projections followed by algebraic
      // interface reinitialization steps
      this->vof_subequations_interface->solve();

      // Overwrite VOF previous solution with the reinitialized result
      FETools::interpolate(
        this->vof_subequations_interface->get_dof_handler(
          VOFSubequationsID::algebraic_interface_reinitialization),
        this->vof_subequations_interface->get_solution(
          VOFSubequationsID::algebraic_interface_reinitialization),
        this->dof_handler,
        this->nonzero_constraints,
        previous_reinitialized_solution_owned);
      this->previous_solutions[0] = previous_reinitialized_solution_owned;
    }

  // Apply filter to solution and set VOF information in the subequation
  // interface
  apply_phase_filter(this->present_solution, this->filtered_solution);
  this->vof_subequations_interface->set_vof_filtered_solution_and_dof_handler(
    this->filtered_solution, this->dof_handler);
  this->vof_subequations_interface->set_vof_solution(this->present_solution);

  // Solve phase gradient and curvature projections followed by algebraic
  // interface reinitialization steps
  this->vof_subequations_interface->solve();

  // Overwrite the VOF solution with the algebraic interface reinitialization
  FETools::interpolate(
    this->vof_subequations_interface->get_dof_handler(
      VOFSubequationsID::algebraic_interface_reinitialization),
    this->vof_subequations_interface->get_solution(
      VOFSubequationsID::algebraic_interface_reinitialization),
    this->dof_handler,
    this->nonzero_constraints,
    this->local_evaluation_point);
  this->present_solution = this->local_evaluation_point;
}

template <int dim>
void
VolumeOfFluid<dim>::compute_level_set_from_phase_fraction(
  const GlobalVectorType &solution,
  GlobalVectorType       &level_set_solution)
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  GlobalVectorType level_set_owned(this->locally_owned_dofs, mpi_communicator);

  for (auto p : this->locally_owned_dofs)
    {
      const double phase = solution[p];
      level_set_owned[p] = (-2.0 * phase + 1.0);
    }

  this->nonzero_constraints.distribute(level_set_owned);

  level_set_solution = level_set_owned;
}

template <int dim>
void
VolumeOfFluid<dim>::compute_phase_fraction_from_level_set(
  const GlobalVectorType &level_set_solution,
  GlobalVectorType       &phase_fraction_solution)
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  GlobalVectorType solution_owned(this->locally_owned_dofs, mpi_communicator);

  for (auto p : this->locally_owned_dofs)
    {
      const double signed_dist = level_set_solution[p];
      solution_owned[p] =
        signed_distance_transformation->transform_signed_distance(signed_dist);
    }
  this->nonzero_constraints.distribute(solution_owned);

  phase_fraction_solution = solution_owned;
}

template <int dim>
void
VolumeOfFluid<dim>::reinitialize_interface_with_geometric_method()
{
  TimerOutput::Scope t(this->computing_timer, "Geometric reinitialization");

  if (simulation_parameters.multiphysics.vof_parameters.regularization_method
        .verbosity != Parameters::Verbosity::quiet)
    {
      announce_string(this->pcout, "VOF geometric interface reinitialization");
      this->pcout << "In redistanciation of the previous solution ..."
                  << std::endl;
    }

  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  signed_distance_solver->setup_dofs();

  GlobalVectorType previous_level_set(locally_owned_dofs,
                                      locally_relevant_dofs,
                                      mpi_communicator);

  if (simulation_parameters.multiphysics.vof_parameters.regularization_method
        .frequency != 1)
    {
      signed_distance_solver->set_level_set_from_background_mesh(
        dof_handler, this->previous_solutions[0]);

      signed_distance_solver->solve();

      GlobalVectorType previous_level_set_owned(this->locally_owned_dofs,
                                                mpi_communicator);

      FETools::interpolate(signed_distance_solver->dof_handler,
                           signed_distance_solver->get_signed_distance(),
                           this->dof_handler,
                           this->nonzero_constraints,
                           previous_level_set_owned);


      previous_level_set = previous_level_set_owned;

      compute_phase_fraction_from_level_set(previous_level_set,
                                            this->previous_solutions[0]);
    }

  if (simulation_parameters.multiphysics.vof_parameters.regularization_method
        .verbosity != Parameters::Verbosity::quiet)
    this->pcout << "In redistanciation of the present solution ..."
                << std::endl;

  signed_distance_solver->set_level_set_from_background_mesh(
    dof_handler, this->present_solution);

  signed_distance_solver->solve();

  GlobalVectorType level_set_owned(this->locally_owned_dofs, mpi_communicator);

  FETools::interpolate(signed_distance_solver->dof_handler,
                       signed_distance_solver->get_signed_distance(),
                       this->dof_handler,
                       this->nonzero_constraints,
                       level_set_owned);


  this->level_set = level_set_owned;

  compute_phase_fraction_from_level_set(this->level_set,
                                        this->present_solution);
}


template class VolumeOfFluid<2>;
template class VolumeOfFluid<3>;
