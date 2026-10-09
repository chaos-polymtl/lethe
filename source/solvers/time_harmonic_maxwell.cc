// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/time_harmonic_maxwell.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_system.h>

#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/vector_tools.h>

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
                    TimerOutput::never,
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
  this->pcout << std::setprecision(simulation_control->get_log_precision())
              << std::scientific;

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
      AssertThrow(dim == 3, TimeHarmonicMaxwellDimensionNotSupported(dim));

      AssertThrow(
        simulation_parameters.fem_parameters.electromagnetics_trial_degree <
          simulation_parameters.fem_parameters.electromagnetics_test_degree,
        ExcMessage(
          "The DPG method requires the test space to be of higher order than the trial space."));

      // Usual case, for quad/hex meshes
      fe_trial_interior = std::make_shared<FESystem<dim>>(
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_degree) ^
          dim,
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_degree) ^
          dim,
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_degree) ^
          dim,
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_degree) ^
          dim);
      fe_trial_skeleton = std::make_shared<FESystem<dim>>(
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_degree),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_degree),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_degree),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_degree));
      fe_test = std::make_shared<FESystem<dim>>(
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_degree),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_degree),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_degree),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_degree));
      mapping = std::make_shared<MappingQ<dim>>(fe_trial_interior->degree);
      cell_quadrature = std::make_shared<QGauss<dim>>(fe_test->degree + 1);
      face_quadrature = std::make_shared<QGauss<dim - 1>>(fe_test->degree + 1);
    }

  // Initialize solutions shared_ptr
  present_solution          = std::make_shared<GlobalVectorType>();
  present_solution_skeleton = std::make_shared<GlobalVectorType>();

  // Allocate solution transfer
  solution_transfer = std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(
    *dof_handler_trial_interior);

  // We may need the temperature field for the physical properties. If so,
  // we create the corresponding FEValues object to evaluate the temperature at
  // the quadrature points and its container to store the values.
  const PhysicalPropertiesManager &physical_properties_manager =
    simulation_parameters.physical_properties_manager;
  unsigned int number_of_material_ids =
    physical_properties_manager.get_number_of_fluids() +
    physical_properties_manager.get_number_of_solids() - 1;

  needs_temperature = false;
  for (unsigned int material_id = 0; material_id <= number_of_material_ids;
       ++material_id)
    {
      needs_temperature =
        (physical_properties_manager.get_electric_conductivity(0, material_id)
           ->depends_on(field::temperature) ||
         physical_properties_manager
           .get_electric_permittivity_real(0, material_id)
           ->depends_on(field::temperature) ||
         physical_properties_manager
           .get_electric_permittivity_imag(0, material_id)
           ->depends_on(field::temperature) ||
         physical_properties_manager
           .get_magnetic_permeability_real(0, material_id)
           ->depends_on(field::temperature) ||
         physical_properties_manager
           .get_magnetic_permeability_imag(0, material_id)
           ->depends_on(field::temperature));

      if (needs_temperature)
        {
          // If temperature-dependent electromagnetic properties are requested,
          // heat transfer must be enabled.

          const std::vector<PhysicsID> active_physics_ids =
            multiphysics->get_active_physics();
          AssertThrow(
            std::ranges::find(active_physics_ids, PhysicsID::heat_transfer) !=
              active_physics_ids.end(),
            ExcMessage(
              "User-defined temperature-dependent electromagnetic properties requested without heat transfer physics. Enable heat transfer in the multiphysics subsection."));
          break;
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::print_THM_setup_memory(
  const TrilinosWrappers::SparsityPattern &sparsity_pattern)
{
  auto       mpi_communicator = triangulation->get_mpi_communicator();
  const auto this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);
  constexpr double bytes_to_gb = 1.0 / (1024.0 * 1024.0 * 1024.0);

  // Fetch memory consumption information on each process
  const auto present_solution_memory =
    this->present_solution->memory_consumption() * bytes_to_gb;
  const auto present_solution_skeleton_memory =
    this->present_solution_skeleton->memory_consumption() * bytes_to_gb;
  const auto system_rhs_memory =
    this->system_rhs.memory_consumption() * bytes_to_gb;
  const auto sparsity_pattern_memory =
    sparsity_pattern.n_nonzero_elements() *
    sizeof(TrilinosWrappers::types::int_type) *
    bytes_to_gb; // We use a proxy for the memory consumption of the
                 // sparsity pattern based on the number of non-zero
                 // elements and the size of the integer type used to store
                 // the sparsity pattern, since the
                 // TrilinosWrappers::SparsityPattern class does not have a
                 // memory_consumption() function implemented.
  const auto system_matrix_memory =
    this->system_matrix.memory_consumption() * bytes_to_gb;
  const auto dof_handler_trial_interior_memory =
    this->dof_handler_trial_interior->memory_consumption() * bytes_to_gb;
  const auto dof_handler_trial_skeleton_memory =
    this->dof_handler_trial_skeleton->memory_consumption() * bytes_to_gb;
  const auto dof_handler_test_memory =
    this->dof_handler_test->memory_consumption() * bytes_to_gb;

  // Gather memory consumption information from all ranks to rank 0
  const auto present_solution_memory_by_rank =
    Utilities::MPI::gather(mpi_communicator, present_solution_memory, 0);
  const auto present_solution_skeleton_memory_by_rank =
    Utilities::MPI::gather(mpi_communicator,
                           present_solution_skeleton_memory,
                           0);
  const auto system_rhs_memory_by_rank =
    Utilities::MPI::gather(mpi_communicator, system_rhs_memory, 0);
  const auto sparsity_pattern_memory_by_rank =
    Utilities::MPI::gather(mpi_communicator, sparsity_pattern_memory, 0);
  const auto system_matrix_memory_by_rank =
    Utilities::MPI::gather(mpi_communicator, system_matrix_memory, 0);
  const auto dof_handler_trial_interior_memory_by_rank =
    Utilities::MPI::gather(mpi_communicator,
                           dof_handler_trial_interior_memory,
                           0);
  const auto dof_handler_trial_skeleton_memory_by_rank =
    Utilities::MPI::gather(mpi_communicator,
                           dof_handler_trial_skeleton_memory,
                           0);
  const auto dof_handler_test_memory_by_rank =
    Utilities::MPI::gather(mpi_communicator, dof_handler_test_memory, 0);

  // Sum memory consumption across all ranks to get total memory usage
  const auto present_solution_memory_total =
    Utilities::MPI::sum(present_solution_memory, mpi_communicator);
  const auto present_solution_skeleton_memory_total =
    Utilities::MPI::sum(present_solution_skeleton_memory, mpi_communicator);
  const auto system_rhs_memory_total =
    Utilities::MPI::sum(system_rhs_memory, mpi_communicator);
  const auto sparsity_pattern_memory_total =
    Utilities::MPI::sum(sparsity_pattern_memory, mpi_communicator);
  const auto system_matrix_memory_total =
    Utilities::MPI::sum(system_matrix_memory, mpi_communicator);
  const auto dof_handler_trial_interior_memory_total =
    Utilities::MPI::sum(dof_handler_trial_interior_memory, mpi_communicator);
  const auto dof_handler_trial_skeleton_memory_total =
    Utilities::MPI::sum(dof_handler_trial_skeleton_memory, mpi_communicator);
  const auto dof_handler_test_memory_total =
    Utilities::MPI::sum(dof_handler_test_memory, mpi_communicator);

  // Print memory consumption information on rank 0
  if (this_mpi_process == 0)
    {
      const auto total_memory =
        present_solution_memory_total + present_solution_skeleton_memory_total +
        system_rhs_memory_total + sparsity_pattern_memory_total +
        system_matrix_memory_total + dof_handler_trial_interior_memory_total +
        dof_handler_trial_skeleton_memory_total + dof_handler_test_memory_total;

      announce_string(this->pcout,
                      "Time-Harmonic Maxwell Memory Diagnostics",
                      65,
                      '=');

      // When debugging memory issues if rank is needed turn this parameter
      // to true to have the memory consumption of each rank printed,
      // otherwise only the total memory consumption across all ranks will
      // be printed.
      bool print_memory_by_rank = true;
      if (print_memory_by_rank)
        {
          this->pcout << "  *** By rank memory *** " << std::endl;

          print_memory_consumption(this->pcout,
                                   "present_solution",
                                   present_solution_memory_by_rank.size(),
                                   present_solution_memory_by_rank);
          print_memory_consumption(
            this->pcout,
            "present_solution_skeleton",
            present_solution_skeleton_memory_by_rank.size(),
            present_solution_skeleton_memory_by_rank);
          print_memory_consumption(this->pcout,
                                   "system_rhs",
                                   system_rhs_memory_by_rank.size(),
                                   system_rhs_memory_by_rank);

          print_memory_consumption(this->pcout,
                                   "sparsity_pattern",
                                   sparsity_pattern_memory_by_rank.size(),
                                   sparsity_pattern_memory_by_rank);
          print_memory_consumption(this->pcout,
                                   "system_matrix",
                                   system_matrix_memory_by_rank.size(),
                                   system_matrix_memory_by_rank);
          print_memory_consumption(
            this->pcout,
            "dof_handler_trial_interior",
            dof_handler_trial_interior_memory_by_rank.size(),
            dof_handler_trial_interior_memory_by_rank);
          print_memory_consumption(
            this->pcout,
            "dof_handler_trial_skeleton",
            dof_handler_trial_skeleton_memory_by_rank.size(),
            dof_handler_trial_skeleton_memory_by_rank);
          print_memory_consumption(this->pcout,
                                   "dof_handler_test",
                                   dof_handler_test_memory_by_rank.size(),
                                   dof_handler_test_memory_by_rank);

          this->pcout
            << "=================================================================="
            << std::endl;
        }

      this->pcout << "  *** Total memory consumption (GB) *** " << std::endl;
      this->pcout << "  present_solution : " << present_solution_memory_total
                  << std::endl;
      this->pcout << "  present_solution_skeleton : "
                  << present_solution_skeleton_memory_total << std::endl;
      this->pcout << "  system_rhs : " << system_rhs_memory_total << std::endl;
      this->pcout << "  sparsity_pattern : " << sparsity_pattern_memory_total
                  << std::endl;
      this->pcout << "  system_matrix : " << system_matrix_memory_total
                  << std::endl;
      this->pcout << "  dof_handler_trial_interior : "
                  << dof_handler_trial_interior_memory_total << std::endl;
      this->pcout << "  dof_handler_trial_skeleton : "
                  << dof_handler_trial_skeleton_memory_total << std::endl;
      this->pcout << "  dof_handler_test : " << dof_handler_test_memory_total
                  << std::endl;
      this->pcout << "  Total DPG solver memory consumption : " << total_memory
                  << std::endl;
      this->pcout
        << " =================================================================="
        << std::endl;
    }
}


template <int dim>
void
TimeHarmonicMaxwell<dim>::update_material_properties(
  const PhysicalPropertiesManager &physical_properties_manager,
  const std::map<field, double>   &field_values,
  const unsigned int               material_id,
  std::complex<double>            &effective_electric_permittivity,
  std::complex<double>            &effective_magnetic_permeability)
{
  effective_electric_permittivity = {
    physical_properties_manager.get_electric_permittivity_real(0, material_id)
      ->value(field_values),
    physical_properties_manager.get_electric_permittivity_imag(0, material_id)
        ->value(field_values) +
      physical_properties_manager.get_electric_conductivity(0, material_id)
        ->value(field_values)};

  effective_magnetic_permeability = {
    physical_properties_manager.get_magnetic_permeability_real(0, material_id)
      ->value(field_values),
    physical_properties_manager.get_magnetic_permeability_imag(0, material_id)
      ->value(field_values)};
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::compute_electromagnetic_scaling(
  const PhysicalPropertiesManager &physical_properties_manager)
{
  // Define reusable variables for the computation of the scaling factor.
  auto        mpi_communicator = this->triangulation->get_mpi_communicator();
  const auto &electromagnetic_parameters =
    this->simulation_parameters.multiphysics.time_harmonic_maxwell_parameters;
  constexpr double void_impedance =
    4 * numbers::PI *
    29.9792458; // Void impedance in Ohms, the unit conversion of this value is
                // done when computing the scaling factor. The factor 29.9792458
                // is the speed of light in cm/ns and is used to make it
                // convenient when multiplying by the void permeability to
                // obtain the void impedance in the desired units (e.g., Ohms).

  // Initialize the amplitude. The default value is put to zero.
  // If there are no waveguide inlets, the maximum of amplitude will be taken
  // between 0 and a strictly positive value ensuring that the strictly positive
  // value remains.
  double max_port_electric_amplitude = 0.0;

  // If there are waveguide inlets we need to compute their scaling from the
  // input power provided by the user.
  if (electromagnetic_parameters.number_of_waveguide_inlets > 0)
    {
      const auto &waveguide_powers = electromagnetic_parameters.waveguide_power;

      // The following vector will store the modal power for each waveguide
      // inlet boundary condition computed from integrating the Poynting vector
      // of the incident wave over the inlet face.
      std::vector<double> waveguide_modal_powers(
        electromagnetic_parameters.number_of_waveguide_inlets, 0.);

      // Resize the vector that will store the waveguide port electric field
      // amplitudes. Those are obtain by computing the square root of the ratio
      // between the input power provided by the user and the modal power
      // computed from integrating the Poynting vector of the incident wave over
      // the inlet face.
      this->waveguide_ports_electric_amplitudes = std::vector<double>(
        electromagnetic_parameters.number_of_waveguide_inlets, 0.);

      // We need to define a FEFaceValues object to perform the integral of the
      // Poynting vector on faces
      FEFaceValues<dim>  fe_face_values_trial_skeleton(*this->mapping,
                                                      *this->fe_trial_skeleton,
                                                      *this->face_quadrature,
                                                      update_quadrature_points |
                                                        update_normal_vectors |
                                                        update_JxW_values);
      const unsigned int n_face_q_points =
        fe_face_values_trial_skeleton.n_quadrature_points;
      BoundaryConditions::BoundaryType bc_type;

      // Containers for electromagnetic quantities at the quadrature points
      Tensor<1, dim, std::complex<double>> E_inc;
      Tensor<1, dim, std::complex<double>> H_inc;
      std::vector<std::complex<double>>    effective_electric_permittivities(
        n_face_q_points);
      std::vector<std::complex<double>> effective_magnetic_permeabilities(
        n_face_q_points);
      unsigned int material_id;
      bool         cell_material_needs_temperature = false;

      // The temperature field is only needed when at least one electromagnetic
      // property depends on temperature.

      std::unique_ptr<FEFaceValues<dim>> fe_face_values_temperature;
      const DoFHandler<dim>             *dof_handler_temperature = nullptr;
      const GlobalVectorType            *temperature_solution    = nullptr;
      std::vector<double> temperature_face_values(n_face_q_points);
      std::map<field, std::vector<double>> field_values_vector;

      if (needs_temperature)
        {
          dof_handler_temperature =
            &this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);
          fe_face_values_temperature = std::make_unique<FEFaceValues<dim>>(
            *this->mapping,
            dof_handler_temperature->get_fe(),
            *this->face_quadrature,
            update_values | update_quadrature_points);
          temperature_solution =
            &this->multiphysics->get_solution(PhysicsID::heat_transfer);
        }


      // Here we perform the Poynting vector integration over the waveguide
      // inlet faces.
      for (const auto &cell : triangulation->active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              material_id = cell->material_id();

              // We check if the physical properties depend on the temperature
              // field. If so, we will need to evaluate the temperature field at
              // the quadrature points.
              cell_material_needs_temperature =
                physical_properties_manager
                  .get_electric_conductivity(0, material_id)
                  ->depends_on(field::temperature) ||
                physical_properties_manager
                  .get_electric_permittivity_real(0, material_id)
                  ->depends_on(field::temperature) ||
                physical_properties_manager
                  .get_electric_permittivity_imag(0, material_id)
                  ->depends_on(field::temperature) ||
                physical_properties_manager
                  .get_magnetic_permeability_real(0, material_id)
                  ->depends_on(field::temperature) ||
                physical_properties_manager
                  .get_magnetic_permeability_imag(0, material_id)
                  ->depends_on(field::temperature);

              for (const auto &face : cell->face_iterators())
                {
                  if (!(face->at_boundary()))
                    continue;

                  fe_face_values_trial_skeleton.reinit(cell, face);

                  if (cell_material_needs_temperature)
                    {
                      const typename DoFHandler<dim>::active_cell_iterator
                        cell_temperature = cell->as_dof_handler_iterator(
                          *dof_handler_temperature);
                      fe_face_values_temperature->reinit(cell_temperature,
                                                         face);
                      fe_face_values_temperature->get_function_values(
                        *temperature_solution, temperature_face_values);
                    }
                  else
                    {
                      std::ranges::fill(temperature_face_values, 0.);
                    }

                  bc_type =
                    this->simulation_parameters
                      .boundary_conditions_time_harmonic_electromagnetics.type
                      .at(face->boundary_id());

                  // We only compute the modal power if its a port boundary
                  if (!(bc_type ==
                        BoundaryConditions::BoundaryType::waveguide_port))
                    continue;

                  // Get the index of the waveguide port boundary
                  // condition to apply the correct mode parameters
                  unsigned int boundary_index = std::distance(
                    electromagnetic_parameters.waveguide_boundary_ids.begin(),
                    std::ranges::find(
                      electromagnetic_parameters.waveguide_boundary_ids,
                      face->boundary_id()));

                  // We update the material properties for the current face
                  // quadrature points.
                  field_values_vector[field::temperature] =
                    temperature_face_values;

                  compute_effective_electromagnetic_properties(
                    physical_properties_manager,
                    field_values_vector,
                    material_id,
                    effective_electric_permittivities,
                    effective_magnetic_permeabilities);

                  // Loop over all face quadrature points
                  for (unsigned int q_point = 0; q_point < n_face_q_points;
                       ++q_point)
                    {
                      // Initialize reusable variables
                      const auto &position =
                        fe_face_values_trial_skeleton.quadrature_point(q_point);
                      const auto &normal =
                        fe_face_values_trial_skeleton.normal_vector(q_point);
                      const double JxW_face =
                        fe_face_values_trial_skeleton.JxW(q_point);

                      std::tie(E_inc, H_inc) =
                        compute_waveguide_port_incident_fields(
                          electromagnetic_parameters,
                          position,
                          normal,
                          effective_electric_permittivities[q_point],
                          effective_magnetic_permeabilities[q_point],
                          boundary_index);

                      // The Poynting vector is given by S = 0.5 * Re(E
                      // x H*), where H* is the complex conjugate of H.
                      // The normal points outward from the domain, so
                      // the power entering the domain is given by the
                      // negative of the flux of the Poynting vector
                      // through the face (that's why we have a negative
                      // sign in front of the integral). The cross product
                      // is written explicitly in 3D, so it is only compiled
                      // for dim = 3 (the waveguide port fields are not
                      // defined in 2D).
                      if constexpr (dim == 3)
                        waveguide_modal_powers[boundary_index] -=
                          0.5 *
                          std::real(
                            normal[0] * (E_inc[1] * std::conj(H_inc[2]) -
                                         E_inc[2] * std::conj(H_inc[1])) +
                            normal[1] * (E_inc[2] * std::conj(H_inc[0]) -
                                         E_inc[0] * std::conj(H_inc[2])) +
                            normal[2] * (E_inc[0] * std::conj(H_inc[1]) -
                                         E_inc[1] * std::conj(H_inc[0]))) *
                          JxW_face;
                    }
                }
            }
        }

      // Now we need to aggregate from all the MPI processes to get the global
      // modal power for each waveguide inlet
      for (unsigned int inlets = 0;
           inlets < electromagnetic_parameters.number_of_waveguide_inlets;
           ++inlets)
        {
          // The power amplitude comes from the non-dimensionalization of power:
          // P = 0.5 * int(ExH*)= E_0 * H_0 * P_modal. So we can isolate the
          // electric field amplitude E_0. We don't need to multiply by the
          // dimensionality because we assume user input in W and impedance in
          // Ohms so the resulting units of this integral will always be
          // Volt/length unit (of the mesh).
          this->waveguide_ports_electric_amplitudes[inlets] =
            std::sqrt(waveguide_powers[inlets] * void_impedance /
                      Utilities::MPI::sum(waveguide_modal_powers[inlets],
                                          mpi_communicator));
        }

      max_port_electric_amplitude =
        *std::ranges::max_element(this->waveguide_ports_electric_amplitudes);
    }

  // Now we need to define the electromagnetic scaling factor according to the
  // type chosen by the user.
  if (electromagnetic_parameters.electromagnetic_scaling_type ==
      Parameters::ElectromagneticScalingType::electric_field)
    {
      // We assume that the electric field amplitude provided by the user is MKS
      // (V/m).
      this->electromagnetic_scaling =
        electromagnetic_parameters.electric_field_amplitude;
    }
  else if (electromagnetic_parameters.electromagnetic_scaling_type ==
           Parameters::ElectromagneticScalingType::magnetic_field)
    {
      // We assume that the magnetic field amplitude provided by the user is in
      // MKS (A/m). We convert it to the electric field for convenience to not
      // carry multiple conversion factors when applying the scaling to the
      // solution.
      this->electromagnetic_scaling =
        void_impedance * electromagnetic_parameters.magnetic_field_amplitude;
    }
  else if (electromagnetic_parameters.electromagnetic_scaling_type ==
           Parameters::ElectromagneticScalingType::power)
    { // From the Poynting vector integration and assumptions made in the
      // computation of the waveguide port excitation, the unit here are in V/L,
      // with L being the length unit of the mesh. Therefore, we need to divide
      // by the length unit to get the correct scaling factor in V/m.
      this->electromagnetic_scaling =
        max_port_electric_amplitude /
        this->simulation_parameters.dimensionality.length;
    }
  else if (electromagnetic_parameters.electromagnetic_scaling_type ==
           Parameters::ElectromagneticScalingType::none)
    {
      // If there is at least one waveguide inlet, we need to apply the scaling
      // because some inlet may have different power even if the user chose to
      // not apply the scaling to the solution after computation.
      if (electromagnetic_parameters.number_of_waveguide_inlets > 0)
        {
          this->electromagnetic_scaling = max_port_electric_amplitude;
        }
      else
        {
          this->electromagnetic_scaling = 1.0;
        }
    }
  else
    {
      AssertThrow(false, ExcMessage("Unknown electromagnetic scaling type."));
    }
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::scale_solution_components(
  const DoFHandler<dim> &dof_handler,
  GlobalVectorType      &solution)
{
  const auto &fe = dof_handler.get_fe();

  if (this->simulation_parameters.multiphysics.time_harmonic_maxwell_parameters
        .electromagnetic_scaling_type ==
      Parameters::ElectromagneticScalingType::none)
    return;

  constexpr double void_impedance = 4 * numbers::PI * 29.9792458;
  // Here, everything is in MKS units, so we need to convert it to the user
  // desired units.
  const double electric_scale =
    this->electromagnetic_scaling *
    this->simulation_parameters.dimensionality.electric_amplitude_scaling;
  const double magnetic_scale =
    this->electromagnetic_scaling / void_impedance *
    this->simulation_parameters.dimensionality.magnetic_amplitude_scaling;

  const IndexSet   &locally_owned_dofs = solution.locally_owned_elements();
  std::vector<bool> dof_scaled(locally_owned_dofs.n_elements(), false);
  std::vector<types::global_dof_index> local_dof_indices(fe.n_dofs_per_cell());

  // We loop over all the cells and apply the correct scaling to the electric
  // and magnetic field components of the solution. We need to do so because the
  // solution is stored as a single vector but contains different physical
  // fields (electric and magnetic) that need to be scaled differently and the
  // system_to_base_index function which allows us to identify which
  // component of the solution corresponds to which dof is based on the local
  // dof index and not the global one. Therefore, we need to loop over all the
  // cells and their local dofs to apply the correct scaling to each component
  // of the solution.
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Get the global dof indices for the current cell
          cell->get_dof_indices(local_dof_indices);

          // Loop over the local dofs of the cell
          for (unsigned int local_dof = 0; local_dof < local_dof_indices.size();
               ++local_dof)
            {
              // Get the global dof index corresponding to the local dof index
              const types::global_dof_index global_dof =
                local_dof_indices[local_dof];

              // Check if the global dof index is locally owned. If not, we skip
              // it since we only want to scale the locally owned part of the
              // solution vector.
              if (!locally_owned_dofs.is_element(global_dof))
                continue;

              const unsigned int locally_owned_position =
                locally_owned_dofs.index_within_set(global_dof);

              // Since we loop over all the cells, we will encounter each dof
              // multiple times (once for each cell it belongs to). We only want
              // to apply the scaling once for each dof, so we keep track of
              // which dofs have already been scaled using the dof_scaled
              // vector. If a dof has already been scaled, we skip it.
              if (dof_scaled[locally_owned_position])
                continue;

              dof_scaled[locally_owned_position] = true;

              // Get the base element of the current dof. Remember that the
              // convention for the FE system we are using is [E_real, E_imag,
              // H_real, H_imag].
              const unsigned int base_element =
                fe.system_to_base_index(local_dof).first.first;

              if (base_element == 0 || base_element == 1)
                solution[global_dof] *= electric_scale;
              else if (base_element == 2 || base_element == 3)
                solution[global_dof] *= magnetic_scale;
            }
        }
    }
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
TimeHarmonicMaxwell<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;

  // Interior output setup
  std::vector<std::string> solution_interior_names(dim, "E_real");
  for (int i = 0; i < dim; ++i)
    {
      solution_interior_names.emplace_back("E_imag");
    }
  for (int i = 0; i < dim; ++i)
    {
      solution_interior_names.emplace_back("H_real");
    }
  for (int i = 0; i < dim; ++i)
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

  // We also want to output the DPG estimator as a separate field for
  // visualization and postprocessing purposes when computed.
  if (this->simulation_parameters.mesh_adaptation.var_adaptation_param
        .error_estimator ==
      Parameters::MultipleAdaptationParameters::ErrorEstimator::dpg)
    {
      solution_output_structs.emplace_back(
        std::in_place_type<OutputStructCellVector>,
        this->local_estimated_error_per_cell,
        std::string("dpg_error_norm"));
    }

  // Skeleton output setup
  // TODO:  it will need its own writer probably because we need to use
  // DataOutFaces object, add the skeleton output to all physics when required?

  return solution_output_structs;
}

template <int dim>
std::vector<double>
TimeHarmonicMaxwell<dim>::calculate_L2_error()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  // Interior L2 error
  FEValues<dim> fe_values_trial_interior(*this->mapping,
                                         *this->fe_trial_interior,
                                         *this->cell_quadrature,
                                         update_values |
                                           update_quadrature_points |
                                           update_JxW_values);

  const unsigned int n_q_points = this->cell_quadrature->size();

  // The exact solution will be defined by the user but will need to be on all
  // possible fields of the ultraweak formulation so we need 4*dim components.
  std::vector<Vector<double>> exact_solution_values(n_q_points,
                                                    Vector<double>(4 * dim));
  auto                       &exact_solution_function =
    simulation_parameters.analytical_solution->electromagnetics;

  // When looping on each cell we will extract the different field
  // solution obtained numerically. The containers used to store the
  // interpolated solution at the quadrature points are declared below.
  std::vector<Tensor<1, dim>> local_E_values_real(n_q_points);
  std::vector<Tensor<1, dim>> local_E_values_imag(n_q_points);
  std::vector<Tensor<1, dim>> local_H_values_real(n_q_points);
  std::vector<Tensor<1, dim>> local_H_values_imag(n_q_points);

  // We create variables that will store all the integration results we are
  // interested in.
  double L2_error_E_real = 0;
  double L2_error_E_imag = 0;
  double L2_error_H_real = 0;
  double L2_error_H_imag = 0;

  for (const auto &cell : dof_handler_trial_interior->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_trial_interior.reinit(cell);

          // Get the simulated solution at quadrature points
          fe_values_trial_interior[extractor_E_real].get_function_values(
            *present_solution, local_E_values_real);
          fe_values_trial_interior[extractor_E_imag].get_function_values(
            *present_solution, local_E_values_imag);
          fe_values_trial_interior[extractor_H_real].get_function_values(
            *present_solution, local_H_values_real);
          fe_values_trial_interior[extractor_H_imag].get_function_values(
            *present_solution, local_H_values_imag);

          // Get the exact solution at quadrature points
          exact_solution_function.vector_value_list(
            fe_values_trial_interior.get_quadrature_points(),
            exact_solution_values);

          // Loop on quadrature points to compute the L2 error contributions
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double JxW = fe_values_trial_interior.JxW(q);

              // Loop on dimensions to compute the squared error
              for (int d = 0; d < dim; ++d)
                {
                  // E real part
                  L2_error_E_real +=
                    Utilities::fixed_power<2>(local_E_values_real[q][d] -
                                              exact_solution_values[q][d]) *
                    JxW;

                  // E imag part
                  L2_error_E_imag += Utilities::fixed_power<2>(
                                       local_E_values_imag[q][d] -
                                       exact_solution_values[q][d + dim]) *
                                     JxW;

                  // H real part
                  L2_error_H_real += Utilities::fixed_power<2>(
                                       local_H_values_real[q][d] -
                                       exact_solution_values[q][d + 2 * dim]) *
                                     JxW;

                  // H imag part
                  L2_error_H_imag += Utilities::fixed_power<2>(
                                       local_H_values_imag[q][d] -
                                       exact_solution_values[q][d + 3 * dim]) *
                                     JxW;
                }
            }
        }
    }
  // Skeleton L2 error
  // TODO

  L2_error_E_real = Utilities::MPI::sum(L2_error_E_real, mpi_communicator);
  L2_error_E_imag = Utilities::MPI::sum(L2_error_E_imag, mpi_communicator);
  L2_error_H_real = Utilities::MPI::sum(L2_error_H_real, mpi_communicator);
  L2_error_H_imag = Utilities::MPI::sum(L2_error_H_imag, mpi_communicator);

  return {std::sqrt(L2_error_E_real),
          std::sqrt(L2_error_E_imag),
          std::sqrt(L2_error_H_real),
          std::sqrt(L2_error_H_imag)};
}


template <int dim>
void
TimeHarmonicMaxwell<dim>::finish_simulation()
{
  auto         mpi_communicator = this->triangulation->get_mpi_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (this_mpi_process == 0 &&
      simulation_parameters.analytical_solution->verbosity !=
        Parameters::Verbosity::quiet)
    {
      ConvergenceTable &error_table = this->error_table;

      error_table.omit_column_from_convergence_rate_evaluation("cells");

      error_table.evaluate_all_convergence_rates(
        ConvergenceTable::reduction_rate_log2);

      error_table.set_scientific("error_E_real", true);
      error_table.set_scientific("error_E_imag", true);
      error_table.set_scientific("error_H_real", true);
      error_table.set_scientific("error_H_imag", true);
      error_table.set_precision("error_E_real",
                                this->simulation_control->get_log_precision());
      error_table.set_precision("error_E_imag",
                                this->simulation_control->get_log_precision());
      error_table.set_precision("error_H_real",
                                this->simulation_control->get_log_precision());
      error_table.set_precision("error_H_imag",
                                this->simulation_control->get_log_precision());
      error_table.write_text(std::cout);
    }

  if (this->simulation_parameters.timer.type == Parameters::Timer::Type::end)
    {
      announce_string(this->pcout, "Time Harmonic Electromagnetics");
      this->pcout << std::defaultfloat;
      this->computing_timer.print_summary();
      this->pcout << std::scientific;
    }
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::update_material_properties_dependencies()
{
  if (needs_temperature)
    {
      this->temperature_last_solved_solution =
        this->multiphysics->get_solution(PhysicsID::heat_transfer);
    }
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
  // No modification of the solution is required at the moment
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::update_boundary_conditions()
{
  if (!this->simulation_parameters
         .boundary_conditions_time_harmonic_electromagnetics.time_dependent)
    return;

  AssertThrow(
    false,
    ExcMessage(
      "Time-dependent boundary conditions not yet implemented for TimeHarmonicMaxwell."));
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::postprocess(bool first_iteration)
{
  if (simulation_parameters.analytical_solution->calculate_error() == true &&
      !first_iteration)
    {
      std::vector<double> errors       = calculate_L2_error();
      double              E_real_error = errors[0];
      double              E_imag_error = errors[1];
      double              H_real_error = errors[2];
      double              H_imag_error = errors[3];

      this->error_table.add_value("cells",
                                  this->triangulation->n_global_active_cells());
      this->error_table.add_value("error_E_real", E_real_error);
      this->error_table.add_value("error_E_imag", E_imag_error);
      this->error_table.add_value("error_H_real", H_real_error);
      this->error_table.add_value("error_H_imag", H_imag_error);

      if (simulation_parameters.analytical_solution->verbosity !=
          Parameters::Verbosity::quiet)
        {
          this->pcout << "L2 error E real: " << E_real_error << std::endl;
          this->pcout << "L2 error E imag: " << E_imag_error << std::endl;
          this->pcout << "L2 error H real: " << H_real_error << std::endl;
          this->pcout << "L2 error H imag: " << H_imag_error << std::endl;
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::output_per_iteration_timer()
{
  announce_string(this->pcout, "Time Harmonic Electromagnetics");
  this->pcout << std::defaultfloat;
  this->computing_timer.print_summary();
  this->pcout << std::scientific;
  this->computing_timer.reset();
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::pre_mesh_adaptation()
{
  this->solution_transfer->prepare_for_coarsening_and_refinement(
    *this->present_solution);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::post_mesh_adaptation()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  // Set up the vectors for the transfer
  GlobalVectorType tmp(this->locally_owned_dofs_trial_interior,
                       mpi_communicator);

  // Interpolate the solution at time and previous time
  this->solution_transfer->interpolate(tmp);

  // Distribute constraints does not need to be call for the DPG method transfer
  // as they live on the skeleton only, but we transfer the interior solution
  // only for multiphysics purposes.

  // Fix on the new mesh
  *this->present_solution = tmp;
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::write_checkpoint()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  // Checkpoint interior
  std::vector<const GlobalVectorType *> sol_set_transfer;
  solution_transfer = std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(
    *dof_handler_trial_interior);

  sol_set_transfer.emplace_back(&(*present_solution));

  solution_transfer->prepare_for_serialization(sol_set_transfer);

  // When temperature-dependent electromagnetic properties are used by the
  // simulation, we need to checkpoint the last temperature solution vector
  if (needs_temperature)
    {
      std::vector<const GlobalVectorType *> temperature_sol_set_transfer;

      temperature_last_solved_solution_transfer =
        std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(
          this->multiphysics->get_dof_handler(PhysicsID::heat_transfer));

      temperature_sol_set_transfer.emplace_back(
        &temperature_last_solved_solution);

      temperature_last_solved_solution_transfer->prepare_for_serialization(
        temperature_sol_set_transfer);
    }

  // Serialize all post-processing tables that are currently used with the
  // TimeHarmonicMaxwell solver
  const std::vector<OutputStructTableHandler> &table_output_structs =
    this->gather_tables();
  serialize_tables_vector(table_output_structs, mpi_communicator);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::read_checkpoint()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();
  this->pcout << "Reading time-harmonic Maxwell checkpoint" << std::endl;

  // Time-harmonic Maxwell solution interior
  std::vector<GlobalVectorType *> input_vectors(1);
  GlobalVectorType distributed_system(locally_owned_dofs_trial_interior,
                                      mpi_communicator);
  input_vectors[0] = &distributed_system;

  solution_transfer = std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(
    *dof_handler_trial_interior);
  solution_transfer->deserialize(input_vectors);

  *present_solution = distributed_system;

  // Temperature solution
  if (needs_temperature)
    {
      const auto &temperature_dof_handler =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);

      const IndexSet locally_relevant_temperature_dofs =
        DoFTools::extract_locally_relevant_dofs(temperature_dof_handler);
      temperature_last_solved_solution.reinit(
        temperature_dof_handler.locally_owned_dofs(),
        locally_relevant_temperature_dofs,
        mpi_communicator);

      GlobalVectorType distributed_temperature_last_solved(
        temperature_dof_handler.locally_owned_dofs(), mpi_communicator);

      temperature_last_solved_solution_transfer =
        std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(
          temperature_dof_handler);

      std::vector<GlobalVectorType *> temperature_input_vectors(1);

      temperature_input_vectors[0] = &distributed_temperature_last_solved;

      temperature_last_solved_solution_transfer->deserialize(
        temperature_input_vectors);

      temperature_last_solved_solution = distributed_temperature_last_solved;
    }

  // Deserialize all post-processing tables that are currently used with the
  // TimeHarmonicMaxwell solver
  std::vector<OutputStructTableHandler> table_output_structs =
    this->gather_tables();
  deserialize_tables_vector(table_output_structs, mpi_communicator);
}

template <int dim>
std::vector<OutputStructTableHandler>
TimeHarmonicMaxwell<dim>::gather_tables()
{
  std::vector<OutputStructTableHandler> table_output_structs;

  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder;
  std::string suffix = ".checkpoint";

  if (this->simulation_parameters.analytical_solution->calculate_error())
    table_output_structs.emplace_back(
      this->error_table,
      prefix + this->simulation_parameters.analytical_solution->get_filename() +
        "_THM" + suffix);


  return table_output_structs;
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::compute_error_estimate(
  const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                        &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  if (ivar.first == Variable::electric_field)
    {
      AssertThrow(
        ivar.second.error_estimator ==
          Parameters::MultipleAdaptationParameters::ErrorEstimator::kelly,
        ExcMessage(
          "Only the Kelly error estimator is currently implemented for the "
          "<electric field> field."));

      ComponentMask electric_field_mask =
        this->fe_trial_interior->component_mask(extractor_E_real) |
        this->fe_trial_interior->component_mask(extractor_E_imag);
      compute_kelly(estimated_error_per_cell, electric_field_mask);
    }
  else if (ivar.first == Variable::magnetic_field)
    {
      AssertThrow(
        ivar.second.error_estimator ==
          Parameters::MultipleAdaptationParameters::ErrorEstimator::kelly,
        ExcMessage(
          "Only the Kelly error estimator is currently implemented for the "
          "<magnetic field> field."));

      ComponentMask magnetic_field_mask =
        this->fe_trial_interior->component_mask(extractor_H_real) |
        this->fe_trial_interior->component_mask(extractor_H_imag);
      compute_kelly(estimated_error_per_cell, magnetic_field_mask);
    }
  else if (ivar.first == Variable::electromagnetic_fields)
    {
      if (ivar.second.error_estimator ==
          Parameters::MultipleAdaptationParameters::ErrorEstimator::kelly)
        {
          ComponentMask electromagnetics_mask =
            this->fe_trial_interior->component_mask(extractor_E_real) |
            this->fe_trial_interior->component_mask(extractor_E_imag) |
            this->fe_trial_interior->component_mask(extractor_H_real) |
            this->fe_trial_interior->component_mask(extractor_H_imag);
          compute_kelly(estimated_error_per_cell, electromagnetics_mask);
        }

      else if (ivar.second.error_estimator ==
               Parameters::MultipleAdaptationParameters::ErrorEstimator::dpg)
        {
          compute_dpg_error(estimated_error_per_cell);
        }
      else
        {
          AssertThrow(
            false,
            ExcMessage(
              "Unknown error estimator type for variable <electromagnetic fields>."));
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::compute_kelly(
  dealii::Vector<float> &estimated_error_per_cell,
  const ComponentMask   &component_mask)
{
  KellyErrorEstimator<dim>::estimate(
    *this->mapping,
    *this->dof_handler_trial_interior,
    *this->face_quadrature,
    typename std::map<types::boundary_id, const Function<dim, double> *>(),
    *this->present_solution,
    estimated_error_per_cell,
    component_mask);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::compute_dpg_error(
  dealii::Vector<float> &estimated_error_per_cell)
{
  // For efficiency, the DPG error estimator is computed in the same loop as the
  // reconstruction of the interior solution, so we do not implement it here as
  // a separate loop. The estimated error per cell is computed and stored in the
  // member variable local_estimated_error_per_cell during the reconstruction if
  // the DPG error estimator is activated, and then used for marking the cells
  // for refinement.
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  // Here we add a flag to be sure if the local_estimated_error_per_cell has
  // been computed before calling this function. If not, all entry of the
  // local_estimated_error_per_cell vector will still be 0.
  bool dpg_error_computed = false;

  for (const auto &cell :
       this->dof_handler_trial_interior->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          if (this->local_estimated_error_per_cell[cell->active_cell_index()] !=
              0)
            dpg_error_computed = true;

          const unsigned int cell_index = cell->active_cell_index();
          estimated_error_per_cell[cell_index] =
            this->local_estimated_error_per_cell[cell_index];
        }
    }

  // Reduce the flag across all MPI ranks so that processes that don't own
  // any cells don't falsely triggers the assertion.
  const bool is_global_dpg_error_computed =
    Utilities::MPI::max(static_cast<unsigned int>(dpg_error_computed),
                        mpi_communicator) != 0;

  AssertThrow(
    is_global_dpg_error_computed,
    ExcMessage(
      "The DPG error estimator has not been computed before calling compute_dpg_error. Please make sure that the system assembly has been performed and the local_estimated_error_per_cell vector has been filled with the DPG error estimates for each cell before calling this function."));
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::setup_dofs()
{
  verify_consistency_of_boundary_conditions();

  auto mpi_communicator = triangulation->get_mpi_communicator();

  // Setup each dof handlers
  this->dof_handler_trial_interior->distribute_dofs(*this->fe_trial_interior);
  DoFRenumbering::Cuthill_McKee(*this->dof_handler_trial_interior);
  this->dof_handler_trial_skeleton->distribute_dofs(*this->fe_trial_skeleton);
  DoFRenumbering::Cuthill_McKee(*this->dof_handler_trial_skeleton);
  this->dof_handler_test->distribute_dofs(*this->fe_test);
  DoFRenumbering::Cuthill_McKee(*this->dof_handler_test);

  // Get the locally owned dofs
  this->locally_owned_dofs_trial_interior =
    this->dof_handler_trial_interior->locally_owned_dofs();
  this->locally_owned_dofs_trial_skeleton =
    this->dof_handler_trial_skeleton->locally_owned_dofs();

  // Get the locally relevant dofs
  this->locally_relevant_dofs_trial_interior =
    DoFTools::extract_locally_relevant_dofs(*this->dof_handler_trial_interior);
  this->locally_relevant_dofs_trial_skeleton =
    DoFTools::extract_locally_relevant_dofs(*this->dof_handler_trial_skeleton);

  // Initialize the solution vectors and the error estimate per cell
  this->present_solution->reinit(this->locally_owned_dofs_trial_interior,
                                 this->locally_relevant_dofs_trial_interior,
                                 mpi_communicator);
  this->present_solution_skeleton->reinit(
    this->locally_owned_dofs_trial_skeleton,
    this->locally_relevant_dofs_trial_skeleton,
    mpi_communicator);
  this->local_estimated_error_per_cell.reinit(triangulation->n_active_cells());

  // We reinitialize the system rhs with the skeleton dofs because we have
  // performed a static condensation of the interior dofs using the Schur
  // complement.
  this->system_rhs.reinit(this->locally_owned_dofs_trial_skeleton,
                          mpi_communicator);

  // We reinitialize the initial guess vector
  this->initial_guess_iterative_solver.reinit(
    this->locally_owned_dofs_trial_skeleton, mpi_communicator);

  // Define constraints
  define_constraints();

  // Sparse matrices initialization
  // In DPG, the sparse matrix and the dynamic sparsity pattern are very
  // expensive so we recast the dynamic sparsity pattern to a
  // sparsity pattern before initializing the system matrix to save memory
  // because the dynamic sparsity pattern is more expensive in terms of
  // memory consumption than the static sparsity pattern.
  TrilinosWrappers::SparsityPattern
    sparsity_pattern; // This needs to be defined outside the following block
                      // because it is used in extra_verbose to report the
                      // memory consumption of the sparsity pattern.

  {
    DynamicSparsityPattern dsp(this->locally_relevant_dofs_trial_skeleton);
    DoFTools::make_sparsity_pattern(*this->dof_handler_trial_skeleton,
                                    dsp,
                                    this->nonzero_constraints,
                                    /*keep_constrained_dofs = */ false);
    SparsityTools::distribute_sparsity_pattern(
      dsp,
      this->locally_owned_dofs_trial_skeleton,
      mpi_communicator,
      this->locally_relevant_dofs_trial_skeleton);

    sparsity_pattern.reinit(this->locally_owned_dofs_trial_skeleton,
                            this->locally_owned_dofs_trial_skeleton,
                            dsp,
                            mpi_communicator);
    sparsity_pattern.compress();
  }

  this->system_matrix.reinit(sparsity_pattern);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
        .verbosity == Parameters::Verbosity::extra_verbose)
    {
      print_THM_setup_memory(sparsity_pattern);
    }
  this->pcout << "  DPG system for Time-Harmonic Maxwell Equations:"
              << std::endl;
  this->pcout
    << "   Number of skeleton degrees of freedom for Time-Harmonic Maxwell: "
    << this->dof_handler_trial_skeleton->n_dofs() << std::endl;
  this->pcout
    << "   Number of interior degrees of freedom for Time-Harmonic Maxwell: "
    << this->dof_handler_trial_interior->n_dofs() << std::endl;

  // Update the multiphysics interface with the TimeHarmonicMaxwell dof_handler,
  // mapping, and solution. Note that we provide the interior dof_handler and
  // interior solution as it is the one useful for the other physics (we do not
  // provide the skeleton dof_handler and skeleton solution).
  multiphysics->set_dof_handler(PhysicsID::electromagnetics,
                                this->dof_handler_trial_interior);
  multiphysics->set_mapping(PhysicsID::electromagnetics, this->mapping);
  multiphysics->set_solution(PhysicsID::electromagnetics,
                             this->present_solution);
}


template <int dim>
void
TimeHarmonicMaxwell<dim>::set_initial_conditions()
{
  // This tmp vector is used instead of the newton update vector as they don't
  // exist for this solver.
  GlobalVectorType tmp(this->locally_owned_dofs_trial_interior,
                       this->triangulation->get_mpi_communicator());

  VectorTools::interpolate(
    *this->mapping,
    *this->dof_handler_trial_interior,
    simulation_parameters.initial_condition->electromagnetics,
    tmp,
    fe_trial_interior->component_mask(extractor_E_real));

  VectorTools::interpolate(
    *this->mapping,
    *this->dof_handler_trial_interior,
    simulation_parameters.initial_condition->electromagnetics,
    tmp,
    fe_trial_interior->component_mask(extractor_E_imag));

  VectorTools::interpolate(
    *this->mapping,
    *this->dof_handler_trial_interior,
    simulation_parameters.initial_condition->electromagnetics,
    tmp,
    fe_trial_interior->component_mask(extractor_H_real));

  VectorTools::interpolate(
    *this->mapping,
    *this->dof_handler_trial_interior,
    simulation_parameters.initial_condition->electromagnetics,
    tmp,
    fe_trial_interior->component_mask(extractor_H_imag));

  // Note that we don't apply the constraints as the initial condition only
  // gives values on the interior dofs and not on the skeleton dofs.

  // No percolation of time vectors is needed as this is always a steady-state
  // solver.

  *this->present_solution = tmp;
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::setup_preconditioner()
{
  preconditioner = std::make_shared<TrilinosWrappers::PreconditionIdentity>();
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::define_constraints()
{
  // Clear previous constraints
  this->nonzero_constraints.clear();
  this->nonzero_constraints.reinit(this->locally_owned_dofs_trial_skeleton,
                                   this->locally_relevant_dofs_trial_skeleton);

  DoFTools::make_hanging_node_constraints(*this->dof_handler_trial_skeleton,
                                          this->nonzero_constraints);

  // Loop over all defined boundary conditions
  for (const auto &[id, type] :
       this->simulation_parameters
         .boundary_conditions_time_harmonic_electromagnetics.type)
    {
      if (type == BoundaryConditions::BoundaryType::pec)
        {
          // Perfect electric conductor (PEC) boundary condition
          // Real
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            0,
            Functions::ZeroFunction<dim>(4 * dim),
            id,
            this->nonzero_constraints);

          // Imaginary
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            dim,
            Functions::ZeroFunction<dim>(4 * dim),
            id,
            this->nonzero_constraints);
        }
      if (type == BoundaryConditions::BoundaryType::pmc)
        {
          // Perfect magnetic conductor (PMC) boundary condition
          // Real
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            2 * dim,
            Functions::ZeroFunction<dim>(4 * dim),
            id,
            this->nonzero_constraints);

          // Imaginary
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            3 * dim,
            Functions::ZeroFunction<dim>(4 * dim),
            id,
            this->nonzero_constraints);
        }
      if (type == BoundaryConditions::BoundaryType::electric_field)
        {
          // Imposed electric field boundary condition

          // Real
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            0,
            TimeHarmonicMaxwellElectricFieldDefined<dim>(
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_x_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_y_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_z_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_x_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_y_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_z_imag),
            id,
            this->nonzero_constraints);

          // Imaginary
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            dim,
            TimeHarmonicMaxwellElectricFieldDefined<dim>(
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_x_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_y_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_z_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_x_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_y_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_z_imag),
            id,
            this->nonzero_constraints);
        }

      if (type == BoundaryConditions::BoundaryType::magnetic_field)
        {
          // Imposed magnetic field boundary condition

          // Real
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            2 * dim,
            TimeHarmonicMaxwellMagneticFieldDefined<dim>(
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_x_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_y_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_z_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_x_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_y_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_z_imag),
            id,
            this->nonzero_constraints);

          // Imaginary
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            3 * dim,
            TimeHarmonicMaxwellMagneticFieldDefined<dim>(
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_x_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_y_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_z_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_x_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_y_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_z_imag),
            id,
            this->nonzero_constraints);
        }
    }

  // The DPG method requires the use of skeleton elements and shape functions.
  // Because we want to use high-order Nedelec elements for the face, but dealii
  // does not support a FE_FaceNedelec, we hack our way around this by using
  // full cell elements, but freezing the interior dofs using the function
  // constrain_dof_to_zero. To do so, we create a container for all dof indices
  // and a container for the face dof indices and we loop on all the faces of
  // each cell to find which dofs are on the face and flag them. The remaining
  // dofs are then constrained to zero.

  std::vector<types::global_dof_index> cell_dof_indices(
    this->fe_trial_skeleton->n_dofs_per_cell());
  std::vector<types::global_dof_index> face_dof_indices(
    this->fe_trial_skeleton->n_dofs_per_face());
  std::vector<bool> is_dof_on_face(this->fe_trial_skeleton->n_dofs_per_cell(),
                                   false);

  // Loop on all skeleton dofs and set the interior constraints to zero
  for (const auto &cell :
       this->dof_handler_trial_skeleton->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Get all dof indices on the cell
          cell->get_dof_indices(cell_dof_indices);

          // Loop on all the faces of the cell
          for (const auto &face : cell->face_iterators())
            {
              face->get_dof_indices(face_dof_indices);

              // Loop on all dofs on the face
              for (const auto &face_dof : face_dof_indices)
                {
                  // Find the first iterator in the cell dof indices that
                  // matches the face dof (find where the face dof is in the
                  // cell dof indices)
                  const auto it = std::ranges::find(cell_dof_indices, face_dof);

                  // If the dof is on the face (find returns the second
                  // iterator if no match is found), set the corresponding
                  // flag to true
                  if (it != cell_dof_indices.end())
                    {
                      is_dof_on_face[std::distance(cell_dof_indices.begin(),
                                                   it)] = true;
                    }
                }
            }

          // Loop on all dofs on the cell and constrain the interior ones to
          // zero
          for (unsigned int index = 0; index < cell_dof_indices.size(); ++index)
            {
              // If the dof is not on a face, then it is an interior dof
              if (!is_dof_on_face[index])
                {
                  this->nonzero_constraints.constrain_dof_to_zero(
                    cell_dof_indices[index]);
                }
            }
        }
    }
  this->nonzero_constraints.make_consistent_in_parallel(
    this->locally_owned_dofs_trial_skeleton,
    this->locally_relevant_dofs_trial_skeleton,
    this->triangulation->get_mpi_communicator());
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::solve_linear_system()
{
  TimerOutput::Scope t(this->computing_timer, "Solve linear system");

  AssertThrow(
    simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::none,
    ExcMessage(
      "The time harmonic electromagnetism physics does not support a preconditioner. Please set: subsection linear solver -> subsection electromagnetics -> set preconditioner = none."));

  // Define the linear solver tolerance
  const double absolute_residual =
    simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
      .minimum_residual;
  const double relative_residual =
    simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
      .relative_residual;

  const double rescale_metric    = this->get_residual_rescale_metric();
  const double rescaled_residual = this->system_rhs.l2_norm() / rescale_metric;
  const double linear_solver_tolerance =
    std::max(relative_residual * rescaled_residual, absolute_residual);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  // Set up the solver control
  SolverControl solver_control(this->dof_handler_trial_skeleton->n_dofs(),
                               linear_solver_tolerance);

  // Solve
  TrilinosWrappers::SolverCG solver(solver_control);

  solver.solve(this->system_matrix,
               this->initial_guess_iterative_solver,
               this->system_rhs,
               *this->preconditioner);

  if (simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -CG iterative solver took : "
                  << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() / rescale_metric << std::endl;
    }

  // Update the solution vector for the skeleton. We keep the solution in the
  // initial_guess_iterative_solver vector as it is used for the next iteration
  // of the solver and the present_solution_skeleton vector is used for
  // everything else (output, multiphysics coupling, etc.). The major difference
  // is that the present_solution_skeleton vector may be renormalized which
  // would break the iterative solver if it was used for the next iteration.
  nonzero_constraints.distribute(this->initial_guess_iterative_solver);
  *this->present_solution_skeleton = this->initial_guess_iterative_solver;

  // Reconstruct the interior solution from the skeleton solution
  reconstruct_interior_solution();

  // We need to apply the scaling to the interior solution as it is the one
  // used for the multiphysics coupling and output.
  scale_solution_components(*this->dof_handler_trial_interior,
                            *this->present_solution);

  // We also apply the scaling to the skeleton solution for
  // consistency, even if it is not used for the multiphysics coupling nor
  // outputted at the moment.
  scale_solution_components(*this->dof_handler_trial_skeleton,
                            *this->present_solution_skeleton);

  update_material_properties_dependencies();
}

template <int dim>
bool
TimeHarmonicMaxwell<dim>::should_solve_auxiliary_physics()
{
  const Parameters::TimeHarmonicMaxwell<dim> &thm_parameters =
    this->simulation_parameters.multiphysics.time_harmonic_maxwell_parameters;
  // For steady simulations, each outer iteration (e.g., mesh adaptation
  // cycles) should solve electromagnetics regardless of time-coupling settings.
  if (this->simulation_control->get_assembly_method() ==
      Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      return true;
    }

  // Always solve at the first step of the simulation (simulation start as 0
  // when the set up is performed, and then it is incremented before solving the
  // physics for the first time, so we want this function to return true in both
  // cases.
  if (this->simulation_control->get_iteration_number() <= 1)
    {
      return true;
    }
  else
    {
      switch (thm_parameters.time_coupling_strategy)
        {
          case Parameters::TimeHarmonicMaxwellCouplingStrategy::none:
            return false;


          case Parameters::TimeHarmonicMaxwellCouplingStrategy::iteration:
            // Solve only if we are at a multiple of the specified iteration
            // frequency. We subtract 1 since the step number starts at 1.
            if ((this->simulation_control->get_iteration_number() - 1) %
                  thm_parameters.coupling_iteration ==
                0)
              return true;
            else
              return false;
          case Parameters::TimeHarmonicMaxwellCouplingStrategy::time:
            {
              // Solve only if the current time has passed a multiple of
              // the specified time frequency since the last time the
              // electromagnetics were solved. This is done by comparing the
              // floor of the current time divided by the time coupling
              // parameter to the floor of the previous time divided by the time
              // coupling parameter. If they are different, it means we have
              // passed a multiple of the time coupling parameter and we should
              // solve the electromagnetics.
              const int difference_in_time_steps = static_cast<int>(
                std::floor(this->simulation_control->get_current_time() /
                           thm_parameters.coupling_time) -
                std::floor(this->simulation_control->get_previous_time() /
                           thm_parameters.coupling_time));
              if (difference_in_time_steps == 1)
                return true;
              else if (difference_in_time_steps > 1)
                {
                  AssertThrow(
                    false,
                    ExcMessage(
                      "The time coupling strategy for the time-harmonic Maxwell solver is set to 'time' with a coupling time of " +
                      std::to_string(thm_parameters.coupling_time) +
                      ", but the simulation time has advanced by " +
                      std::to_string(
                        this->simulation_control->get_current_time() -
                        this->simulation_control->get_previous_time()) +
                      " since the last time the electromagnetics were solved. Please either reduce the simulation time step or increase the coupling time to avoid missing coupling points."));
                  return true;
                }
              else
                return false;
            }
          case Parameters::TimeHarmonicMaxwellCouplingStrategy::threshold:
            {
              auto mpi_communicator = triangulation->get_mpi_communicator();
              // Initialize the temporary variables to determine if the physical
              // properties have changed by more than the specified threshold
              // since the last time the electromagnetics were solved.
              unsigned int                     material_id;
              const PhysicalPropertiesManager &physical_properties_manager =
                this->simulation_parameters.physical_properties_manager;
              std::vector<std::complex<double>>
                effective_electric_permittivities_last_solved;
              std::vector<std::complex<double>>
                effective_magnetic_permeabilities_last_solved;
              std::map<field, std::vector<double>>
                                  field_values_vector_last_solved;
              std::vector<double> temperature_last_solved_values;
              std::vector<std::complex<double>>
                effective_electric_permittivities_current;
              std::vector<std::complex<double>>
                effective_magnetic_permeabilities_current;
              std::map<field, std::vector<double>> field_values_vector_current;
              std::vector<double>                  temperature_current_values;
              double                               max_relative_change = 0.0;
              unsigned int                         n_q_points          = 0;

              // The temperature field is only needed when at least one
              // electromagnetic property depends on temperature.
              std::unique_ptr<FEValues<dim>> fe_values_temperature;
              const DoFHandler<dim>         *dof_handler_temperature = nullptr;
              const GlobalVectorType *temperature_current_solution   = nullptr;


              if (needs_temperature)
                {
                  dof_handler_temperature =
                    &this->multiphysics->get_dof_handler(
                      PhysicsID::heat_transfer);
                  const QGauss<dim> temperature_quadrature(
                    dof_handler_temperature->get_fe().degree + 1);
                  n_q_points = temperature_quadrature.size();
                  temperature_current_solution =
                    &this->multiphysics->get_solution(PhysicsID::heat_transfer);

                  fe_values_temperature = std::make_unique<FEValues<dim>>(
                    *this->mapping,
                    dof_handler_temperature->get_fe(),
                    temperature_quadrature,
                    update_values | update_quadrature_points);

                  temperature_last_solved_values.resize(n_q_points);
                  temperature_current_values.resize(n_q_points);

                  // Loop over all cell quadrature points and check if the
                  // physical properties have changed by more than the specified
                  // threshold since the last time the electromagnetics were
                  // solved according to the temperature field. If so, we need
                  // to solve the electromagnetics again.
                  for (const auto &cell :
                       dof_handler_temperature->active_cell_iterators())
                    {
                      if (cell->is_locally_owned())
                        {
                          fe_values_temperature->reinit(cell);
                          fe_values_temperature->get_function_values(
                            *temperature_current_solution,
                            temperature_current_values);

                          fe_values_temperature->get_function_values(
                            this->temperature_last_solved_solution,
                            temperature_last_solved_values);

                          field_values_vector_last_solved[field::temperature] =
                            temperature_last_solved_values;
                          field_values_vector_current[field::temperature] =
                            temperature_current_values;

                          material_id = cell->material_id();

                          compute_effective_electromagnetic_properties(
                            physical_properties_manager,
                            field_values_vector_last_solved,
                            material_id,
                            effective_electric_permittivities_last_solved,
                            effective_magnetic_permeabilities_last_solved);

                          compute_effective_electromagnetic_properties(
                            physical_properties_manager,
                            field_values_vector_current,
                            material_id,
                            effective_electric_permittivities_current,
                            effective_magnetic_permeabilities_current);

                          for (unsigned int q = 0; q < n_q_points; ++q)
                            {
                              // Here we check if the magnitude of the complex
                              // property difference normalized by the magnitude
                              // of the last solved property is greater than the
                              // threshold. We start by computing only the
                              // absolute square to do the comparison, and then
                              // we perform the square root to get the actual
                              // threshold value to have a faster computation.
                              double relative_change_electric_permittivity =
                                std::norm(
                                  effective_electric_permittivities_current[q] -
                                  effective_electric_permittivities_last_solved
                                    [q]) /
                                std::norm(
                                  effective_electric_permittivities_last_solved
                                    [q]);
                              double relative_change_magnetic_permeability =
                                std::norm(
                                  effective_magnetic_permeabilities_current[q] -
                                  effective_magnetic_permeabilities_last_solved
                                    [q]) /
                                std::norm(
                                  effective_magnetic_permeabilities_last_solved
                                    [q]);
                              max_relative_change =
                                std::max({relative_change_electric_permittivity,
                                          relative_change_magnetic_permeability,
                                          max_relative_change});
                            }
                        }
                    }
                }



              // Reduce the maximum relative change across all MPI ranks
              max_relative_change =
                Utilities::MPI::max(max_relative_change, mpi_communicator);

              if (simulation_parameters.linear_solver
                    .at(PhysicsID::electromagnetics)
                    .verbosity != Parameters::Verbosity::quiet)
                {
                  this->pcout
                    << "  - The maximum relative change in the electromagnetics physical properties since the last time the time-harmonic Maxwell equations were solved is "
                    << std::sqrt(max_relative_change) << std::endl;
                }

              // Return the threshold comparison result. We take the square
              // of the coupling_threshold to get the actual relative change,
              // because max_relative_change is already square and it is more
              // efficient to elevate to the power 2 rather than computing a
              // square root.
              if (max_relative_change >
                  Utilities::fixed_power<2>(thm_parameters.coupling_threshold))
                {
                  return true;
                }
              else
                {
                  return false;
                }
            }


          default:
            AssertThrow(false, ExcMessage("Unknown time coupling strategy."));
            return false;
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::setup_assemblers()
{
  this->assemblers.clear();
  this->face_assemblers.clear();

  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters =
    this->simulation_parameters.multiphysics.time_harmonic_maxwell_parameters;
  const BoundaryConditions::TimeHarmonicMaxwellBoundaryConditions<dim> &
    boundary_conditions = this->simulation_parameters
                            .boundary_conditions_time_harmonic_electromagnetics;

  // Cell terms: Gram matrix, interior bilinear form and imposed current
  // density.
  this->assemblers.emplace_back(
    std::make_shared<TimeHarmonicMaxwellAssemblerCore<dim>>(
      time_harmonic_maxwell_parameters));

  // Skeleton terms, which are present on every face of the mesh.
  this->face_assemblers.emplace_back(
    std::make_shared<TimeHarmonicMaxwellAssemblerSkeleton<dim>>(
      boundary_conditions));

  // Robin boundary conditions. The waveguide port amplitudes must have been
  // computed by compute_electromagnetic_scaling before this call.
  bool has_robin_boundary = false;
  for (const auto &bc : boundary_conditions.type)
    if (is_robin_boundary_type(bc.second))
      has_robin_boundary = true;
  if (has_robin_boundary)
    this->face_assemblers.emplace_back(
      std::make_shared<TimeHarmonicMaxwellAssemblerRobinBC<dim>>(
        time_harmonic_maxwell_parameters,
        boundary_conditions,
        this->waveguide_ports_electric_amplitudes));
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::assemble_system_matrix()
{
  // We always need to compute the scaling for the electromagnetic fields even
  // if the user does not want to apply it to the solution because it is used to
  // normalize the waveguide inlets relative input power
  compute_electromagnetic_scaling(
    this->simulation_parameters.physical_properties_manager);

  TimerOutput::Scope t(this->computing_timer, "Assemble matrix and RHS");

  // The time-harmonic Maxwell system is assembled again every time it is
  // solved (i.e., at each coupling step), so the global system needs to be
  // reset before the assembly.
  this->system_matrix = 0;
  this->system_rhs    = 0;

  setup_assemblers();

  auto scratch_data = TimeHarmonicMaxwellScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe_trial_interior,
    *this->fe_trial_skeleton,
    *this->fe_test,
    *this->cell_quadrature,
    *this->face_quadrature,
    *this->mapping);

  // We may need the temperature field for the physical properties.
  if (needs_temperature)
    scratch_data.enable_temperature(
      this->multiphysics->get_dof_handler(PhysicsID::heat_transfer).get_fe(),
      *this->cell_quadrature,
      *this->face_quadrature,
      *this->mapping);

  // As it is standard, we loop over the cells of the triangulation. We use
  // the DoFHandler associated with the interior trial space for this loop.
  WorkStream::run(this->dof_handler_trial_interior->begin_active(),
                  this->dof_handler_trial_interior->end(),
                  *this,
                  &TimeHarmonicMaxwell::assemble_local_system_matrix,
                  &TimeHarmonicMaxwell::copy_local_matrix_to_global_matrix,
                  scratch_data,
                  DPGCopyData(this->fe_test->n_dofs_per_cell(),
                              this->fe_trial_interior->n_dofs_per_cell(),
                              this->fe_trial_skeleton->n_dofs_per_cell()));

  // After the loop over the cells, we finalize the assembly by compressing
  // the vectors because of the MPI parallelization.
  this->system_matrix.compress(VectorOperation::add);
  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::assemble_local_dpg_system(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  TimeHarmonicMaxwellScratchData<dim>                  &scratch_data,
  DPGCopyData                                          &copy_data)
{
  // We get the same cell for the test and skeleton trial spaces to make sure
  // that all the FEValues objects are reinitialized on the same physical cell.
  const typename DoFHandler<dim>::active_cell_iterator cell_test =
    cell->as_dof_handler_iterator(*this->dof_handler_test);
  const typename DoFHandler<dim>::active_cell_iterator cell_skeleton =
    cell->as_dof_handler_iterator(*this->dof_handler_trial_skeleton);

  scratch_data.reinit(cell, cell_test);

  // The temperature field is only evaluated at the quadrature points if the
  // physical properties of the cell material depend on it.
  typename DoFHandler<dim>::active_cell_iterator cell_temperature;
  const GlobalVectorType                        *temperature_solution = nullptr;
  if (scratch_data.cell_material_needs_temperature)
    {
      cell_temperature = cell->as_dof_handler_iterator(
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer));
      temperature_solution =
        &this->multiphysics->get_solution(PhysicsID::heat_transfer);
      scratch_data.reinit_temperature(cell_temperature, *temperature_solution);
    }
  scratch_data.calculate_physical_properties();

  // We reset the local matrices and vector where we are aggregating the
  // information of the current cell.
  copy_data.reset();

  // Cell contributions to the local DPG system
  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_matrix(scratch_data, copy_data);
      assembler->assemble_rhs(scratch_data, copy_data);
    }

  // Skeleton contributions to the local DPG system. All the faces of the cell
  // belong to the skeleton.
  for (const unsigned int face_no : cell->face_indices())
    {
      scratch_data.reinit_face(cell_skeleton, cell_test, face_no);
      if (scratch_data.cell_material_needs_temperature)
        scratch_data.reinit_face_temperature(cell_temperature,
                                             face_no,
                                             *temperature_solution);
      scratch_data.calculate_face_physical_properties();

      for (auto &assembler : this->face_assemblers)
        {
          assembler->assemble_matrix(scratch_data, copy_data);
          assembler->assemble_rhs(scratch_data, copy_data);
        }
    }

  // After having assembled all the matrices and vectors, we compute the
  // condensation operators that are required both for the assembly of the
  // skeleton system and for the reconstruction of the interior solution.

  // We only need the inverse of the Gram matrix $G$, so we invert it.
  copy_data.G_matrix.invert();

  // We construct $M_4 = B^\dagger G^{-1}$ with it.
  copy_data.B_matrix.Tmmult(copy_data.M4_matrix, copy_data.G_matrix);

  // Then using $M_4$ we compute the condensed matrices $M_1 = B^\dagger G^{-1}
  // B$ and $M_2 = B^\dagger G^{-1} \hat{B}$.
  copy_data.M4_matrix.mmult(copy_data.M1_matrix, copy_data.B_matrix);
  copy_data.M4_matrix.mmult(copy_data.M2_matrix, copy_data.B_hat_matrix);

  // Finally, as for the $G$ matrix, we invert the $M_1$ matrix.
  copy_data.M1_matrix.invert();
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  TimeHarmonicMaxwellScratchData<dim>                  &scratch_data,
  DPGCopyData                                          &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  assemble_local_dpg_system(cell, scratch_data, copy_data);

  // We construct $M_5 = \hat{B}^\dagger G^{-1}$ and the matrix
  // $M_3 = \hat{B}^\dagger G^{-1} \hat{B}$.
  copy_data.B_hat_matrix.Tmmult(copy_data.M5_matrix, copy_data.G_matrix);
  copy_data.M5_matrix.mmult(copy_data.M3_matrix, copy_data.B_hat_matrix);

  // Now, we have to compute the local matrix and the local RHS for the
  // condensed system.

  // The cell matrix is obtained with the formula $(M_3 -
  // M_2^\dagger M_1^{-1} M_2)$:
  copy_data.M2_matrix.Tmmult(copy_data.tmp_matrix_M2M1, copy_data.M1_matrix);
  copy_data.tmp_matrix_M2M1.mmult(copy_data.tmp_matrix_M2M1M2,
                                  copy_data.M2_matrix);
  copy_data.tmp_matrix_M2M1M2.add(-1.0, copy_data.M3_matrix);
  copy_data.tmp_matrix_M2M1M2 *= -1.0;
  // This line is used to convert the LAPACK matrix to a full matrix so we can
  // perform the distribution to the global system.
  copy_data.local_matrix = copy_data.tmp_matrix_M2M1M2;

  // Then we compute the cell RHS using $(M_5 - M_2^\dagger M_1^{-1} M_4)l$.
  copy_data.tmp_matrix_M2M1.mmult(copy_data.tmp_matrix_M2M1M4,
                                  copy_data.M4_matrix);
  copy_data.M5_matrix.add(-1.0, copy_data.tmp_matrix_M2M1M4);
  copy_data.M5_matrix.vmult(copy_data.local_rhs, copy_data.l_vector);

  // Get the local dof indices for the skeleton trial space to be able to
  // distribute the local matrix and RHS to the global system. Cannot use
  // cell->get_dof_indices() because the skeleton trial space is not the same as
  // the interior trial space so we need to use the as_dof_handler_iterator()
  // function to change the cell iterator to the skeleton trial space and then
  // get the dof indices from that.
  cell->as_dof_handler_iterator(*this->dof_handler_trial_skeleton)
    ->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::copy_local_matrix_to_global_matrix(
  const DPGCopyData &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  this->nonzero_constraints.distribute_local_to_global(
    copy_data.local_matrix,
    copy_data.local_rhs,
    copy_data.local_dof_indices,
    this->system_matrix,
    this->system_rhs);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::assemble_system_rhs()
{
  // Because of the static condensation, the condensed right-hand side
  // requires the full local DPG system. It is therefore assembled with the
  // matrix in assemble_system_matrix().
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::reconstruct_local_interior_solution(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  TimeHarmonicMaxwellScratchData<dim>                  &scratch_data,
  DPGCopyData                                          &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  assemble_local_dpg_system(cell, scratch_data, copy_data);

  // Now, we already have the solution on the skeleton and only need to
  // perform $u_h = M_1^{-1} (M_4 l - M_2 \hat{u}_h)$ on each cell. When this
  // is obtained, we can compute at the same time the error indicator
  // $\Psi = G^{-1}(l - B u_h - \hat{B}\hat{u}_h)$ if the dpg error estimator
  // is activated.

  // We first get the skeleton solution vector for this cell.
  cell->as_dof_handler_iterator(*this->dof_handler_trial_skeleton)
    ->get_dof_values(*this->present_solution_skeleton,
                     copy_data.local_skeleton_solution);

  // Then we do the matrix-vector products to obtain the interior unknowns.
  copy_data.M2_matrix.vmult(copy_data.tmp_vector_interior,
                            copy_data.local_skeleton_solution);
  copy_data.M4_matrix.vmult(copy_data.local_interior_rhs, copy_data.l_vector);
  copy_data.local_interior_rhs -= copy_data.tmp_vector_interior;
  copy_data.M1_matrix.vmult(copy_data.local_interior_solution,
                            copy_data.local_interior_rhs);

  // We can also compute the error indicator on this cell if the dpg error
  // estimator is activated. The residual R = l - B u_h - \hat{B}\hat{u}_h is
  // stored in l_vector, its Riesz representation Psi = G^{-1} R in
  // local_residual, and the squared energy norm of the residual is
  // ||R||^2_V = R^T G^-1 R = R^T Psi.
  if (this->simulation_parameters.mesh_adaptation.var_adaptation_param
        .error_estimator ==
      Parameters::MultipleAdaptationParameters::ErrorEstimator::dpg)
    {
      copy_data.B_matrix.vmult(copy_data.tmp_vector_error_indicator,
                               copy_data.local_interior_solution);
      copy_data.B_hat_matrix.vmult_add(copy_data.tmp_vector_error_indicator,
                                       copy_data.local_skeleton_solution);
      copy_data.l_vector -= copy_data.tmp_vector_error_indicator;
      copy_data.G_matrix.vmult(copy_data.local_residual, copy_data.l_vector);

      copy_data.local_residual_norm_squared =
        copy_data.l_vector * copy_data.local_residual;
    }
  copy_data.active_cell_index = cell->active_cell_index();

  cell->get_dof_indices(copy_data.local_dof_indices_trial_interior);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::copy_local_interior_solution_to_global(
  const DPGCopyData &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  // We map the cell interior solution to the global interior solution. Adding
  // the values at the dof indices of the cell is what
  // cell->distribute_local_to_global(local_vector, global_vector) does, which
  // cannot be used here since the copier does not have access to the cell.
  // Since the interior trial space is discontinuous, each dof belongs to a
  // single cell and adding the value is equivalent to setting it.
  this->locally_owned_solution_interior.add(
    copy_data.local_dof_indices_trial_interior,
    copy_data.local_interior_solution);

  // Store the error indicator of the cell if the dpg error estimator is
  // activated
  if (this->simulation_parameters.mesh_adaptation.var_adaptation_param
        .error_estimator ==
      Parameters::MultipleAdaptationParameters::ErrorEstimator::dpg)
    {
      this->local_estimated_error_per_cell(copy_data.active_cell_index) =
        std::sqrt(copy_data.local_residual_norm_squared);
      this->squared_residual_L2_norm += copy_data.local_residual_norm_squared;
    }
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::reconstruct_interior_solution()
{
  MPI_Comm mpi_communicator = this->triangulation->get_mpi_communicator();

  const bool dpg_error_estimator_enabled =
    this->simulation_parameters.mesh_adaptation.var_adaptation_param
      .error_estimator ==
    Parameters::MultipleAdaptationParameters::ErrorEstimator::dpg;

  // The interior solution is assembled cell by cell. However,
  // present_solution is a ghosted vector (it also stores the locally relevant
  // dofs owned by the other processes, which are needed by the other physics
  // and the output) and a ghosted vector can only be read. Therefore, the
  // contributions of the locally owned cells are first added to the following
  // vector, which only stores the locally owned dofs. After the loop on the
  // cells, compress(VectorOperation::add) finalizes it, and the assignment to
  // the ghosted vector updates its ghost values. This vector is only needed
  // during the reconstruction, so it is allocated here and released at the end
  // of the function to avoid keeping it in memory between the solves.
  this->locally_owned_solution_interior.reinit(
    this->locally_owned_dofs_trial_interior, mpi_communicator);

  // Squared L2 norm of the residual for the dpg error indicator, accumulated
  // over the locally owned cells
  this->squared_residual_L2_norm = 0.0;

  auto scratch_data = TimeHarmonicMaxwellScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe_trial_interior,
    *this->fe_trial_skeleton,
    *this->fe_test,
    *this->cell_quadrature,
    *this->face_quadrature,
    *this->mapping);

  if (needs_temperature)
    scratch_data.enable_temperature(
      this->multiphysics->get_dof_handler(PhysicsID::heat_transfer).get_fe(),
      *this->cell_quadrature,
      *this->face_quadrature,
      *this->mapping);


  WorkStream::run(this->dof_handler_trial_interior->begin_active(),
                  this->dof_handler_trial_interior->end(),
                  *this,
                  &TimeHarmonicMaxwell::reconstruct_local_interior_solution,
                  &TimeHarmonicMaxwell::copy_local_interior_solution_to_global,
                  scratch_data,
                  DPGCopyData(this->fe_test->n_dofs_per_cell(),
                              this->fe_trial_interior->n_dofs_per_cell(),
                              this->fe_trial_skeleton->n_dofs_per_cell()));

  // After the loop over the cells, we finalize the assembly by compressing
  // the vector because of the MPI parallelization.
  this->locally_owned_solution_interior.compress(VectorOperation::add);

  *this->present_solution = this->locally_owned_solution_interior;

  // The non-ghosted vector is not needed anymore, so we release its memory by
  // swapping it with an empty vector, which is deallocated at the end of the
  // scope. This is used instead of clear(), which is not available for all the
  // types of GlobalVectorType.
  {
    GlobalVectorType empty_vector;
    this->locally_owned_solution_interior.swap(empty_vector);
  }

  // We also output the global error indicator if the dpg error estimator is
  // activated and in verbose mode
  if (dpg_error_estimator_enabled &&
      (this->simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
         .verbosity != Parameters::Verbosity::quiet))
    {
      this->pcout << "   Time-Harmonic Maxwell DPG residual: "
                  << std::sqrt(
                       Utilities::MPI::sum(this->squared_residual_L2_norm,
                                           mpi_communicator))
                  << std::endl;
    }
}


template class TimeHarmonicMaxwell<2>;
template class TimeHarmonicMaxwell<3>;
