// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/time_harmonic_maxwell.h>

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
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order));
      fe_test = std::make_shared<FESystem<dim>>(
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order),
        FE_Nedelec<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order));
      mapping = std::make_shared<MappingQ<dim>>(fe_trial_interior->degree);
      cell_quadrature = std::make_shared<QGauss<dim>>(fe_test->degree + 1);
      face_quadrature = std::make_shared<QGauss<dim - 1>>(fe_test->degree + 1);
    }

  // Initialize solutions and DPG error indicator shared_ptr
  present_solution            = std::make_shared<GlobalVectorType>();
  present_solution_skeleton   = std::make_shared<GlobalVectorType>();
  present_DPG_error_indicator = std::make_shared<GlobalVectorType>();

  // Allocate solution transfer
  solution_transfer = std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(
    *dof_handler_trial_interior);
}


template <>
std::pair<Tensor<1, 2, std::complex<double>>,
          Tensor<1, 2, std::complex<double>>>
TimeHarmonicMaxwell<2>::compute_waveguide_port_incident_fields(
  const Point<2> & /*p*/,
  const Tensor<1, 2> & /*normal*/,
  const std::complex<double> & /*effective_electric_permittivity*/,
  const std::complex<double> & /*effective_magnetic_permeability*/,
  const types::boundary_id /*boundary_id_index*/)
{
  // The waveguide port excitation would be completely different in 2D as curls
  // and cross products are not defined the same way as in 3D. Therefore, we do
  // not implement it for now and throw an error if someone tries to use the 2D
  // version.
  AssertThrow(false, TimeHarmonicMaxwellDimensionNotSupported(2));
  return std::pair<Tensor<1, 2, std::complex<double>>,
                   Tensor<1, 2, std::complex<double>>>();
}

template <>
std::pair<Tensor<1, 3, std::complex<double>>,
          Tensor<1, 3, std::complex<double>>>
TimeHarmonicMaxwell<3>::compute_waveguide_port_incident_fields(
  const Point<3>             &p,
  const Tensor<1, 3>         &normal,
  const std::complex<double> &effective_electric_permittivity,
  const std::complex<double> &effective_magnetic_permeability,
  const types::boundary_id    boundary_id_index)
{
  // Define some constexpr values for the computation
  static constexpr std::complex<double> imag{0., 1.};
  static constexpr double               PI  = numbers::PI;
  static constexpr unsigned int         dim = 3;

  // Gather the relevant waveguide parameters
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters =
    this->simulation_parameters.multiphysics.time_harmonic_maxwell_parameters;
  const double omega =
    2.0 * PI * time_harmonic_maxwell_parameters.electromagnetic_frequency;
  const auto &waveguide_corners =
    time_harmonic_maxwell_parameters.waveguide_corners[boundary_id_index];
  const Parameters::WaveguideMode mode =
    time_harmonic_maxwell_parameters.waveguide_mode[boundary_id_index];
  unsigned int m =
    time_harmonic_maxwell_parameters.mode_order_m[boundary_id_index];
  unsigned int n =
    time_harmonic_maxwell_parameters.mode_order_n[boundary_id_index];

  // We first define the transverse face vectors of the waveguide in the global
  // system
  Tensor<1, dim> transverse_vector_1 =
    waveguide_corners[1] - waveguide_corners[0];
  Tensor<1, dim> transverse_vector_2 =
    waveguide_corners[2] - waveguide_corners[0];
  double         length_t1 = transverse_vector_1.norm();
  double         length_t2 = transverse_vector_2.norm();
  Tensor<1, dim> e_t1      = transverse_vector_1 / length_t1;
  Tensor<1, dim> e_t2      = transverse_vector_2 / length_t2;

  // Check if the transverse vectors form a perfect rectangle (orthogonal and
  // aligned with the axes)
  AssertThrow(
    std::abs(e_t1 * e_t2) < 1e-12,
    ExcMessage(
      "The transverse plane defined by the waveguide corners for the waveguide port at boundary ID " +
      std::to_string(time_harmonic_maxwell_parameters
                       .waveguide_boundary_ids[boundary_id_index]) +
      " is not a perfect rectangle (i.e., the vector created by the waveguide corners are not orthogonal). Please check the waveguide corners definition in the input prm file."));


  // Also verify that those transverse vectors are perpendicular to the normal
  // of the
  // face and form an orthogonal basis.
  if ((std::abs(normal * e_t1) > 1e-12) || (std::abs(normal * e_t2) > 1e-12))
    AssertThrow(
      false,
      ExcMessage(
        "The transverse plane defined by the waveguide corners for the waveguide port at boundary ID " +
        std::to_string(time_harmonic_maxwell_parameters
                         .waveguide_boundary_ids[boundary_id_index]) +
        " is not orthogonal to the boundary face normal. Please check the waveguide corners definition in the input prm file."));

  // Create a third vector to complete the right-handed coordinate system
  Tensor<1, dim> e_t3 = cross_product_3d(e_t1, e_t2);

  // Determine if the system needs to be flipped so the e_t3 vector points in
  // the direction opposite to the outward normal of the face boundary. For an
  // incident wave at the inlet, the propagation direction should point into the
  // domain (opposite to the outward normal). We use the sign of (normal · t3)
  // to determine this:
  //   - If normal · t3 > 0: t3 points outward, need to flip entire system
  //   - If normal · t3 < 0: t3 points inward (correct for incident wave)
  // Note that by swapping the basis vectors e_t1 and e_t2, we change the parity
  // of the system since it is equivalent to a reflection. This means that
  // pseudo vectors (like the magnetic field) will not change sign while regular
  // vectors (like the electric field) will. Even though this does not affect
  // the physics of the solution, it is important to be consistent with the
  // definition of the mode profiles that we use, which assume a specific parity
  // for the system. Therefore, we keep track of the status of the parity and
  // apply it later on to the magnetic field to make it consistent with the
  // switch of sign that has been applied to the electric field when we compute
  // the excitation.
  double parity_factor = 1.0;
  if ((normal * e_t3) > 0)
    {
      std::swap(e_t1, e_t2);
      std::swap(length_t1, length_t2);
      std::swap(m, n);
      e_t3 =
        cross_product_3d(e_t1, e_t2); // Recompute t3 after swapping t1 and t2
                                      // to be sure the system is right-handed
      parity_factor = -1.0;
    }

  // Compute the various wavenumbers k in using this global coordinate system
  double k_t1 = m * PI / length_t1;                   // Transverse wavenumber 1
  double k_t2 = n * PI / length_t2;                   // Transverse wavenumber 2
  double k_c  = std::sqrt(k_t1 * k_t1 + k_t2 * k_t2); // Cutoff wavenumber
  std::complex<double> k =
    omega *
    std::sqrt(effective_electric_permittivity *
              effective_magnetic_permeability); // Wavenumber in the medium
  std::complex<double> k_l = std::sqrt(
    k * k - std::complex<double>(k_c * k_c, 0)); // Longitudinal wavenumber

  // Verify that the mode is not evanescent, i.e. k_l is not purely imaginary
  // (k_c^2 < k^2). std::norm computes the squared magnitude of a complex
  // number.
  AssertThrow(std::norm(k) > (k_c * k_c),
              ExcMessage(
                "The chosen mode for the waveguide port at boundary ID " +
                std::to_string(time_harmonic_maxwell_parameters
                                 .waveguide_boundary_ids[boundary_id_index]) +
                " is evanescent at the given frequency. Please "
                "choose another mode or increase the frequency."));

  // Now we want to compute the electromagnetic field and surface admittance at
  // point p as if the waveguide center was at the origin. So we will perform a
  // change of basis to a local coordinate system where the waveguide center is
  // at the origin.
  Tensor<1, dim> origin_local =
    0.25 * (waveguide_corners[0] + waveguide_corners[1] + waveguide_corners[2] +
            waveguide_corners[3]);

  Tensor<1, dim> p_local =
    p - origin_local; // Coordinates of point p in this local system.

  double x_local = p_local * e_t1; // Coordinate along e_t1 in the local system
  double y_local = p_local * e_t2; // Coordinate along e_t2 in the local system
  // We assume that the z_local coordinate is 0 since we are on the face.

  // Compute the E and H field components for the TE mode in the local
  // coordinate system {t1, t2, t3} = {x', y', z'}. We assume z' = 0 at the
  // boundary.
  Tensor<1, dim, std::complex<double>> E_inc_local;
  Tensor<1, dim, std::complex<double>> H_inc_local;

  if (mode == Parameters::WaveguideMode::TE)
    {
      std::complex<double> factor =
        imag * omega * effective_magnetic_permeability / (k_c * k_c);

      E_inc_local[0] = -factor * k_t2 *
                       std::cos(k_t1 * (x_local + length_t1 / 2)) *
                       std::sin(k_t2 * (y_local + length_t2 / 2));
      E_inc_local[1] = factor * k_t1 *
                       std::sin(k_t1 * (x_local + length_t1 / 2)) *
                       std::cos(k_t2 * (y_local + length_t2 / 2));
      E_inc_local[2] = 0.0;

      H_inc_local[0] = -imag * k_l * k_t1 / (k_c * k_c) *
                       std::sin(k_t1 * (x_local + length_t1 / 2)) *
                       std::cos(k_t2 * (y_local + length_t2 / 2));
      H_inc_local[1] = -imag * k_l * k_t2 / (k_c * k_c) *
                       std::cos(k_t1 * (x_local + length_t1 / 2)) *
                       std::sin(k_t2 * (y_local + length_t2 / 2));
      H_inc_local[2] = std::cos(k_t1 * (x_local + length_t1 / 2)) *
                       std::cos(k_t2 * (y_local + length_t2 / 2));
    }
  else if (mode == Parameters::WaveguideMode::TM)
    {
      std::complex<double> factor =
        imag * omega * effective_electric_permittivity / (k_c * k_c);

      H_inc_local[0] = factor * k_t2 *
                       std::sin(k_t1 * (x_local + length_t1 / 2)) *
                       std::cos(k_t2 * (y_local + length_t2 / 2));
      H_inc_local[1] = -factor * k_t1 *
                       std::cos(k_t1 * (x_local + length_t1 / 2)) *
                       std::sin(k_t2 * (y_local + length_t2 / 2));
      H_inc_local[2] = 0.0;

      E_inc_local[0] = imag * k_l * k_t1 / (k_c * k_c) *
                       std::cos(k_t1 * (x_local + length_t1 / 2)) *
                       std::sin(k_t2 * (y_local + length_t2 / 2));
      E_inc_local[1] = imag * k_l * k_t2 / (k_c * k_c) *
                       std::sin(k_t1 * (x_local + length_t1 / 2)) *
                       std::cos(k_t2 * (y_local + length_t2 / 2));
      E_inc_local[2] = std::sin(k_t1 * (x_local + length_t1 / 2)) *
                       std::sin(k_t2 * (y_local + length_t2 / 2));
    }
  else
    {
      AssertThrow(false, ExcMessage("Unknown waveguide mode type."));
    }

  // Convert the E and H field components from the local coordinate system back
  // to the global coordinate system using the basis vectors e_t1, e_t2, e_t3
  Tensor<1, dim, std::complex<double>> E_inc =
    E_inc_local[0] * e_t1 + E_inc_local[1] * e_t2 + E_inc_local[2] * e_t3;
  Tensor<1, dim, std::complex<double>> H_inc =
    parity_factor *
    (H_inc_local[0] * e_t1 + H_inc_local[1] * e_t2 + H_inc_local[2] * e_t3);

  return std::make_pair(E_inc, H_inc);
}


template <>
std::pair<Tensor<1, 2, std::complex<double>>, std::complex<double>>
TimeHarmonicMaxwell<2>::compute_waveguide_port_excitation(
  const Point<2> & /*p*/,
  const Tensor<1, 2> & /*normal*/,
  const std::complex<double> & /*effective_electric_permittivity*/,
  const std::complex<double> & /*effective_magnetic_permeability*/,
  const types::boundary_id /*boundary_id_index*/)
{
  // The waveguide port excitation would be completely different in 2D as curls
  // and cross products are not defined the same way as in 3D. Therefore, we do
  // not implement it for now and throw an error if someone tries to use the 2D
  // version.
  AssertThrow(false, TimeHarmonicMaxwellDimensionNotSupported(2));
  return std::pair<Tensor<1, 2, std::complex<double>>, std::complex<double>>();
}

template <>
std::pair<Tensor<1, 3, std::complex<double>>, std::complex<double>>
TimeHarmonicMaxwell<3>::compute_waveguide_port_excitation(
  const Point<3>             &p,
  const Tensor<1, 3>         &normal,
  const std::complex<double> &effective_electric_permittivity,
  const std::complex<double> &effective_magnetic_permeability,
  const types::boundary_id    boundary_id_index)
{
  static constexpr unsigned int dim = 3;
  static constexpr double       PI  = numbers::PI;
  const auto                    incident_fields =
    compute_waveguide_port_incident_fields(p,
                                           normal,
                                           effective_electric_permittivity,
                                           effective_magnetic_permeability,
                                           boundary_id_index);

  // Gather the relevant waveguide parameters
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters =
    this->simulation_parameters.multiphysics.time_harmonic_maxwell_parameters;
  const double omega =
    2.0 * PI * time_harmonic_maxwell_parameters.electromagnetic_frequency;
  const auto &waveguide_corners =
    time_harmonic_maxwell_parameters.waveguide_corners[boundary_id_index];
  const Parameters::WaveguideMode mode =
    time_harmonic_maxwell_parameters.waveguide_mode[boundary_id_index];
  unsigned int m =
    time_harmonic_maxwell_parameters.mode_order_m[boundary_id_index];
  unsigned int n =
    time_harmonic_maxwell_parameters.mode_order_n[boundary_id_index];

  const double k_t1 = m * PI /
                      (waveguide_corners[1] - waveguide_corners[0])
                        .norm(); // Transverse wavenumber 1
  const double k_t2 = n * PI /
                      (waveguide_corners[2] - waveguide_corners[0])
                        .norm(); // Transverse wavenumber 2
  const std::complex<double> k_l =
    std::sqrt(omega * omega * effective_electric_permittivity *
                effective_magnetic_permeability -
              (k_t1 * k_t1 + k_t2 * k_t2)); // Longitudinal wavenumber

  std::complex<double> surface_admittance =
    (mode == Parameters::WaveguideMode::TE) ?
      k_l / (omega * effective_magnetic_permeability) :
      omega * effective_electric_permittivity / k_l;

  const Tensor<1, dim, std::complex<double>> &E_inc = incident_fields.first;
  const Tensor<1, dim, std::complex<double>> &H_inc = incident_fields.second;

  // We normalize the excitation by the electromagnetic scaling factor
  const Tensor<1, dim, std::complex<double>> excitation =
    (this->waveguide_ports_electric_amplitudes[boundary_id_index] /
     this->electromagnetic_scaling) *
    (cross_product_3d(normal, H_inc) +
     map_H12(surface_admittance * E_inc, normal));

  return std::make_pair(excitation, surface_admittance);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::update_material_properties(
  const PhysicalPropertiesManager &properties_manager,
  std::complex<double>            &effective_electric_permittivity,
  std::complex<double>            &effective_magnetic_permeability,
  const unsigned int               material_id)
{
  std::map<field, double>
    field_values; // Empty map since no field dependence for now
  effective_electric_permittivity = {
    properties_manager.get_electric_permittivity_real(0, material_id)
      ->value(field_values),
    properties_manager.get_electric_permittivity_imag(0, material_id)
        ->value(field_values) +
      properties_manager.get_electric_conductivity(0, material_id)
        ->value(field_values)};

  effective_magnetic_permeability = {
    properties_manager.get_magnetic_permeability_real(0, material_id)
      ->value(field_values),
    properties_manager.get_magnetic_permeability_imag(0, material_id)
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
      const unsigned int n_q_points_face =
        fe_face_values_trial_skeleton.n_quadrature_points;
      BoundaryConditions::BoundaryType bc_type;

      // Containers for electromagnetic quantities at the quadrature points
      Tensor<1, dim, std::complex<double>> E_inc;
      Tensor<1, dim, std::complex<double>> H_inc;
      std::complex<double>                 effective_electric_permittivity;
      std::complex<double>                 effective_magnetic_permeability;
      unsigned int                         material_id;

      // Here we perform the Poynting vector integration over the waveguide
      // inlet faces.
      for (const auto &cell : triangulation->active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              material_id = cell->material_id();
              update_material_properties(physical_properties_manager,
                                         effective_electric_permittivity,
                                         effective_magnetic_permeability,
                                         material_id);

              for (const auto &face : cell->face_iterators())
                {
                  fe_face_values_trial_skeleton.reinit(cell, face);

                  if (!(face->at_boundary()))
                    continue;
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

                  // Loop over all face quadrature points
                  for (unsigned int q_point = 0; q_point < n_q_points_face;
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
                          position,
                          normal,
                          effective_electric_permittivity,
                          effective_magnetic_permeability,
                          boundary_index);

                      // The Poynting vector is given by S = 0.5 * Re(E
                      // x H*), where H* is the complex conjugate of H.
                      // The normal points outward from the domain, so
                      // the power entering the domain is given by the
                      // negative of the flux of the Poynting vector
                      // through the face (that's why we have a negative
                      // sign in front of the integral).
                      waveguide_modal_powers[boundary_index] -=
                        0.5 *
                        std::real(normal[0] * (E_inc[1] * std::conj(H_inc[2]) -
                                               E_inc[2] * std::conj(H_inc[1])) +
                                  normal[1] * (E_inc[2] * std::conj(H_inc[0]) -
                                               E_inc[0] * std::conj(H_inc[2])) +
                                  normal[2] *
                                    (E_inc[0] * std::conj(H_inc[1]) -
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
          // electric field amplitude E_0.
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
      // We assume that the electric field amplitude provided by the user is
      // already scaled by the unit conversion factor, so we can directly use it
      // as the scaling factor.
      this->electromagnetic_scaling =
        electromagnetic_parameters.electric_field_amplitude;
    }
  else if (electromagnetic_parameters.electromagnetic_scaling_type ==
           Parameters::ElectromagneticScalingType::magnetic_field)
    {
      // We assume that the magnetic field amplitude provided by the user is
      // already scaled by the unit conversion factor, but the void impedance
      // needs to be rescaled from the dimensionality.
      this->electromagnetic_scaling =
        (void_impedance *
         electromagnetic_parameters.electric_field_dimensionality /
         electromagnetic_parameters.magnetic_field_dimensionality) *
        electromagnetic_parameters.magnetic_field_amplitude;
    }
  else if (electromagnetic_parameters.electromagnetic_scaling_type ==
           Parameters::ElectromagneticScalingType::power)
    {
      this->electromagnetic_scaling =
        electromagnetic_parameters.electric_field_dimensionality *
        max_port_electric_amplitude;
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
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters =
    this->simulation_parameters.multiphysics.time_harmonic_maxwell_parameters;

  const auto &fe = dof_handler.get_fe();

  if (time_harmonic_maxwell_parameters.electromagnetic_scaling_type ==
      Parameters::ElectromagneticScalingType::none)
    return;

  constexpr double void_impedance = 4 * numbers::PI * 29.9792458;
  const double     electric_scale = this->electromagnetic_scaling;
  const double     magnetic_scale =
    this->electromagnetic_scaling *
    time_harmonic_maxwell_parameters.magnetic_field_dimensionality /
    (void_impedance *
     time_harmonic_maxwell_parameters.electric_field_dimensionality);

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
  if (!this->simulation_parameters.boundary_conditions.time_dependent)
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

  if (this->simulation_parameters.timer.type ==
      Parameters::Timer::Type::iteration)
    {
      announce_string(this->pcout, "Time Harmonic Electromagnetics");
      this->pcout << std::defaultfloat;
      this->computing_timer.print_summary();
      this->pcout << std::scientific;
      this->computing_timer.reset();
    }
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
  std::vector<const GlobalVectorType *> sol_set_transfer;

  solution_transfer = std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(
    *dof_handler_trial_interior);

  sol_set_transfer.emplace_back(&(*present_solution));

  std::string checkpoint_file_prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;

  solution_transfer->prepare_for_serialization(sol_set_transfer);

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

  std::vector<GlobalVectorType *> input_vectors(1);
  GlobalVectorType distributed_system(locally_owned_dofs_trial_interior,
                                      mpi_communicator);
  input_vectors[0] = &distributed_system;

  SolutionTransfer<dim, GlobalVectorType> system_trans_vectors(
    *this->dof_handler_trial_interior);

  std::string checkpoint_file_prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;

  solution_transfer->deserialize(input_vectors);

  *present_solution = distributed_system;

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
  // assembly of the system matrix and rhs, so we do not implement it here as a
  // separate loop. The estimated error per cell is computed and stored in the
  // member variable local_estimated_error_per_cell during the assembly, and
  // then used for marking the cells for refinement.
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
  // any cells don't falsely trigger the assertion.
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
  this->locally_owned_dofs_test = this->dof_handler_test->locally_owned_dofs();

  // Get the locally relevant dofs
  this->locally_relevant_dofs_trial_interior =
    DoFTools::extract_locally_relevant_dofs(*this->dof_handler_trial_interior);
  this->locally_relevant_dofs_trial_skeleton =
    DoFTools::extract_locally_relevant_dofs(*this->dof_handler_trial_skeleton);
  this->locally_relevant_dofs_test =
    DoFTools::extract_locally_relevant_dofs(*this->dof_handler_test);

  // Initialize the solution vectors and error indicator
  this->present_solution->reinit(this->locally_owned_dofs_trial_interior,
                                 this->locally_relevant_dofs_trial_interior,
                                 mpi_communicator);
  this->present_solution_skeleton->reinit(
    this->locally_owned_dofs_trial_skeleton,
    this->locally_relevant_dofs_trial_skeleton,
    mpi_communicator);
  this->present_DPG_error_indicator->reinit(this->locally_owned_dofs_test,
                                            this->locally_relevant_dofs_test,
                                            mpi_communicator);
  this->local_estimated_error_per_cell.reinit(triangulation->n_active_cells());

  // We reinitialize the system rhs with the skeleton dofs because we have
  // performed a static condensation of the interior dofs using the Schur
  // complement.
  this->system_rhs.reinit(this->locally_owned_dofs_trial_skeleton,
                          mpi_communicator);

  // Define constraints
  define_constraints();

  // Sparse matrices initialization
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

  this->system_matrix.reinit(this->locally_owned_dofs_trial_skeleton,
                             this->locally_owned_dofs_trial_skeleton,
                             dsp,
                             mpi_communicator);
  this->pcout << "  DPG system for Time-Harmonic Maxwell Equations:"
              << std::endl;
  this->pcout
    << "   Number of skeleton degrees of freedom for Time-Harmonic Maxwell: "
    << this->dof_handler_trial_skeleton->n_dofs() << std::endl;
  this->pcout
    << "   Number of interior degrees of freedom for Time-Harmonic Maxwell: "
    << this->dof_handler_trial_interior->n_dofs() << std::endl;

  // Update the multiphysics interface with the TimeHarmonicMaxwell dof_handler
  // and solution. Note that we provide the interior dof_handler and interior
  // solution as it is the one useful for the other physics (we do not provide
  // the skeleton dof_handler and skeleton solution).
  multiphysics->set_dof_handler(PhysicsID::electromagnetics,
                                this->dof_handler_trial_interior);
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
  this->nonzero_constraints.close();
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::solve_linear_system()
{
  TimerOutput::Scope t(this->computing_timer, "Solve linear system");

  auto mpi_communicator = this->triangulation->get_mpi_communicator();

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
        .verbosity == Parameters::Verbosity::extra_verbose)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  // Set up the solver control
  SolverControl solver_control(this->dof_handler_trial_skeleton->n_dofs(),
                               linear_solver_tolerance);

  // Solve
  GlobalVectorType completely_distributed_solution(
    this->locally_owned_dofs_trial_skeleton, mpi_communicator);
  TrilinosWrappers::SolverCG solver(solver_control);

  solver.solve(this->system_matrix,
               completely_distributed_solution,
               this->system_rhs,
               *this->preconditioner);

  if (simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
        .verbosity == Parameters::Verbosity::extra_verbose)
    {
      this->pcout << "  -CG iterative solver took : "
                  << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() / rescale_metric << std::endl;
    }

  // Update the solution vector for the skeleton
  nonzero_constraints.distribute(completely_distributed_solution);
  *this->present_solution_skeleton = completely_distributed_solution;

  // Reconstruct the interior solution from the skeleton solution
  reconstruct_interior_solution();


  // We need to apply the scaling to the interior solution as it is the one
  // used for the multiphysics coupling and output.
  scale_solution_components(*this->dof_handler_trial_interior,
                            *this->present_solution);

  // We also apply the scaling to the skeleton solution for
  // consistency, even if it is not used for the multiphysics coupling nor
  // outputed at the moment.
  scale_solution_components(*this->dof_handler_trial_skeleton,
                            *this->present_solution_skeleton);
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

  // Always solve at the first step of the simulation (simulation start as 0 but
  // it is then incremented before solving the physics for the first time, so
  // the first time this function is called, the step number is 1).
  if (this->simulation_control->get_step_number() == 1)
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
            if ((this->simulation_control->get_step_number() - 1) %
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
            AssertThrow(
              false,
              ExcMessage(
                "Time coupling strategy 'threshold' is not yet implemented "
                "for the time-harmonic Maxwell solver since the physical "
                "properties only support a constant model."));
            return true;
          default:
            AssertThrow(false, ExcMessage("Unknown time coupling strategy."));
            return false;
        }
    }
}

template <>
void
TimeHarmonicMaxwell<2>::assemble_system_matrix()
{
  // We need a specific definition of the function for the 2D so the compiler
  // doesn't try to build curl and cross operations at compile time that will
  // never be used anyway. Indeed, curl and cross operations behave very
  // differently in 2D than in 3D. This physics only support 3D problem at the
  // moment. So even though the class is templated in dim, we only want to
  // compile the whole function when dim=3.
  AssertThrow(false, TimeHarmonicMaxwellDimensionNotSupported(2));
}

template <>
void
TimeHarmonicMaxwell<3>::assemble_system_matrix()
{
  // Constexpr values that are used in the assembly. Since we are in the
  // specialized 3D function we define dim=3 here. This makes easier to read the
  // code below to see what are templated in dim and what are not.
  static constexpr double               PI = numbers::PI;
  static constexpr std::complex<double> imag{0., 1.};
  static constexpr int                  dim = 3;

  // Get properties manager and define model physical properties
  const auto &properties_manager =
    this->simulation_parameters.physical_properties_manager;
  std::complex<double> effective_electric_permittivity;
  std::complex<double> effective_magnetic_permeability;
  unsigned int         material_id;

  /// Excitation properties
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters =
    this->simulation_parameters.multiphysics.time_harmonic_maxwell_parameters;
  const double omega =
    2.0 * PI * time_harmonic_maxwell_parameters.electromagnetic_frequency;

  // We always need to compute the scaling for the electromagnetic fields even
  // if the user does not want to apply it to the solution because it is used to
  // normalize the waveguide inlets relative input power
  compute_electromagnetic_scaling(properties_manager);

  TimerOutput::Scope t(this->computing_timer, "Assemble matrix and RHS");


  // We then create the corresponding FEValues and FEFaceValues objects. Note
  // that only the test space needs gradients because of the ultraweak
  // formulation. Similarly, because everything is on the same triangulation,
  // we only need to update the quadrature points and JxW values in one of the
  // spaces. Here we choose the trial space.
  const unsigned int n_q_points      = this->cell_quadrature->size();
  const unsigned int n_face_q_points = this->face_quadrature->size();

  FEValues<dim> fe_values_trial_interior(*this->mapping,
                                         *this->fe_trial_interior,
                                         *this->cell_quadrature,
                                         update_values |
                                           update_quadrature_points |
                                           update_JxW_values);

  FEValues<dim> fe_values_test(*this->mapping,
                               *this->fe_test,
                               *this->cell_quadrature,
                               update_values | update_gradients);

  FEFaceValues<dim> fe_face_values_trial_skeleton(*this->mapping,
                                                  *this->fe_trial_skeleton,
                                                  *this->face_quadrature,
                                                  update_values |
                                                    update_quadrature_points |
                                                    update_normal_vectors |
                                                    update_JxW_values);

  FEFaceValues<dim> fe_face_values_test(*this->mapping,
                                        *this->fe_test,
                                        *this->face_quadrature,
                                        update_values);

  // We also create all the relevant matrices and vector to build the DPG
  // system. To do so we first need the number of dofs per cell for each of
  // the finite element spaces.
  const unsigned int dofs_per_cell_test = this->fe_test->n_dofs_per_cell();
  const unsigned int dofs_per_cell_trial_interior =
    this->fe_trial_interior->n_dofs_per_cell();
  const unsigned int dofs_per_cell_trial_skeleton =
    this->fe_trial_skeleton->n_dofs_per_cell();

  // To avoid unecessary computations, we will precompute the shape functions
  // of all our spaces once and then use those precomputed values to assemble
  // the local matrices. Consequently, we need to create containers to store
  // these values. Since those are complex-valued functions, we use two
  // different containers for the real and imaginary parts because the
  // conjugate of a complex tensor is not implemented in deal.II.
  std::vector<Tensor<1, dim, std::complex<double>>> F(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> F_conj(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> I(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> I_conj(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> curl_F(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> curl_F_conj(
    dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> curl_I(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> curl_I_conj(
    dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> F_face(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> F_face_conj(
    dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> I_face_conj(
    dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> n_cross_I_face(
    dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> n_cross_I_face_conj(
    dofs_per_cell_test);

  std::vector<Tensor<1, dim, std::complex<double>>> E(
    dofs_per_cell_trial_interior);
  std::vector<Tensor<1, dim, std::complex<double>>> H(
    dofs_per_cell_trial_interior);
  std::vector<Tensor<1, dim, std::complex<double>>> E_hat(
    dofs_per_cell_trial_skeleton);
  std::vector<Tensor<1, dim, std::complex<double>>> n_cross_E_hat(
    dofs_per_cell_trial_skeleton);
  std::vector<Tensor<1, dim, std::complex<double>>> n_cross_H_hat(
    dofs_per_cell_trial_skeleton);

  // Also, to avoid multiple "if" calls during the assembly to understand where
  // each term needs to be assembled in the DPG global matrix, we will use
  // containers to store dofs relationships for each of the local matrices.

  // G matrix stands for the Riesz map and needs the relationships between test
  // functions F and I
  std::vector<std::pair<unsigned int, unsigned int>> G_FF;
  std::vector<std::pair<unsigned int, unsigned int>> G_FI;
  std::vector<std::pair<unsigned int, unsigned int>> G_IF;
  std::vector<std::pair<unsigned int, unsigned int>> G_II;

  // B matrix stands for the bilinear form and needs the relationships between
  // interior trial space (E and H) and the test (F and I)
  std::vector<std::pair<unsigned int, unsigned int>> B_FE;
  std::vector<std::pair<unsigned int, unsigned int>> B_IE;
  std::vector<std::pair<unsigned int, unsigned int>> B_FH;
  std::vector<std::pair<unsigned int, unsigned int>> B_IH;

  // B_hat matrix stands for the bilinear form, but on the skeleton, and needs
  // the relationships between skeleton trial space (E_hat, H_hat) and the test
  // (F and I). Note that the B_hat_FE is required for the Robin boundary
  // condition that we chose to apply on the electric field instead of the
  // magnetic field. This choice is arbitrary but it cannot be applied on both
  // fields at the same time.
  std::vector<std::pair<unsigned int, unsigned int>> B_hat_IE;
  std::vector<std::pair<unsigned int, unsigned int>> B_hat_FH;
  std::vector<std::pair<unsigned int, unsigned int>> B_hat_FE;

  // l vector stands for the linear form of the problem and needs the
  // relationships between test functions. Note that in its simplest form, the
  // time-harmonic Maxwell equation does not have a magnetic source term, so
  // there is no contribution from the I test functions.
  std::vector<unsigned int> l_F;

  // Reserve memory to avoid reallocations of each relationship vectors
  G_FF.reserve(dofs_per_cell_test * dofs_per_cell_test);
  G_FI.reserve(dofs_per_cell_test * dofs_per_cell_test);
  G_IF.reserve(dofs_per_cell_test * dofs_per_cell_test);
  G_II.reserve(dofs_per_cell_test * dofs_per_cell_test);

  B_FE.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);
  B_IE.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);
  B_FH.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);
  B_IH.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);

  B_hat_IE.reserve(dofs_per_cell_trial_skeleton * dofs_per_cell_test);
  B_hat_FH.reserve(dofs_per_cell_trial_skeleton * dofs_per_cell_test);
  B_hat_FE.reserve(dofs_per_cell_trial_skeleton * dofs_per_cell_test);

  l_F.reserve(dofs_per_cell_test);

  // Here we create the DPG local matrices and vector used for the assembly
  // before condensation.
  LAPACKFullMatrix<double> G_matrix(dofs_per_cell_test, dofs_per_cell_test);

  LAPACKFullMatrix<double> B_matrix(dofs_per_cell_test,
                                    dofs_per_cell_trial_interior);

  LAPACKFullMatrix<double> B_hat_matrix(dofs_per_cell_test,
                                        dofs_per_cell_trial_skeleton);

  Vector<double> l_vector(dofs_per_cell_test);

  // We create the condensation matrices which are defined as :
  // $M_1 = B^\dagger G^{-1}B$;
  // $M_2 = B^\dagger G^{-1}\hat{B}$;
  // $M_3 = \hat{B}^\dagger G^{-1}\hat{B}$;
  // $M_4 = B^\dagger G^{-1}$;
  // $M_5 = \hat{B}^\dagger G^{-1}$.
  LAPACKFullMatrix<double> M1_matrix(dofs_per_cell_trial_interior,
                                     dofs_per_cell_trial_interior);
  LAPACKFullMatrix<double> M2_matrix(dofs_per_cell_trial_interior,
                                     dofs_per_cell_trial_skeleton);
  LAPACKFullMatrix<double> M3_matrix(dofs_per_cell_trial_skeleton,
                                     dofs_per_cell_trial_skeleton);
  LAPACKFullMatrix<double> M4_matrix(dofs_per_cell_trial_interior,
                                     dofs_per_cell_test);
  LAPACKFullMatrix<double> M5_matrix(dofs_per_cell_trial_skeleton,
                                     dofs_per_cell_test);

  // During the calculation of matrix vector product, we require intermediary
  // matrices and vector that we also allocate here. The temporary matrices
  // are defined as :
  // $tmp_matrix_M2M1  = M_2^\dagger M_1^{-1}$;
  // $tmp_matrix_M2M1M2 = M_2^\dagger M_1^{-1} M_2$;
  // $tmp_matrix_M2M1M4 = M_2^\dagger M_1^{-1} M_4$.

  LAPACKFullMatrix<double> tmp_matrix_M2M1(dofs_per_cell_trial_skeleton,
                                           dofs_per_cell_trial_interior);

  LAPACKFullMatrix<double> tmp_matrix_M2M1M2(dofs_per_cell_trial_skeleton,
                                             dofs_per_cell_trial_skeleton);

  LAPACKFullMatrix<double> tmp_matrix_M2M1M4(dofs_per_cell_trial_skeleton,
                                             dofs_per_cell_test);

  // We create the cell matrix and the RHS that will be distributed in the
  // full system after the assembly along with the index’s mapping.
  FullMatrix<double> cell_matrix(dofs_per_cell_trial_skeleton,
                                 dofs_per_cell_trial_skeleton);
  Vector<double>     cell_skeleton_rhs(dofs_per_cell_trial_skeleton);

  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell_trial_skeleton);

  // We also create objects used for the various Robin boundary conditions.
  BoundaryConditions::BoundaryType bc_type(
    BoundaryConditions::BoundaryType::none);
  Tensor<1, dim, std::complex<double>> g_inc;
  std::complex<double>                 boundary_surface_admittance;
  std::complex<double>                 conj_boundary_surface_admittance;

  // As it is standard we first loop over the cells of the triangulation. Here
  // we have the choice of the DoFHandler to perform this loop. We use the
  // DoFHandler associated with the trial space.
  for (const auto &cell :
       this->dof_handler_trial_interior->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // We update the material properties for the current cell and define
          // some constants that will be used during the assembly.
          material_id = cell->material_id();
          update_material_properties(properties_manager,
                                     effective_electric_permittivity,
                                     effective_magnetic_permeability,
                                     material_id);
          const std::complex<double> iweffective_magnetic_permeability =
            imag * omega * effective_magnetic_permeability;
          const std::complex<double> conj_iweffective_magnetic_permeability =
            std::conj(iweffective_magnetic_permeability);
          const std::complex<double> iweps_r =
            imag * omega * effective_electric_permittivity;
          const std::complex<double> conj_iweps_r = std::conj(iweps_r);

          // We reinitialize the FEValues objects to the current cell.
          fe_values_trial_interior.reinit(cell);

          // We will also need to reinitialize the FEValues for the test
          // space and make sure that is the same cell as the one used for the
          // trial space.
          const typename DoFHandler<dim>::active_cell_iterator cell_test =
            cell->as_dof_handler_iterator(*this->dof_handler_test);
          fe_values_test.reinit(cell_test);

          // Similarly, we reinitialize the FEValues for the trial space on
          // the skeleton, but this will not be used before we also loop on
          // the cells faces.
          const typename DoFHandler<dim>::active_cell_iterator cell_skeleton =
            cell->as_dof_handler_iterator(*this->dof_handler_trial_skeleton);

          // We then reinitialize all the matrices where we are aggregating
          // information for the current cell.
          G_matrix     = 0;
          B_matrix     = 0;
          B_hat_matrix = 0;
          l_vector     = 0;

          // We also need to reinitialize the $M_1$ condensation matrix
          // between each iteration on cell to get rid of its inverse status.
          M1_matrix = 0;

          // Here we reset the dofs relationships containers. These containers
          // are used to store all the relationships between the dof i
          // dof j according to which space they belong to. Indeed, when
          // looping over all the test and trial dofs, the terms  we need to
          // compute depend on the specific combination of spaces
          // involved. The following containers are therefore used to avoid
          // multiple "if" statements inside the dofs loops and are filled with
          // all the relevant dofs pairs that we need for each of the
          // different terms.

          // For example, in the ultraweak form of Maxwell equation we have a
          // term (E, curl(I)) in the interior so the matrix B has a term
          // (curl(I_i), E_j), so we want to only compute this term for the
          // dofs i that are in the test space I and the dofs j that are in
          // the trial space for E. Therefore, when looping on all the
          // interior and test dofs in a cell we add the pairs of dofs that
          // are in those spaces to the container B_IE. In the end, the
          // container B_IE has all the required pairs of dofs indices for which
          // we need to compute the term. We do this for all the different terms
          // of the formulation.

          // Finally, we can loop on all the pairs that we have assigned in
          // each vector container to compute the desired terms.
          G_FF.clear();
          G_FI.clear();
          G_IF.clear();
          G_II.clear();

          B_FE.clear();
          B_IE.clear();
          B_FH.clear();
          B_IH.clear();

          l_F.clear();

          // We fill the dofs relationship containers at the cell level. To do
          // so, we first loop on the test space dofs.
          for (unsigned int i : fe_values_test.dof_indices())
            {
              // Get the information on which element the dof is
              const unsigned int current_element_test_i =
                this->fe_test->system_to_base_index(i).first.first;

              // Fill the load vector relationship
              if ((current_element_test_i == 0) ||
                  (current_element_test_i == 1))
                {
                  l_F.emplace_back(i);
                }

              // Loop over the test dofs a second time to fill the dofs
              // relationship for the G matrix (Riesz map)
              for (unsigned int j : fe_values_test.dof_indices())
                {
                  const unsigned int current_element_test_j =
                    this->fe_test->system_to_base_index(j).first.first;
                  if (((current_element_test_i == 0) ||
                       (current_element_test_i == 1)) &&
                      ((current_element_test_j == 0) ||
                       (current_element_test_j == 1)))
                    {
                      G_FF.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 0) ||
                       (current_element_test_i == 1)) &&
                      ((current_element_test_j == 2) ||
                       (current_element_test_j == 3)))
                    {
                      G_FI.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 2) ||
                       (current_element_test_i == 3)) &&
                      ((current_element_test_j == 0) ||
                       (current_element_test_j == 1)))
                    {
                      G_IF.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 2) ||
                       (current_element_test_i == 3)) &&
                      ((current_element_test_j == 2) ||
                       (current_element_test_j == 3)))
                    {
                      G_II.emplace_back(i, j);
                    }
                }

              // Then we loop over the trial dofs space to fill the dofs
              // relationship for the B matrix (bilinear form)
              for (unsigned int j : fe_values_trial_interior.dof_indices())
                {
                  const unsigned int current_element_trial_j =
                    this->fe_trial_interior->system_to_base_index(j)
                      .first.first;

                  if (((current_element_test_i == 0) ||
                       (current_element_test_i == 1)) &&
                      ((current_element_trial_j == 0) ||
                       (current_element_trial_j == 1)))
                    {
                      B_FE.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 0) ||
                       (current_element_test_i == 1)) &&
                      ((current_element_trial_j == 2) ||
                       (current_element_trial_j == 3)))
                    {
                      B_FH.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 2) ||
                       (current_element_test_i == 3)) &&
                      ((current_element_trial_j == 0) ||
                       (current_element_trial_j == 1)))
                    {
                      B_IE.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 2) ||
                       (current_element_test_i == 3)) &&
                      ((current_element_trial_j == 2) ||
                       (current_element_trial_j == 3)))
                    {
                      B_IH.emplace_back(i, j);
                    }
                }
            }

          // Now we loop over all quadrature points of the cell
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              // To avoid unnecessary computation, we fill the shape values
              // containers for the real and imaginary parts of the electric
              // and magnetic fields and the dofs relationship at the current
              // quadrature point.
              const double &JxW = fe_values_trial_interior.JxW(q_point);

              for (unsigned int i : fe_values_test.dof_indices())
                {
                  F[i] =
                    fe_values_test[extractor_E_real].value(i, q_point) +
                    imag * fe_values_test[extractor_E_imag].value(i, q_point);
                  F_conj[i] =
                    fe_values_test[extractor_E_real].value(i, q_point) -
                    imag * fe_values_test[extractor_E_imag].value(i, q_point);

                  curl_F[i] =
                    fe_values_test[extractor_E_real].curl(i, q_point) +
                    imag * fe_values_test[extractor_E_imag].curl(i, q_point);
                  curl_F_conj[i] =
                    fe_values_test[extractor_E_real].curl(i, q_point) -
                    imag * fe_values_test[extractor_E_imag].curl(i, q_point);

                  I[i] =
                    fe_values_test[extractor_H_real].value(i, q_point) +
                    imag * fe_values_test[extractor_H_imag].value(i, q_point);
                  I_conj[i] =
                    fe_values_test[extractor_H_real].value(i, q_point) -
                    imag * fe_values_test[extractor_H_imag].value(i, q_point);

                  curl_I[i] =
                    fe_values_test[extractor_H_real].curl(i, q_point) +
                    imag * fe_values_test[extractor_H_imag].curl(i, q_point);
                  curl_I_conj[i] =
                    fe_values_test[extractor_H_real].curl(i, q_point) -
                    imag * fe_values_test[extractor_H_imag].curl(i, q_point);
                }

              for (unsigned int i : fe_values_trial_interior.dof_indices())
                {
                  E[i] =
                    fe_values_trial_interior[extractor_E_real].value(i,
                                                                     q_point) +
                    imag *
                      fe_values_trial_interior[extractor_E_imag].value(i,
                                                                       q_point);
                  H[i] =
                    fe_values_trial_interior[extractor_H_real].value(i,
                                                                     q_point) +
                    imag *
                      fe_values_trial_interior[extractor_H_imag].value(i,
                                                                       q_point);
                }

              // Now we loop on each relationship container to assemble the
              // relevant matrices
              for (const auto &[i, j] : G_FF)
                {
                  G_matrix(i, j) +=
                    (((F[j] * F_conj[i]) + (curl_F[j] * curl_F_conj[i]) +
                      (conj_iweps_r * F[j] * iweps_r * F_conj[i])) *
                     JxW)
                      .real();
                }

              for (const auto &[i, j] : G_FI)
                {
                  G_matrix(i, j) += (((curl_I[j] * iweps_r * F_conj[i]) -
                                      (conj_iweffective_magnetic_permeability *
                                       I[j] * curl_F_conj[i])) *
                                     JxW)
                                      .real();
                }

              for (const auto &[i, j] : G_IF)
                {
                  G_matrix(i, j) +=
                    (((conj_iweps_r * F[j] * curl_I_conj[i]) -
                      (curl_F[j] * iweffective_magnetic_permeability *
                       I_conj[i])) *
                     JxW)
                      .real();
                }

              for (const auto &[i, j] : G_II)
                {
                  G_matrix(i, j) +=
                    (((I[j] * I_conj[i]) + (curl_I[j] * curl_I_conj[i]) +
                      (conj_iweffective_magnetic_permeability * I[j] *
                       iweffective_magnetic_permeability * I_conj[i])) *
                     JxW)
                      .real();
                }

              for (const auto &[i, j] : B_FE)
                {
                  B_matrix(i, j) += (iweps_r * E[j] * F_conj[i] * JxW).real();
                }

              for (const auto &[i, j] : B_FH)
                {
                  B_matrix(i, j) += (H[j] * curl_F_conj[i] * JxW).real();
                }

              for (const auto &[i, j] : B_IE)
                {
                  B_matrix(i, j) += (E[j] * curl_I_conj[i] * JxW).real();
                }

              for (const auto &[i, j] : B_IH)
                {
                  B_matrix(i, j) -=
                    (iweffective_magnetic_permeability * H[j] * I_conj[i] * JxW)
                      .real();
                }

              for (const auto &i : l_F)
                {
                  l_vector[i] += 0.0;
                }
            }

          // We now build the skeleton terms. Similarly, we choose to loop on
          // the skeleton trial space faces.
          for (const auto &face : cell_skeleton->face_iterators())
            {
              // We reinitialize the FEFaceValues objects to the current
              // faces.
              fe_face_values_test.reinit(cell_test, face);
              fe_face_values_trial_skeleton.reinit(cell_skeleton, face);

              // Get the boundary condition type on the current face
              bc_type = BoundaryConditions::BoundaryType::none;

              if (face->at_boundary())
                bc_type =
                  this->simulation_parameters
                    .boundary_conditions_time_harmonic_electromagnetics.type.at(
                      face->boundary_id());

              // Reset the face dofs relationships
              G_FF.clear();
              G_FI.clear();
              G_IF.clear();
              G_II.clear();

              B_hat_FH.clear();
              B_hat_IE.clear();
              B_hat_FE.clear();

              l_F.clear();

              // We fill the dofs relationship containers at the face level.
              // To do so, we first loop on the test space dofs.
              for (unsigned int i : fe_face_values_test.dof_indices())
                {
                  // Get the information on which element the dof is
                  const unsigned int current_element_test_i =
                    this->fe_test->system_to_base_index(i).first.first;

                  // Apply the different Robin boundary conditions if needed.
                  // The load vector and the B_hat matrix will have a
                  // contribution in addition to a modification of the Riesz
                  // map (G matrix) because of the energy norm that we want to
                  // minimize there.
                  if ((bc_type ==
                       BoundaryConditions::BoundaryType::silver_muller) ||
                      (bc_type ==
                       BoundaryConditions::BoundaryType::impedance_boundary) ||
                      (bc_type ==
                       BoundaryConditions::BoundaryType::waveguide_port))
                    {
                      if ((current_element_test_i == 0) ||
                          (current_element_test_i == 1))
                        {
                          l_F.emplace_back(i);
                        }

                      // Loop over the dofs test to fill the G_matrix dofs
                      // relationship
                      for (unsigned int j : fe_face_values_test.dof_indices())
                        {
                          const unsigned int current_element_test_j =
                            this->fe_test->system_to_base_index(j).first.first;

                          if (((current_element_test_i == 0) ||
                               (current_element_test_i == 1)) &&
                              ((current_element_test_j == 0) ||
                               (current_element_test_j == 1)))
                            {
                              G_FF.emplace_back(i, j);
                            }
                          if (((current_element_test_i == 0) ||
                               (current_element_test_i == 1)) &&
                              ((current_element_test_j == 2) ||
                               (current_element_test_j == 3)))
                            {
                              G_FI.emplace_back(i, j);
                            }
                          if (((current_element_test_i == 2) ||
                               (current_element_test_i == 3)) &&
                              ((current_element_test_j == 0) ||
                               (current_element_test_j == 1)))
                            {
                              G_IF.emplace_back(i, j);
                            }
                          if (((current_element_test_i == 2) ||
                               (current_element_test_i == 3)) &&
                              ((current_element_test_j == 2) ||
                               (current_element_test_j == 3)))
                            {
                              G_II.emplace_back(i, j);
                            }
                        }
                      // Loop over the dofs trial space to fill the B_hat
                      // matrix dofs relationship for the Robin boundary
                      // condition
                      for (unsigned int j :
                           fe_face_values_trial_skeleton.dof_indices())
                        {
                          const unsigned int current_element_trial_j =
                            this->fe_trial_skeleton->system_to_base_index(j)
                              .first.first;

                          if (((current_element_test_i == 0) ||
                               (current_element_test_i == 1)) &&
                              ((current_element_trial_j == 0) ||
                               (current_element_trial_j == 1)))
                            {
                              B_hat_FE.emplace_back(i, j);
                            }
                          if (((current_element_test_i == 2) ||
                               (current_element_test_i == 3)) &&
                              ((current_element_trial_j == 0) ||
                               (current_element_trial_j == 1)))
                            {
                              B_hat_IE.emplace_back(i, j);
                            }
                        }
                    }
                  else
                    {
                      // If not on a Robin B.C., assemble all the other
                      // relevant skeleton terms
                      for (unsigned int j :
                           fe_face_values_trial_skeleton.dof_indices())
                        {
                          const unsigned int current_element_trial_j =
                            this->fe_trial_skeleton->system_to_base_index(j)
                              .first.first;

                          if (((current_element_test_i == 0) ||
                               (current_element_test_i == 1)) &&
                              ((current_element_trial_j == 2) ||
                               (current_element_trial_j == 3)))
                            {
                              B_hat_FH.emplace_back(i, j);
                            }
                          if (((current_element_test_i == 2) ||
                               (current_element_test_i == 3)) &&
                              ((current_element_trial_j == 0) ||
                               (current_element_trial_j == 1)))
                            {
                              B_hat_IE.emplace_back(i, j);
                            }
                        }
                    }
                }

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

                  // As for the cell, we first loop over the test dofs to fill
                  // the face values containers
                  for (unsigned int i : fe_face_values_test.dof_indices())
                    {
                      F_face[i] =
                        fe_face_values_test[extractor_E_real].value(i,
                                                                    q_point) +
                        imag *
                          fe_face_values_test[extractor_E_imag].value(i,
                                                                      q_point);
                      F_face_conj[i] =
                        fe_face_values_test[extractor_E_real].value(i,
                                                                    q_point) -
                        imag *
                          fe_face_values_test[extractor_E_imag].value(i,
                                                                      q_point);

                      I_face_conj[i] =
                        fe_face_values_test[extractor_H_real].value(i,
                                                                    q_point) -
                        imag *
                          fe_face_values_test[extractor_H_imag].value(i,
                                                                      q_point);

                      n_cross_I_face[i] = cross_product_3d(
                        normal,
                        fe_face_values_test[extractor_H_real].value(i,
                                                                    q_point) +
                          imag * fe_face_values_test[extractor_H_imag].value(
                                   i, q_point));
                      n_cross_I_face_conj[i] = cross_product_3d(
                        normal,
                        fe_face_values_test[extractor_H_real].value(i,
                                                                    q_point) -
                          imag * fe_face_values_test[extractor_H_imag].value(
                                   i, q_point));
                    }

                  // Then, similarly we loop over the trial dofs to fill the
                  // face values containers. Note that to be in
                  // H^-1/2(curl), the fields needs to have the tangential
                  // property mapping (n x (E x n)) which effectively extract
                  // the tangential component of the field at the face. So
                  // here we apply this operation using the map_H12 function
                  // that we defined earlier. Stricly speeking, nx(E_parallel)
                  // = n x E, and we would not need to use the map_H12
                  // function, but we keep it for consistency.
                  for (unsigned int i :
                       fe_face_values_trial_skeleton.dof_indices())
                    {
                      E_hat[i] = map_H12(
                        fe_face_values_trial_skeleton[extractor_E_real].value(
                          i, q_point) +
                          imag * fe_face_values_trial_skeleton[extractor_E_imag]
                                   .value(i, q_point),
                        normal);

                      n_cross_E_hat[i] = cross_product_3d(
                        normal,
                        map_H12(
                          fe_face_values_trial_skeleton[extractor_E_real].value(
                            i, q_point) +
                            imag *
                              fe_face_values_trial_skeleton[extractor_E_imag]
                                .value(i, q_point),
                          normal));

                      n_cross_H_hat[i] = cross_product_3d(
                        normal,
                        map_H12(
                          fe_face_values_trial_skeleton[extractor_H_real].value(
                            i, q_point) +
                            imag *
                              fe_face_values_trial_skeleton[extractor_H_imag]
                                .value(i, q_point),
                          normal));
                    }

                  // Here we apply the excitation at the relevant boundary.
                  if (bc_type ==
                      BoundaryConditions::BoundaryType::silver_muller)
                    {
                      boundary_surface_admittance =
                        sqrt(effective_electric_permittivity /
                             effective_magnetic_permeability);
                      conj_boundary_surface_admittance =
                        std::conj(boundary_surface_admittance);
                      g_inc = 0.;
                    }
                  if (bc_type ==
                      BoundaryConditions::BoundaryType::impedance_boundary)
                    {
                      unsigned int face_id = face->boundary_id();

                      boundary_surface_admittance =
                        this->simulation_parameters
                          .boundary_conditions_time_harmonic_electromagnetics
                          .surface_admittance_real.at(face_id)
                          ->value(position) +
                        imag *
                          this->simulation_parameters
                            .boundary_conditions_time_harmonic_electromagnetics
                            .surface_admittance_imag.at(face_id)
                            ->value(position);

                      conj_boundary_surface_admittance =
                        std::conj(boundary_surface_admittance);

                      // Get the incident electromagnetic field at this face
                      g_inc[0] =
                        this->simulation_parameters
                          .boundary_conditions_time_harmonic_electromagnetics
                          .excitation_x_real.at(face_id)
                          ->value(position) +
                        imag *
                          this->simulation_parameters
                            .boundary_conditions_time_harmonic_electromagnetics
                            .excitation_x_imag.at(face_id)
                            ->value(position);

                      g_inc[1] =
                        this->simulation_parameters
                          .boundary_conditions_time_harmonic_electromagnetics
                          .excitation_y_real.at(face_id)
                          ->value(position) +
                        imag *
                          this->simulation_parameters
                            .boundary_conditions_time_harmonic_electromagnetics
                            .excitation_y_imag.at(face_id)
                            ->value(position);

                      g_inc[2] =
                        this->simulation_parameters
                          .boundary_conditions_time_harmonic_electromagnetics
                          .excitation_z_real.at(face_id)
                          ->value(position) +
                        imag *
                          this->simulation_parameters
                            .boundary_conditions_time_harmonic_electromagnetics
                            .excitation_z_imag.at(face_id)
                            ->value(position);
                    }
                  if (bc_type ==
                      BoundaryConditions::BoundaryType::waveguide_port)
                    {
                      unsigned int boundary_index = std::distance(
                        time_harmonic_maxwell_parameters.waveguide_boundary_ids
                          .begin(),
                        std::ranges::find(time_harmonic_maxwell_parameters
                                            .waveguide_boundary_ids,
                                          face->boundary_id()));

                      std::tie(g_inc, boundary_surface_admittance) =
                        compute_waveguide_port_excitation(
                          position,
                          normal,
                          effective_electric_permittivity,
                          effective_magnetic_permeability,
                          boundary_index);

                      conj_boundary_surface_admittance =
                        std::conj(boundary_surface_admittance);
                    }

                  // Now we loop on each relationship container to assemble
                  // the relevant matrices.
                  for (const auto &[i, j] : G_FF)
                    {
                      G_matrix(i, j) +=
                        (conj_boundary_surface_admittance * F_face[j] *
                         boundary_surface_admittance * F_face_conj[i] *
                         JxW_face)
                          .real();
                    }

                  for (const auto &[i, j] : G_FI)
                    {
                      G_matrix(i, j) +=
                        (n_cross_I_face[j] * boundary_surface_admittance *
                         F_face_conj[i] * JxW_face)
                          .real();
                    }

                  for (const auto &[i, j] : G_IF)
                    {
                      G_matrix(i, j) +=
                        (conj_boundary_surface_admittance * F_face[j] *
                         n_cross_I_face_conj[i] * JxW_face)
                          .real();
                    }

                  for (const auto &[i, j] : G_II)
                    {
                      G_matrix(i, j) +=
                        (n_cross_I_face[j] * n_cross_I_face_conj[i] * JxW_face)
                          .real();
                    }

                  for (const auto &[i, j] : B_hat_FH)
                    {
                      B_hat_matrix(i, j) +=
                        (n_cross_H_hat[j] * F_face_conj[i] * JxW_face).real();
                    }

                  for (const auto &[i, j] : B_hat_IE)
                    {
                      B_hat_matrix(i, j) +=
                        (n_cross_E_hat[j] * I_face_conj[i] * JxW_face).real();
                    }

                  for (const auto &[i, j] : B_hat_FE)
                    {
                      B_hat_matrix(i, j) -=
                        (boundary_surface_admittance * E_hat[j] *
                         F_face_conj[i] * JxW_face)
                          .real();
                    }

                  for (const auto &i : l_F)
                    {
                      l_vector[i] -= (g_inc * F_face_conj[i] * JxW_face).real();
                    }
                }
            } // End of face loop

          // Finally, after having assembled all the matrices and vectors, we
          // build the condensed version of the system.

          // We only need the inverse of the Gram matrix $G$, so we
          // invert it.
          G_matrix.invert();

          // We construct $M_4 = B^\dagger G^{-1}$ and $M_5 = \hat{B}^\dagger
          // G^{-1}$ with it:
          B_matrix.Tmmult(M4_matrix, G_matrix);
          B_hat_matrix.Tmmult(M5_matrix, G_matrix);

          // Then using $M_4$ we compute the condensed matrix $M_1 = B^\dagger
          // G^{-1} B$ and $M_2 = B^\dagger G^{-1} \hat{B}$:
          M4_matrix.mmult(M1_matrix, B_matrix);
          M4_matrix.mmult(M2_matrix, B_hat_matrix);

          // We also compute the matrix $M_3 = \hat{B}^\dagger G^{-1} \hat{B}$
          M5_matrix.mmult(M3_matrix, B_hat_matrix);

          // Finally, as for the $G$ matrix, we invert the $M_1$
          // matrix:
          M1_matrix.invert();

          // Now, we have to compute the local matrix and the local RHS for the
          // condensed system.

          // The cell matrix is obtained with the formula $(M_3 -
          // M_2^\dagger M_1^{-1} M_2)$:
          M2_matrix.Tmmult(tmp_matrix_M2M1, M1_matrix);
          tmp_matrix_M2M1.mmult(tmp_matrix_M2M1M2, M2_matrix);
          tmp_matrix_M2M1M2.add(-1.0, M3_matrix);
          tmp_matrix_M2M1M2 *= -1.0;
          // This line is used to convert the LAPACK matrix to a full
          // matrix so we can perform the distribution to the global
          // system below.
          cell_matrix = tmp_matrix_M2M1M2;

          // Then we compute the cell RHS using $(M_5 -
          // M_2^\dagger M_1^{-1} M_4)l -
          // G$.
          tmp_matrix_M2M1.mmult(tmp_matrix_M2M1M4, M4_matrix);
          M5_matrix.add(-1.0, tmp_matrix_M2M1M4);
          M5_matrix.vmult(cell_skeleton_rhs, l_vector);

          // Map to global matrix
          cell_skeleton->get_dof_indices(local_dof_indices);
          this->nonzero_constraints.distribute_local_to_global(
            cell_matrix,
            cell_skeleton_rhs,
            local_dof_indices,
            this->system_matrix,
            this->system_rhs);
        }
    }

  // After the loop over the cells, we finalize the assembly by compressing
  // the vectors because of the MPI parallelization.
  this->system_matrix.compress(VectorOperation::add);
  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::assemble_system_rhs()
{
  // At the moment everything is done in assemble_system_matrix(). In the
  // future, we want to add the possibility of assembling both the matrix and
  // rhs in the same function and loop for efficiency.
}

template <>
void
TimeHarmonicMaxwell<2>::reconstruct_interior_solution()
{
  // As for the assembly, we need to perform operations like curl and cross
  // product, that are completely different in 2D in comparison to 3D. At the
  // moment, we do not support 2D time-harmonic Maxwell so this function throws
  // an error.
  AssertThrow(false, TimeHarmonicMaxwellDimensionNotSupported(2));
}


template <>
void
TimeHarmonicMaxwell<3>::reconstruct_interior_solution()
{
  // Constexpr values and used in the assembly. Since we are in the specialized
  // 3D function we define dim = 3 here. This makes easier to read the code
  // below to see what are templated in dim and what are not.
  static constexpr double               PI = numbers::PI;
  static constexpr std::complex<double> imag{0., 1.};
  static constexpr int                  dim = 3;

  // Get properties manager and define model physical properties
  const auto &properties_manager =
    this->simulation_parameters.physical_properties_manager;
  std::complex<double> effective_electric_permittivity;
  std::complex<double> effective_magnetic_permeability;
  unsigned int         material_id;

  /// Excitation properties
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters =
    this->simulation_parameters.multiphysics.time_harmonic_maxwell_parameters;
  const double omega =
    2.0 * PI * time_harmonic_maxwell_parameters.electromagnetic_frequency;

  // We then create the corresponding FEValues and FEFaceValues objects. Note
  // that only the test space needs gradients because of the ultraweak
  // formulation. Similarly, because everything is on the same triangulation,
  // we only need to update the quadrature points and JxW values in one of the
  // spaces. Here we choose the trial space.
  const unsigned int n_q_points      = this->cell_quadrature->size();
  const unsigned int n_face_q_points = this->face_quadrature->size();

  FEValues<dim> fe_values_trial_interior(*this->mapping,
                                         *this->fe_trial_interior,
                                         *this->cell_quadrature,
                                         update_values |
                                           update_quadrature_points |
                                           update_JxW_values);

  FEValues<dim> fe_values_test(*this->mapping,
                               *this->fe_test,
                               *this->cell_quadrature,
                               update_values | update_gradients);

  FEFaceValues<dim> fe_face_values_trial_skeleton(*this->mapping,
                                                  *this->fe_trial_skeleton,
                                                  *this->face_quadrature,
                                                  update_values |
                                                    update_quadrature_points |
                                                    update_normal_vectors |
                                                    update_JxW_values);

  FEFaceValues<dim> fe_face_values_test(*this->mapping,
                                        *this->fe_test,
                                        *this->face_quadrature,
                                        update_values);

  // We also create all the relevant matrices and vector to build the DPG
  // system. To do so we first need the number of dofs per cell for each of
  // the finite element spaces.
  const unsigned int dofs_per_cell_test = this->fe_test->n_dofs_per_cell();
  const unsigned int dofs_per_cell_trial_interior =
    this->fe_trial_interior->n_dofs_per_cell();
  const unsigned int dofs_per_cell_trial_skeleton =
    this->fe_trial_skeleton->n_dofs_per_cell();

  // To avoid unecessary computations, we will precompute the shape functions
  // of all our spaces once and then use those precomputed values to assemble
  // the local matrices. Consequently, we need to create containers to store
  // these values. Since those are complex-valued functions, we use two
  // different containers for the real and imaginary parts because the
  // conjugate of a complex tensor is not implemented in deal.II.
  std::vector<Tensor<1, dim, std::complex<double>>> F(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> F_conj(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> I(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> I_conj(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> curl_F(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> curl_F_conj(
    dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> curl_I(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> curl_I_conj(
    dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> F_face(dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> F_face_conj(
    dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> I_face_conj(
    dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> n_cross_I_face(
    dofs_per_cell_test);
  std::vector<Tensor<1, dim, std::complex<double>>> n_cross_I_face_conj(
    dofs_per_cell_test);

  std::vector<Tensor<1, dim, std::complex<double>>> E(
    dofs_per_cell_trial_interior);
  std::vector<Tensor<1, dim, std::complex<double>>> H(
    dofs_per_cell_trial_interior);
  std::vector<Tensor<1, dim, std::complex<double>>> E_hat(
    dofs_per_cell_trial_skeleton);
  std::vector<Tensor<1, dim, std::complex<double>>> n_cross_E_hat(
    dofs_per_cell_trial_skeleton);
  std::vector<Tensor<1, dim, std::complex<double>>> n_cross_H_hat(
    dofs_per_cell_trial_skeleton);

  // Also, to avoid multiple "if" calls during the assembly to understand where
  // each term needs to be assembled in the DPG global matrix, we will use
  // containers to store dofs relationships for each of the local matrices.

  // G matrix stands for the Riesz map and needs the relationships between test
  // functions F and I
  std::vector<std::pair<unsigned int, unsigned int>> G_FF;
  std::vector<std::pair<unsigned int, unsigned int>> G_FI;
  std::vector<std::pair<unsigned int, unsigned int>> G_IF;
  std::vector<std::pair<unsigned int, unsigned int>> G_II;

  // B matrix stands for the bilinear form and needs the relationships between
  // interior trial space (E and H) and the test (F and I)
  std::vector<std::pair<unsigned int, unsigned int>> B_FE;
  std::vector<std::pair<unsigned int, unsigned int>> B_IE;
  std::vector<std::pair<unsigned int, unsigned int>> B_FH;
  std::vector<std::pair<unsigned int, unsigned int>> B_IH;

  // B_hat matrix stands for the bilinear form, but on the skeleton, and needs
  // the relationships between skeleton trial space (E_hat, H_hat) and the test
  // (F and I). Note that the B_hat_FE is required for the Robin boundary
  // condition that we chose to apply on the electric field instead of the
  // magnetic field. This choice is arbitrary but it cannot be applied on both
  // fields at the same time.
  std::vector<std::pair<unsigned int, unsigned int>> B_hat_IE;
  std::vector<std::pair<unsigned int, unsigned int>> B_hat_FH;
  std::vector<std::pair<unsigned int, unsigned int>> B_hat_FE;

  // l vector stands for the linear form of the problem and needs the
  // relationships between test functions. Note that in its simplest form, the
  // time-harmonic Maxwell equation does not have a magnetic source term, so
  // there is no contribution from the I test functions.
  std::vector<unsigned int> l_F;

  // Reserve memory to avoid reallocations of each relationship vectors
  G_FF.reserve(dofs_per_cell_test * dofs_per_cell_test);
  G_FI.reserve(dofs_per_cell_test * dofs_per_cell_test);
  G_IF.reserve(dofs_per_cell_test * dofs_per_cell_test);
  G_II.reserve(dofs_per_cell_test * dofs_per_cell_test);

  B_FE.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);
  B_IE.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);
  B_FH.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);
  B_IH.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);

  B_hat_IE.reserve(dofs_per_cell_trial_skeleton * dofs_per_cell_test);
  B_hat_FH.reserve(dofs_per_cell_trial_skeleton * dofs_per_cell_test);
  B_hat_FE.reserve(dofs_per_cell_trial_skeleton * dofs_per_cell_test);

  l_F.reserve(dofs_per_cell_test);

  // Here we create the DPG local matrices and vector used for the assembly
  // before condensation.
  LAPACKFullMatrix<double> G_matrix(dofs_per_cell_test, dofs_per_cell_test);

  LAPACKFullMatrix<double> B_matrix(dofs_per_cell_test,
                                    dofs_per_cell_trial_interior);

  LAPACKFullMatrix<double> B_hat_matrix(dofs_per_cell_test,
                                        dofs_per_cell_trial_skeleton);

  Vector<double> l_vector(dofs_per_cell_test);

  // We create the condensation matrices which are defined as :
  // $M_1 = B^\dagger G^{-1}B$;
  // $M_2 = B^\dagger G^{-1}\hat{B}$;
  // $M_3 = \hat{B}^\dagger G^{-1}\hat{B}$;
  // $M_4 = B^\dagger G^{-1}$;
  // $M_5 = \hat{B}^\dagger G^{-1}$.
  LAPACKFullMatrix<double> M1_matrix(dofs_per_cell_trial_interior,
                                     dofs_per_cell_trial_interior);
  LAPACKFullMatrix<double> M2_matrix(dofs_per_cell_trial_interior,
                                     dofs_per_cell_trial_skeleton);
  LAPACKFullMatrix<double> M3_matrix(dofs_per_cell_trial_skeleton,
                                     dofs_per_cell_trial_skeleton);
  LAPACKFullMatrix<double> M4_matrix(dofs_per_cell_trial_interior,
                                     dofs_per_cell_test);
  LAPACKFullMatrix<double> M5_matrix(dofs_per_cell_trial_skeleton,
                                     dofs_per_cell_test);

  // We create the cell matrix and the RHS that will be distributed in the
  // full system after the assembly along with the index’s mapping.
  FullMatrix<double> cell_matrix(dofs_per_cell_trial_skeleton,
                                 dofs_per_cell_trial_skeleton);
  Vector<double>     cell_skeleton_rhs(dofs_per_cell_trial_skeleton);

  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell_trial_skeleton);

  // We also create objects used for the various Robin boundary conditions.
  BoundaryConditions::BoundaryType     bc_type;
  Tensor<1, dim, std::complex<double>> g_inc;
  std::complex<double>                 boundary_surface_admittance;
  std::complex<double>                 conj_boundary_surface_admittance;

  // L2 norm of the residual for the dpg error indicator
  double residual_L2_norm = 0.0;

  // The above declarations are the same as the ones in the assembly of the
  // skeleton system. Below we add new ones that are specific to the
  // reconstruction of the interior solution.

  MPI_Comm mpi_communicator = this->triangulation->get_mpi_communicator();

  // We initialize vectors to store the locally owned solution and the error
  // indicator.
  GlobalVectorType locally_owned_solution_interior(
    this->locally_owned_dofs_trial_interior, mpi_communicator);
  GlobalVectorType locally_owned_error_indicator(this->locally_owned_dofs_test,
                                                 mpi_communicator);

  // Temporary vector used when reconstructing the interior solution :
  // $tmp_{interior} = M_2 * x_{skeleton}$
  Vector<double> tmp_vector_interior(dofs_per_cell_trial_interior);

  // Temporary vector used when computing the error indicator :
  // $tmp_{error\_indicator} = B * x_{interior}$
  Vector<double> tmp_vector_error_indicator(dofs_per_cell_test);

  // Finally, when reconstructing the interior solution from the skeleton, we
  // require additional vectors that we allocate here.
  Vector<double> cell_interior_rhs(dofs_per_cell_trial_interior);
  Vector<double> cell_interior_solution(dofs_per_cell_trial_interior);
  Vector<double> cell_skeleton_solution(dofs_per_cell_trial_skeleton);
  Vector<double> cell_residual(dofs_per_cell_test);

  // The assembly on each cell is the same as previously done, so we first loop
  // over the cells of the triangulation.
  for (const auto &cell :
       this->dof_handler_trial_interior->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // We update the material properties for the current cell and define
          // some constants that will be used during the assembly.
          material_id = cell->material_id();
          update_material_properties(properties_manager,
                                     effective_electric_permittivity,
                                     effective_magnetic_permeability,
                                     material_id);
          const std::complex<double> iweffective_magnetic_permeability =
            imag * omega * effective_magnetic_permeability;
          const std::complex<double> conj_iweffective_magnetic_permeability =
            std::conj(iweffective_magnetic_permeability);
          const std::complex<double> iweps_r =
            imag * omega * effective_electric_permittivity;
          const std::complex<double> conj_iweps_r = std::conj(iweps_r);

          // We reinitialize the FEValues objects to the current cell.
          fe_values_trial_interior.reinit(cell);

          // We will also need to reinitialize the FEValues for the test
          // space and make sure that is the same cell as the one used for the
          // trial space.
          const typename DoFHandler<dim>::active_cell_iterator cell_test =
            cell->as_dof_handler_iterator(*this->dof_handler_test);
          fe_values_test.reinit(cell_test);

          // Similarly, we reinitialize the FEValues for the trial space on
          // the skeleton, but this will not be used before we also loop on
          // the cells faces.
          const typename DoFHandler<dim>::active_cell_iterator cell_skeleton =
            cell->as_dof_handler_iterator(*this->dof_handler_trial_skeleton);

          // We then reinitialize all the matrices that we are aggregating
          // information for each cell.
          G_matrix     = 0;
          B_matrix     = 0;
          B_hat_matrix = 0;
          l_vector     = 0;

          // We also need to reinitialize the $M_1$ condensation matrix
          // between each iteration on cell to get rid of its inverse status.
          M1_matrix = 0;

          // Here we reset the dofs relationships containers. These containers
          // are used to store all the relationships between the dof i
          // dof j according to which space they belong to. Indeed, when
          // looping over all the test and trial dofs, the terms  we need to
          // compute depend on the specific combination of spaces
          // involved. The following containers are therefore used to avoid
          // multiple "if" statements inside the dofs loops and are filled with
          // all the relevant dofs pairs that we need for each of the
          // different terms.

          // For example, in the ultraweak form of Maxwell equation we have a
          // term (E, curl(I)) in the interior so the matrix B has a term
          // (curl(I_i), E_j), so we want to only compute this term for the
          // dofs i that are in the test space I and the dofs j that are in
          // the trial space for E. Therefore, when looping on all the
          // interior and test dofs in a cell we add the pairs of dofs that
          // are in those spaces to the container B_IE. We do this for all the
          // different terms of the formulation.

          // Finally, we can loop on all the pairs that we have assigned in
          // each vector container to compute the desired terms.
          G_FF.clear();
          G_FI.clear();
          G_IF.clear();
          G_II.clear();

          B_FE.clear();
          B_IE.clear();
          B_FH.clear();
          B_IH.clear();

          l_F.clear();

          // We fill the dofs relationship containers at the cell level. To do
          // so, we first loop on the test space dofs.
          for (unsigned int i : fe_values_test.dof_indices())
            {
              // Get the information on which element the dof is
              const unsigned int current_element_test_i =
                this->fe_test->system_to_base_index(i).first.first;

              // Fill the load vector relationship
              if ((current_element_test_i == 0) ||
                  (current_element_test_i == 1))
                {
                  l_F.emplace_back(i);
                }

              // Loop over the test dofs a second time to fill the dofs
              // relationship for the G matrix (Riesz map)
              for (unsigned int j : fe_values_test.dof_indices())
                {
                  const unsigned int current_element_test_j =
                    this->fe_test->system_to_base_index(j).first.first;
                  if (((current_element_test_i == 0) ||
                       (current_element_test_i == 1)) &&
                      ((current_element_test_j == 0) ||
                       (current_element_test_j == 1)))
                    {
                      G_FF.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 0) ||
                       (current_element_test_i == 1)) &&
                      ((current_element_test_j == 2) ||
                       (current_element_test_j == 3)))
                    {
                      G_FI.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 2) ||
                       (current_element_test_i == 3)) &&
                      ((current_element_test_j == 0) ||
                       (current_element_test_j == 1)))
                    {
                      G_IF.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 2) ||
                       (current_element_test_i == 3)) &&
                      ((current_element_test_j == 2) ||
                       (current_element_test_j == 3)))
                    {
                      G_II.emplace_back(i, j);
                    }
                }

              // Then we loop over the trial dofs space to fill the dofs
              // relationship for the B matrix (bilinear form)
              for (unsigned int j : fe_values_trial_interior.dof_indices())
                {
                  const unsigned int current_element_trial_j =
                    this->fe_trial_interior->system_to_base_index(j)
                      .first.first;

                  if (((current_element_test_i == 0) ||
                       (current_element_test_i == 1)) &&
                      ((current_element_trial_j == 0) ||
                       (current_element_trial_j == 1)))
                    {
                      B_FE.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 0) ||
                       (current_element_test_i == 1)) &&
                      ((current_element_trial_j == 2) ||
                       (current_element_trial_j == 3)))
                    {
                      B_FH.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 2) ||
                       (current_element_test_i == 3)) &&
                      ((current_element_trial_j == 0) ||
                       (current_element_trial_j == 1)))
                    {
                      B_IE.emplace_back(i, j);
                    }
                  if (((current_element_test_i == 2) ||
                       (current_element_test_i == 3)) &&
                      ((current_element_trial_j == 2) ||
                       (current_element_trial_j == 3)))
                    {
                      B_IH.emplace_back(i, j);
                    }
                }
            }

          // Now we loop over all quadrature points of the cell
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              // To avoid unnecessary computation, we fill the shape values
              // containers for the real and imaginary parts of the electric
              // and magnetic fields and the dofs relationship at the current
              // quadrature point.
              const double &JxW = fe_values_trial_interior.JxW(q_point);

              for (unsigned int i : fe_values_test.dof_indices())
                {
                  F[i] =
                    fe_values_test[extractor_E_real].value(i, q_point) +
                    imag * fe_values_test[extractor_E_imag].value(i, q_point);
                  F_conj[i] =
                    fe_values_test[extractor_E_real].value(i, q_point) -
                    imag * fe_values_test[extractor_E_imag].value(i, q_point);

                  curl_F[i] =
                    fe_values_test[extractor_E_real].curl(i, q_point) +
                    imag * fe_values_test[extractor_E_imag].curl(i, q_point);
                  curl_F_conj[i] =
                    fe_values_test[extractor_E_real].curl(i, q_point) -
                    imag * fe_values_test[extractor_E_imag].curl(i, q_point);

                  I[i] =
                    fe_values_test[extractor_H_real].value(i, q_point) +
                    imag * fe_values_test[extractor_H_imag].value(i, q_point);
                  I_conj[i] =
                    fe_values_test[extractor_H_real].value(i, q_point) -
                    imag * fe_values_test[extractor_H_imag].value(i, q_point);

                  curl_I[i] =
                    fe_values_test[extractor_H_real].curl(i, q_point) +
                    imag * fe_values_test[extractor_H_imag].curl(i, q_point);
                  curl_I_conj[i] =
                    fe_values_test[extractor_H_real].curl(i, q_point) -
                    imag * fe_values_test[extractor_H_imag].curl(i, q_point);
                }

              for (unsigned int i : fe_values_trial_interior.dof_indices())
                {
                  E[i] =
                    fe_values_trial_interior[extractor_E_real].value(i,
                                                                     q_point) +
                    imag *
                      fe_values_trial_interior[extractor_E_imag].value(i,
                                                                       q_point);
                  H[i] =
                    fe_values_trial_interior[extractor_H_real].value(i,
                                                                     q_point) +
                    imag *
                      fe_values_trial_interior[extractor_H_imag].value(i,
                                                                       q_point);
                }

              // Now we loop on each relationship container to assemble the
              // relevant matrices
              for (const auto &[i, j] : G_FF)
                {
                  G_matrix(i, j) +=
                    (((F[j] * F_conj[i]) + (curl_F[j] * curl_F_conj[i]) +
                      (conj_iweps_r * F[j] * iweps_r * F_conj[i])) *
                     JxW)
                      .real();
                }

              for (const auto &[i, j] : G_FI)
                {
                  G_matrix(i, j) += (((curl_I[j] * iweps_r * F_conj[i]) -
                                      (conj_iweffective_magnetic_permeability *
                                       I[j] * curl_F_conj[i])) *
                                     JxW)
                                      .real();
                }

              for (const auto &[i, j] : G_IF)
                {
                  G_matrix(i, j) +=
                    (((conj_iweps_r * F[j] * curl_I_conj[i]) -
                      (curl_F[j] * iweffective_magnetic_permeability *
                       I_conj[i])) *
                     JxW)
                      .real();
                }

              for (const auto &[i, j] : G_II)
                {
                  G_matrix(i, j) +=
                    (((I[j] * I_conj[i]) + (curl_I[j] * curl_I_conj[i]) +
                      (conj_iweffective_magnetic_permeability * I[j] *
                       iweffective_magnetic_permeability * I_conj[i])) *
                     JxW)
                      .real();
                }

              for (const auto &[i, j] : B_FE)
                {
                  B_matrix(i, j) += (iweps_r * E[j] * F_conj[i] * JxW).real();
                }

              for (const auto &[i, j] : B_FH)
                {
                  B_matrix(i, j) += (H[j] * curl_F_conj[i] * JxW).real();
                }

              for (const auto &[i, j] : B_IE)
                {
                  B_matrix(i, j) += (E[j] * curl_I_conj[i] * JxW).real();
                }

              for (const auto &[i, j] : B_IH)
                {
                  B_matrix(i, j) -=
                    (iweffective_magnetic_permeability * H[j] * I_conj[i] * JxW)
                      .real();
                }

              for (const auto &i : l_F)
                {
                  l_vector[i] += 0.0;
                }
            }

          // We now build the skeleton terms. Similarly, we choose to loop on
          // the skeleton trial space faces.
          for (const auto &face : cell_skeleton->face_iterators())
            {
              // We reinitialize the FEFaceValues objects to the current
              // faces.
              fe_face_values_test.reinit(cell_test, face);
              fe_face_values_trial_skeleton.reinit(cell_skeleton, face);

              // Get the boundary condition type on the current face
              bc_type = BoundaryConditions::BoundaryType::none;

              if (face->at_boundary())
                bc_type =
                  this->simulation_parameters
                    .boundary_conditions_time_harmonic_electromagnetics.type.at(
                      face->boundary_id());

              // Reset the face dofs relationships
              G_FF.clear();
              G_FI.clear();
              G_IF.clear();
              G_II.clear();

              B_hat_FH.clear();
              B_hat_IE.clear();
              B_hat_FE.clear();

              l_F.clear();

              // We fill the dofs relationship containers at the face level.
              // To do so, we first loop on the test space dofs.
              for (unsigned int i : fe_face_values_test.dof_indices())
                {
                  // Get the information on which element the dof is
                  const unsigned int current_element_test_i =
                    this->fe_test->system_to_base_index(i).first.first;

                  // Apply the different Robin boundary conditions if needed.
                  // The load vector and the B_hat matrix will have a
                  // contribution in addition to a modification of the Riesz
                  // map (G matrix) because of the energy norm that we want to
                  // minimize there.
                  if ((bc_type ==
                       BoundaryConditions::BoundaryType::silver_muller) ||
                      (bc_type ==
                       BoundaryConditions::BoundaryType::impedance_boundary) ||
                      (bc_type ==
                       BoundaryConditions::BoundaryType::waveguide_port))
                    {
                      if ((current_element_test_i == 0) ||
                          (current_element_test_i == 1))
                        {
                          l_F.emplace_back(i);
                        }

                      // Loop over the dofs test to fill the G_matrix dofs
                      // relationship
                      for (unsigned int j : fe_face_values_test.dof_indices())
                        {
                          const unsigned int current_element_test_j =
                            this->fe_test->system_to_base_index(j).first.first;

                          if (((current_element_test_i == 0) ||
                               (current_element_test_i == 1)) &&
                              ((current_element_test_j == 0) ||
                               (current_element_test_j == 1)))
                            {
                              G_FF.emplace_back(i, j);
                            }
                          if (((current_element_test_i == 0) ||
                               (current_element_test_i == 1)) &&
                              ((current_element_test_j == 2) ||
                               (current_element_test_j == 3)))
                            {
                              G_FI.emplace_back(i, j);
                            }
                          if (((current_element_test_i == 2) ||
                               (current_element_test_i == 3)) &&
                              ((current_element_test_j == 0) ||
                               (current_element_test_j == 1)))
                            {
                              G_IF.emplace_back(i, j);
                            }
                          if (((current_element_test_i == 2) ||
                               (current_element_test_i == 3)) &&
                              ((current_element_test_j == 2) ||
                               (current_element_test_j == 3)))
                            {
                              G_II.emplace_back(i, j);
                            }
                        }
                      // Loop over the dofs trial space to fill the B_hat
                      // matrix dofs relationship for the Robin boundary
                      // condition
                      for (unsigned int j :
                           fe_face_values_trial_skeleton.dof_indices())
                        {
                          const unsigned int current_element_trial_j =
                            this->fe_trial_skeleton->system_to_base_index(j)
                              .first.first;

                          if (((current_element_test_i == 0) ||
                               (current_element_test_i == 1)) &&
                              ((current_element_trial_j == 0) ||
                               (current_element_trial_j == 1)))
                            {
                              B_hat_FE.emplace_back(i, j);
                            }
                          if (((current_element_test_i == 2) ||
                               (current_element_test_i == 3)) &&
                              ((current_element_trial_j == 0) ||
                               (current_element_trial_j == 1)))
                            {
                              B_hat_IE.emplace_back(i, j);
                            }
                        }
                    }
                  else
                    {
                      // If not on a Robin B.C., assemble all the other
                      // relevant skeleton terms
                      for (unsigned int j :
                           fe_face_values_trial_skeleton.dof_indices())
                        {
                          const unsigned int current_element_trial_j =
                            this->fe_trial_skeleton->system_to_base_index(j)
                              .first.first;

                          if (((current_element_test_i == 0) ||
                               (current_element_test_i == 1)) &&
                              ((current_element_trial_j == 2) ||
                               (current_element_trial_j == 3)))
                            {
                              B_hat_FH.emplace_back(i, j);
                            }
                          if (((current_element_test_i == 2) ||
                               (current_element_test_i == 3)) &&
                              ((current_element_trial_j == 0) ||
                               (current_element_trial_j == 1)))
                            {
                              B_hat_IE.emplace_back(i, j);
                            }
                        }
                    }
                }

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

                  // As for the cell, we first loop over the test dofs to fill
                  // the face values containers
                  for (unsigned int i : fe_face_values_test.dof_indices())
                    {
                      F_face[i] =
                        fe_face_values_test[extractor_E_real].value(i,
                                                                    q_point) +
                        imag *
                          fe_face_values_test[extractor_E_imag].value(i,
                                                                      q_point);
                      F_face_conj[i] =
                        fe_face_values_test[extractor_E_real].value(i,
                                                                    q_point) -
                        imag *
                          fe_face_values_test[extractor_E_imag].value(i,
                                                                      q_point);

                      I_face_conj[i] =
                        fe_face_values_test[extractor_H_real].value(i,
                                                                    q_point) -
                        imag *
                          fe_face_values_test[extractor_H_imag].value(i,
                                                                      q_point);

                      n_cross_I_face[i] = cross_product_3d(
                        normal,
                        fe_face_values_test[extractor_H_real].value(i,
                                                                    q_point) +
                          imag * fe_face_values_test[extractor_H_imag].value(
                                   i, q_point));
                      n_cross_I_face_conj[i] = cross_product_3d(
                        normal,
                        fe_face_values_test[extractor_H_real].value(i,
                                                                    q_point) -
                          imag * fe_face_values_test[extractor_H_imag].value(
                                   i, q_point));
                    }

                  // Then, similarly we loop over the trial dofs to fill the
                  // face values containers. Note that to be in
                  // H^-1/2(curl), the fields needs to have the tangential
                  // property mapping (n x (E x n)) which effectively extracts
                  // the tangential component of the field at the face. So
                  // here we apply this operation using the map_H12 function
                  // that we defined earlier. Stricly speeking, nx(E_parallel)
                  // = n x E, and we would not need to use the map_H12
                  // function, but we keep it for consistency.
                  for (unsigned int i :
                       fe_face_values_trial_skeleton.dof_indices())
                    {
                      E_hat[i] = map_H12(
                        fe_face_values_trial_skeleton[extractor_E_real].value(
                          i, q_point) +
                          imag * fe_face_values_trial_skeleton[extractor_E_imag]
                                   .value(i, q_point),
                        normal);

                      n_cross_E_hat[i] = cross_product_3d(
                        normal,
                        map_H12(
                          fe_face_values_trial_skeleton[extractor_E_real].value(
                            i, q_point) +
                            imag *
                              fe_face_values_trial_skeleton[extractor_E_imag]
                                .value(i, q_point),
                          normal));

                      n_cross_H_hat[i] = cross_product_3d(
                        normal,
                        map_H12(
                          fe_face_values_trial_skeleton[extractor_H_real].value(
                            i, q_point) +
                            imag *
                              fe_face_values_trial_skeleton[extractor_H_imag]
                                .value(i, q_point),
                          normal));
                    }

                  // Here we apply the excitation at the relevant boundary.
                  if (bc_type ==
                      BoundaryConditions::BoundaryType::silver_muller)
                    {
                      boundary_surface_admittance =
                        sqrt(effective_electric_permittivity /
                             effective_magnetic_permeability);
                      conj_boundary_surface_admittance =
                        std::conj(boundary_surface_admittance);
                      g_inc = 0.;
                    }
                  if (bc_type ==
                      BoundaryConditions::BoundaryType::impedance_boundary)
                    {
                      unsigned int face_id = face->boundary_id();

                      boundary_surface_admittance =
                        this->simulation_parameters
                          .boundary_conditions_time_harmonic_electromagnetics
                          .surface_admittance_real.at(face_id)
                          ->value(position) +
                        imag *
                          this->simulation_parameters
                            .boundary_conditions_time_harmonic_electromagnetics
                            .surface_admittance_imag.at(face_id)
                            ->value(position);

                      conj_boundary_surface_admittance =
                        std::conj(boundary_surface_admittance);

                      // Get the incident electromagnetic field at this face
                      g_inc[0] =
                        this->simulation_parameters
                          .boundary_conditions_time_harmonic_electromagnetics
                          .excitation_x_real.at(face_id)
                          ->value(position) +
                        imag *
                          this->simulation_parameters
                            .boundary_conditions_time_harmonic_electromagnetics
                            .excitation_x_imag.at(face_id)
                            ->value(position);

                      g_inc[1] =
                        this->simulation_parameters
                          .boundary_conditions_time_harmonic_electromagnetics
                          .excitation_y_real.at(face_id)
                          ->value(position) +
                        imag *
                          this->simulation_parameters
                            .boundary_conditions_time_harmonic_electromagnetics
                            .excitation_y_imag.at(face_id)
                            ->value(position);

                      g_inc[2] =
                        this->simulation_parameters
                          .boundary_conditions_time_harmonic_electromagnetics
                          .excitation_z_real.at(face_id)
                          ->value(position) +
                        imag *
                          this->simulation_parameters
                            .boundary_conditions_time_harmonic_electromagnetics
                            .excitation_z_imag.at(face_id)
                            ->value(position);
                    }
                  if (bc_type ==
                      BoundaryConditions::BoundaryType::waveguide_port)
                    {
                      unsigned int boundary_index = std::distance(
                        time_harmonic_maxwell_parameters.waveguide_boundary_ids
                          .begin(),
                        std::ranges::find(time_harmonic_maxwell_parameters
                                            .waveguide_boundary_ids,
                                          face->boundary_id()));

                      std::tie(g_inc, boundary_surface_admittance) =
                        compute_waveguide_port_excitation(
                          position,
                          normal,
                          effective_electric_permittivity,
                          effective_magnetic_permeability,
                          boundary_index);

                      conj_boundary_surface_admittance =
                        std::conj(boundary_surface_admittance);
                    }

                  // Now we loop on each relationship container to assemble
                  // the relevant matrices.
                  for (const auto &[i, j] : G_FF)
                    {
                      G_matrix(i, j) +=
                        (conj_boundary_surface_admittance * F_face[j] *
                         boundary_surface_admittance * F_face_conj[i] *
                         JxW_face)
                          .real();
                    }

                  for (const auto &[i, j] : G_FI)
                    {
                      G_matrix(i, j) +=
                        (n_cross_I_face[j] * boundary_surface_admittance *
                         F_face_conj[i] * JxW_face)
                          .real();
                    }

                  for (const auto &[i, j] : G_IF)
                    {
                      G_matrix(i, j) +=
                        (conj_boundary_surface_admittance * F_face[j] *
                         n_cross_I_face_conj[i] * JxW_face)
                          .real();
                    }

                  for (const auto &[i, j] : G_II)
                    {
                      G_matrix(i, j) +=
                        (n_cross_I_face[j] * n_cross_I_face_conj[i] * JxW_face)
                          .real();
                    }

                  for (const auto &[i, j] : B_hat_FH)
                    {
                      B_hat_matrix(i, j) +=
                        (n_cross_H_hat[j] * F_face_conj[i] * JxW_face).real();
                    }

                  for (const auto &[i, j] : B_hat_IE)
                    {
                      B_hat_matrix(i, j) +=
                        (n_cross_E_hat[j] * I_face_conj[i] * JxW_face).real();
                    }

                  for (const auto &[i, j] : B_hat_FE)
                    {
                      B_hat_matrix(i, j) -=
                        (boundary_surface_admittance * E_hat[j] *
                         F_face_conj[i] * JxW_face)
                          .real();
                    }

                  for (const auto &i : l_F)
                    {
                      l_vector[i] -= (g_inc * F_face_conj[i] * JxW_face).real();
                    }
                }
            } // End of face loop

          // Finally, after having assembled all the matrices and vectors, we
          // build the condensed version of the system.

          // We only need the inverse of the Gram matrix $G$, so we
          // invert it.
          G_matrix.invert();

          // We construct $M_4 = B^\dagger G^{-1}$ and $M_5 = \hat{B}^\dagger
          // G^{-1}$ with it:
          B_matrix.Tmmult(M4_matrix, G_matrix);
          B_hat_matrix.Tmmult(M5_matrix, G_matrix);

          // Then using $M_4$ we compute the condensed matrix $M_1 = B^\dagger
          // G^{-1} B$ and $M_2 = B^\dagger G^{-1} \hat{B}$:
          M4_matrix.mmult(M1_matrix, B_matrix);
          M4_matrix.mmult(M2_matrix, B_hat_matrix);

          // We also compute the matrix $M_3 = \hat{B}^\dagger G^{-1} \hat{B}$
          M5_matrix.mmult(M3_matrix, B_hat_matrix);

          // Finally, as for the $G$ matrix, we invert the $M_1$
          // matrix:
          M1_matrix.invert();

          // Now,  we have already the
          // solution on the skeleton and only need to perform $u_h = M_1^{-1}
          // (M_4 l - M_2 \hat{u}_h)$ on each cell. When this is obtained, we
          // can perform at the same time the error indicator (\Psi =
          // G^{-1}(l-B u_h
          // - \hat{B}\hat{u}_h)).

          // We first get the solution vector for this cell.
          cell_skeleton->get_dof_values(*present_solution_skeleton,
                                        cell_skeleton_solution);

          // Then we do the matrix-vector products to obtain the interior
          // unknowns.
          M2_matrix.vmult(tmp_vector_interior, cell_skeleton_solution);
          M4_matrix.vmult(cell_interior_rhs, l_vector);
          cell_interior_rhs -= tmp_vector_interior;
          M1_matrix.vmult(cell_interior_solution, cell_interior_rhs);

          // Finally, we map the cell interior solution to the global
          // interior solution.
          cell->distribute_local_to_global(cell_interior_solution,
                                           locally_owned_solution_interior);

          // We can also compute the error indicator on this cell.
          B_matrix.vmult(tmp_vector_error_indicator, cell_interior_solution);
          B_hat_matrix.vmult_add(tmp_vector_error_indicator,
                                 cell_skeleton_solution);
          l_vector -= tmp_vector_error_indicator;
          G_matrix.vmult(cell_residual, l_vector);

          // Compute the error indicator on the cell if the dpg error_indicator
          // is activated
          if (this->simulation_parameters.mesh_adaptation.var_adaptation_param
                .error_estimator ==
              Parameters::MultipleAdaptationParameters::ErrorEstimator::dpg)
            {
              double error_sqared =
                l_vector * cell_residual; // ||R||^2_V = R^T G^-1 R = R^T Psi
              local_estimated_error_per_cell(cell->active_cell_index()) =
                std::sqrt(error_sqared); // ||R||^2_V = R^T G^-1 R = R^T Psi
              residual_L2_norm += error_sqared;
            }

          cell_test->distribute_local_to_global(cell_residual,
                                                locally_owned_error_indicator);
        }
    }

  // After the loop over the cells, we finalize the assembly by compressing
  // the vectors because of the MPI parallelization.
  locally_owned_solution_interior.compress(VectorOperation::add);
  locally_owned_error_indicator.compress(VectorOperation::add);

  *this->present_solution            = locally_owned_solution_interior;
  *this->present_DPG_error_indicator = locally_owned_error_indicator;

  // We also output the global error indicator if the dpg error estimator is
  // activated and in verbose mode
  if ((this->simulation_parameters.mesh_adaptation.var_adaptation_param
         .error_estimator ==
       Parameters::MultipleAdaptationParameters::ErrorEstimator::dpg) &&
      (this->simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
         .verbosity != Parameters::Verbosity::quiet))
    {
      this->pcout << "   Time-Harmonic Maxwell DPG residual: "
                  << std::sqrt(
                       Utilities::MPI::sum(residual_L2_norm, mpi_communicator))
                  << std::endl;
    }
}


template class TimeHarmonicMaxwell<2>;
template class TimeHarmonicMaxwell<3>;
