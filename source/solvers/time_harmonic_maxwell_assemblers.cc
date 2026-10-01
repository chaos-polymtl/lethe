// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/time_harmonic_maxwell_assemblers.h>

#include <deal.II/base/numbers.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <tuple>

template <int dim>
std::pair<Tensor<1, dim, std::complex<double>>,
          Tensor<1, dim, std::complex<double>>>
compute_waveguide_port_incident_fields(
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters,
  const Point<dim>                           &p,
  const Tensor<1, dim>                       &normal,
  const std::complex<double>                 &effective_electric_permittivity,
  const std::complex<double>                 &effective_magnetic_permeability,
  const types::boundary_id                    boundary_id_index)
{
  // The waveguide port excitation would be completely different in 2D as curls
  // and cross products are not defined the same way as in 3D. Therefore, we do
  // not implement it for now and throw an error if someone tries to use the 2D
  // version.
  if constexpr (dim != 3)
    {
      (void)time_harmonic_maxwell_parameters;
      (void)p;
      (void)normal;
      (void)effective_electric_permittivity;
      (void)effective_magnetic_permeability;
      (void)boundary_id_index;
      AssertThrow(false, TimeHarmonicMaxwellDimensionNotSupported(dim));
      return std::pair<Tensor<1, dim, std::complex<double>>,
                       Tensor<1, dim, std::complex<double>>>();
    }
  else
    {
      // Define some constexpr values for the computation
      static constexpr std::complex<double> imag{0., 1.};
      static constexpr double               PI = numbers::PI;

      // Gather the relevant waveguide parameters
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

      // We first define the transverse face vectors of the waveguide in the
      // global system
      Tensor<1, dim> transverse_vector_1 =
        waveguide_corners[1] - waveguide_corners[0];
      Tensor<1, dim> transverse_vector_2 =
        waveguide_corners[2] - waveguide_corners[0];
      double         length_t1 = transverse_vector_1.norm();
      double         length_t2 = transverse_vector_2.norm();
      Tensor<1, dim> e_t1      = transverse_vector_1 / length_t1;
      Tensor<1, dim> e_t2      = transverse_vector_2 / length_t2;

      // Check if the transverse vectors form a perfect rectangle (orthogonal
      // and aligned with the axes)
      AssertThrow(
        std::abs(e_t1 * e_t2) < 1e-12,
        ExcMessage(
          "The transverse plane defined by the waveguide corners for the waveguide port at boundary ID " +
          std::to_string(time_harmonic_maxwell_parameters
                           .waveguide_boundary_ids[boundary_id_index]) +
          " is not a perfect rectangle (i.e., the vector created by the waveguide corners are not orthogonal). Please check the waveguide corners definition in the input prm file."));


      // Also verify that those transverse vectors are perpendicular to the
      // normal of the face and form an orthogonal basis.
      if ((std::abs(normal * e_t1) > 1e-12) ||
          (std::abs(normal * e_t2) > 1e-12))
        AssertThrow(
          false,
          ExcMessage(
            "The transverse plane defined by the waveguide corners for the waveguide port at boundary ID " +
            std::to_string(time_harmonic_maxwell_parameters
                             .waveguide_boundary_ids[boundary_id_index]) +
            " is not orthogonal to the boundary face normal. Please check the waveguide corners definition in the input prm file."));

      // Create a third vector to complete the right-handed coordinate system
      Tensor<1, dim> e_t3 = cross_product_3d(e_t1, e_t2);

      // Determine if the system needs to be flipped so the e_t3 vector points
      // in the direction opposite to the outward normal of the face boundary.
      // For an incident wave at the inlet, the propagation direction should
      // point into the domain (opposite to the outward normal). We use the
      // sign of (normal · t3) to determine this:
      //   - If normal · t3 > 0: t3 points outward, need to flip entire system
      //   - If normal · t3 < 0: t3 points inward (correct for incident wave)
      // Note that by swapping the basis vectors e_t1 and e_t2, we change the
      // parity of the system since it is equivalent to a reflection. This
      // means that pseudo vectors (like the magnetic field) will not change
      // sign while regular vectors (like the electric field) will. Even though
      // this does not affect the physics of the solution, it is important to
      // be consistent with the definition of the mode profiles that we use,
      // which assume a specific parity for the system. Therefore, we keep
      // track of the status of the parity and apply it later on to the
      // magnetic field to make it consistent with the switch of sign that has
      // been applied to the electric field when we compute the excitation.
      double parity_factor = 1.0;
      if ((normal * e_t3) > 0)
        {
          std::swap(e_t1, e_t2);
          std::swap(length_t1, length_t2);
          std::swap(m, n);
          e_t3 = cross_product_3d(e_t1,
                                  e_t2); // Recompute t3 after swapping t1 and
                                         // t2 to be sure the system is
                                         // right-handed
          parity_factor = -1.0;
        }

      // Compute the various wavenumbers k in using this global coordinate
      // system
      double k_t1 = m * PI / length_t1; // Transverse wavenumber 1
      double k_t2 = n * PI / length_t2; // Transverse wavenumber 2
      double k_c  = std::sqrt(k_t1 * k_t1 + k_t2 * k_t2); // Cutoff wavenumber
      std::complex<double> k =
        omega *
        std::sqrt(effective_electric_permittivity *
                  effective_magnetic_permeability); // Wavenumber in the medium
      std::complex<double> k_l = std::sqrt(
        k * k - std::complex<double>(k_c * k_c, 0)); // Longitudinal wavenumber

      // Verify that the mode is not evanescent, i.e. k_l is not purely
      // imaginary (k_c^2 < k^2). std::norm computes the squared magnitude of a
      // complex number.
      AssertThrow(std::norm(k) > (k_c * k_c),
                  ExcMessage(
                    "The chosen mode for the waveguide port at boundary ID " +
                    std::to_string(
                      time_harmonic_maxwell_parameters
                        .waveguide_boundary_ids[boundary_id_index]) +
                    " is evanescent at the given frequency. Please "
                    "choose another mode or increase the frequency."));

      // Now we want to compute the electromagnetic field and surface
      // admittance at point p as if the waveguide center was at the origin. So
      // we will perform a change of basis to a local coordinate system where
      // the waveguide center is at the origin.
      Tensor<1, dim> origin_local =
        0.25 * (waveguide_corners[0] + waveguide_corners[1] +
                waveguide_corners[2] + waveguide_corners[3]);

      Tensor<1, dim> p_local =
        p - origin_local; // Coordinates of point p in this local system.

      double x_local =
        p_local * e_t1; // Coordinate along e_t1 in the local system
      double y_local =
        p_local * e_t2; // Coordinate along e_t2 in the local system
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

      // Convert the E and H field components from the local coordinate system
      // back to the global coordinate system using the basis vectors e_t1,
      // e_t2, e_t3
      Tensor<1, dim, std::complex<double>> E_inc =
        E_inc_local[0] * e_t1 + E_inc_local[1] * e_t2 + E_inc_local[2] * e_t3;
      Tensor<1, dim, std::complex<double>> H_inc =
        parity_factor *
        (H_inc_local[0] * e_t1 + H_inc_local[1] * e_t2 + H_inc_local[2] * e_t3);

      return std::make_pair(E_inc, H_inc);
    }
}

template <int dim>
std::pair<Tensor<1, dim, std::complex<double>>, std::complex<double>>
compute_waveguide_port_excitation(
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters,
  const std::vector<double>  &waveguide_ports_electric_amplitudes,
  const Point<dim>           &p,
  const Tensor<1, dim>       &normal,
  const std::complex<double> &effective_electric_permittivity,
  const std::complex<double> &effective_magnetic_permeability,
  const types::boundary_id    boundary_id_index)
{
  // The waveguide port excitation would be completely different in 2D as curls
  // and cross products are not defined the same way as in 3D. Therefore, we do
  // not implement it for now and throw an error if someone tries to use the 2D
  // version.
  if constexpr (dim != 3)
    {
      (void)time_harmonic_maxwell_parameters;
      (void)waveguide_ports_electric_amplitudes;
      (void)p;
      (void)normal;
      (void)effective_electric_permittivity;
      (void)effective_magnetic_permeability;
      (void)boundary_id_index;
      AssertThrow(false, TimeHarmonicMaxwellDimensionNotSupported(dim));
      return std::pair<Tensor<1, dim, std::complex<double>>,
                       std::complex<double>>();
    }
  else
    {
      static constexpr double PI = numbers::PI;
      const auto              incident_fields =
        compute_waveguide_port_incident_fields(time_harmonic_maxwell_parameters,
                                               p,
                                               normal,
                                               effective_electric_permittivity,
                                               effective_magnetic_permeability,
                                               boundary_id_index);

      // Gather the relevant waveguide parameters
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
      const Tensor<1, dim, std::complex<double>> &H_inc =
        incident_fields.second;

      // We normalize the excitation by the maximum amplitude accross all
      // the waveguide ports to ensure that everything is normalized
      double scaling_factor =
        waveguide_ports_electric_amplitudes[boundary_id_index] /
        *std::ranges::max_element(waveguide_ports_electric_amplitudes);
      const Tensor<1, dim, std::complex<double>> excitation =
        scaling_factor * (cross_product_3d(normal, H_inc) +
                          map_H12(surface_admittance * E_inc, normal));

      return std::make_pair(excitation, surface_admittance);
    }
}


template <int dim>
TimeHarmonicMaxwellAssemblerCore<dim>::TimeHarmonicMaxwellAssemblerCore(
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters)
  : omega(2.0 * numbers::PI *
          time_harmonic_maxwell_parameters.electromagnetic_frequency)
{}

template <int dim>
void
TimeHarmonicMaxwellAssemblerCore<dim>::assemble_matrix(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
  THMCopyData                               &copy_data)
{
  static constexpr std::complex<double> imag{0., 1.};

  const std::vector<unsigned int> &test_dofs_F =
    scratch_data.test_dofs_electric;
  const std::vector<unsigned int> &test_dofs_I =
    scratch_data.test_dofs_magnetic;
  const std::vector<unsigned int> &trial_dofs_E =
    scratch_data.trial_interior_dofs_electric;
  const std::vector<unsigned int> &trial_dofs_H =
    scratch_data.trial_interior_dofs_magnetic;

  LAPACKFullMatrix<double> &G_matrix = copy_data.G_matrix;
  LAPACKFullMatrix<double> &B_matrix = copy_data.B_matrix;

  // Loop over the quadrature points of the cell
  for (unsigned int q = 0; q < scratch_data.n_q_points; ++q)
    {
      const std::complex<double> iweffective_magnetic_permeability =
        imag * omega * scratch_data.effective_magnetic_permeabilities[q];
      const std::complex<double> conj_iweffective_magnetic_permeability =
        std::conj(iweffective_magnetic_permeability);
      const std::complex<double> iweps_r =
        imag * omega * scratch_data.effective_electric_permittivities[q];
      const std::complex<double> conj_iweps_r = std::conj(iweps_r);

      const double &JxW = scratch_data.JxW[q];

      const auto F           = scratch_data.phi_F[q];
      const auto F_conj      = scratch_data.phi_F_conj[q];
      const auto curl_F      = scratch_data.curl_phi_F[q];
      const auto curl_F_conj = scratch_data.curl_phi_F_conj[q];
      const auto I           = scratch_data.phi_I[q];
      const auto I_conj      = scratch_data.phi_I_conj[q];
      const auto curl_I      = scratch_data.curl_phi_I[q];
      const auto curl_I_conj = scratch_data.curl_phi_I_conj[q];
      const auto E           = scratch_data.phi_E[q];
      const auto H           = scratch_data.phi_H[q];

      // Gram matrix G (Riesz map of the test space norm)
      for (const unsigned int i : test_dofs_F)
        {
          for (const unsigned int j : test_dofs_F)
            G_matrix(i, j) +=
              (((F[j] * F_conj[i]) + (curl_F[j] * curl_F_conj[i]) +
                (conj_iweps_r * F[j] * iweps_r * F_conj[i])) *
               JxW)
                .real();

          for (const unsigned int j : test_dofs_I)
            G_matrix(i, j) += (((curl_I[j] * iweps_r * F_conj[i]) -
                                (conj_iweffective_magnetic_permeability * I[j] *
                                 curl_F_conj[i])) *
                               JxW)
                                .real();
        }

      for (const unsigned int i : test_dofs_I)
        {
          for (const unsigned int j : test_dofs_F)
            G_matrix(i, j) +=
              (((conj_iweps_r * F[j] * curl_I_conj[i]) -
                (curl_F[j] * iweffective_magnetic_permeability * I_conj[i])) *
               JxW)
                .real();

          for (const unsigned int j : test_dofs_I)
            G_matrix(i, j) +=
              (((I[j] * I_conj[i]) + (curl_I[j] * curl_I_conj[i]) +
                (conj_iweffective_magnetic_permeability * I[j] *
                 iweffective_magnetic_permeability * I_conj[i])) *
               JxW)
                .real();
        }

      // Interior bilinear form B
      for (const unsigned int i : test_dofs_F)
        {
          for (const unsigned int j : trial_dofs_E)
            B_matrix(i, j) += (iweps_r * E[j] * F_conj[i] * JxW).real();

          for (const unsigned int j : trial_dofs_H)
            B_matrix(i, j) += (H[j] * curl_F_conj[i] * JxW).real();
        }

      for (const unsigned int i : test_dofs_I)
        {
          for (const unsigned int j : trial_dofs_E)
            B_matrix(i, j) += (E[j] * curl_I_conj[i] * JxW).real();

          for (const unsigned int j : trial_dofs_H)
            B_matrix(i, j) -=
              (iweffective_magnetic_permeability * H[j] * I_conj[i] * JxW)
                .real();
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwellAssemblerCore<dim>::assemble_rhs(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
  THMCopyData                               &copy_data)
{
  Vector<double> &l_vector = copy_data.l_vector;

  // Loop over the quadrature points of the cell. Note that in its simplest
  // form, the time-harmonic Maxwell equations do not have a magnetic source
  // term, so there is no contribution from the I test functions.
  for (unsigned int q = 0; q < scratch_data.n_q_points; ++q)
    {
      const double &JxW = scratch_data.JxW[q];

      const Tensor<1, dim, std::complex<double>> &J =
        scratch_data.current_density_values[q];
      const auto F_conj = scratch_data.phi_F_conj[q];

      for (const unsigned int i : scratch_data.test_dofs_electric)
        l_vector[i] += (J * F_conj[i] * JxW).real();
    }
}


template <int dim>
void
TimeHarmonicMaxwellAssemblerSkeleton<dim>::assemble_matrix(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
  THMCopyData                               &copy_data)
{
  // On the faces where a Robin boundary condition is applied, the magnetic
  // trace term is replaced by the Robin boundary condition.
  const bool robin_face = scratch_data.face_at_boundary &&
                          is_robin_boundary_type(boundary_conditions.type.at(
                            scratch_data.face_boundary_id));

  const std::vector<unsigned int> &test_dofs_F =
    scratch_data.test_dofs_electric;
  const std::vector<unsigned int> &test_dofs_I =
    scratch_data.test_dofs_magnetic;
  const std::vector<unsigned int> &trial_dofs_E_hat =
    scratch_data.trial_skeleton_dofs_electric;
  const std::vector<unsigned int> &trial_dofs_H_hat =
    scratch_data.trial_skeleton_dofs_magnetic;

  LAPACKFullMatrix<double> &B_hat_matrix = copy_data.B_hat_matrix;

  // Loop over the quadrature points of the face
  for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
    {
      const double JxW_face = scratch_data.face_JxW[q];

      const auto F_face_conj   = scratch_data.phi_F_face_conj[q];
      const auto I_face_conj   = scratch_data.phi_I_face_conj[q];
      const auto n_cross_E_hat = scratch_data.n_cross_phi_E_hat[q];
      const auto n_cross_H_hat = scratch_data.n_cross_phi_H_hat[q];

      if (!robin_face)
        {
          for (const unsigned int i : test_dofs_F)
            for (const unsigned int j : trial_dofs_H_hat)
              B_hat_matrix(i, j) +=
                (n_cross_H_hat[j] * F_face_conj[i] * JxW_face).real();
        }

      for (const unsigned int i : test_dofs_I)
        for (const unsigned int j : trial_dofs_E_hat)
          B_hat_matrix(i, j) +=
            (n_cross_E_hat[j] * I_face_conj[i] * JxW_face).real();
    }
}

template <int dim>
void
TimeHarmonicMaxwellAssemblerSkeleton<dim>::assemble_rhs(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
  THMCopyData                               &copy_data)
{
  // The surface current density is only imposed on the interior faces
  if (scratch_data.face_at_boundary)
    return;

  Vector<double> &l_vector = copy_data.l_vector;

  // Loop over the quadrature points of the face. Each of the two cells sharing
  // the face receives half of the surface current density (see the class
  // documentation).
  for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
    {
      const double JxW_face = scratch_data.face_JxW[q];

      const Tensor<1, dim, std::complex<double>> &K_t =
        scratch_data.face_surface_current_density_values[q];
      const auto F_face_conj = scratch_data.phi_F_face_conj[q];

      for (const unsigned int i : scratch_data.test_dofs_electric)
        l_vector[i] += 0.5 * (K_t * F_face_conj[i] * JxW_face).real();
    }
}


template <int dim>
unsigned int
TimeHarmonicMaxwellAssemblerRobinBC<dim>::get_waveguide_port_index(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data) const
{
  return std::distance(
    time_harmonic_maxwell_parameters.waveguide_boundary_ids.begin(),
    std::ranges::find(time_harmonic_maxwell_parameters.waveguide_boundary_ids,
                      scratch_data.face_boundary_id));
}

template <int dim>
std::pair<Tensor<1, dim, std::complex<double>>, std::complex<double>>
TimeHarmonicMaxwellAssemblerRobinBC<dim>::compute_robin_boundary_data(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
  const unsigned int                         q,
  const BoundaryConditions::BoundaryType     bc_type,
  const unsigned int                         waveguide_port_index) const
{
  static constexpr std::complex<double> imag{0., 1.};

  Tensor<1, dim, std::complex<double>> g_inc;
  std::complex<double>                 boundary_surface_admittance;

  const Point<dim>        &position = scratch_data.face_quadrature_points[q];
  const types::boundary_id face_id  = scratch_data.face_boundary_id;

  if (bc_type == BoundaryConditions::BoundaryType::silver_muller)
    {
      boundary_surface_admittance =
        std::sqrt(scratch_data.face_effective_electric_permittivities[q] /
                  scratch_data.face_effective_magnetic_permeabilities[q]);
      g_inc = 0.;
    }
  if (bc_type == BoundaryConditions::BoundaryType::impedance_boundary)
    {
      boundary_surface_admittance =
        boundary_conditions.surface_admittance_real.at(face_id)->value(
          position) +
        imag * boundary_conditions.surface_admittance_imag.at(face_id)->value(
                 position);

      // Get the incident electromagnetic field at this face. The excitation
      // is only defined for 3D problems.
      if constexpr (dim == 3)
        {
          g_inc[0] =
            boundary_conditions.excitation_x_real.at(face_id)->value(position) +
            imag * boundary_conditions.excitation_x_imag.at(face_id)->value(
                     position);

          g_inc[1] =
            boundary_conditions.excitation_y_real.at(face_id)->value(position) +
            imag * boundary_conditions.excitation_y_imag.at(face_id)->value(
                     position);

          g_inc[2] =
            boundary_conditions.excitation_z_real.at(face_id)->value(position) +
            imag * boundary_conditions.excitation_z_imag.at(face_id)->value(
                     position);
        }
      else
        AssertThrow(false, TimeHarmonicMaxwellDimensionNotSupported(dim));
    }
  if (bc_type == BoundaryConditions::BoundaryType::waveguide_port)
    {
      std::tie(g_inc, boundary_surface_admittance) =
        compute_waveguide_port_excitation(
          time_harmonic_maxwell_parameters,
          waveguide_ports_electric_amplitudes,
          position,
          scratch_data.face_normals[q],
          scratch_data.face_effective_electric_permittivities[q],
          scratch_data.face_effective_magnetic_permeabilities[q],
          waveguide_port_index);
    }

  return std::make_pair(g_inc, boundary_surface_admittance);
}

template <int dim>
void
TimeHarmonicMaxwellAssemblerRobinBC<dim>::assemble_matrix(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
  THMCopyData                               &copy_data)
{
  if (!scratch_data.face_at_boundary)
    return;

  const BoundaryConditions::BoundaryType bc_type =
    boundary_conditions.type.at(scratch_data.face_boundary_id);
  if (!is_robin_boundary_type(bc_type))
    return;

  const unsigned int waveguide_port_index =
    (bc_type == BoundaryConditions::BoundaryType::waveguide_port) ?
      get_waveguide_port_index(scratch_data) :
      0;

  const std::vector<unsigned int> &test_dofs_F =
    scratch_data.test_dofs_electric;
  const std::vector<unsigned int> &test_dofs_I =
    scratch_data.test_dofs_magnetic;
  const std::vector<unsigned int> &trial_dofs_E_hat =
    scratch_data.trial_skeleton_dofs_electric;

  LAPACKFullMatrix<double> &G_matrix     = copy_data.G_matrix;
  LAPACKFullMatrix<double> &B_hat_matrix = copy_data.B_hat_matrix;

  // Loop over the quadrature points of the face
  for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
    {
      const double JxW_face = scratch_data.face_JxW[q];

      const std::complex<double> boundary_surface_admittance =
        compute_robin_boundary_data(scratch_data,
                                    q,
                                    bc_type,
                                    waveguide_port_index)
          .second;
      const std::complex<double> conj_boundary_surface_admittance =
        std::conj(boundary_surface_admittance);

      const auto F_face              = scratch_data.phi_F_face[q];
      const auto F_face_conj         = scratch_data.phi_F_face_conj[q];
      const auto n_cross_I_face      = scratch_data.n_cross_phi_I_face[q];
      const auto n_cross_I_face_conj = scratch_data.n_cross_phi_I_face_conj[q];
      const auto E_hat               = scratch_data.phi_E_hat[q];

      // Boundary terms of the Gram matrix G coming from the energy norm that
      // is minimized on the Robin boundaries
      for (const unsigned int i : test_dofs_F)
        {
          for (const unsigned int j : test_dofs_F)
            G_matrix(i, j) +=
              (conj_boundary_surface_admittance * F_face[j] *
               boundary_surface_admittance * F_face_conj[i] * JxW_face)
                .real();

          for (const unsigned int j : test_dofs_I)
            G_matrix(i, j) += (n_cross_I_face[j] * boundary_surface_admittance *
                               F_face_conj[i] * JxW_face)
                                .real();
        }

      for (const unsigned int i : test_dofs_I)
        {
          for (const unsigned int j : test_dofs_F)
            G_matrix(i, j) += (conj_boundary_surface_admittance * F_face[j] *
                               n_cross_I_face_conj[i] * JxW_face)
                                .real();

          for (const unsigned int j : test_dofs_I)
            G_matrix(i, j) +=
              (n_cross_I_face[j] * n_cross_I_face_conj[i] * JxW_face).real();
        }

      // Robin term of the skeleton bilinear form B_hat, which replaces the
      // magnetic trace term
      for (const unsigned int i : test_dofs_F)
        for (const unsigned int j : trial_dofs_E_hat)
          B_hat_matrix(i, j) -=
            (boundary_surface_admittance * E_hat[j] * F_face_conj[i] * JxW_face)
              .real();
    }
}

template <int dim>
void
TimeHarmonicMaxwellAssemblerRobinBC<dim>::assemble_rhs(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
  THMCopyData                               &copy_data)
{
  if (!scratch_data.face_at_boundary)
    return;

  const BoundaryConditions::BoundaryType bc_type =
    boundary_conditions.type.at(scratch_data.face_boundary_id);
  if (!is_robin_boundary_type(bc_type))
    return;

  const unsigned int waveguide_port_index =
    (bc_type == BoundaryConditions::BoundaryType::waveguide_port) ?
      get_waveguide_port_index(scratch_data) :
      0;

  Vector<double> &l_vector = copy_data.l_vector;

  // Loop over the quadrature points of the face
  for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
    {
      const double JxW_face = scratch_data.face_JxW[q];

      const Tensor<1, dim, std::complex<double>> g_inc =
        compute_robin_boundary_data(scratch_data,
                                    q,
                                    bc_type,
                                    waveguide_port_index)
          .first;

      const auto F_face_conj = scratch_data.phi_F_face_conj[q];

      for (const unsigned int i : scratch_data.test_dofs_electric)
        l_vector[i] -= (g_inc * F_face_conj[i] * JxW_face).real();
    }
}


template std::pair<Tensor<1, 2, std::complex<double>>,
                   Tensor<1, 2, std::complex<double>>>
compute_waveguide_port_incident_fields<2>(
  const Parameters::TimeHarmonicMaxwell<2> &time_harmonic_maxwell_parameters,
  const Point<2>                           &p,
  const Tensor<1, 2>                       &normal,
  const std::complex<double>               &effective_electric_permittivity,
  const std::complex<double>               &effective_magnetic_permeability,
  const types::boundary_id                  boundary_id_index);
template std::pair<Tensor<1, 3, std::complex<double>>,
                   Tensor<1, 3, std::complex<double>>>
compute_waveguide_port_incident_fields<3>(
  const Parameters::TimeHarmonicMaxwell<3> &time_harmonic_maxwell_parameters,
  const Point<3>                           &p,
  const Tensor<1, 3>                       &normal,
  const std::complex<double>               &effective_electric_permittivity,
  const std::complex<double>               &effective_magnetic_permeability,
  const types::boundary_id                  boundary_id_index);

template std::pair<Tensor<1, 2, std::complex<double>>, std::complex<double>>
compute_waveguide_port_excitation<2>(
  const Parameters::TimeHarmonicMaxwell<2> &time_harmonic_maxwell_parameters,
  const std::vector<double>                &waveguide_ports_electric_amplitudes,
  const Point<2>                           &p,
  const Tensor<1, 2>                       &normal,
  const std::complex<double>               &effective_electric_permittivity,
  const std::complex<double>               &effective_magnetic_permeability,
  const types::boundary_id                  boundary_id_index);
template std::pair<Tensor<1, 3, std::complex<double>>, std::complex<double>>
compute_waveguide_port_excitation<3>(
  const Parameters::TimeHarmonicMaxwell<3> &time_harmonic_maxwell_parameters,
  const std::vector<double>                &waveguide_ports_electric_amplitudes,
  const Point<3>                           &p,
  const Tensor<1, 3>                       &normal,
  const std::complex<double>               &effective_electric_permittivity,
  const std::complex<double>               &effective_magnetic_permeability,
  const types::boundary_id                  boundary_id_index);

template class TimeHarmonicMaxwellAssemblerCore<2>;
template class TimeHarmonicMaxwellAssemblerCore<3>;
template class TimeHarmonicMaxwellAssemblerSkeleton<2>;
template class TimeHarmonicMaxwellAssemblerSkeleton<3>;
template class TimeHarmonicMaxwellAssemblerRobinBC<2>;
template class TimeHarmonicMaxwellAssemblerRobinBC<3>;
