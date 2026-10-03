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
  const Tensor<1, dim, std::complex<double>> &H_inc = incident_fields.second;

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


namespace
{
  /**
   * @brief Add a complex value to a local matrix of the real DPG system. The
   * rows and columns of the real system are associated with the real and
   * imaginary parts of the shape functions. Since a shape function is a real
   * base function multiplied by 1 (real part) or i (imaginary part), the
   * entries associated with the base functions a (rows) and b (columns) are
   * Re(s_b conj(s_a) z), which gives the 2x2 block
   * [[Re z, -Im z], [Im z, Re z]].
   *
   * @param[in,out] matrix Local matrix.
   *
   * @param[in] row_real Local dof index of the real part of the row base
   * function.
   *
   * @param[in] row_imag Local dof index of the imaginary part of the row base
   * function.
   *
   * @param[in] column_real Local dof index of the real part of the column base
   * function.
   *
   * @param[in] column_imag Local dof index of the imaginary part of the column
   * base function.
   *
   * @param[in] z Complex value associated with the pair of base functions.
   */
  inline void
  add_complex_block(LAPACKFullMatrix<double>   &matrix,
                    const unsigned int          row_real,
                    const unsigned int          row_imag,
                    const unsigned int          column_real,
                    const unsigned int          column_imag,
                    const std::complex<double> &z)
  {
    matrix(row_real, column_real) += z.real();
    matrix(row_real, column_imag) -= z.imag();
    matrix(row_imag, column_real) += z.imag();
    matrix(row_imag, column_imag) += z.real();
  }

  /**
   * @brief Add a real value to a local matrix of the real DPG system. This is
   * the 2x2 block of add_complex_block() for a real z: the real and imaginary
   * parts are not coupled.
   *
   * @param[in,out] matrix Local matrix.
   *
   * @param[in] row_real Local dof index of the real part of the row base
   * function.
   *
   * @param[in] row_imag Local dof index of the imaginary part of the row base
   * function.
   *
   * @param[in] column_real Local dof index of the real part of the column base
   * function.
   *
   * @param[in] column_imag Local dof index of the imaginary part of the column
   * base function.
   *
   * @param[in] z Real value associated with the pair of base functions.
   */
  inline void
  add_real_block(LAPACKFullMatrix<double> &matrix,
                 const unsigned int        row_real,
                 const unsigned int        row_imag,
                 const unsigned int        column_real,
                 const unsigned int        column_imag,
                 const double              z)
  {
    matrix(row_real, column_real) += z;
    matrix(row_imag, column_imag) += z;
  }
} // namespace


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
  DPGCopyData                               &copy_data)
{
  using ScratchData = TimeHarmonicMaxwellScratchData<dim>;
  static constexpr std::complex<double> imag{0., 1.};

  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_test     = scratch_data.n_base_functions_test;
  const unsigned int n_interior = scratch_data.n_base_functions_trial_interior;

  // Local dof indices of the real and imaginary parts of the base functions:
  // electric (F) and magnetic (I) test functions, and electric (E) and
  // magnetic (H) interior trial functions
  const auto &F_real = scratch_data.test_base_dofs[ScratchData::E_real_field];
  const auto &F_imag = scratch_data.test_base_dofs[ScratchData::E_imag_field];
  const auto &I_real = scratch_data.test_base_dofs[ScratchData::H_real_field];
  const auto &I_imag = scratch_data.test_base_dofs[ScratchData::H_imag_field];
  const auto &E_real =
    scratch_data.trial_interior_base_dofs[ScratchData::E_real_field];
  const auto &E_imag =
    scratch_data.trial_interior_base_dofs[ScratchData::E_imag_field];
  const auto &H_real =
    scratch_data.trial_interior_base_dofs[ScratchData::H_real_field];
  const auto &H_imag =
    scratch_data.trial_interior_base_dofs[ScratchData::H_imag_field];

  LAPACKFullMatrix<double> &G_matrix = copy_data.G_matrix;
  LAPACKFullMatrix<double> &B_matrix = copy_data.B_matrix;

  // Coefficients of the terms at the quadrature points, multiplied by JxW:
  // (1 + |i omega eps|^2) JxW, (1 + |i omega mu|^2) JxW, i omega eps JxW,
  // conj(i omega mu) JxW and i omega mu JxW
  std::vector<double>               JxW_FF(n_q_points);
  std::vector<double>               JxW_II(n_q_points);
  std::vector<std::complex<double>> iweps_JxW(n_q_points);
  std::vector<std::complex<double>> conj_iwmu_JxW(n_q_points);
  std::vector<std::complex<double>> iwmu_JxW(n_q_points);
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const std::complex<double> iweps_r =
        imag * omega * scratch_data.effective_electric_permittivities[q];
      const std::complex<double> iwmu_r =
        imag * omega * scratch_data.effective_magnetic_permeabilities[q];
      const double JxW = scratch_data.JxW[q];

      JxW_FF[q]        = (1. + std::norm(iweps_r)) * JxW;
      JxW_II[q]        = (1. + std::norm(iwmu_r)) * JxW;
      iweps_JxW[q]     = iweps_r * JxW;
      conj_iwmu_JxW[q] = std::conj(iwmu_r) * JxW;
      iwmu_JxW[q]      = iwmu_r * JxW;
    }

  // Gram matrix G (Riesz map of the test space norm). G is symmetric, so we
  // only compute the pairs of base functions a <= b and add the transposed
  // blocks. For the base functions a (rows) and b (columns), the terms are:
  // - FF: (1 + |iweps|^2) (phi_b, phi_a) + (curl phi_b, curl phi_a)
  // - II: (1 + |iwmu|^2) (phi_b, phi_a) + (curl phi_b, curl phi_a)
  // - FI: iweps (curl phi_b, phi_a) - conj(iwmu) (phi_b, curl phi_a)
  // - IF: the transpose of FI, which is the block of the complex conjugate,
  // with iweps = i omega eps and iwmu = i omega mu.
  for (unsigned int a = 0; a < n_test; ++a)
    {
      const Tensor<1, dim> *phi_a  = &scratch_data.phi_test(a, 0);
      const Tensor<1, dim> *curl_a = &scratch_data.curl_phi_test(a, 0);

      for (unsigned int b = a; b < n_test; ++b)
        {
          const Tensor<1, dim> *phi_b  = &scratch_data.phi_test(b, 0);
          const Tensor<1, dim> *curl_b = &scratch_data.curl_phi_test(b, 0);

          double               z_FF = 0.;
          double               z_II = 0.;
          std::complex<double> z_FI_ab(0.);
          std::complex<double> z_FI_ba(0.);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double phi_phi      = phi_b[q] * phi_a[q];
              const double curl_curl    = curl_b[q] * curl_a[q];
              const double curl_b_phi_a = curl_b[q] * phi_a[q];
              const double phi_b_curl_a = phi_b[q] * curl_a[q];

              z_FF += JxW_FF[q] * phi_phi + scratch_data.JxW[q] * curl_curl;
              z_II += JxW_II[q] * phi_phi + scratch_data.JxW[q] * curl_curl;
              z_FI_ab +=
                iweps_JxW[q] * curl_b_phi_a - conj_iwmu_JxW[q] * phi_b_curl_a;
              z_FI_ba +=
                iweps_JxW[q] * phi_b_curl_a - conj_iwmu_JxW[q] * curl_b_phi_a;
            }

          add_real_block(
            G_matrix, F_real[a], F_imag[a], F_real[b], F_imag[b], z_FF);
          add_real_block(
            G_matrix, I_real[a], I_imag[a], I_real[b], I_imag[b], z_II);
          add_complex_block(
            G_matrix, F_real[a], F_imag[a], I_real[b], I_imag[b], z_FI_ab);
          add_complex_block(G_matrix,
                            I_real[b],
                            I_imag[b],
                            F_real[a],
                            F_imag[a],
                            std::conj(z_FI_ab));

          if (b != a)
            {
              add_real_block(
                G_matrix, F_real[b], F_imag[b], F_real[a], F_imag[a], z_FF);
              add_real_block(
                G_matrix, I_real[b], I_imag[b], I_real[a], I_imag[a], z_II);
              add_complex_block(
                G_matrix, F_real[b], F_imag[b], I_real[a], I_imag[a], z_FI_ba);
              add_complex_block(G_matrix,
                                I_real[a],
                                I_imag[a],
                                F_real[b],
                                F_imag[b],
                                std::conj(z_FI_ba));
            }
        }
    }

  // Interior bilinear form B. For the test base function a (rows) and the
  // interior trial base function c (columns), the terms are:
  // - FE: i omega eps (phi_c, phi_a)
  // - FH and IE: (phi_c, curl phi_a)
  // - IH: -i omega mu (phi_c, phi_a)
  for (unsigned int a = 0; a < n_test; ++a)
    {
      const Tensor<1, dim> *phi_a  = &scratch_data.phi_test(a, 0);
      const Tensor<1, dim> *curl_a = &scratch_data.curl_phi_test(a, 0);

      for (unsigned int c = 0; c < n_interior; ++c)
        {
          const Tensor<1, dim> *phi_c = &scratch_data.phi_trial_interior(c, 0);

          std::complex<double> z_FE(0.);
          std::complex<double> z_IH(0.);
          double               z_FH = 0.;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double phi_c_phi_a  = phi_c[q] * phi_a[q];
              const double phi_c_curl_a = phi_c[q] * curl_a[q];

              z_FE += iweps_JxW[q] * phi_c_phi_a;
              z_FH += scratch_data.JxW[q] * phi_c_curl_a;
              z_IH -= iwmu_JxW[q] * phi_c_phi_a;
            }

          add_complex_block(
            B_matrix, F_real[a], F_imag[a], E_real[c], E_imag[c], z_FE);
          add_real_block(
            B_matrix, F_real[a], F_imag[a], H_real[c], H_imag[c], z_FH);
          add_real_block(
            B_matrix, I_real[a], I_imag[a], E_real[c], E_imag[c], z_FH);
          add_complex_block(
            B_matrix, I_real[a], I_imag[a], H_real[c], H_imag[c], z_IH);
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwellAssemblerCore<dim>::assemble_rhs(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
  DPGCopyData                               &copy_data)
{
  using ScratchData = TimeHarmonicMaxwellScratchData<dim>;

  Vector<double> &l_vector = copy_data.l_vector;

  // Loop over the quadrature points of the cell. No volume source term is
  // considered at the moment, so the interior load vector is zero. Note that
  // in its simplest form, the time-harmonic Maxwell equations do not have a
  // magnetic source term, so there is no contribution from the I test
  // functions.
  for (unsigned int q = 0; q < scratch_data.n_q_points; ++q)
    {
      for (unsigned int a = 0; a < scratch_data.n_base_functions_test; ++a)
        {
          l_vector[scratch_data.test_base_dofs[ScratchData::E_real_field][a]] +=
            0.0;
          l_vector[scratch_data.test_base_dofs[ScratchData::E_imag_field][a]] +=
            0.0;
        }
    }
}


template <int dim>
void
TimeHarmonicMaxwellAssemblerSkeleton<dim>::assemble_matrix(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
  DPGCopyData                               &copy_data)
{
  using ScratchData = TimeHarmonicMaxwellScratchData<dim>;

  // On the faces where a Robin boundary condition is applied, the magnetic
  // trace term is replaced by the Robin boundary condition.
  const bool robin_face = scratch_data.face_at_boundary &&
                          is_robin_boundary_type(boundary_conditions.type.at(
                            scratch_data.face_boundary_id));

  const auto &F_real = scratch_data.test_base_dofs[ScratchData::E_real_field];
  const auto &F_imag = scratch_data.test_base_dofs[ScratchData::E_imag_field];
  const auto &I_real = scratch_data.test_base_dofs[ScratchData::H_real_field];
  const auto &I_imag = scratch_data.test_base_dofs[ScratchData::H_imag_field];
  const auto &E_hat_real =
    scratch_data.trial_skeleton_base_dofs[ScratchData::E_real_field];
  const auto &E_hat_imag =
    scratch_data.trial_skeleton_base_dofs[ScratchData::E_imag_field];
  const auto &H_hat_real =
    scratch_data.trial_skeleton_base_dofs[ScratchData::H_real_field];
  const auto &H_hat_imag =
    scratch_data.trial_skeleton_base_dofs[ScratchData::H_imag_field];

  LAPACKFullMatrix<double> &B_hat_matrix = copy_data.B_hat_matrix;

  // The skeleton terms only involve the tangential traces of the test and
  // skeleton trial functions, so only the base functions associated with the
  // face contribute. Since the electric and magnetic skeleton fields share the
  // same base functions, the IE and FH terms have the same value
  // <n x e_hat_d, phi_a> for the test base function a and the skeleton base
  // function d.
  for (const unsigned int a :
       scratch_data.face_test_base_functions[scratch_data.face_no])
    {
      const Tensor<1, dim> *phi_a = &scratch_data.phi_test_face(a, 0);

      for (const unsigned int d :
           scratch_data
             .face_trial_skeleton_base_functions[scratch_data.face_no])
        {
          const Tensor<1, dim> *n_cross_e_hat_d =
            &scratch_data.n_cross_tangential_phi_trial_skeleton(d, 0);

          double z = 0.;
          for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
            z += scratch_data.face_JxW[q] * (n_cross_e_hat_d[q] * phi_a[q]);

          add_real_block(B_hat_matrix,
                         I_real[a],
                         I_imag[a],
                         E_hat_real[d],
                         E_hat_imag[d],
                         z);
          if (!robin_face)
            add_real_block(B_hat_matrix,
                           F_real[a],
                           F_imag[a],
                           H_hat_real[d],
                           H_hat_imag[d],
                           z);
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwellAssemblerSkeleton<dim>::assemble_rhs(
  const TimeHarmonicMaxwellScratchData<dim> & /*scratch_data*/,
  DPGCopyData & /*copy_data*/)
{}


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
  DPGCopyData                               &copy_data)
{
  using ScratchData = TimeHarmonicMaxwellScratchData<dim>;

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

  const unsigned int n_face_q_points = scratch_data.n_face_q_points;
  const unsigned int n_test          = scratch_data.n_base_functions_test;

  const auto &F_real = scratch_data.test_base_dofs[ScratchData::E_real_field];
  const auto &F_imag = scratch_data.test_base_dofs[ScratchData::E_imag_field];
  const auto &I_real = scratch_data.test_base_dofs[ScratchData::H_real_field];
  const auto &I_imag = scratch_data.test_base_dofs[ScratchData::H_imag_field];
  const auto &E_hat_real =
    scratch_data.trial_skeleton_base_dofs[ScratchData::E_real_field];
  const auto &E_hat_imag =
    scratch_data.trial_skeleton_base_dofs[ScratchData::E_imag_field];

  LAPACKFullMatrix<double> &G_matrix     = copy_data.G_matrix;
  LAPACKFullMatrix<double> &B_hat_matrix = copy_data.B_hat_matrix;

  // Surface admittance at the face quadrature points, multiplied by JxW
  std::vector<double>               JxW_YY(n_face_q_points);
  std::vector<std::complex<double>> Y_JxW(n_face_q_points);
  for (unsigned int q = 0; q < n_face_q_points; ++q)
    {
      const std::complex<double> boundary_surface_admittance =
        compute_robin_boundary_data(scratch_data,
                                    q,
                                    bc_type,
                                    waveguide_port_index)
          .second;
      JxW_YY[q] =
        std::norm(boundary_surface_admittance) * scratch_data.face_JxW[q];
      Y_JxW[q] = boundary_surface_admittance * scratch_data.face_JxW[q];
    }

  // Boundary terms of the Gram matrix G coming from the energy norm that is
  // minimized on the Robin boundaries. They involve the full trace of the test
  // functions, so all the test base functions contribute. As for the cell
  // terms, G is symmetric and the IF block is the transpose of the FI block.
  // For the base functions a (rows) and b (columns), the terms are:
  // - FF: |Y|^2 <phi_b, phi_a>
  // - II: <n x phi_b, n x phi_a>
  // - FI: Y <n x phi_b, phi_a>
  for (unsigned int a = 0; a < n_test; ++a)
    {
      const Tensor<1, dim> *phi_a = &scratch_data.phi_test_face(a, 0);
      const Tensor<1, dim> *n_cross_a =
        &scratch_data.n_cross_phi_test_face(a, 0);

      for (unsigned int b = a; b < n_test; ++b)
        {
          const Tensor<1, dim> *phi_b = &scratch_data.phi_test_face(b, 0);
          const Tensor<1, dim> *n_cross_b =
            &scratch_data.n_cross_phi_test_face(b, 0);

          double               z_FF = 0.;
          double               z_II = 0.;
          std::complex<double> z_FI_ab(0.);
          std::complex<double> z_FI_ba(0.);
          for (unsigned int q = 0; q < n_face_q_points; ++q)
            {
              z_FF += JxW_YY[q] * (phi_b[q] * phi_a[q]);
              z_II += scratch_data.face_JxW[q] * (n_cross_b[q] * n_cross_a[q]);
              z_FI_ab += Y_JxW[q] * (n_cross_b[q] * phi_a[q]);
              z_FI_ba += Y_JxW[q] * (n_cross_a[q] * phi_b[q]);
            }

          add_real_block(
            G_matrix, F_real[a], F_imag[a], F_real[b], F_imag[b], z_FF);
          add_real_block(
            G_matrix, I_real[a], I_imag[a], I_real[b], I_imag[b], z_II);
          add_complex_block(
            G_matrix, F_real[a], F_imag[a], I_real[b], I_imag[b], z_FI_ab);
          add_complex_block(G_matrix,
                            I_real[b],
                            I_imag[b],
                            F_real[a],
                            F_imag[a],
                            std::conj(z_FI_ab));

          if (b != a)
            {
              add_real_block(
                G_matrix, F_real[b], F_imag[b], F_real[a], F_imag[a], z_FF);
              add_real_block(
                G_matrix, I_real[b], I_imag[b], I_real[a], I_imag[a], z_II);
              add_complex_block(
                G_matrix, F_real[b], F_imag[b], I_real[a], I_imag[a], z_FI_ba);
              add_complex_block(G_matrix,
                                I_real[a],
                                I_imag[a],
                                F_real[b],
                                F_imag[b],
                                std::conj(z_FI_ba));
            }
        }
    }

  // Robin term of the skeleton bilinear form B_hat, which replaces the
  // magnetic trace term: -Y <e_hat_d, phi_a>. It only involves the tangential
  // trace of the test functions, so only the base functions associated with
  // the face contribute.
  for (const unsigned int a :
       scratch_data.face_test_base_functions[scratch_data.face_no])
    {
      const Tensor<1, dim> *phi_a = &scratch_data.phi_test_face(a, 0);

      for (const unsigned int d :
           scratch_data
             .face_trial_skeleton_base_functions[scratch_data.face_no])
        {
          const Tensor<1, dim> *e_hat_d =
            &scratch_data.tangential_phi_trial_skeleton(d, 0);

          std::complex<double> z(0.);
          for (unsigned int q = 0; q < n_face_q_points; ++q)
            z -= Y_JxW[q] * (e_hat_d[q] * phi_a[q]);

          add_complex_block(B_hat_matrix,
                            F_real[a],
                            F_imag[a],
                            E_hat_real[d],
                            E_hat_imag[d],
                            z);
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwellAssemblerRobinBC<dim>::assemble_rhs(
  const TimeHarmonicMaxwellScratchData<dim> &scratch_data,
  DPGCopyData                               &copy_data)
{
  using ScratchData = TimeHarmonicMaxwellScratchData<dim>;

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

  const unsigned int n_face_q_points = scratch_data.n_face_q_points;

  const auto &F_real = scratch_data.test_base_dofs[ScratchData::E_real_field];
  const auto &F_imag = scratch_data.test_base_dofs[ScratchData::E_imag_field];

  Vector<double> &l_vector = copy_data.l_vector;

  // Excitation at the face quadrature points, multiplied by JxW
  std::vector<Tensor<1, dim, std::complex<double>>> g_inc_JxW(n_face_q_points);
  for (unsigned int q = 0; q < n_face_q_points; ++q)
    g_inc_JxW[q] = compute_robin_boundary_data(scratch_data,
                                               q,
                                               bc_type,
                                               waveguide_port_index)
                     .first *
                   scratch_data.face_JxW[q];

  // Excitation term -<g, F>: for the test base function a, the real and
  // imaginary parts of w = <g, phi_a> go to the real and imaginary parts of
  // the electric test function.
  for (unsigned int a = 0; a < scratch_data.n_base_functions_test; ++a)
    {
      const Tensor<1, dim> *phi_a = &scratch_data.phi_test_face(a, 0);

      std::complex<double> w(0.);
      for (unsigned int q = 0; q < n_face_q_points; ++q)
        w += g_inc_JxW[q] * phi_a[q];

      l_vector[F_real[a]] -= w.real();
      l_vector[F_imag[a]] -= w.imag();
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
