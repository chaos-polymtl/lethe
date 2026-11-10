// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/tensors_and_points_dimension_manipulation.h>

#include <solvers/time_harmonic_maxwell.h>
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
  const unsigned int                          waveguide_port_index)
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
      (void)waveguide_port_index;
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
    time_harmonic_maxwell_parameters.waveguide_corners[waveguide_port_index];
  const Parameters::WaveguideMode mode =
    time_harmonic_maxwell_parameters.waveguide_mode[waveguide_port_index];
  unsigned int m =
    time_harmonic_maxwell_parameters.mode_order_m[waveguide_port_index];
  unsigned int n =
    time_harmonic_maxwell_parameters.mode_order_n[waveguide_port_index];

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
                       .waveguide_boundary_ids[waveguide_port_index]) +
      " is not a perfect rectangle (i.e., the vector created by the waveguide corners are not orthogonal). Please check the waveguide corners definition in the input prm file."));


  // Also verify that those transverse vectors are perpendicular to the normal
  // of the face and form an orthogonal basis.
  if ((std::abs(normal * e_t1) > 1e-12) || (std::abs(normal * e_t2) > 1e-12))
    AssertThrow(
      false,
      ExcMessage(
        "The transverse plane defined by the waveguide corners for the waveguide port at boundary ID " +
        std::to_string(time_harmonic_maxwell_parameters
                         .waveguide_boundary_ids[waveguide_port_index]) +
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
  // The flip swaps e_t1 and e_t2 (and the associated lengths and mode orders,
  // so the mode keeps the same transverse profile) and recomputes e_t3 as their
  // cross product so the resulting system is right-handed.
  if ((normal * e_t3) > 0)
    {
      std::swap(e_t1, e_t2);
      std::swap(length_t1, length_t2);
      std::swap(m, n);
      e_t3 =
        cross_product_3d(e_t1, e_t2); // Recompute t3 after swapping t1 and t2
                                      // to be sure the system is right-handed
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
  AssertThrow(
    std::norm(k) > (k_c * k_c),
    ExcMessage("The chosen mode for the waveguide port at boundary ID " +
               std::to_string(time_harmonic_maxwell_parameters
                                .waveguide_boundary_ids[waveguide_port_index]) +
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

  const double x_local =
    p_local * e_t1; // Coordinate along e_t1 in the local system
  const double y_local =
    p_local * e_t2; // Coordinate along e_t2 in the local system
  // We assume that the z_local coordinate is 0 since we are on the face.

  // Compute the E and H field components for the TE mode in the local
  // coordinate system {t1, t2, t3} = {x', y', z'}. We assume z' = 0 at the
  // boundary.
  Tensor<1, dim, std::complex<double>> E_inc_local;
  Tensor<1, dim, std::complex<double>> H_inc_local;

  // sin and cos functions that are repeatedly used in the computation of the E
  // and H fields.
  const double sin_k_t1_x = std::sin(k_t1 * (x_local + 0.5 * length_t1));
  const double cos_k_t1_x = std::cos(k_t1 * (x_local + 0.5 * length_t1));
  const double sin_k_t2_y = std::sin(k_t2 * (y_local + 0.5 * length_t2));
  const double cos_k_t2_y = std::cos(k_t2 * (y_local + 0.5 * length_t2));

  if (mode == Parameters::WaveguideMode::TE)
    {
      std::complex<double> factor =
        imag * omega * effective_magnetic_permeability / (k_c * k_c);

      E_inc_local[0] = -factor * k_t2 * cos_k_t1_x * sin_k_t2_y;
      E_inc_local[1] = factor * k_t1 * sin_k_t1_x * cos_k_t2_y;
      E_inc_local[2] = 0.0;

      H_inc_local[0] =
        -imag * k_l * k_t1 / (k_c * k_c) * sin_k_t1_x * cos_k_t2_y;
      H_inc_local[1] =
        -imag * k_l * k_t2 / (k_c * k_c) * cos_k_t1_x * sin_k_t2_y;
      H_inc_local[2] = cos_k_t1_x * cos_k_t2_y;
    }
  else if (mode == Parameters::WaveguideMode::TM)
    {
      std::complex<double> factor =
        imag * omega * effective_electric_permittivity / (k_c * k_c);

      H_inc_local[0] = -factor * k_t2 * sin_k_t1_x * cos_k_t2_y;
      H_inc_local[1] = factor * k_t1 * cos_k_t1_x * sin_k_t2_y;
      H_inc_local[2] = 0.0;

      E_inc_local[0] =
        imag * k_l * k_t1 / (k_c * k_c) * cos_k_t1_x * sin_k_t2_y;
      E_inc_local[1] =
        imag * k_l * k_t2 / (k_c * k_c) * sin_k_t1_x * cos_k_t2_y;
      E_inc_local[2] = sin_k_t1_x * sin_k_t2_y;
    }
  else
    {
      AssertThrow(false, ExcMessage("Unknown waveguide mode type."));
    }

  // Convert the E and H field components from the local coordinate system back
  // to the global coordinate system using the basis vectors e_t1, e_t2, e_t3
  const Tensor<1, dim, std::complex<double>> E_inc =
    E_inc_local[0] * e_t1 + E_inc_local[1] * e_t2 + E_inc_local[2] * e_t3;
  const Tensor<1, dim, std::complex<double>> H_inc =
    H_inc_local[0] * e_t1 + H_inc_local[1] * e_t2 + H_inc_local[2] * e_t3;

  return std::make_pair(E_inc, H_inc);
}

template <int dim>
unsigned int
get_waveguide_port_index(
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters,
  const types::boundary_id                    boundary_id)
{
  // The waveguide port index is the position of the boundary id in the list of
  // the port boundary ids of the waveguide inlets
  const std::vector<int> &waveguide_boundary_ids =
    time_harmonic_maxwell_parameters.waveguide_boundary_ids;
  const unsigned int waveguide_port_index = static_cast<unsigned int>(
    std::distance(waveguide_boundary_ids.begin(),
                  std::ranges::find(waveguide_boundary_ids,
                                    static_cast<int>(boundary_id))));

  return waveguide_port_index;
}

namespace
{
  /**
   * @brief Add a complex value to a 2x2 block of a local matrix of the real DPG
   * system. Each shape function is a real base function multiplied by s = 1
   * (real part) or by the imaginary unit (imaginary part). This gives the block
   * [[Re z, -Im z], [Im z, Re z]] for the (real, imaginary) rows of i and the
   * (real, imaginary) columns of j when multiplying shape functions.
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
  void
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
   * @brief Add a real value to a 2x2 block of a local matrix of the real DPG
   * system. This is the block of add_complex_block() for a real z: the real
   * and imaginary parts are not coupled, so only the diagonal of the block is
   * filled.
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
  void
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

  const unsigned int n_q_points            = scratch_data.n_q_points;
  const unsigned int n_base_functions_test = scratch_data.n_base_functions_test;
  const unsigned int n_base_functions_trial_interior =
    scratch_data.n_base_functions_trial_interior;

  // Local dof indices of the real and imaginary parts of each base function:
  // electric (F) and magnetic (I) test functions, and electric (E) and magnetic
  // (H) interior trial functions
  const std::vector<unsigned int> &test_dofs_F_real =
    scratch_data.test_local_dof_index_map[ScratchData::E_real];
  const std::vector<unsigned int> &test_dofs_F_imag =
    scratch_data.test_local_dof_index_map[ScratchData::E_imag];
  const std::vector<unsigned int> &test_dofs_I_real =
    scratch_data.test_local_dof_index_map[ScratchData::H_real];
  const std::vector<unsigned int> &test_dofs_I_imag =
    scratch_data.test_local_dof_index_map[ScratchData::H_imag];
  const std::vector<unsigned int> &trial_dofs_E_real =
    scratch_data.trial_interior_local_dof_index_map[ScratchData::E_real];
  const std::vector<unsigned int> &trial_dofs_E_imag =
    scratch_data.trial_interior_local_dof_index_map[ScratchData::E_imag];
  const std::vector<unsigned int> &trial_dofs_H_real =
    scratch_data.trial_interior_local_dof_index_map[ScratchData::H_real];
  const std::vector<unsigned int> &trial_dofs_H_imag =
    scratch_data.trial_interior_local_dof_index_map[ScratchData::H_imag];

  // Coefficients of the terms at the quadrature points, multiplied by JxW and
  // computed by the scratch data, with iweps = i omega eps and
  // iwmu = i omega mu. Since JxW is real, conj(iwmu) JxW is the complex
  // conjugate of iwmu JxW.
  const std::vector<double> &JxW              = scratch_data.JxW;
  const std::vector<double> &JxW_FF_container = scratch_data.JxW_FF_container;
  const std::vector<double> &JxW_II_container = scratch_data.JxW_II_container;
  const std::vector<std::complex<double>> &iweps_JxW_container =
    scratch_data.iweps_JxW_container;
  const std::vector<std::complex<double>> &iwmu_JxW_container =
    scratch_data.iwmu_JxW_container;

  LAPACKFullMatrix<double> &G_matrix = copy_data.G_matrix;
  LAPACKFullMatrix<double> &B_matrix = copy_data.B_matrix;

  // Gram matrix G (Riesz map of the test space norm). For the test base
  // functions i (rows) and j (columns), the blocks are:
  // - FF: (1 + |iweps|^2) (phi_j, phi_i) + (curl phi_j, curl phi_i)
  // - II: (1 + |iwmu|^2) (phi_j, phi_i) + (curl phi_j, curl phi_i)
  // - FI: iweps (curl phi_j, phi_i) - conj(iwmu) (phi_j, curl phi_i)
  // - IF: since G is Hermitian, the block IF(j, i) is conj(FI(i, j)).
  // FF and II are real and symmetric in i and j, so we only loop over the pairs
  // i <= j and add the blocks of the pair (j, i) at the same time.
  for (unsigned int i = 0; i < n_base_functions_test; ++i)
    {
      const auto phi_i      = scratch_data.phi_test[i];
      const auto curl_phi_i = scratch_data.curl_phi_test[i];

      for (unsigned int j = i; j < n_base_functions_test; ++j)
        {
          const auto phi_j      = scratch_data.phi_test[j];
          const auto curl_phi_j = scratch_data.curl_phi_test[j];

          double               z_FF = 0.;
          double               z_II = 0.;
          std::complex<double> z_FI_ij(0.);
          std::complex<double> z_FI_ji(0.);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // The four dot products are shared by the four blocks
              const double phi_j_phi_i   = phi_j[q] * phi_i[q];
              const double curl_j_curl_i = curl_phi_j[q] * curl_phi_i[q];
              const double curl_j_phi_i  = curl_phi_j[q] * phi_i[q];
              const double phi_j_curl_i  = phi_j[q] * curl_phi_i[q];

              z_FF +=
                JxW_FF_container[q] * phi_j_phi_i + JxW[q] * curl_j_curl_i;
              z_II +=
                JxW_II_container[q] * phi_j_phi_i + JxW[q] * curl_j_curl_i;
              z_FI_ij += iweps_JxW_container[q] * curl_j_phi_i -
                         std::conj(iwmu_JxW_container[q]) * phi_j_curl_i;
              z_FI_ji += iweps_JxW_container[q] * phi_j_curl_i -
                         std::conj(iwmu_JxW_container[q]) * curl_j_phi_i;
            }

          // Blocks of the pair (i, j)
          add_real_block(G_matrix,
                         test_dofs_F_real[i],
                         test_dofs_F_imag[i],
                         test_dofs_F_real[j],
                         test_dofs_F_imag[j],
                         z_FF);
          add_real_block(G_matrix,
                         test_dofs_I_real[i],
                         test_dofs_I_imag[i],
                         test_dofs_I_real[j],
                         test_dofs_I_imag[j],
                         z_II);
          add_complex_block(G_matrix,
                            test_dofs_F_real[i],
                            test_dofs_F_imag[i],
                            test_dofs_I_real[j],
                            test_dofs_I_imag[j],
                            z_FI_ij);
          add_complex_block(G_matrix,
                            test_dofs_I_real[j],
                            test_dofs_I_imag[j],
                            test_dofs_F_real[i],
                            test_dofs_F_imag[i],
                            std::conj(z_FI_ij));

          // Blocks of the pair (j, i), i != j because the diagonal blocks are
          // already filled by the pair (i, j)
          if (j != i)
            {
              add_real_block(G_matrix,
                             test_dofs_F_real[j],
                             test_dofs_F_imag[j],
                             test_dofs_F_real[i],
                             test_dofs_F_imag[i],
                             z_FF);
              add_real_block(G_matrix,
                             test_dofs_I_real[j],
                             test_dofs_I_imag[j],
                             test_dofs_I_real[i],
                             test_dofs_I_imag[i],
                             z_II);
              add_complex_block(G_matrix,
                                test_dofs_F_real[j],
                                test_dofs_F_imag[j],
                                test_dofs_I_real[i],
                                test_dofs_I_imag[i],
                                z_FI_ji);
              add_complex_block(G_matrix,
                                test_dofs_I_real[i],
                                test_dofs_I_imag[i],
                                test_dofs_F_real[j],
                                test_dofs_F_imag[j],
                                std::conj(z_FI_ji));
            }
        }
    }

  // Interior bilinear form B. For the test base function i (rows) and the
  // interior trial base function j (columns), the blocks are:
  // - FE: iweps (phi_j, phi_i)
  // - FH and IE: (phi_j, curl phi_i)
  // - IH: -iwmu (phi_j, phi_i)
  for (unsigned int i = 0; i < n_base_functions_test; ++i)
    {
      const auto phi_i      = scratch_data.phi_test[i];
      const auto curl_phi_i = scratch_data.curl_phi_test[i];

      for (unsigned int j = 0; j < n_base_functions_trial_interior; ++j)
        {
          const auto phi_j = scratch_data.phi_trial_interior[j];

          std::complex<double> z_FE(0.);
          double               z_FH = 0.;
          std::complex<double> z_IH(0.);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double phi_j_phi_i  = phi_j[q] * phi_i[q];
              const double phi_j_curl_i = phi_j[q] * curl_phi_i[q];

              z_FE += iweps_JxW_container[q] * phi_j_phi_i;
              z_FH += JxW[q] * phi_j_curl_i;
              z_IH -= iwmu_JxW_container[q] * phi_j_phi_i;
            }

          add_complex_block(B_matrix,
                            test_dofs_F_real[i],
                            test_dofs_F_imag[i],
                            trial_dofs_E_real[j],
                            trial_dofs_E_imag[j],
                            z_FE);
          add_real_block(B_matrix,
                         test_dofs_F_real[i],
                         test_dofs_F_imag[i],
                         trial_dofs_H_real[j],
                         trial_dofs_H_imag[j],
                         z_FH);
          add_real_block(B_matrix,
                         test_dofs_I_real[i],
                         test_dofs_I_imag[i],
                         trial_dofs_E_real[j],
                         trial_dofs_E_imag[j],
                         z_FH);
          add_complex_block(B_matrix,
                            test_dofs_I_real[i],
                            test_dofs_I_imag[i],
                            trial_dofs_H_real[j],
                            trial_dofs_H_imag[j],
                            z_IH);
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

  // We return early because there are no volume source terms in the
  // time-harmonic Maxwell equations at the moment.
  return;

  // Loop over the quadrature points of the cell. No volume source term is
  // considered at the moment, so the interior load vector is zero. The
  // electric test functions of the base function i are its E_real and E_imag
  // shape functions.
  for (unsigned int q = 0; q < scratch_data.n_q_points; ++q)
    {
      for (unsigned int i = 0; i < scratch_data.n_base_functions_test; ++i)
        {
          l_vector[scratch_data
                     .test_local_dof_index_map[ScratchData::E_real][i]] += 0.0;
          l_vector[scratch_data
                     .test_local_dof_index_map[ScratchData::E_imag][i]] += 0.0;
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

  // Local dof indices of the real and imaginary parts of each base function:
  // electric (F) and magnetic (I) test functions, and electric (E_hat) and
  // magnetic (H_hat) skeleton trial functions
  const std::vector<unsigned int> &test_dofs_F_real =
    scratch_data.test_local_dof_index_map[ScratchData::E_real];
  const std::vector<unsigned int> &test_dofs_F_imag =
    scratch_data.test_local_dof_index_map[ScratchData::E_imag];
  const std::vector<unsigned int> &test_dofs_I_real =
    scratch_data.test_local_dof_index_map[ScratchData::H_real];
  const std::vector<unsigned int> &test_dofs_I_imag =
    scratch_data.test_local_dof_index_map[ScratchData::H_imag];
  const std::vector<unsigned int> &trial_dofs_E_hat_real =
    scratch_data.trial_skeleton_local_dof_index_map[ScratchData::E_real];
  const std::vector<unsigned int> &trial_dofs_E_hat_imag =
    scratch_data.trial_skeleton_local_dof_index_map[ScratchData::E_imag];
  const std::vector<unsigned int> &trial_dofs_H_hat_real =
    scratch_data.trial_skeleton_local_dof_index_map[ScratchData::H_real];
  const std::vector<unsigned int> &trial_dofs_H_hat_imag =
    scratch_data.trial_skeleton_local_dof_index_map[ScratchData::H_imag];

  const std::vector<double> &face_JxW = scratch_data.face_JxW;

  LAPACKFullMatrix<double> &B_hat_matrix = copy_data.B_hat_matrix;

  // The skeleton terms only involve the tangential traces of the test and
  // skeleton trial functions, so only the base functions associated with the
  // face contribute. For the test base function i (rows) and the skeleton
  // trial base function j (columns), the blocks are:
  // - IE_hat: <phi_i, n x e_hat_j>
  // - FH_hat: <phi_i, n x h_hat_j>, the same value since the electric and
  //   magnetic skeleton fields share the same base functions
  for (const unsigned int i :
       scratch_data.test_tangential_base_functions[scratch_data.face_no])
    {
      const auto phi_i = scratch_data.phi_test_face[i];

      for (const unsigned int j :
           scratch_data
             .trial_skeleton_tangential_base_functions[scratch_data.face_no])
        {
          const auto n_cross_phi_j =
            scratch_data.n_cross_tangential_phi_trial_skeleton[j];

          double z = 0.;
          for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
            z += face_JxW[q] * (n_cross_phi_j[q] * phi_i[q]);

          add_real_block(B_hat_matrix,
                         test_dofs_I_real[i],
                         test_dofs_I_imag[i],
                         trial_dofs_E_hat_real[j],
                         trial_dofs_E_hat_imag[j],
                         z);

          // The magnetic trace term is replaced by the Robin boundary condition
          // on the faces where it is applied.
          if (!robin_face)
            add_real_block(B_hat_matrix,
                           test_dofs_F_real[i],
                           test_dofs_F_imag[i],
                           trial_dofs_H_hat_real[j],
                           trial_dofs_H_hat_imag[j],
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
TimeHarmonicMaxwellAssemblerRobinBC<dim>::TimeHarmonicMaxwellAssemblerRobinBC(
  const Parameters::TimeHarmonicMaxwell<dim> &time_harmonic_maxwell_parameters,
  const BoundaryConditions::TimeHarmonicMaxwellBoundaryConditions<dim>
                            &boundary_conditions,
  const std::vector<double> &waveguide_ports_electric_amplitudes)
  : time_harmonic_maxwell_parameters(time_harmonic_maxwell_parameters)
  , boundary_conditions(boundary_conditions)
  , waveguide_ports_relative_amplitudes(waveguide_ports_electric_amplitudes)
{
  // In the case of multiple waveguide ports, we normalize their amplitudes by
  // their maximum here before assembling the Robin boundary condition to avoid
  // extra computations as this factor is constant across all quadrature points.
  // This normalization is done to ensure that the waveguide port with the
  // highest amplitude is set to 1, and the other ports are scaled accordingly.
  if (waveguide_ports_relative_amplitudes.empty())
    return;

  const double max_port_electric_amplitude =
    *std::ranges::max_element(waveguide_ports_relative_amplitudes);
  for (double &amplitude : waveguide_ports_relative_amplitudes)
    amplitude /= max_port_electric_amplitude;
}

template <int dim>
std::pair<Tensor<1, dim, std::complex<double>>, std::complex<double>>
TimeHarmonicMaxwellAssemblerRobinBC<dim>::compute_waveguide_port_excitation(
  const Point<dim>           &p,
  const Tensor<1, dim>       &normal,
  const std::complex<double> &effective_electric_permittivity,
  const std::complex<double> &effective_magnetic_permeability,
  const unsigned int          waveguide_port_index) const
{
  // The waveguide port excitation would be completely different in 2D as curls
  // and cross products are not defined the same way as in 3D. Therefore, we do
  // not implement it for now and throw an error if someone tries to use the 2D
  // version.
  if constexpr (dim != 3)
    {
      (void)p;
      (void)normal;
      (void)effective_electric_permittivity;
      (void)effective_magnetic_permeability;
      (void)waveguide_port_index;
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
                                           waveguide_port_index);

  // Gather the relevant waveguide parameters
  const double omega =
    2.0 * PI * time_harmonic_maxwell_parameters.electromagnetic_frequency;
  const auto &waveguide_corners =
    time_harmonic_maxwell_parameters.waveguide_corners[waveguide_port_index];
  const Parameters::WaveguideMode mode =
    time_harmonic_maxwell_parameters.waveguide_mode[waveguide_port_index];
  unsigned int m =
    time_harmonic_maxwell_parameters.mode_order_m[waveguide_port_index];
  unsigned int n =
    time_harmonic_maxwell_parameters.mode_order_n[waveguide_port_index];

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

  const Tensor<1, dim, std::complex<double>> excitation =
    waveguide_ports_relative_amplitudes[waveguide_port_index] *
    (cross_product_3d(normal, H_inc) +
     tangential_component(surface_admittance * E_inc, normal));

  return std::make_pair(excitation, surface_admittance);
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
  if (!scratch_data.face_at_boundary)
    return;

  const BoundaryConditions::BoundaryType bc_type =
    boundary_conditions.type.at(scratch_data.face_boundary_id);
  if (!is_robin_boundary_type(bc_type))
    return;

  const unsigned int waveguide_port_index =
    (bc_type == BoundaryConditions::BoundaryType::waveguide_port) ?
      get_waveguide_port_index(time_harmonic_maxwell_parameters,
                               scratch_data.face_boundary_id) :
      0;

  using ScratchData = TimeHarmonicMaxwellScratchData<dim>;

  const unsigned int n_face_q_points = scratch_data.n_face_q_points;

  // Base functions of the test space associated with the face
  const std::vector<unsigned int> &test_tangential_base_functions =
    scratch_data.test_tangential_base_functions[scratch_data.face_no];
  const unsigned int n_test_tangential_base_functions =
    test_tangential_base_functions.size();

  // Local dof indices of the real and imaginary parts of each base function:
  // electric (F) and magnetic (I) test functions, and electric (E_hat)
  // skeleton trial functions
  const std::vector<unsigned int> &test_dofs_F_real =
    scratch_data.test_local_dof_index_map[ScratchData::E_real];
  const std::vector<unsigned int> &test_dofs_F_imag =
    scratch_data.test_local_dof_index_map[ScratchData::E_imag];
  const std::vector<unsigned int> &test_dofs_I_real =
    scratch_data.test_local_dof_index_map[ScratchData::H_real];
  const std::vector<unsigned int> &test_dofs_I_imag =
    scratch_data.test_local_dof_index_map[ScratchData::H_imag];
  const std::vector<unsigned int> &trial_dofs_E_hat_real =
    scratch_data.trial_skeleton_local_dof_index_map[ScratchData::E_real];
  const std::vector<unsigned int> &trial_dofs_E_hat_imag =
    scratch_data.trial_skeleton_local_dof_index_map[ScratchData::E_imag];

  const std::vector<double> &face_JxW = scratch_data.face_JxW;

  LAPACKFullMatrix<double> &G_matrix     = copy_data.G_matrix;
  LAPACKFullMatrix<double> &B_hat_matrix = copy_data.B_hat_matrix;

  // Surface admittance Y at the face quadrature points, multiplied by JxW:
  // |Y|^2 JxW and Y JxW. It depends on the boundary condition, so it is
  // computed here and not in the scratch data.
  std::vector<double>               JxW_YY_container(n_face_q_points);
  std::vector<std::complex<double>> Y_JxW_container(n_face_q_points);
  for (unsigned int q = 0; q < n_face_q_points; ++q)
    {
      const std::complex<double> boundary_surface_admittance =
        compute_robin_boundary_data(scratch_data,
                                    q,
                                    bc_type,
                                    waveguide_port_index)
          .second;
      JxW_YY_container[q] =
        std::norm(boundary_surface_admittance) * face_JxW[q];
      Y_JxW_container[q] = boundary_surface_admittance * face_JxW[q];
    }

  // Boundary terms of the Gram matrix G coming from the test norm on the Robin
  // boundaries:
  // (n x I + Y^* F_t, n x I~ + Y^* F~_t)_{Gamma_R},
  // where F_t = n x (F x n) is the tangential trace of F and the inner product
  // is conjugated in its first argument. This term only involves the
  // tangential traces of the test functions, which are in L^2(Gamma_R) for the
  // test space. Since |n| = 1, F_t . F~_t = (n x F) . (n x F~), so for the test
  // base functions i (rows) and j (columns), the blocks are:
  // - FF: |Y|^2 <n x phi_j, n x phi_i>
  // - II: <n x phi_j, n x phi_i>
  // - FI: Y <n x phi_j, phi_i>
  // - IF: since G is Hermitian, the block IF(j, i) is conj(FI(i, j)).
  // Only the base functions associated with the face have a tangential trace
  // on the face, so only them contribute. As for the cell terms, we only loop
  // over the pairs i <= j and add the blocks of the pair (j, i) at the same
  // time.
  for (unsigned int position_i = 0;
       position_i < n_test_tangential_base_functions;
       ++position_i)
    {
      const unsigned int i     = test_tangential_base_functions[position_i];
      const auto         phi_i = scratch_data.phi_test_face[i];
      const auto         n_cross_phi_i = scratch_data.n_cross_phi_test_face[i];

      for (unsigned int position_j = position_i;
           position_j < n_test_tangential_base_functions;
           ++position_j)
        {
          const unsigned int j     = test_tangential_base_functions[position_j];
          const auto         phi_j = scratch_data.phi_test_face[j];
          const auto n_cross_phi_j = scratch_data.n_cross_phi_test_face[j];

          double               z_FF = 0.;
          double               z_II = 0.;
          std::complex<double> z_FI_ij(0.);
          std::complex<double> z_FI_ji(0.);
          for (unsigned int q = 0; q < n_face_q_points; ++q)
            {
              // The FF and II blocks share the dot product of the tangential
              // traces
              const double n_cross_j_n_cross_i =
                n_cross_phi_j[q] * n_cross_phi_i[q];

              z_FF += JxW_YY_container[q] * n_cross_j_n_cross_i;
              z_II += face_JxW[q] * n_cross_j_n_cross_i;
              z_FI_ij += Y_JxW_container[q] * (n_cross_phi_j[q] * phi_i[q]);
              z_FI_ji += Y_JxW_container[q] * (n_cross_phi_i[q] * phi_j[q]);
            }

          // Blocks of the pair (i, j)
          add_real_block(G_matrix,
                         test_dofs_F_real[i],
                         test_dofs_F_imag[i],
                         test_dofs_F_real[j],
                         test_dofs_F_imag[j],
                         z_FF);
          add_real_block(G_matrix,
                         test_dofs_I_real[i],
                         test_dofs_I_imag[i],
                         test_dofs_I_real[j],
                         test_dofs_I_imag[j],
                         z_II);
          add_complex_block(G_matrix,
                            test_dofs_F_real[i],
                            test_dofs_F_imag[i],
                            test_dofs_I_real[j],
                            test_dofs_I_imag[j],
                            z_FI_ij);
          add_complex_block(G_matrix,
                            test_dofs_I_real[j],
                            test_dofs_I_imag[j],
                            test_dofs_F_real[i],
                            test_dofs_F_imag[i],
                            std::conj(z_FI_ij));

          // Blocks of the pair (j, i), i != j because the diagonal blocks are
          // already filled by the pair (i, j)
          if (j != i)
            {
              add_real_block(G_matrix,
                             test_dofs_F_real[j],
                             test_dofs_F_imag[j],
                             test_dofs_F_real[i],
                             test_dofs_F_imag[i],
                             z_FF);
              add_real_block(G_matrix,
                             test_dofs_I_real[j],
                             test_dofs_I_imag[j],
                             test_dofs_I_real[i],
                             test_dofs_I_imag[i],
                             z_II);
              add_complex_block(G_matrix,
                                test_dofs_F_real[j],
                                test_dofs_F_imag[j],
                                test_dofs_I_real[i],
                                test_dofs_I_imag[i],
                                z_FI_ji);
              add_complex_block(G_matrix,
                                test_dofs_I_real[i],
                                test_dofs_I_imag[i],
                                test_dofs_F_real[j],
                                test_dofs_F_imag[j],
                                std::conj(z_FI_ji));
            }
        }
    }

  // Robin term of the skeleton bilinear form B_hat, which replaces the
  // magnetic trace term. For the test base function i (rows) and the skeleton
  // trial base function j (columns), the FE_hat block is
  // -Y < phi_i, n x (e_hat_j x n) >. It only involves the tangential trace of
  // the test functions, so only the base functions associated with the face
  // contribute.
  for (const unsigned int i : test_tangential_base_functions)
    {
      const auto phi_i = scratch_data.phi_test_face[i];

      for (const unsigned int j :
           scratch_data
             .trial_skeleton_tangential_base_functions[scratch_data.face_no])
        {
          const auto tangential_phi_j =
            scratch_data.tangential_phi_trial_skeleton[j];

          std::complex<double> z(0.);
          for (unsigned int q = 0; q < n_face_q_points; ++q)
            z -= Y_JxW_container[q] * (tangential_phi_j[q] * phi_i[q]);

          add_complex_block(B_hat_matrix,
                            test_dofs_F_real[i],
                            test_dofs_F_imag[i],
                            trial_dofs_E_hat_real[j],
                            trial_dofs_E_hat_imag[j],
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
  if (!scratch_data.face_at_boundary)
    return;

  const BoundaryConditions::BoundaryType bc_type =
    boundary_conditions.type.at(scratch_data.face_boundary_id);
  if (!is_robin_boundary_type(bc_type))
    return;

  const unsigned int waveguide_port_index =
    (bc_type == BoundaryConditions::BoundaryType::waveguide_port) ?
      get_waveguide_port_index(time_harmonic_maxwell_parameters,
                               scratch_data.face_boundary_id) :
      0;

  using ScratchData = TimeHarmonicMaxwellScratchData<dim>;

  const unsigned int n_face_q_points = scratch_data.n_face_q_points;

  // Local dof indices of the real and imaginary parts of each base function of
  // the electric test functions (F)
  const std::vector<unsigned int> &test_dofs_F_real =
    scratch_data.test_local_dof_index_map[ScratchData::E_real];
  const std::vector<unsigned int> &test_dofs_F_imag =
    scratch_data.test_local_dof_index_map[ScratchData::E_imag];

  Vector<double> &l_vector = copy_data.l_vector;

  // Tangential component of the excitation g at the face quadrature points,
  // multiplied by JxW. The left-hand side of the Robin boundary condition
  // n x H + Y n x (E x n) = g is tangential to the boundary, so only the
  // tangential component of g is meaningful. It is also the only component
  // for which -<F, g> is bounded on the test space, which only controls the
  // tangential trace of F on the Robin boundaries.
  std::vector<Tensor<1, dim, std::complex<double>>> g_inc_JxW_container(
    n_face_q_points);
  for (unsigned int q = 0; q < n_face_q_points; ++q)
    g_inc_JxW_container[q] =
      tangential_component(compute_robin_boundary_data(
                             scratch_data, q, bc_type, waveguide_port_index)
                             .first,
                           scratch_data.face_normals[q]) *
      scratch_data.face_JxW[q];

  // Excitation term -<F,g>. For the test base function i, the entries of the
  // real and imaginary parts of the electric test function are the real and
  // imaginary parts of w = <phi_i, g>, since Re(w conj(s)) is Re(w) for s = 1
  // and Im(w) for s equal to the imaginary unit. Since g is tangential, only
  // the base functions associated with the face contribute.
  for (const unsigned int i :
       scratch_data.test_tangential_base_functions[scratch_data.face_no])
    {
      const auto phi_i = scratch_data.phi_test_face[i];

      std::complex<double> w(0.);
      for (unsigned int q = 0; q < n_face_q_points; ++q)
        w += g_inc_JxW_container[q] * phi_i[q];

      l_vector[test_dofs_F_real[i]] -= w.real();
      l_vector[test_dofs_F_imag[i]] -= w.imag();
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
  const unsigned int                        waveguide_port_index);
template std::pair<Tensor<1, 3, std::complex<double>>,
                   Tensor<1, 3, std::complex<double>>>
compute_waveguide_port_incident_fields<3>(
  const Parameters::TimeHarmonicMaxwell<3> &time_harmonic_maxwell_parameters,
  const Point<3>                           &p,
  const Tensor<1, 3>                       &normal,
  const std::complex<double>               &effective_electric_permittivity,
  const std::complex<double>               &effective_magnetic_permeability,
  const unsigned int                        waveguide_port_index);

template unsigned int
get_waveguide_port_index<2>(
  const Parameters::TimeHarmonicMaxwell<2> &time_harmonic_maxwell_parameters,
  const types::boundary_id                  boundary_id);
template unsigned int
get_waveguide_port_index<3>(
  const Parameters::TimeHarmonicMaxwell<3> &time_harmonic_maxwell_parameters,
  const types::boundary_id                  boundary_id);

template class TimeHarmonicMaxwellAssemblerCore<2>;
template class TimeHarmonicMaxwellAssemblerCore<3>;
template class TimeHarmonicMaxwellAssemblerSkeleton<2>;
template class TimeHarmonicMaxwellAssemblerSkeleton<3>;
template class TimeHarmonicMaxwellAssemblerRobinBC<2>;
template class TimeHarmonicMaxwellAssemblerRobinBC<3>;
