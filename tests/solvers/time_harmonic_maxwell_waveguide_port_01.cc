// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the incident fields of a rectangular waveguide port,
 * for a TE and a TM mode in a lossy medium (complex permittivity and
 * permeability) and for both orientations of the outward normal of the port
 * with respect to the frame defined by the waveguide corners. It checks that:
 * - the fields form a wave entering the domain, i.e. n x H_inc = Y E_inc,t on
 *   the port, where n is the outward normal, E_inc,t is the tangential
 *   component of the incident electric field and Y is the modal admittance
 *   (k_l / (omega mu) for TE modes and omega epsilon / k_l for TM modes);
 * - the tangential electric field on the port is real, with a maximum
 *   amplitude of 1;
 * - the fields do not depend on the ordering of corner 1 and corner 2, when
 *   the mode orders m and n are swapped accordingly.
 *
 * For each case, it outputs the tangential electric field and n x H_inc / Y at
 * two points, which must be equal, the maximum over the port of
 * |n x H_inc - Y E_inc,t| (0), of |Im(E_inc,t)| (0) and of |E_inc,t| (1), and
 * the maximum change of the fields when corner 1 and corner 2 are swapped (0).
 */

// Lethe
#include <core/parameters_multiphysics.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <solvers/time_harmonic_maxwell_assemblers.h>

// Tests
#include <../tests/tests.h>

void
test()
{
  constexpr int dim = 3;

  // Lossy medium filling the waveguide
  const std::complex<double> permittivity(4.0, 1.0);
  const std::complex<double> permeability(1.0, 0.2);

  // Rectangular waveguide port of size 1 x 0.5 in the plane x = 0. The frame
  // defined by the corners is e_t1 = e_y, e_t2 = e_z and e_t1 x e_t2 = e_x.
  const double width  = 1.0;
  const double height = 0.5;

  Parameters::TimeHarmonicMaxwell<dim> parameters;
  parameters.electromagnetic_frequency  = 1.0;
  parameters.number_of_waveguide_inlets = 1;
  parameters.waveguide_boundary_ids     = {0};
  parameters.waveguide_corners.resize(1);

  parameters.waveguide_corners[0][0] = Point<dim>(0., 0., 0.);
  parameters.waveguide_corners[0][1] = Point<dim>(0., width, 0.);
  parameters.waveguide_corners[0][2] = Point<dim>(0., 0., height);
  parameters.waveguide_corners[0][3] = Point<dim>(0., width, height);

  // Same port with corner 1 and corner 2 swapped (so e_t1 = e_z, e_t2 = e_y and
  // e_t1 x e_t2 = -e_x).
  Parameters::TimeHarmonicMaxwell<dim> swapped_parameters = parameters;

  swapped_parameters.waveguide_corners[0][0] = Point<dim>(0., 0., 0.);
  swapped_parameters.waveguide_corners[0][1] = Point<dim>(0., 0., height);
  swapped_parameters.waveguide_corners[0][2] = Point<dim>(0., width, 0.);
  swapped_parameters.waveguide_corners[0][3] = Point<dim>(0., width, height);

  const double omega = 2. * numbers::PI * parameters.electromagnetic_frequency;

  // The real and imaginary parts at the level of the roundoff errors are
  // replaced by 0, so the output does not depend on the platform
  const double tolerance = 1e-12;
  const auto   filter    = [tolerance](const std::complex<double> value) {
    return std::complex<double>(
      filter_out_small_numbers(value.real(), tolerance),
      filter_out_small_numbers(value.imag(), tolerance));
  };
  const auto print_vector =
    [&filter](const Tensor<1, dim, std::complex<double>> &vector) {
      deallog << filter(vector[0]) << " " << filter(vector[1]) << " "
              << filter(vector[2]) << std::endl;
    };

  // Points where the fields are printed: a generic point and the point where
  // the tangential electric field of both modes reaches its maximum amplitude
  const std::vector<Point<dim>> printed_points = {Point<dim>(0., 0.3, 0.1),
                                                  Point<dim>(0., 0.5, 0.)};

  struct Mode
  {
    std::string               name;
    Parameters::WaveguideMode type;
    unsigned int              m;
    unsigned int              n;
  };
  const std::vector<Mode> modes = {
    {"TE10", Parameters::WaveguideMode::TE, 1, 0},
    {"TM11", Parameters::WaveguideMode::TM, 1, 1}};

  for (const Mode &mode : modes)
    {
      parameters.waveguide_mode         = {mode.type};
      parameters.mode_order_m           = {mode.m};
      parameters.mode_order_n           = {mode.n};
      swapped_parameters.waveguide_mode = {mode.type};
      swapped_parameters.mode_order_m   = {mode.n};
      swapped_parameters.mode_order_n   = {mode.m};

      // Modal admittance of the mode in the lossy medium
      const double k_c_squared =
        Utilities::fixed_power<2>(mode.m * numbers::PI / width) +
        Utilities::fixed_power<2>(mode.n * numbers::PI / height);
      const std::complex<double> k_l =
        std::sqrt(omega * omega * permittivity * permeability - k_c_squared);
      const std::complex<double> admittance =
        (mode.type == Parameters::WaveguideMode::TE) ?
          k_l / (omega * permeability) :
          omega * permittivity / k_l;

      deallog << mode.name << " mode, modal admittance Y = " << admittance
              << std::endl;

      // The outward normal is either opposite to e_t1 x e_t2 (no flip of the
      // corner frame) or along it (flip of the corner frame)
      for (const double normal_direction : {-1., 1.})
        {
          Tensor<1, dim> normal;
          normal[0] = normal_direction;

          // The sampling points include the points where the tangential
          // electric field of both modes reaches its maximum amplitude
          double max_relation_error  = 0.;
          double max_imaginary_part  = 0.;
          double max_amplitude       = 0.;
          double max_ordering_change = 0.;
          for (unsigned int i = 0; i <= 10; ++i)
            for (unsigned int j = 0; j <= 10; ++j)
              {
                const Point<dim> point(0., 0.1 * i * width, 0.1 * j * height);
                const auto [E_inc, H_inc] =
                  compute_waveguide_port_incident_fields(
                    parameters, point, normal, permittivity, permeability, 0);
                const auto [E_inc_swapped, H_inc_swapped] =
                  compute_waveguide_port_incident_fields(swapped_parameters,
                                                         point,
                                                         normal,
                                                         permittivity,
                                                         permeability,
                                                         0);

                const Tensor<1, dim, std::complex<double>> E_inc_t =
                  tangential_component(E_inc, normal);
                const Tensor<1, dim, std::complex<double>> n_cross_H_inc =
                  cross_product_3d(normal, H_inc);

                max_relation_error =
                  std::max(max_relation_error,
                           (n_cross_H_inc - admittance * E_inc_t).norm());
                for (unsigned int d = 0; d < dim; ++d)
                  max_imaginary_part =
                    std::max(max_imaginary_part, std::abs(E_inc_t[d].imag()));
                max_amplitude = std::max(max_amplitude, E_inc_t.norm());
                max_ordering_change =
                  std::max(max_ordering_change,
                           (E_inc - E_inc_swapped).norm() +
                             (H_inc - H_inc_swapped).norm());
              }

          deallog << "  Outward normal "
                  << ((normal_direction > 0.) ? "+e_x" : "-e_x") << std::endl;

          // Incident tangential electric field and n x H_inc / Y, which must
          // be equal for a wave entering the domain
          for (const Point<dim> &point : printed_points)
            {
              const auto [E_inc, H_inc] =
                compute_waveguide_port_incident_fields(
                  parameters, point, normal, permittivity, permeability, 0);
              deallog << "    Point " << point << std::endl;
              deallog << "      E_inc,t       : ";
              print_vector(tangential_component(E_inc, normal));
              deallog << "      n x H_inc / Y : ";
              print_vector(cross_product_3d(normal, H_inc) / admittance);
            }

          deallog << "    Max |n x H_inc - Y E_inc,t| over the port : "
                  << filter_out_small_numbers(max_relation_error, tolerance)
                  << std::endl;
          deallog << "    Max |Im(E_inc,t)| over the port           : "
                  << filter_out_small_numbers(max_imaginary_part, tolerance)
                  << std::endl;
          deallog << "    Max |E_inc,t| over the port               : "
                  << max_amplitude << std::endl;
          deallog << "    Max change with the corner ordering       : "
                  << filter_out_small_numbers(max_ordering_change, tolerance)
                  << std::endl;
        }
    }
}

int
main()
{
  try
    {
      initlog();
      test();
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
