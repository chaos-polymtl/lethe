// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests that the incident fields of a rectangular waveguide
 * port form a wave entering the domain, i.e. that they satisfy n x H_inc = Y
 * E_inc,t on the port, where n is the outward normal of the port, E_inc,t is
 * the tangential component of the incident electric field and Y is the modal
 * admittance (k_l / (omega mu) for TE modes and omega epsilon / k_l for TM
 * modes). The relation is checked for a TE and a TM mode in a lossy medium
 * (complex permittivity and permeability) and for both orientations of the
 * outward normal with respect to the frame defined by the waveguide corners.
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

  const double omega = 2. * numbers::PI * parameters.electromagnetic_frequency;

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
      parameters.waveguide_mode = {mode.type};
      parameters.mode_order_m   = {mode.m};
      parameters.mode_order_n   = {mode.n};

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

      // The outward normal is either opposite to e_t1 x e_t2 (no swap of the
      // corner frame) or along it (swap of the corner frame)
      for (const double normal_direction : {-1., 1.})
        {
          Tensor<1, dim> normal;
          normal[0] = normal_direction;

          double max_error     = 0.;
          double max_reference = 0.;
          for (const double y : {0.2, 0.5, 0.8})
            for (const double z : {0.1, 0.25, 0.4})
              {
                const Point<dim> point(0., y, z);
                const auto [E_inc, H_inc] =
                  compute_waveguide_port_incident_fields(
                    parameters, point, normal, permittivity, permeability, 0);

                const Tensor<1, dim, std::complex<double>> n_cross_H_inc =
                  cross_product_3d(normal, H_inc);
                const Tensor<1, dim, std::complex<double>> Y_E_inc_t =
                  tangential_component(admittance * E_inc, normal);

                max_error =
                  std::max(max_error, (n_cross_H_inc - Y_E_inc_t).norm());
                max_reference = std::max(max_reference, Y_E_inc_t.norm());
              }

          deallog << mode.name << ", outward normal "
                  << (normal_direction > 0. ? "+" : "-")
                  << "e_x: n x H_inc = Y E_inc,t is "
                  << (max_error < 1e-12 * max_reference ? "satisfied" :
                                                          "NOT satisfied")
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
