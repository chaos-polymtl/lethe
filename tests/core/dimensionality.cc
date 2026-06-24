// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the Dimensionality scaling factors produced by
 * Dimensionality::define_all_scales().
 *
 * Each test case sets (L, M, T, I, theta) directly on the struct and calls
 * define_all_scales(), then checks every scaling member against the expected
 * analytical value.
 *
 * Convention used in dimensionality.cc:
 *   scaling = (SI unit expressed in new base units)^{-1}
 * i.e. a physical quantity Q_SI is converted to Q_new via
 *   Q_new = Q_SI * scaling.
 *
 * Three scenarios are tested:
 *  1. SI base (all scales = 1)          → all scalings must be 1
 *  2. L = 0.01, M = 0.001               → Verify the CGS system
 *  3. Custom: L = 3.7e-2, M = 4.2e1, T = 0.8, I = 1.3, theta = 2.1 → Verify a
 * custom set of units
 */

// Lethe
#include <core/dimensionality.h>

// Tests (with common definitions)
#include <../tests/tests.h>

// ── helpers ──────────────────────────────────────────────────────────────────

/** Relative tolerance for floating-point comparisons. */
constexpr double tol = 1e-14;

bool
approx_equal(double a, double b)
{
  return std::abs(a - b) <= tol * std::max(1.0, std::abs(b));
}

/**
 * Build a Dimensionality object with the given base-unit values, call
 * define_all_scales(), and return it.
 */
Parameters::Dimensionality
make_dim(double L     = 1.0,
         double M     = 1.0,
         double T     = 1.0,
         double I     = 1.0,
         double theta = 1.0)
{
  Parameters::Dimensionality d;
  d.length           = L;
  d.mass             = M;
  d.time             = T;
  d.electric_current = I;
  d.temperature      = theta;
  d.define_all_scales();
  return d;
}

// ── individual checks
// ─────────────────────────────────────────────────────────

void
check(const std::string &label, double computed, double expected)
{
  if (approx_equal(computed, expected))
    {
      deallog << "  OK  " << label << " = " << computed << std::endl;
    }
  else
    {
      deallog << "  FAIL " << label << ": got " << computed << ", expected "
              << expected << std::endl;
    }
}

// ── test functions
// ────────────────────────────────────────────────────────────

/**
 * Scenario 1 – SI base units (L=M=T=I=theta=1).
 * Every scaling factor must equal 1.
 */
void
test_si_base()
{
  deallog << "--- Scenario 1: SI base units (all = 1) ---" << std::endl;

  auto d = make_dim();

  check("density_scaling", d.density_scaling, 1.0);
  check("specific_gas_constant_scaling", d.specific_gas_constant_scaling, 1.0);
  check("viscosity_scaling", d.viscosity_scaling, 1.0);
  check("specific_heat_scaling", d.specific_heat_scaling, 1.0);
  check("thermal_conductivity_scaling", d.thermal_conductivity_scaling, 1.0);
  check("enthalpy_scaling", d.enthalpy_scaling, 1.0);
  check("diffusivity_scaling", d.diffusivity_scaling, 1.0);
  check("thermal_expansion_scaling", d.thermal_expansion_scaling, 1.0);
  check("surface_tension_scaling", d.surface_tension_scaling, 1.0);
  check("surface_tension_gradient_scaling",
        d.surface_tension_gradient_scaling,
        1.0);
  check("electric_amplitude_scaling", d.electric_amplitude_scaling, 1.0);
  check("magnetic_amplitude_scaling", d.magnetic_amplitude_scaling, 1.0);
  check("vacuum_permittivity_scaling", d.vacuum_permittivity_scaling, 1.0);
  check("vacuum_permeability_scaling", d.vacuum_permeability_scaling, 1.0);
}

/**
 * Scenario 2 – Length in centimetres: L = 0.01, M = 0.001, T=I=theta=1.
 * Note that this is not exactly the CGS system, since in CGS, the electric
 * current is removed as a fundamental dimension using the definition of the
 * statcoulomb. This would imply to rewrite the electromagnetic solvers in terms
 * of the statcoulomb, which is not the desired behavior. So we simply have an
 * SI-like system with L in cm and M in g.
 */
void
test_CGS()
{
  deallog << "--- Scenario 2: CGS ---" << std::endl;

  const double I     = 1.0;
  const double L     = 0.01;
  const double M     = 0.001;
  const double T     = 1.0;
  const double theta = 1.0;
  auto         d     = make_dim(L, M, T, I, theta);

  // density: L^3/M
  check("density_scaling", d.density_scaling, L * L * L / M);

  // specific gas constant: 1/L^2 * T^2 * theta  (T=theta=1)
  check("specific_gas_constant_scaling",
        d.specific_gas_constant_scaling,
        1.0 / (L * L));

  // viscosity: 1/L^2 * T
  check("viscosity_scaling", d.viscosity_scaling, 1.0 / (L * L));

  // specific heat: same formula as specific gas constant
  check("specific_heat_scaling", d.specific_heat_scaling, 1.0 / (L * L));

  // thermal conductivity: 1/(M*L) * T^3 * theta
  check("thermal_conductivity_scaling",
        d.thermal_conductivity_scaling,
        1.0 / (M * L));

  // enthalpy: 1/(M*L^2) * T^2
  check("enthalpy_scaling", d.enthalpy_scaling, 1.0 / (M * L * L));

  // diffusivity: 1/L^2 * T
  check("diffusivity_scaling", d.diffusivity_scaling, 1.0 / (L * L));

  // thermal expansion: theta
  check("thermal_expansion_scaling", d.thermal_expansion_scaling, 1.0);

  // surface tension: T^2/M
  check("surface_tension_scaling", d.surface_tension_scaling, 1.0 / M);

  // surface tension gradient: theta*T^2/M
  check("surface_tension_gradient_scaling",
        d.surface_tension_gradient_scaling,
        1.0 / M);

  // electric field E : = I*T^3/(M*L)
  check("electric_amplitude_scaling",
        d.electric_amplitude_scaling,
        1.0 / (M * L));

  // magnetic field H : L/I
  check("magnetic_amplitude_scaling", d.magnetic_amplitude_scaling, L);

  // vacuum permittivity : M*L^3/(I^2*T^4)
  check("vacuum_permittivity_scaling",
        d.vacuum_permittivity_scaling,
        1.0 * M * L * L * L / 1.0 / 1.0 / 1.0 / 1.0);

  // vacuum permeability : T^2*I^2/(M*L)
  check("vacuum_permeability_scaling",
        d.vacuum_permeability_scaling,
        1.0 / (M * L));

  // electromagnetic frequency scaling = L/c_SI
  const double c_SI = 299792458.0;
  check("electromagnetic_frequency_scaling",
        d.electromagnetic_frequency_scaling,
        L / c_SI);
}

/**
 * Scenario 3 – Randomized scaling & outputs: L = 3.7e-2, M = 4.2e1, T = 0.8, I
 * = 1.3, theta = 2.1.
 */
void
test_random_scaled_outputs()
{
  deallog << "--- Scenario 3b: Randomized scaling & outputs ---" << std::endl;

  const double L     = 3.7e-2;
  const double M     = 4.2e1;
  const double T     = 0.8;
  const double I     = 1.3;
  const double theta = 2.1;

  auto d = make_dim(L, M, T, I, theta);

  // density: L^3/M
  check("density_scaling", d.density_scaling, L * L * L / M);

  // specific gas constant: 1/L^2 * T^2 * theta
  check("specific_gas_constant_scaling",
        d.specific_gas_constant_scaling,
        1. / L / L * T * T * theta);

  // viscosity: 1/L^2 * T
  check("viscosity_scaling", d.viscosity_scaling, 1. / L / L * T);

  // specific heat: same formula as specific gas constant
  check("specific_heat_scaling",
        d.specific_heat_scaling,
        1. / L / L * T * T * theta);

  // thermal conductivity: 1/(M*L) * T^3 * theta
  check("thermal_conductivity_scaling",
        d.thermal_conductivity_scaling,
        1. / M / L * T * T * T * theta);

  // enthalpy: 1/(M*L^2) * T^2
  check("enthalpy_scaling", d.enthalpy_scaling, 1. / M / L / L * T * T);

  // diffusivity: 1/L^2 * T
  check("diffusivity_scaling", d.diffusivity_scaling, 1. / L / L * T);

  // thermal expansion: theta
  check("thermal_expansion_scaling", d.thermal_expansion_scaling, 1. * theta);

  // surface tension: T^2/M
  check("surface_tension_scaling", d.surface_tension_scaling, 1. * T * T / M);

  // surface tension gradient: theta*T^2/M
  check("surface_tension_gradient_scaling",
        d.surface_tension_gradient_scaling,
        1. * theta * T * T / M);

  // cahn-hilliard mobility: 1. * M / L / L / L / T
  check("cahn_hilliard_mobility_scaling",
        d.cahn_hilliard_mobility_scaling,
        1. * M / L / L / L / T);

  // cahn-hilliard epsilon: 1. / L
  check("cahn_hilliard_epsilon_scaling",
        d.cahn_hilliard_epsilon_scaling,
        1. / L);


  // electric field E : scaling = I*T^3/(M*L)
  check("electric_amplitude_scaling",
        d.electric_amplitude_scaling,
        1 * I * T * T * T / M / L);

  // magnetic field H: scaling = L/I
  check("magnetic_amplitude_scaling", d.magnetic_amplitude_scaling, 1 * L / I);

  // vacuum permittivity :
  check("vacuum_permittivity_scaling",
        d.vacuum_permittivity_scaling,
        1 * M * L * L * L / I / I / T / T / T / T);

  // vacuum permeability : T^2*I^2/(M*L)
  check("vacuum_permeability_scaling",
        d.vacuum_permeability_scaling,
        1 * T * T * I * I / M / L);

  // electromagnetic frequency scaling = L/c_SI
  const double c_SI = 299792458.0;
  check("electromagnetic_frequency_scaling",
        d.electromagnetic_frequency_scaling,
        L / c_SI);
}

// ── main
// ──────────────────────────────────────────────────────────────────────

int
main()
{
  try
    {
      initlog();

      deallog << "Beginning dimensionality scaling tests" << std::endl;

      test_si_base();
      test_CGS();
      test_random_scaled_outputs();

      deallog << "OK" << std::endl;
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}
