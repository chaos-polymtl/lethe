// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the parsing of the QMOM subsection, including its "phase
 * evolution", "pharma" and "physical constants" subsections, and of the
 * functions prescribing the QMOM fields. It checks that:
 * - an empty QMOM subsection gives the default values;
 * - every entry of the subsections is parsed into its own member;
 * - an odd number of moments uses one more quadrature node;
 * - a number of scalar names that differs from the number of scalars and an
 *   invalid scalar name throw;
 * - the number of expressions of the scalars and of the particle size
 *   distribution is parsed, and a number of distribution expressions that
 *   does not match the distribution type throws.
 */

#include <core/parameters_multiphysics.h>
#include <core/utilities.h>

#include <deal.II/base/parameter_handler.h>

#include <../tests/tests.h>

#include <exception>
#include <string>

/**
 * @brief Declare and parse a QMOM subsection.
 *
 * @param[in] input Content of the QMOM subsection.
 *
 * @return The parsed QMOM parameters.
 */
Parameters::QMOM
parse_qmom(const std::string &input)
{
  ParameterHandler       prm;
  const Parameters::QMOM declarations;
  declarations.declare_parameters(prm);

  prm.parse_input_from_string("subsection QMOM\n" + input + "end\n");

  Parameters::QMOM qmom;
  qmom.parse_parameters(prm);
  return qmom;
}

/**
 * @brief Print all the QMOM parameters.
 *
 * @param[in] qmom QMOM parameters.
 */
void
print_qmom(const Parameters::QMOM &qmom)
{
  deallog << "number of scalars: " << qmom.number_of_scalars << std::endl;
  deallog << "scalar names:";
  for (const std::string &scalar_name : qmom.scalar_names)
    deallog << " " << scalar_name;
  deallog << std::endl;
  deallog << "number of moments: " << qmom.number_of_moments << std::endl;
  deallog << "number of nodes: " << qmom.get_number_of_nodes() << std::endl;
  deallog << "splitting scheme: " << enum_to_string(qmom.splitting_scheme)
          << std::endl;
  deallog << "enable nucleation: " << qmom.enable_nucleation << std::endl;
  deallog << "enable coalescence: " << qmom.enable_coalescence << std::endl;
  deallog << "enable growth: " << qmom.enable_growth << std::endl;
  deallog << "enable adaptive node count: "
          << qmom.phase_evolution.enable_adaptive_node_count << std::endl;

  const Parameters::QMOMPharma &pharma = qmom.pharma;
  deallog << "operating temperature: " << pharma.operating_temperature
          << std::endl;
  deallog << "flow rate ratio: " << pharma.flow_rate_ratio << std::endl;
  deallog << "total flow rate: " << pharma.total_flow_rate << std::endl;
  deallog << "mixing time scale: " << pharma.mixing_time_scale << std::endl;
  deallog << "initial concentration organic: "
          << pharma.initial_concentration_organic << std::endl;
  deallog << "ion concentration mixture: " << pharma.ion_concentration_mixture
          << std::endl;
  deallog << "density aqueous: " << pharma.density_aqueous << std::endl;
  deallog << "molar mass aqueous: " << pharma.molar_mass_aqueous << std::endl;
  deallog << "mole fraction solubility aqueous: "
          << pharma.mole_fraction_solubility_aqueous << std::endl;
  deallog << "permittivity aqueous: " << pharma.permittivity_aqueous
          << std::endl;
  deallog << "density organic: " << pharma.density_organic << std::endl;
  deallog << "molar mass organic: " << pharma.molar_mass_organic << std::endl;
  deallog << "mole fraction solubility organic: "
          << pharma.mole_fraction_solubility_organic << std::endl;
  deallog << "permittivity organic: " << pharma.permittivity_organic
          << std::endl;
  deallog << "viscosity mixture: " << pharma.viscosity_mixture << std::endl;
  deallog << "shape factor: " << pharma.shape_factor << std::endl;

  const Parameters::QMOMPhysicalConstants &constants = qmom.physical_constants;
  deallog << "boltzmann constant: " << constants.boltzmann_constant
          << std::endl;
  deallog << "avogadro number: " << constants.avogadro_number << std::endl;
  deallog << "elementary charge: " << constants.elementary_charge << std::endl;
  deallog << "permittivity vacuum: " << constants.permittivity_vacuum
          << std::endl;
}

/**
 * @brief Parse a QMOM subsection and report whether the parsing threw.
 *
 * @param[in] label Description of the scenario printed in the log.
 * @param[in] input Content of the QMOM subsection.
 */
void
check_qmom_throws(const std::string &label, const std::string &input)
{
  try
    {
      parse_qmom(input);
      deallog << label << ": did not throw" << std::endl;
    }
  catch (const std::exception &)
    {
      deallog << label << ": threw" << std::endl;
    }
}

/**
 * @brief Declare and parse the functions prescribing the QMOM fields.
 *
 * @param[in] input Content of the subsection holding the functions.
 *
 * @return The parsed functions.
 */
Parameters::QMOMFunctions<2>
parse_functions(const std::string &input)
{
  ParameterHandler prm;
  prm.enter_subsection("functions");
  Parameters::QMOMFunctions<2>::declare_parameters(prm);
  prm.leave_subsection();

  prm.parse_input_from_string("subsection functions\n" + input + "end\n");

  Parameters::QMOMFunctions<2> functions;
  prm.enter_subsection("functions");
  functions.parse_parameters(prm);
  prm.leave_subsection();
  return functions;
}

/**
 * @brief Print the distribution type and the number of expressions of the
 * functions prescribing the QMOM fields.
 *
 * @param[in] functions Functions to print.
 */
void
print_functions(const Parameters::QMOMFunctions<2> &functions)
{
  deallog << "scalar expressions: " << functions.scalars->n_components
          << std::endl;
  deallog << "distribution type: "
          << enum_to_string(functions.distribution_type) << std::endl;
  deallog << "distribution expressions: "
          << functions.distribution->n_components << std::endl;
}

void
test()
{
  deallog << "Default values" << std::endl;
  print_qmom(parse_qmom(""));

  deallog << "Parsed values" << std::endl;
  print_qmom(parse_qmom("  set number of scalars  = 3\n"
                        "  set scalar names       = ethanol, lipids, ions\n"
                        "  set number of moments  = 5\n"
                        "  set splitting scheme   = strang\n"
                        "  set enable nucleation  = true\n"
                        "  set enable coalescence = true\n"
                        "  set enable growth      = true\n"
                        "  subsection phase evolution\n"
                        "    set enable adaptive node count = true\n"
                        "  end\n"
                        "  subsection pharma\n"
                        "    set operating temperature            = 300\n"
                        "    set flow rate ratio                  = 2.5\n"
                        "    set total flow rate                  = 1e-7\n"
                        "    set mixing time scale                = 0.01\n"
                        "    set initial concentration organic    = 12.5\n"
                        "    set ion concentration mixture        = 150\n"
                        "    set density aqueous                  = 1000\n"
                        "    set molar mass aqueous               = 0.018\n"
                        "    set mole fraction solubility aqueous = 1e-8\n"
                        "    set permittivity aqueous             = 80\n"
                        "    set density organic                  = 790\n"
                        "    set molar mass organic               = 0.046\n"
                        "    set mole fraction solubility organic = 3e-3\n"
                        "    set permittivity organic             = 25\n"
                        "    set viscosity mixture                = 1e-3\n"
                        "    set shape factor                     = 0.5\n"
                        "  end\n"
                        "  subsection physical constants\n"
                        "    set boltzmann constant  = 1.4e-23\n"
                        "    set avogadro number     = 6e23\n"
                        "    set elementary charge   = 1.6e-19\n"
                        "    set permittivity vacuum = 8.9e-12\n"
                        "  end\n"));

  check_qmom_throws("Too few scalar names",
                    "  set number of scalars = 2\n"
                    "  set scalar names      = ethanol\n");
  check_qmom_throws("Invalid scalar name",
                    "  set number of scalars = 1\n"
                    "  set scalar names      = free lipids\n");

  deallog << "Default functions" << std::endl;
  print_functions(parse_functions(""));

  deallog << "Lognormal distribution" << std::endl;
  print_functions(
    parse_functions("  set Function expression = 1; x\n"
                    "  subsection distribution\n"
                    "    set type                = lognormal\n"
                    "    set Function expression = 1e18; log(1e-7); 0.3\n"
                    "  end\n"));

  deallog << "Chi squared distribution" << std::endl;
  print_functions(parse_functions("  subsection distribution\n"
                                  "    set type                = chi_squared\n"
                                  "    set Function expression = 1e18; 3\n"
                                  "  end\n"));

  try
    {
      parse_functions("  subsection distribution\n"
                      "    set type                = gamma\n"
                      "    set Function expression = 1e18; 2\n"
                      "  end\n");
      deallog << "Wrong number of distribution expressions: did not throw"
              << std::endl;
    }
  catch (const std::exception &)
    {
      deallog << "Wrong number of distribution expressions: threw" << std::endl;
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
