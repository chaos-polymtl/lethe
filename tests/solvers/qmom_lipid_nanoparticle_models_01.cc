// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the models of the lipid nanoparticle production of the QMOM
 * physics (solubility, growth, nucleation and coalescence) and the consumption
 * of the free lipids by the particles.
 *
 * The reference values were computed with the 0D implementation of the same
 * models (QMOM_pharma_arena_0D of the DG-QMOM repository of the CHAOS
 * laboratory), with its default parameters: cholesterol, a flow rate ratio of
 * 3 with instantaneous mixing, and a free lipid concentration of 11.73 / 4
 * kg/m^3. The flow rate ratio of 3 corresponds to an ethanol concentration of
 * rho_o / 4. The 0D implementation uses an elementary charge of 1.6021766e-19
 * C, a mole fraction solubility in ethanol of 2.6e-4 + 22 * 1.2e-4 and an ion
 * concentration in the mixture of 25 * 3 / 4 mol/m^3, which are set here as
 * well.
 */

// Lethe
#include <core/parameters_multiphysics.h>

#include <solvers/qmom_source_terms.h>

// Deal.II
#include <deal.II/base/exceptions.h>

// magic_enum
#include <magic_enum.hpp>

// Tests
#include <../tests/tests.h>

// Standard
#include <array>
#include <cmath>
#include <exception>
#include <string>
#include <vector>

namespace
{
  /// Relative tolerance of the comparisons with the 0D implementation. The
  /// models are the same, so only round-off separates them.
  constexpr double reference_tolerance = 1e-10;

  /// Number of moments of the process
  constexpr unsigned int n_moments = 8;

  /**
   * @brief Throw if a verification fails.
   *
   * @param[in] condition Result of the verification.
   * @param[in] description Description of the verification.
   */
  void
  verify(const bool condition, const std::string &description)
  {
    AssertThrow(condition, ExcMessage("Verification failed: " + description));
  }

  /**
   * @brief Verify that a value matches its reference within the relative
   * tolerance.
   *
   * @param[in] value Value computed by Lethe.
   * @param[in] reference Value computed by the 0D implementation.
   * @param[in] description Description of the value.
   */
  void
  verify_close(const double       value,
               const double       reference,
               const std::string &description)
  {
    verify(std::abs(value - reference) <=
             reference_tolerance * std::abs(reference),
           description + ": " + std::to_string(value) + " instead of " +
             std::to_string(reference));
  }

  /**
   * @brief Return the QMOM parameters of the lipid nanoparticle production,
   * with the constants of the 0D implementation and all the mechanisms
   * enabled.
   *
   * @return QMOM parameters.
   */
  Parameters::QMOM
  get_qmom_parameters()
  {
    Parameters::QMOM qmom_parameters;
    qmom_parameters.process =
      Parameters::QMOMProcess::lipid_nanoparticle_production;
    qmom_parameters.number_of_scalars  = 2;
    qmom_parameters.number_of_moments  = n_moments;
    qmom_parameters.enable_growth      = true;
    qmom_parameters.enable_nucleation  = true;
    qmom_parameters.enable_coalescence = true;
    qmom_parameters.lipid_nanoparticle_production
      .mole_fraction_solubility_organic = 2.6e-4 + (295. - 273.) * 1.2e-4;
    qmom_parameters.lipid_nanoparticle_production.ion_concentration_mixture =
      25. * 0.75;
    qmom_parameters.physical_constants.elementary_charge = 1.6021766e-19;
    return qmom_parameters;
  }

  /**
   * @brief Return the local scalars of the reference state: ethanol
   * concentration of rho_o / 4 and free lipid concentration of 11.73 / 4
   * kg/m^3.
   *
   * @param[in] mixture Properties of the mixture.
   *
   * @return Local scalars.
   */
  std::vector<double>
  get_reference_scalars(const QMOMLipidNanoparticleMixture &mixture)
  {
    std::vector<double> local_scalars(
      magic_enum::enum_count<Parameters::QMOMLipidNanoparticleScalar>());
    local_scalars[Parameters::get_scalar_index(
      Parameters::QMOMLipidNanoparticleScalar::ethanol_concentration)] =
      0.25 * mixture.parameters.density_organic;
    local_scalars[Parameters::get_scalar_index(
      Parameters::QMOMLipidNanoparticleScalar::free_lipids_concentration)] =
      11.73 / 4.;
    return local_scalars;
  }

  /**
   * @brief Test the composition, solubility and supersaturation of the
   * mixture.
   */
  void
  test_local_state()
  {
    const Parameters::QMOM             qmom_parameters = get_qmom_parameters();
    const QMOMLipidNanoparticleMixture mixture(
      qmom_parameters.lipid_nanoparticle_production,
      qmom_parameters.physical_constants);

    const QMOMLipidNanoparticleLocalState state =
      mixture.compute_local_state(get_reference_scalars(mixture));
    verify_close(state.aqueous_volume_fraction, 0.75, "aqueous fraction");
    verify_close(state.solubility, 0.0069095295633899165, "solubility");
    verify_close(state.supersaturation, 424.41384367726369, "supersaturation");
    verify_close(mixture.capillary_length,
                 5.9095764727006115e-09,
                 "capillary length");
    verify_close(mixture.compute_critical_size(state.supersaturation),
                 9.7667503885595939e-10,
                 "critical size");
    verify_close(mixture.compute_inverse_debye_length(
                   state.aqueous_volume_fraction),
                 498354547.54421246,
                 "inverse Debye length");

    deallog << "Local state: passed" << std::endl;
  }

  /**
   * @brief Test the growth rates, including the dissolution of a particle
   * smaller than its equilibrium size.
   */
  void
  test_growth()
  {
    const Parameters::QMOM             qmom_parameters = get_qmom_parameters();
    const QMOMLipidNanoparticleMixture mixture(
      qmom_parameters.lipid_nanoparticle_production,
      qmom_parameters.physical_constants);
    const QMOMLipidNanoparticleGrowth growth(mixture);

    const std::vector<double> abscissas{5e-9, 5e-8, 5e-10, 0.};
    std::vector<double>       growth_rates(abscissas.size());
    growth.compute_growth_rates(0.,
                                get_reference_scalars(mixture),
                                abscissas,
                                growth_rates);
    verify_close(growth_rates[0], 1.7285225935895732e-07, "growth at 5 nm");
    verify_close(growth_rates[1], 1.7372858203483693e-06, "growth at 50 nm");
    // S*(0.5 nm) = exp(11.8) is larger than S = 424: Ostwald ripening
    verify(growth_rates[2] < 0., "dissolution below the equilibrium size");
    verify(growth_rates[3] == 0., "no growth of a node of zero size");

    deallog << "Growth: passed" << std::endl;
  }

  /**
   * @brief Test the nucleation rate and the moments of the nuclei, and the
   * absence of nucleation in an undersaturated mixture.
   */
  void
  test_nucleation()
  {
    const Parameters::QMOM             qmom_parameters = get_qmom_parameters();
    const QMOMLipidNanoparticleMixture mixture(
      qmom_parameters.lipid_nanoparticle_production,
      qmom_parameters.physical_constants);
    const QMOMLipidNanoparticleNucleation nucleation(mixture);

    constexpr double reference_nucleation_rate = 1.7214216687897362e+21;
    constexpr double reference_critical_size   = 9.7667503885595939e-10;

    std::vector<double> local_scalars = get_reference_scalars(mixture);
    std::vector<double> moment_rates(n_moments);
    nucleation.compute_moment_rates(0., local_scalars, moment_rates);

    // dM_0/dt = J and dM_1/dt = J L_c, since the mean size of the nuclei is
    // the critical size
    verify_close(moment_rates[0], reference_nucleation_rate, "nucleation rate");
    verify_close(moment_rates[1],
                 reference_nucleation_rate * reference_critical_size,
                 "nucleation rate of the first moment");

    // Below the solubility, there is no nucleation
    local_scalars[Parameters::get_scalar_index(
      Parameters::QMOMLipidNanoparticleScalar::free_lipids_concentration)] =
      1e-3;
    nucleation.compute_moment_rates(0., local_scalars, moment_rates);
    for (const double moment_rate : moment_rates)
      verify(moment_rate == 0., "no nucleation in an undersaturated mixture");

    deallog << "Nucleation: passed" << std::endl;
  }

  /**
   * @brief Test the coalescence kernel of three pairs of particles.
   */
  void
  test_coalescence()
  {
    const Parameters::QMOM             qmom_parameters = get_qmom_parameters();
    const QMOMLipidNanoparticleMixture mixture(
      qmom_parameters.lipid_nanoparticle_production,
      qmom_parameters.physical_constants);
    const QMOMLipidNanoparticleCoalescence coalescence(mixture);

    const std::vector<double> local_scalars = get_reference_scalars(mixture);

    constexpr std::array<std::array<double, 3>, 3> reference_kernels{
      {{{5e-9, 5e-8, 1.0414431150314881e-22}},
       {{2e-8, 2e-8, 3.261528950613222e-23}},
       {{2e-9, 3e-9, 3.5613426715303653e-23}}}};
    for (const auto &[first_size, second_size, reference_kernel] :
         reference_kernels)
      {
        verify_close(coalescence.compute_kernel(
                       0., local_scalars, first_size, second_size),
                     reference_kernel,
                     "coalescence kernel");
        // The kernel is symmetric
        verify_close(coalescence.compute_kernel(
                       0., local_scalars, second_size, first_size),
                     reference_kernel,
                     "symmetric coalescence kernel");
      }

    deallog << "Coalescence: passed" << std::endl;
  }

  /**
   * @brief Test the consumption of the free lipids by the particles: the
   * lipid mass is conserved, a negative concentration within round-off is set
   * to zero and a larger one is reported.
   */
  void
  test_lipid_consumption()
  {
    const Parameters::QMOM qmom_parameters = get_qmom_parameters();
    const QMOMSourceTerms  source_terms(qmom_parameters);
    verify(source_terms.scalars_are_coupled(), "the scalars are coupled");

    const Parameters::QMOMLipidNanoparticleProduction &parameters =
      qmom_parameters.lipid_nanoparticle_production;
    const double lipid_mass_per_volume_moment =
      parameters.lipid_density * parameters.shape_factor;
    constexpr unsigned int free_lipids_index = Parameters::get_scalar_index(
      Parameters::QMOMLipidNanoparticleScalar::free_lipids_concentration);
    constexpr unsigned int ethanol_index = Parameters::get_scalar_index(
      Parameters::QMOMLipidNanoparticleScalar::ethanol_concentration);

    // The particles gain a volume of 1e-4 per m^3: the free lipids lose the
    // corresponding mass, while the ethanol is unchanged
    const std::vector<double> moments_before(n_moments, 0.);
    std::vector<double>       moments_after(n_moments, 0.);
    moments_after[3] = 1e-4;
    std::vector<double> scalars{200., 1.};
    verify(source_terms.update_local_scalars(moments_before,
                                             moments_after,
                                             scalars),
           "consumption of part of the free lipids");
    verify_close(scalars[free_lipids_index],
                 1. - lipid_mass_per_volume_moment * 1e-4,
                 "free lipids after consumption");
    verify(scalars[ethanol_index] == 200., "ethanol is not consumed");

    // Consuming exactly the free lipids leaves at most a round-off negative
    // concentration, which is set to zero
    moments_after[3] = 1. / lipid_mass_per_volume_moment;
    scalars          = {200., 1.};
    verify(source_terms.update_local_scalars(moments_before,
                                             moments_after,
                                             scalars),
           "consumption of all the free lipids");
    verify(scalars[free_lipids_index] >= 0. &&
             scalars[free_lipids_index] <= 1e-12,
           "free lipids after consumption of all of them");

    // Consuming twice the free lipids is a failure
    moments_after[3] = 2. / lipid_mass_per_volume_moment;
    scalars          = {200., 1.};
    verify(!source_terms.update_local_scalars(moments_before,
                                              moments_after,
                                              scalars),
           "consumption of more free lipids than available");

    // The generic process does not couple the scalars
    const QMOMSourceTerms generic_source_terms{Parameters::QMOM()};
    verify(!generic_source_terms.scalars_are_coupled(),
           "the scalars of the generic process are not coupled");

    deallog << "Lipid consumption: passed" << std::endl;
  }

} // namespace

int
main()
{
  try
    {
      initlog();
      test_local_state();
      test_growth();
      test_nucleation();
      test_coalescence();
      test_lipid_consumption();
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
