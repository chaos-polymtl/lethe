// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file qmom_source_terms.h
 * @brief Mechanisms that evolve the particle size distribution of the QMOM
 * physics in phase space: growth, nucleation and coalescence.
 *
 * Each mechanism is a family of models: an abstract base class, one derived
 * class per model of the parameter file and a factory, model_cast(), that
 * builds the model selected in the parameters. The models only use standard
 * vectors and do not depend on QMOMKit. QMOMSourceTerms gathers the models of
 * the three mechanisms and wraps them into the callbacks of QMOMKit.
 *
 * The phase space has a single dimension, the internal coordinate xi of the
 * particles. A quadrature of the distribution is a set of nodes, each with a
 * weight (number of particles per m^3) and an abscissa (value of xi). The
 * growth and nucleation models apply to any internal coordinate; the
 * coalescence model, which conserves the volume of the particles, requires xi
 * to be the particle size, in m. The units below are given for that case.
 *
 * Local scalars: the phase-space evolution is done independently at every
 * degree of freedom, and every model receives the values of the QMOM scalars
 * at the degree of freedom being evolved, the "local scalars". They are
 * ordered as the scalars of the QMOM subsection. They are constant during the
 * phase-space step: they are the values at the end of the transport in
 * physical space. None of the models of this file uses them: they are part of
 * the interface so that the models of a process, whose rates depend on the
 * local composition (e.g. the supersaturation of the lipid nanoparticle
 * production), can be added without changing it.
 */

#ifndef lethe_qmom_source_terms_h
#define lethe_qmom_source_terms_h

#include <core/parameters_multiphysics.h>

#include <deal.II/base/function_parser.h>

#ifdef LETHE_WITH_QMOMKIT
#  include <qmomkit/moment_source.hpp>
#  include <qmomkit/node_advection.hpp>
#endif

#include <memory>
#include <string>
#include <vector>

using namespace dealii;

/**
 * @brief Base class of the growth models.
 *
 * The growth rate G is the velocity of the particle size xi in phase space, in
 * m/s: the abscissa of every node of the quadrature moves with dxi/dt = G,
 * while its weight is unchanged.
 */
class QMOMGrowthModelBase
{
public:
  virtual ~QMOMGrowthModelBase() = default;

  /**
   * @brief Build the growth model selected in the "growth" subsection.
   *
   * @param[in] growth_parameters Parameters of the "growth" subsection.
   *
   * @return Growth model.
   */
  static std::unique_ptr<QMOMGrowthModelBase>
  model_cast(const Parameters::QMOMGrowth &growth_parameters);

  /**
   * @brief Compute the growth rate of every node of the quadrature.
   *
   * @param[in] time Time, in s.
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (see the file description).
   * @param[in] abscissas Sizes of the nodes, in m.
   * @param[out] growth_rates Growth rates of the nodes, in m/s, of the size of
   * abscissas.
   */
  virtual void
  compute_growth_rates(const double               time,
                       const std::vector<double> &local_scalars,
                       const std::vector<double> &abscissas,
                       std::vector<double>       &growth_rates) const = 0;
};

/**
 * @brief Constant growth: G = a. All the sizes increase by the same amount.
 */
class QMOMConstantGrowth : public QMOMGrowthModelBase
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] growth_coefficient Growth rate a, in m/s.
   */
  explicit QMOMConstantGrowth(const double growth_coefficient);

  /**
   * @brief Compute the growth rate of every node: G = a.
   *
   * @param[in] time Time, in s (unused).
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (unused).
   * @param[in] abscissas Sizes of the nodes, in m.
   * @param[out] growth_rates Growth rates of the nodes, in m/s.
   */
  void
  compute_growth_rates(const double               time,
                       const std::vector<double> &local_scalars,
                       const std::vector<double> &abscissas,
                       std::vector<double>       &growth_rates) const override;

private:
  /// Growth rate a [m/s]
  const double growth_coefficient;
};

/**
 * @brief Linear growth: G = a xi. All the sizes are multiplied by exp(a t).
 */
class QMOMLinearGrowth : public QMOMGrowthModelBase
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] growth_coefficient Coefficient a, in 1/s.
   */
  explicit QMOMLinearGrowth(const double growth_coefficient);

  /**
   * @brief Compute the growth rate of every node: G = a xi.
   *
   * @param[in] time Time, in s (unused).
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (unused).
   * @param[in] abscissas Sizes of the nodes, in m.
   * @param[out] growth_rates Growth rates of the nodes, in m/s.
   */
  void
  compute_growth_rates(const double               time,
                       const std::vector<double> &local_scalars,
                       const std::vector<double> &abscissas,
                       std::vector<double>       &growth_rates) const override;

private:
  /// Coefficient a [1/s]
  const double growth_coefficient;
};

/**
 * @brief Growth rate given by an expression of the particle size xi and of the
 * time t.
 *
 * The expression is evaluated by a parser for every node, which is slower than
 * the other models.
 */
class QMOMFunctionGrowth : public QMOMGrowthModelBase
{
public:
  /**
   * @brief Constructor. Parses the expression.
   *
   * @param[in] function_expression Expression of the growth rate, in m/s, in
   * terms of xi, in m, and t, in s.
   */
  explicit QMOMFunctionGrowth(const std::string &function_expression);

  /**
   * @brief Compute the growth rate of every node by evaluating the expression
   * at its size.
   *
   * @param[in] time Time, in s.
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (unused).
   * @param[in] abscissas Sizes of the nodes, in m.
   * @param[out] growth_rates Growth rates of the nodes, in m/s.
   */
  void
  compute_growth_rates(const double               time,
                       const std::vector<double> &local_scalars,
                       const std::vector<double> &abscissas,
                       std::vector<double>       &growth_rates) const override;

private:
  /// Growth rate as a function of xi (the "spatial" variable) and t. It is
  /// mutable, since its time is set at every evaluation.
  mutable FunctionParser<1> growth_function;
};

/**
 * @brief Base class of the nucleation models.
 *
 * The nucleation creates particles at the rate J, in 1/(m^3 s), whose sizes
 * follow a log-normal distribution: the rate of change of the moment M_k is J
 * times the moment k of the distribution of the nuclei. It does not depend on
 * the particles that already exist, so it also creates particles where there
 * are none.
 */
class QMOMNucleationModelBase
{
public:
  virtual ~QMOMNucleationModelBase() = default;

  /**
   * @brief Build the nucleation model selected in the "nucleation"
   * subsection.
   *
   * @param[in] nucleation_parameters Parameters of the "nucleation"
   * subsection.
   *
   * @return Nucleation model.
   */
  static std::unique_ptr<QMOMNucleationModelBase>
  model_cast(const Parameters::QMOMNucleation &nucleation_parameters);

  /**
   * @brief Compute the rate of change of the moments due to the nucleation.
   *
   * @param[in] time Time, in s.
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (see the file description).
   * @param[out] moment_rates Rates of change dM_k/dt of the moments, of the
   * size of the number of moments.
   */
  virtual void
  compute_moment_rates(const double               time,
                       const std::vector<double> &local_scalars,
                       std::vector<double>       &moment_rates) const = 0;
};

/**
 * @brief Log-normal distribution of the sizes of the nuclei, given by its mean
 * size and standard deviation.
 *
 * Its moments are E[xi^k] = exp(k mu + k^2 s^2 / 2), with
 * s^2 = ln(1 + sigma^2 / m^2) and mu = ln(m) - s^2 / 2, where m is the mean
 * size and sigma the standard deviation.
 */
class QMOMLogNormalNuclei
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] mean_size Mean size m of the nuclei, in m.
   * @param[in] standard_deviation Standard deviation sigma of the size of the
   * nuclei, in m.
   */
  QMOMLogNormalNuclei(const double mean_size, const double standard_deviation);

  /**
   * @brief Fill the rates of change of the moments due to a nucleation at the
   * rate J: dM_k/dt = J E[xi^k].
   *
   * @param[in] nucleation_rate Nucleation rate J, in 1/(m^3 s).
   * @param[out] moment_rates Rates of change dM_k/dt of the moments.
   */
  void
  fill_moment_rates(const double         nucleation_rate,
                    std::vector<double> &moment_rates) const;

private:
  /// Mean mu of the logarithm of the size
  const double log_mean;

  /// Variance s^2 of the logarithm of the size
  const double log_variance;
};

/**
 * @brief Constant nucleation rate J = a: dM_k/dt = a E[xi^k], with E[xi^k]
 * the moments of the log-normal distribution of the nuclei.
 */
class QMOMConstantNucleation : public QMOMNucleationModelBase
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] nucleation_coefficient Nucleation rate a, in 1/(m^3 s).
   * @param[in] nuclei_mean_size Mean size of the nuclei, in m.
   * @param[in] nuclei_standard_deviation Standard deviation of the size of the
   * nuclei, in m.
   */
  QMOMConstantNucleation(const double nucleation_coefficient,
                         const double nuclei_mean_size,
                         const double nuclei_standard_deviation);

  /**
   * @brief Compute the rate of change of the moments: dM_k/dt = a E[xi^k].
   *
   * @param[in] time Time, in s (unused).
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (unused).
   * @param[out] moment_rates Rates of change dM_k/dt of the moments.
   */
  void
  compute_moment_rates(const double               time,
                       const std::vector<double> &local_scalars,
                       std::vector<double>       &moment_rates) const override;

private:
  /// Nucleation rate a [1/(m^3 s)]
  const double nucleation_coefficient;

  /// Size distribution of the nuclei
  const QMOMLogNormalNuclei nuclei;
};

/**
 * @brief Nucleation rate J given by an expression of the time t.
 *
 * The moments change at the rate dM_k/dt = J(t) E[xi^k], with E[xi^k] the
 * moments of the log-normal distribution of the nuclei.
 */
class QMOMFunctionNucleation : public QMOMNucleationModelBase
{
public:
  /**
   * @brief Constructor. Parses the expression.
   *
   * @param[in] function_expression Expression of the nucleation rate, in
   * 1/(m^3 s), in terms of t, in s.
   * @param[in] nuclei_mean_size Mean size of the nuclei, in m.
   * @param[in] nuclei_standard_deviation Standard deviation of the size of the
   * nuclei, in m.
   */
  QMOMFunctionNucleation(const std::string &function_expression,
                         const double       nuclei_mean_size,
                         const double       nuclei_standard_deviation);

  /**
   * @brief Compute the rate of change of the moments: dM_k/dt =
   * J(t) E[xi^k].
   *
   * @param[in] time Time, in s.
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (unused).
   * @param[out] moment_rates Rates of change dM_k/dt of the moments.
   */
  void
  compute_moment_rates(const double               time,
                       const std::vector<double> &local_scalars,
                       std::vector<double>       &moment_rates) const override;

private:
  /// Nucleation rate as a function of t. The time is its only variable, so it
  /// takes the place of the spatial coordinate of the parser.
  FunctionParser<1> nucleation_function;

  /// Size distribution of the nuclei
  const QMOMLogNormalNuclei nuclei;
};

/**
 * @brief Base class of the coalescence models.
 *
 * Two particles of sizes xi_1 and xi_2 merge at the rate given by the kernel
 * beta(xi_1, xi_2), in m^3/s, into a particle of the same total volume, of size
 * (xi_1^3 + xi_2^3)^(1/3). With a quadrature of weights w_i and abscissas xi_i,
 * the rate of change of the moment M_k is:
 *
 * dM_k/dt = 1/2 sum_i sum_j w_i w_j beta(xi_i, xi_j)
 *           [(xi_i^3 + xi_j^3)^(k/3) - xi_i^k - xi_j^k]
 *
 * so that M_3, proportional to the volume of the particles, is conserved. The
 * derived classes only provide the kernel.
 */
class QMOMCoalescenceModelBase
{
public:
  virtual ~QMOMCoalescenceModelBase() = default;

  /**
   * @brief Build the coalescence model selected in the "coalescence"
   * subsection.
   *
   * @param[in] coalescence_parameters Parameters of the "coalescence"
   * subsection.
   *
   * @return Coalescence model.
   */
  static std::unique_ptr<QMOMCoalescenceModelBase>
  model_cast(const Parameters::QMOMCoalescence &coalescence_parameters);

  /**
   * @brief Compute the coalescence kernel of two particles.
   *
   * @param[in] time Time, in s.
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (see the file description).
   * @param[in] first_size Size xi_1 of the first particle, in m.
   * @param[in] second_size Size xi_2 of the second particle, in m.
   *
   * @return Kernel beta(xi_1, xi_2), in m^3/s.
   */
  virtual double
  compute_kernel(const double               time,
                 const std::vector<double> &local_scalars,
                 const double               first_size,
                 const double               second_size) const = 0;

  /**
   * @brief Compute the rate of change of the moments due to the coalescence
   * of the nodes of a quadrature (see the class description).
   *
   * The pairs that include a node of non-positive size are skipped, since these
   * nodes hold no particle that can merge.
   *
   * @param[in] time Time, in s.
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (see the file description).
   * @param[in] weights Weights of the nodes, in 1/m^3.
   * @param[in] abscissas Sizes of the nodes, in m.
   * @param[out] moment_rates Rates of change dM_k/dt of the moments, of the
   * size of the number of moments.
   */
  void
  compute_moment_rates(const double               time,
                       const std::vector<double> &local_scalars,
                       const std::vector<double> &weights,
                       const std::vector<double> &abscissas,
                       std::vector<double>       &moment_rates) const;
};

/**
 * @brief Constant coalescence kernel: beta = a. Then dM_0/dt = -a M_0^2 / 2.
 */
class QMOMConstantCoalescence : public QMOMCoalescenceModelBase
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] coalescence_coefficient Kernel a, in m^3/s.
   */
  explicit QMOMConstantCoalescence(const double coalescence_coefficient);

  /**
   * @brief Compute the coalescence kernel: beta = a.
   *
   * @param[in] time Time, in s (unused).
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (unused).
   * @param[in] first_size Size of the first particle, in m (unused).
   * @param[in] second_size Size of the second particle, in m (unused).
   *
   * @return Kernel a, in m^3/s.
   */
  double
  compute_kernel(const double               time,
                 const std::vector<double> &local_scalars,
                 const double               first_size,
                 const double               second_size) const override;

private:
  /// Kernel a [m^3/s]
  const double coalescence_coefficient;
};

/**
 * @brief Coalescence kernel given by an expression of the sizes xi_1 and xi_2
 * of the two particles and of the time t.
 *
 * The expression is evaluated by a parser for every pair of nodes, which is
 * slower than the constant model.
 */
class QMOMFunctionCoalescence : public QMOMCoalescenceModelBase
{
public:
  /**
   * @brief Constructor. Parses the expression.
   *
   * @param[in] function_expression Expression of the kernel, in m^3/s, in
   * terms of xi_1 and xi_2, in m, and t, in s.
   */
  explicit QMOMFunctionCoalescence(const std::string &function_expression);

  /**
   * @brief Compute the coalescence kernel by evaluating the expression.
   *
   * @param[in] time Time, in s.
   * @param[in] local_scalars Values of the scalars at the degree of freedom
   * (unused).
   * @param[in] first_size Size xi_1 of the first particle, in m.
   * @param[in] second_size Size xi_2 of the second particle, in m.
   *
   * @return Kernel beta(xi_1, xi_2), in m^3/s.
   */
  double
  compute_kernel(const double               time,
                 const std::vector<double> &local_scalars,
                 const double               first_size,
                 const double               second_size) const override;

private:
  /// Kernel as a function of xi_1 and xi_2 (the "spatial" variables) and t.
  /// It is mutable, since its time is set at every evaluation.
  mutable FunctionParser<2> coalescence_function;
};

/**
 * @brief Phase-space mechanisms of the QMOM physics.
 *
 * They are built from the parameters of the "QMOM" subsection: the models of
 * the enabled mechanisms, wrapped into the callbacks of QMOMKit.
 *
 * The callbacks are built once, at construction, and read the local scalars
 * stored in this object: set_local_scalars() must be called at every degree of
 * freedom, before its moments are evolved. Since the callbacks refer to this
 * object, it can be neither copied nor moved.
 */
class QMOMSourceTerms
{
public:
  /**
   * @brief Constructor. Builds the model of every enabled mechanism and the
   * callbacks of QMOMKit.
   *
   * @param[in] qmom_parameters Parameters of the "QMOM" subsection.
   */
  explicit QMOMSourceTerms(const Parameters::QMOM &qmom_parameters);

  QMOMSourceTerms(const QMOMSourceTerms &) = delete;
  QMOMSourceTerms &
  operator=(const QMOMSourceTerms &) = delete;

  /**
   * @brief Indicate whether the nucleation is enabled.
   *
   * @return True if the nucleation is enabled.
   */
  bool
  nucleation_is_enabled() const
  {
    return nucleation_model != nullptr;
  }

  /**
   * @brief Set the values of the scalars at the degree of freedom whose
   * moments are evolved next (see the file description).
   *
   * @param[in] new_local_scalars Values of the scalars, ordered as the scalars
   * of the QMOM subsection.
   */
  void
  set_local_scalars(const std::vector<double> &new_local_scalars);

  /**
   * @brief Compute the rate of change of the moments due to the nucleation
   * alone, which is zero if the nucleation is disabled.
   *
   * It is used where there are no particles, so that there is no quadrature to
   * grow or coalesce.
   *
   * @param[in] time Time, in s.
   * @param[out] moment_rates Rates of change dM_k/dt of the moments, of the
   * size of the number of moments.
   */
  void
  compute_nucleation_moment_rates(const double         time,
                                  std::vector<double> &moment_rates) const;

#ifdef LETHE_WITH_QMOMKIT
  /**
   * @brief Get the phase-space velocity of the quadrature nodes: the growth
   * rate, or zero if the growth is disabled.
   *
   * @return Reference to the phase-space velocity.
   */
  const qmomkit::PhaseSpaceVelocity &
  get_phase_space_velocity() const
  {
    return phase_space_velocity;
  }

  /**
   * @brief Get the moment source: the coalescence, or zero if it is disabled.
   *
   * @return Reference to the moment source.
   */
  const qmomkit::QmomMomentSource &
  get_moment_source() const
  {
    return moment_source;
  }

  /**
   * @brief Get the phase-space boundary contribution: the nucleation, or zero
   * if it is disabled.
   *
   * QMOMKit calls a source of particles that does not depend on the existing
   * particles a boundary contribution.
   *
   * @return Reference to the boundary contribution.
   */
  const qmomkit::QmomMomentBoundaryContribution &
  get_boundary_contribution() const
  {
    return boundary_contribution;
  }
#endif

private:
  /// Growth model, or nullptr if the growth is disabled
  std::unique_ptr<QMOMGrowthModelBase> growth_model;

  /// Nucleation model, or nullptr if the nucleation is disabled
  std::unique_ptr<QMOMNucleationModelBase> nucleation_model;

  /// Coalescence model, or nullptr if the coalescence is disabled
  std::unique_ptr<QMOMCoalescenceModelBase> coalescence_model;

  /// Values of the scalars at the degree of freedom being evolved
  std::vector<double> local_scalars;

  /// Weights and abscissas of the quadrature nodes, gathered for the
  /// coalescence. They are members to avoid an allocation at every call.
  mutable std::vector<double> node_weights;
  mutable std::vector<double> node_abscissas;

#ifdef LETHE_WITH_QMOMKIT
  /// Phase-space velocity of the nodes, built from the growth model
  qmomkit::PhaseSpaceVelocity phase_space_velocity;

  /// Moment source, built from the coalescence model
  qmomkit::QmomMomentSource moment_source;

  /// Boundary contribution, built from the nucleation model
  qmomkit::QmomMomentBoundaryContribution boundary_contribution;
#endif
};

#endif
