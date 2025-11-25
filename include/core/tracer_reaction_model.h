// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tracer_reaction_prefactor_model_h
#define lethe_tracer_reaction_prefactor_model_h

#include <core/physical_property_model.h>

#include <algorithm>
#include <cmath>

inline double
Ceffective(const double C, const double epsilon)
{
  return std::sqrt(C * C + epsilon * epsilon);
}

inline double
dCeffective_dC(const double C, const double epsilon)
{
  const double Ce = Ceffective(C, epsilon);
  return (Ce > 0.0) ? (C / Ce) : 0.0;
}

/**
 * @brief Abstract class to calculate the tracer reaction prefactor.
 *
 * This model computes the prefactor
 * \f$ k = \alpha \, C^{(n-1)} \f$
 * such that the full reaction term is assembled as
 * \f$ r = -k \, C \f$.
 *
 * Here, \f$ \alpha \f$ represents the reaction constant.
 * This parent class could have children classes that are not based on the power
 * law.
 */
class TracerReactionPrefactorModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a TracerReactionPrefactorModel object by casting it to
   * the proper child class.
   *
   * @param material_properties Parameters for a material.
   */
  static std::shared_ptr<TracerReactionPrefactorModel>
  model_cast(const Parameters::Material &material_properties);
};

/**
 * @brief Dummy tracer reaction prefactor model.
 *
 * This class is a dummy implementation and always returns zero.
 */
class NoneTracerReactionPrefactor : public TracerReactionPrefactorModel
{
public:
  /**
   * @brief Default constructor.
   */
  NoneTracerReactionPrefactor()
  {}

  /**
   * @brief Returns a dummy value.
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return 0.;
  }

  /**
   * @brief Dummy vector value computation.
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> & /*property_vector*/) override
  {}

  /**
   * @brief Returns a dummy jacobian value.
   */
  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  }

  /**
   * @brief Dummy vector jacobian computation.
   */
  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> & /*jacobian_vector*/) override
  {}
};

/**
 * @brief Constant tracer reaction prefactor.
 *
 *This model computes the prefactor as
 * \f[
 * k = alpha * C_eff^{n-1}, with a smooth floor C_eff = sqrt(C^2 + eps^2)
 * \f]
 * where \f$ \alpha \f$ (the reaction constant) and the reaction order \f$ n \f$
 * are provided.
 */
class ConstantTracerReactionPrefactor : public TracerReactionPrefactorModel
{
public:
  /**
   * @brief Constructor.
   *
   * @param p_tracer_reaction_constant The base reaction constant \f$ \alpha \f$.
   * @param p_tracer_reaction_order The reaction order \f$ n \f$.
   */
  ConstantTracerReactionPrefactor(const double p_tracer_reaction_constant,
                                  const double p_tracer_reaction_order,
                                  const double p_tracer_reaction_epsilon)
    : tracer_reaction_constant(p_tracer_reaction_constant)
    , tracer_reaction_order(p_tracer_reaction_order)
    , tracer_reaction_epsilon(p_tracer_reaction_epsilon)
  {
    this->model_depends_on[field::tracer_concentration] = true;
  }

  /**
   * @brief Compute the tracer reaction prefactor.
   *
   * @param fields_value Map of field values.
   * @return The computed prefactor value.
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    Assert(fields_value.find(field::tracer_concentration) != fields_value.end(),
           PhysicialPropertyModelFieldUndefined(
             "ConstantTracerReactionPrefactor", "tracer_concentration"));

    const double Ceff = Ceffective(fields_value.at(field::tracer_concentration),
                                   tracer_reaction_epsilon);
    return tracer_reaction_constant *
           std::pow(Ceff, tracer_reaction_order - 1.);
  }

  /**
   * @brief Compute the vector of tracer reaction prefactor values.
   *
   * @param field_vectors Map of field vectors.
   * @param property_vector Vector to be filled with computed prefactor values.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    Assert(
      field_vectors.find(field::tracer_concentration) != field_vectors.end(),
      PhysicialPropertyModelFieldUndefined("ConstantTracerReactionPrefactor",
                                           "tracer_concentration"));
    const auto &concentration_vector =
      field_vectors.at(field::tracer_concentration);
    for (size_t i = 0; i < property_vector.size(); ++i)
      {
        const double Ceff =
          Ceffective(concentration_vector[i], tracer_reaction_epsilon);
        property_vector[i] =
          tracer_reaction_constant * std::pow(Ceff, tracer_reaction_order - 1.);
      }
  }

  /**
   * @brief Compute the jacobian (partial derivative) of the prefactor with respect to a field.
   *
   * @param field_values Map of field values.
   * @param id The field identifier with respect to which the derivative is computed.
   * @return The computed derivative value.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    Assert(field_values.find(field::tracer_concentration) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "ConstantTracerReactionPrefactor", "tracer_concentration"));
    if (id == field::tracer_concentration)
      {
        const double C    = field_values.at(field::tracer_concentration);
        const double Ceff = Ceffective(C, tracer_reaction_epsilon);
        const double dCe  = dCeffective_dC(C, tracer_reaction_epsilon);

        return tracer_reaction_constant * (tracer_reaction_order - 1.) *
               std::pow(Ceff, tracer_reaction_order - 2.) * dCe;
      }
    return 0.0;
  }

  /**
   * @brief Compute the vector jacobian (partial derivatives) of the prefactor with respect to a field.
   *
   * @param field_vectors Map of field vectors.
   * @param id The field identifier with respect to which the derivative is computed.
   * @param jacobian_vector Vector to be filled with computed derivative values.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    Assert(
      field_vectors.find(field::tracer_concentration) != field_vectors.end(),
      PhysicialPropertyModelFieldUndefined("ConstantTracerReactionPrefactor",
                                           "tracer_concentration"));
    if (id != field::tracer_concentration)
      {
        std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0.0);
        return;
      }
    const auto &concentration_vector =
      field_vectors.at(field::tracer_concentration);
    for (size_t i = 0; i < jacobian_vector.size(); ++i)
      {
        const double Ceff =
          Ceffective(concentration_vector[i], tracer_reaction_epsilon);
        const double dCe =
          dCeffective_dC(concentration_vector[i], tracer_reaction_epsilon);
        jacobian_vector[i] = tracer_reaction_constant *
                             (tracer_reaction_order - 1.) *
                             std::pow(Ceff, tracer_reaction_order - 2.) * dCe;
      }
  };

private:
  const double tracer_reaction_constant;
  const double tracer_reaction_order;
  const double tracer_reaction_epsilon;
};

/**
 * @brief Tanh level-set-dependent tracer reaction prefactor.
 *
 * This model computes a prefactor that depends on the level-set function,
 * smoothly transitioning between a reaction constant \(\alpha\) inside and
 * outside the solid.
 */
class TanhLevelsetTracerReactionPrefactor : public TracerReactionPrefactorModel
{
public:
  /**
   * @brief Constructor.
   *
   * @param p_tracer_reaction_constant_outside Reaction constant \(\alpha\) outside the solid.
   * @param p_tracer_reaction_constant_inside Reaction constant \(\alpha\) inside the solid.
   * @param p_thickness Thickness for the tanh smoothing function.
   * @param p_tracer_reaction_order The reaction order.
   */
  TanhLevelsetTracerReactionPrefactor(
    const double p_tracer_reaction_constant_outside,
    const double p_tracer_reaction_constant_inside,
    const double p_thickness,
    const double p_tracer_reaction_order,
    const double p_tracer_reaction_epsilon)
    : tracer_reaction_constant_outside(p_tracer_reaction_constant_outside)
    , tracer_reaction_constant_inside(p_tracer_reaction_constant_inside)
    , thickness(p_thickness)
    , delta_reaction_constant(tracer_reaction_constant_outside -
                              tracer_reaction_constant_inside)
    , tracer_reaction_order(p_tracer_reaction_order)
    , tracer_reaction_epsilon(p_tracer_reaction_epsilon)
  {
    this->model_depends_on[field::levelset]             = true;
    this->model_depends_on[field::tracer_concentration] = true;
  }

  /**
   * @brief Compute the tracer reaction prefactor.
   *
   * @param field_values Map of field values.
   * @return The computed prefactor value.
   */
  double
  value(const std::map<field, double> &field_values) override
  {
    Assert(field_values.find(field::levelset) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_values.find(field::tracer_concentration) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "tracer_concentration"));
    const double levelset = field_values.at(field::levelset);
    const double Ceff = Ceffective(field_values.at(field::tracer_concentration),
                                   tracer_reaction_epsilon);
    const double k =
      tracer_reaction_constant_inside +
      delta_reaction_constant * (0.5 + 0.5 * std::tanh(levelset / thickness));
    return k * std::pow(Ceff, tracer_reaction_order - 1.);
  }

  /**
   * @brief Compute a vector of tracer reaction prefactor values.
   *
   * @param field_vectors Map of field vectors.
   * @param property_vector Vector to be filled with computed prefactor values.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    Assert(field_vectors.find(field::levelset) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_vectors.find(field::tracer_concentration) !=
             field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "tracer_concentration"));

    const std::vector<double> &levelset_vec = field_vectors.at(field::levelset);
    const std::vector<double> &concentration_vector =
      field_vectors.at(field::tracer_concentration);

    const unsigned int n_values = levelset_vec.size();
    Assert(n_values == concentration_vector.size(),
           SizeOfFields(n_values, concentration_vector.size()));
    for (unsigned int i = 0; i < n_values; ++i)
      {
        const double k = tracer_reaction_constant_inside +
                         delta_reaction_constant *
                           (0.5 + 0.5 * std::tanh(levelset_vec[i] / thickness));
        const double Ceff =
          Ceffective(concentration_vector[i], tracer_reaction_epsilon);
        property_vector[i] = k * std::pow(Ceff, tracer_reaction_order - 1.);
      }
  }

  /**
   * @brief Compute the jacobian (partial derivative) of the prefactor with respect to a specified field.
   *
   * @param field_values Map of field values.
   * @param id The field identifier with respect to which the derivative is computed.
   * @return The computed derivative value.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    Assert(field_values.find(field::levelset) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_values.find(field::tracer_concentration) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "tracer_concentration"));
    const double levelset = field_values.at(field::levelset);
    const double concentration_val =
      field_values.at(field::tracer_concentration);
    const double Ceff = Ceffective(concentration_val, tracer_reaction_epsilon);
    if (id == field::levelset)
      {
        // dk/dlambda for tanh profile: Δk * 0.5 * (1 - tanh^2(lambda/sigma)) *
        // (1/sigma)
        const double tanh = std::tanh(levelset / thickness);
        const double dkdlambda =
          delta_reaction_constant * 0.5 * (1.0 - tanh * tanh) / thickness;
        return dkdlambda * std::pow(Ceff, tracer_reaction_order - 1.);
      }
    else if (id == field::tracer_concentration)
      {
        const double k = tracer_reaction_constant_inside +
                         delta_reaction_constant *
                           (0.5 + 0.5 * std::tanh(levelset / thickness));
        const double dCe =
          dCeffective_dC(concentration_val, tracer_reaction_epsilon);
        return k * (tracer_reaction_order - 1.) *
               std::pow(Ceff, tracer_reaction_order - 2.) * dCe;
      }
    else
      return 0.0;
  }

  /**
   * @brief Compute the vector jacobian (partial derivatives) of the prefactor with respect to a field.
   *
   * @param field_vectors Map of field vectors.
   * @param id The field identifier with respect to which the derivative is computed.
   * @param jacobian_vector Vector to be filled with computed derivative values.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    Assert(field_vectors.find(field::levelset) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_vectors.find(field::tracer_concentration) !=
             field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "tracer_concentration"));
    const std::vector<double> &levelset_vec = field_vectors.at(field::levelset);
    const std::vector<double> &concentration_vec =
      field_vectors.at(field::tracer_concentration);
    const unsigned int n_values = levelset_vec.size();
    if (id == field::levelset)
      {
        for (unsigned int i = 0; i < n_values; ++i)
          {
            const double tanh      = std::tanh(levelset_vec[i] / thickness);
            const double dkdlambda = delta_reaction_constant * 0.5 *
                                     (1.0 - std::pow(tanh, 2)) / thickness;
            const double Ceff =
              Ceffective(concentration_vec[i], tracer_reaction_epsilon);
            jacobian_vector[i] =
              dkdlambda * std::pow(Ceff, tracer_reaction_order - 1.);
          }
      }
    else if (id == field::tracer_concentration)
      {
        for (unsigned int i = 0; i < n_values; ++i)
          {
            const double k =
              tracer_reaction_constant_inside +
              delta_reaction_constant *
                (0.5 + 0.5 * std::tanh(levelset_vec[i] / thickness));
            const double Ceff =
              Ceffective(concentration_vec[i], tracer_reaction_epsilon);
            const double dCe =
              dCeffective_dC(concentration_vec[i], tracer_reaction_epsilon);
            jacobian_vector[i] = k * (tracer_reaction_order - 1.) *
                                 std::pow(Ceff, tracer_reaction_order - 2.) *
                                 dCe;
          }
      }
    else
      {
        std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0.0);
      }
  }

private:
  const double tracer_reaction_constant_outside;
  const double tracer_reaction_constant_inside;
  const double thickness;
  const double delta_reaction_constant;
  const double tracer_reaction_order;
  const double tracer_reaction_epsilon;
};

/**
 * @brief Reaction prefactor that depends on the Gaussian of the level set
 */
class GaussianLevelsetTracerReactionPrefactor
  : public TracerReactionPrefactorModel
{
public:
  /**
   * @brief Constructor of the Gaussian level set-dependent reaction prefactor model.
   *
   * @param[in] p_tracer_reaction_constant_interface Reaction constant at the
   * fluid-solid interface the solid
   * @param[in] p_tracer_reaction_constant_bulk Reaction constant in the bulk of
   * the phases
   * @param[in] p_thickness Thickness of the Gaussian function used to smooth
   * the property jump
   * @param[in] p_tracer_reaction_order The reaction order.
   */
  GaussianLevelsetTracerReactionPrefactor(
    const double p_tracer_reaction_constant_interface,
    const double p_tracer_reaction_constant_bulk,
    const double p_thickness,
    const double p_tracer_reaction_order,
    const double p_tracer_reaction_epsilon)
    : tracer_reaction_constant_interface(p_tracer_reaction_constant_interface)
    , tracer_reaction_constant_bulk(p_tracer_reaction_constant_bulk)
    , delta_reaction_constant(tracer_reaction_constant_interface -
                              tracer_reaction_constant_bulk)
    , squared_thickness(std::pow(p_thickness, 2))
    , tracer_reaction_order(p_tracer_reaction_order)
    , tracer_reaction_epsilon(p_tracer_reaction_epsilon)
  {
    this->model_depends_on[field::levelset]             = true;
    this->model_depends_on[field::tracer_concentration] = true;
  }

  /**
   * @brief Compute the reaction constant.
   *
   * @param[in] field_values Values of the various fields on which the property
   * may depend. In this case, the reaction constant depends on the level set.
   * The map stores a single value per field.
   *
   * @return Value of the reaction constant computed with the @p field_values.
   */
  double
  value(const std::map<field, double> &field_values) override
  {
    Assert(field_values.find(field::levelset) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "GaussianLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_values.find(field::tracer_concentration) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "GaussianLevelsetTracerReactionPrefactor",
             "tracer_concentration"));
    const double levelset_val = field_values.at(field::levelset);
    const double Ceff = Ceffective(field_values.at(field::tracer_concentration),
                                   tracer_reaction_epsilon);

    const double k = tracer_reaction_constant_bulk +
                     delta_reaction_constant *
                       std::exp(-std::pow(levelset_val, 2) / squared_thickness);
    return k * std::pow(Ceff, tracer_reaction_order - 1.);
  }

  /**
   * @brief Compute a vector of reaction constant.
   *
   * @param[in] field_vectors Vectors of the fields on which the reaction
   * constant may depend. In this case, the reaction constant depends on the
   * level set. The map stores a vector of values per field.
   *
   * @param[out] property_vector Vectors of computed reaction constants.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    Assert(field_vectors.find(field::levelset) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "GaussianLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_vectors.find(field::tracer_concentration) !=
             field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "GaussianLevelsetTracerReactionPrefactor",
             "tracer_concentration"));

    const std::vector<double> &levelset_vec = field_vectors.at(field::levelset);
    const std::vector<double> &concentration_vec =
      field_vectors.at(field::tracer_concentration);
    const unsigned int n_values = levelset_vec.size();
    Assert(n_values == concentration_vec.size(),
           SizeOfFields(n_values, concentration_vec.size()));
    for (unsigned int i = 0; i < n_values; ++i)
      {
        const double k =
          tracer_reaction_constant_bulk +
          delta_reaction_constant *
            std::exp(-std::pow(levelset_vec[i], 2) / squared_thickness);
        const double Ceff =
          Ceffective(concentration_vec[i], tracer_reaction_epsilon);
        property_vector[i] = k * std::pow(Ceff, tracer_reaction_order - 1.);
      }
  }

  /**
   * @brief Compute the jacobian (the partial derivative) of the reaction constant
   * with respect to a specified field.
   *
   * @param[in] field_values Values of the various fields on which the specific
   * heat may depend. The map stores a single value per field.
   *
   * @param[in] id Indicator of the field with respect to which the jacobian
   * should be computed.
   *
   * @return Value of the derivative of the reaction constant with respect to the
   * specified field.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    Assert(field_values.find(field::levelset) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "GaussianLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_values.find(field::tracer_concentration) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "GaussianLevelsetTracerReactionPrefactor",
             "tracer_concentration"));

    const double levelset_val = field_values.at(field::levelset);
    const double concentration_val =
      field_values.at(field::tracer_concentration);
    const double Ceff = Ceffective(concentration_val, tracer_reaction_epsilon);
    if (id == field::levelset)
      {
        // dk/dlambda = Δk * exp(-lambda^2/sigma^2) * (-2*lambda/sigma^2)
        const double exponential =
          std::exp(-pow(levelset_val, 2) / squared_thickness);
        const double dkdlambda = delta_reaction_constant * exponential *
                                 (-2.0 * levelset_val / squared_thickness);
        return dkdlambda * std::pow(Ceff, tracer_reaction_order - 1.);
      }
    else if (id == field::tracer_concentration)
      {
        const double k =
          tracer_reaction_constant_bulk +
          delta_reaction_constant *
            std::exp(-std::pow(levelset_val, 2) / squared_thickness);
        const double dCe =
          dCeffective_dC(concentration_val, tracer_reaction_epsilon);
        return k * (tracer_reaction_order - 1.) *
               std::pow(Ceff, tracer_reaction_order - 2.) * dCe;
      }
    else
      return 0.0;
  }

  /**
   * @brief Compute the derivative of the reaction constant with respect to a field.
   *
   * @param[in] field_vectors Vector of values of the fields used to evaluate
   * the reaction constant. The map stores a vector of values per field.
   *
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be computed.
   *
   * @param[out] jacobian Vector of computed derivative values of the
   * reaction constant with respect to the field of the specified @p id.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    Assert(field_vectors.find(field::levelset) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "GaussianLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_vectors.find(field::tracer_concentration) !=
             field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "GaussianLevelsetTracerReactionPrefactor",
             "tracer_concentration"));
    const std::vector<double> &levelset_vec = field_vectors.at(field::levelset);
    const std::vector<double> &concentration_vec =
      field_vectors.at(field::tracer_concentration);
    const unsigned int n_values = levelset_vec.size();
    if (id == field::levelset)
      {
        for (unsigned int i = 0; i < n_values; ++i)
          {
            const double exponential =
              std::exp(-pow(levelset_vec[i], 2) / squared_thickness);
            const double dkdlambda =
              delta_reaction_constant * exponential *
              (-2.0 * levelset_vec[i] / squared_thickness);
            const double Ceff =
              Ceffective(concentration_vec[i], tracer_reaction_epsilon);
            jacobian_vector[i] =
              dkdlambda * std::pow(Ceff, tracer_reaction_order - 1.);
          }
      }
    else if (id == field::tracer_concentration)
      {
        for (unsigned int i = 0; i < n_values; ++i)
          {
            const double k =
              tracer_reaction_constant_bulk +
              delta_reaction_constant *
                std::exp(-std::pow(levelset_vec[i], 2) / squared_thickness);
            const double Ceff =
              Ceffective(concentration_vec[i], tracer_reaction_epsilon);
            const double dCe =
              dCeffective_dC(concentration_vec[i], tracer_reaction_epsilon);
            jacobian_vector[i] = k * (tracer_reaction_order - 1.) *
                                 std::pow(Ceff, tracer_reaction_order - 2.) *
                                 dCe;
          }
      }
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0.0);
  }

private:
  const double tracer_reaction_constant_interface;
  const double tracer_reaction_constant_bulk;
  const double delta_reaction_constant;
  const double squared_thickness;
  const double tracer_reaction_order;
  const double tracer_reaction_epsilon;
};

#endif
