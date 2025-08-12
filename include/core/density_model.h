// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_density_model_h
#define lethe_density_model_h

#include <core/physical_property_model.h>

/**
 * @brief Abstract class for calculating the density of materials.
 */
class DensityModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiate and return a pointer to a DensityModel object by
   * casting it to the proper child class.
   *
   * @param[in] material_properties Property parameters of a material (fluid or
   * solid).
   *
   * @return Casted DensityModel object.
   */
  static std::shared_ptr<DensityModel>
  model_cast(const Parameters::Material &material_properties);


  /**
   * @brief Get the value of the compressibility factor used in the
   * density model.
   *
   * @return Value of the compressibility factor.
   */
  virtual double
  get_psi() const = 0;

  /**
   * @brief Get the value of the reference state density
   * used in the density model.
   *
   * @return Value of the reference state density.
   */
  virtual double
  get_density_ref() const = 0;

  /**
   * @brief Check if the model is a constant density model.
   *
   * @return Boolean indicating if the model corresponds to a constant
   * density model.
   */
  bool
  is_constant_density_model() const
  {
    return constant_density_model;
  }

protected:
  /// Boolean indicating if the model corresponds to a constant density model
  bool constant_density_model = false;
};


/**
 * @brief Constant density model.
 */
class DensityConstant : public DensityModel
{
public:
  /**
   * @brief Constructor of the constant density model.
   *
   * @param[in] p_density Constant density value.
   */
  DensityConstant(const double p_density)
    : density(p_density)
  {
    this->constant_density_model = true;
  }

  /**
   * @brief Compute density value.
   *
   * @param[in] field_values Value of the various fields on which the property
   * may depend. The map stores a single value per field.
   *
   * @note Here, the @p field_values parameter is ignored since the density
   * remains constant.
   *
   * @return Density value.
   */
  double
  value(const std::map<field, double> &field_values) override
  {
    (void)field_values;
    return density;
  }

  /**
   * @brief Compute a vector of density values.
   *
   * @param[in] field_vectors Vectors of the fields on which the density may
   * depend. The map stores a vector of values per field.
   *
   * @param[out] property_vector Vector of computed density values.
   *
   * @note Here, the @p field_vectors parameter is ignored since the density
   * remains constant.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    (void)field_vectors;
    std::fill(property_vector.begin(), property_vector.end(), density);
  }

  /**
   * @brief Compute the jacobian (the partial derivative) of the density with
   * respect to a specified field.
   *
   * @param[in] field_values Values of the various fields on which the density
   * may depend. The map stores a single value per field.
   *
   * @param[in] id Indicator of the field with respect to which the jacobian
   * should be computed.
   *
   * @note Here, the @p field_values and @p id parameters are ignored since the
   * density remains constant.
   *
   * @return Value of the derivative of the density with respect to the
   * specified field. Since the density remains constant, this function returns
   * zero.
   */
  double
  jacobian(const std::map<field, double> &field_values, const field id) override
  {
    (void)field_values;
    (void)id;
    return 0.0;
  }

  /**
   * @brief Computes the derivative of the density with respect to specified
   * fields.
   *
   * @param[in] field_vectors Vector of values of the fields used to evaluate
   * the density. The map stores a vector of values per field.
   *
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be computed.
   *
   * @param[out] jacobian_vector Vector of computed derivative values of the
   * density
   * with respect to the field of the specified @p id. In this case, it returns
   * a vector of zeros since the density remains constant.
   *
   * @note Here, the @p field_vectors and @p id parameters are ignored since the
   * density remains constant.
   *
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    (void)field_vectors;
    (void)id;
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0.0);
  }

  /**
   * @brief Get the value of the compressibility factor used in the density
   * model.
   *
   * @return Value of the compressibility factor which in this case is null.
   */
  double
  get_psi() const override
  {
    return 0.0;
  }

  /**
   * @brief Get the value of the reference state density used in the density
   * model.
   *
   * @return Value of the reference state density, here the constant density
   * value.
   */
  double
  get_density_ref() const override
  {
    return density;
  }

private:
  /// Density of the material.
  const double density;
};


/**
 * @brief Isothermal ideal gas density.
 *
 * Isothermal ideal gas density assumes that the fluid's density varies
 * according the following state equation:
 *
 * \f$\rho = \rho_\text{ref} + \psi p = \rho_\text{ref} + \frac{1}{R T} \ p\f$
 *
 * where \f$\rho_\text{ref}\f$ is the density of the fluid at the reference
 * state, \f$\psi = \frac{1}{R T}\f$ is the compressibility factor derived from
 * the ideal gas law with \f$R= \frac{R_u}{M}\f$ the specific gas constant
 * (universal gas constant (\f$R_u\f$) divided by the molar mass of the gas
 * (\f$M\f$)) and \f$T\f$ the temperature of the gas, finally, \f$p\f$ is the
 * differential pressure between the reference state and the current state.
 * The model is used for weakly compressible flows when temperature
 * fluctuations' influence on density can be neglected.
 */
class DensityIsothermalIdealGas : public DensityModel
{
public:
  /**
   * @brief Constructor of the isothermal ideal gas density model.
   *
   * @param[in] p_density_ref Value of the density of the gas at reference
   * state.
   *
   * @param[in] p_R Value of the specific gas constant \f$\left(R=\frac{R_u}{M}
   * \right)\f$ with \f$R_u\f$ the universal gas constant and \f$M\f$ the molar
   * mass of the gas.
   *
   * @param[in] p_T Value of the temperature of the gas at reference state.
   */
  DensityIsothermalIdealGas(const double p_density_ref,
                            const double p_R,
                            const double p_T)
    : density_ref(p_density_ref)
    , psi(1. / (p_R * p_T))
  {
    this->model_depends_on[field::pressure] = true;
  }

  /**
   * @brief Compute the density of the isothermal ideal gas.
   *
   * @param[in] field_values Values of the various fields on which the property
   * may depend. In this case, the density depends on pressure. The map stores a
   * single value per field.
   *
   * @return Value of the density computed with the @p field_values.
   */
  double
  value(const std::map<field, double> &field_values) override
  {
    Assert(field_values.find(field::pressure) != field_values.end(),
           PhysicialPropertyModelFieldUndefined("DensityIsothermalIdealGas",
                                                "pressure"));
    const double pressure = field_values.at(field::pressure);
    return density_ref + psi * pressure;
  }

  /**
   * @brief Compute a vector of density values for an isothermal ideal gas.
   *
   * @param[in] field_vectors Vectors of the fields on which the density may
   * depend. In this case, the density depends on pressure. The map stores a
   * vector of values per field.
   *
   * @param[out] property_vector Vectors of computed density values.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    Assert(field_vectors.find(field::pressure) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined("DensityIsothermalIdealGas",
                                                "pressure"));
    const std::vector<double> &pressure = field_vectors.at(field::pressure);
    for (unsigned int i = 0; i < property_vector.size(); ++i)
      property_vector[i] = density_ref + psi * pressure[i];
  }

  /**
   * @brief Compute the jacobian (the partial derivative) of the density with
   * respect to a specified field.
   *
   * @param[in] field_values Values of the various fields on which the density
   * may depend. The map stores a single value per field.
   *
   * @param[in] id Indicator of the field with respect to which the jacobian
   * should be computed.
   *
   * @return Value of the derivative of the density with respect to the
   * specified field.
   */
  double
  jacobian(const std::map<field, double> &field_values, const field id) override
  {
    (void)field_values;
    if (id == field::pressure)
      return psi;
    else
      return 0;
  }

  /**
   * @brief Compute the derivative of the density with respect to a field for
   * an isothermal ideal gas.
   *
   * @param[in] field_vectors Vector of values of the fields used to evaluate
   * the density. The map stores a vector of values per field.
   *
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be computed.
   *
   * @param[out] jacobian_vector Vector of computed derivative values of the
   * density
   * with respect to the field of the specified @p id.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    (void)field_vectors;
    if (id == field::pressure)
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), psi);
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

  /**
   * @brief Get the value of the compressibility factor used in the density
   * model.
   *
   * @return Value of the isothermal ideal gas compressibility factor.
   */
  double
  get_psi() const override
  {
    return psi;
  }

  /**
   * @brief Get the value of the reference state density used in the density
   * model.
   *
   * @return Value of the density of the isothermal ideal gas at reference state
   */
  double
  get_density_ref() const override
  {
    return density_ref;
  }

private:
  /// Density of the isothermal ideal gas at reference state
  const double density_ref;
  /// Compressibility factor of the isothermal ideal gas at reference state
  const double psi;
};

#endif
