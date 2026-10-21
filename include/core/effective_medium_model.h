// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_effective_medium_model_h
#define lethe_effective_medium_model_h

#include <core/parameters_cfd_dem.h>

#include <algorithm>
#include <cmath>
#include <complex>

DeclException1(InvalidPhysicalPropertyReal,
               double,
               << "Invalid physical property value: " << arg1
               << ". The physical property values must be non-negative.");

DeclException1(
  InvalidPhysicalPropertyComplex,
  std::complex<double>,
  << "Invalid physical property value: " << arg1
  << ". The physical property values must be in the upper-half right quadrant of the complex plane.");

template <typename ValueType>
class EffectiveMediumModel
{
public:
  /**
   * @brief Default virtual destructor. Does not do anything.
   */
  virtual ~EffectiveMediumModel() = default;

  /**
   * @brief Creates an instance of the effective medium model based on the parameters.
   * @param[in] cfd_dem_parameters The CFD-DEM parameters.
   * @return A shared pointer to the created effective medium model.
   */
  static std::shared_ptr<EffectiveMediumModel<ValueType>>
  model_cast(const Parameters::CFDDEM &cfd_dem_parameters);

  /**
   * @brief Calculates the value of a physical property for a mixture of two materials.
   * @param[in] void_fraction The void fraction of the mixture (between 0 and
   * 1).
   * @param[in] fluid_property The value of the physical property for the
   * fluid.
   * @param[in] solid_property The value of the physical property for the
   * solid.
   *
   * @return value of the physical property calculated with the void_fraction, fluid_property and solid_property.
   */
  virtual ValueType
  value(const double    void_fraction,
        const ValueType fluid_property,
        const ValueType solid_property) const = 0;

  /**
   * @brief Calculates the value of a physical property for a mixture of two materials at multiple points.
   * @param[in] void_fraction_vector The vector of void fractions of the mixture
   * (between 0 and 1).
   * @param[in] fluid_property_vector The vector of values of the physical
   * property for the fluid.
   * @param[in] solid_property_vector The vector of values of the physical
   * property for the solid.
   * @param[out] property_vector The vector of values of the physical property
   * calculated with the void_fraction_vector, fluid_property_vector and
   * solid_property_vector.
   */
  virtual void
  vector_value(const std::vector<double>    &void_fraction_vector,
               const std::vector<ValueType> &fluid_property_vector,
               const std::vector<ValueType> &solid_property_vector,
               std::vector<ValueType>       &property_vector) const = 0;
};

template <typename ValueType>
class NoEffectiveMediumModel final : public EffectiveMediumModel<ValueType>
{
public:
  /**
   * @brief Calculates the value of a physical property for a mixture of two materials without any effective medium model. The value is simply the fluid property. It is useful for cases where the solid phase properties are considered negligible.
   * @param[in] void_fraction The void fraction of the mixture (between 0 and
   * 1).
   * @param[in] fluid_property The value of the physical property for the fluid.
   * @param[in] solid_property The value of the physical property for the solid.
   * @return value of the physical property calculated with the void_fraction, fluid_property and solid_property. In this case, it simply returns the fluid_property.
   */
  ValueType
  value([[maybe_unused]] const double    void_fraction,
        const ValueType                  fluid_property,
        [[maybe_unused]] const ValueType solid_property) const override
  {
    return fluid_property;
  }

  /**
   * @brief Calculates the value of a physical property for a mixture of two materials at multiple points without any effective medium model. The value is simply the fluid property. It is useful for cases where the solid phase properties are considered negligible.
   * @param[in] void_fraction_vector The vector of void fractions of the mixture
   * (between 0 and 1).
   * @param[in] fluid_property_vector The vector of values of the physical
   * property for the fluid.
   * @param[in] solid_property_vector The vector of values of the physical
   * property for the solid.
   * @param[in, out] property_vector The vector of values of the physical
   * property calculated with the void_fraction_vector, fluid_property_vector
   * and solid_property_vector. In this case, it simply returns the
   * fluid_property_vector.
   */
  void
  vector_value(const std::vector<double>    &void_fraction_vector,
               const std::vector<ValueType> &fluid_property_vector,
               const std::vector<ValueType> &solid_property_vector,
               std::vector<ValueType>       &property_vector) const override
  {
    property_vector = fluid_property_vector;
  }
};

template <typename ValueType>
class BruggemanEffectiveMediumModel final
  : public EffectiveMediumModel<ValueType>
{
public:
  /**
   * @brief Calculates the value of a physical property for a mixture of two materials using the Bruggeman effective medium approximation.
   * @param[in] void_fraction The void fraction of the mixture (between 0 and
   * 1).
   * @param[in] fluid_property The value of the physical property for the fluid.
   * @param[in] solid_property The value of the physical property for the solid.
   * @return value of the physical property calculated with the void_fraction, fluid_property and solid_property.
   */
  ValueType
  value(const double    void_fraction,
        const ValueType fluid_property,
        const ValueType solid_property) const override;

  /**
   * @brief Calculates the value of a physical property for a mixture of two materials at multiple points using the Bruggeman effective medium approximation.
   * @param[in] void_fraction_vector The vector of void fractions of the mixture
   * (between 0 and 1).
   * @param[in] fluid_property_vector The vector of values of the physical
   * property for the fluid.
   * @param[in] solid_property_vector The vector of values of the physical
   * property for the solid.
   * @param[in, out] property_vector The vector of values of the physical
   * property calculated with the void_fraction_vector, fluid_property_vector
   * and solid_property_vector.
   */
  void
  vector_value(const std::vector<double>    &void_fraction_vector,
               const std::vector<ValueType> &fluid_property_vector,
               const std::vector<ValueType> &solid_property_vector,
               std::vector<ValueType>       &property_vector) const override;
};

#endif
